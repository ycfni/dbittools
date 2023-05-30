// IMPORTS
import org.locationtech.jts.geom.Geometry
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.PolygonROI
import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.common.ColorTools
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.scripting.QP

import org.apache.commons.math3.util.FastMath
import groovy.transform.*
import org.apache.commons.math3.util.MathArrays
import java.util.stream.Collectors
import java.util.stream.IntStream

import qupath.lib.common.Timeit;
import qupath.lib.common.Version;
import qupath.lib.common.GeneralTools;
import qupath.lib.analysis.images.ContourTracing;
import qupath.bioimageio.spec.BioimageIoSpec;
import qupath.lib.objects.classes.PathClass;
import qupath.lib.regions.ImagePlane;
import qupath.lib.scripting.ScriptAttributes;
import qupath.opencv.dnn.DnnModels;
import qupath.lib.regions.ImageRegion;
import qupath.lib.objects.hierarchy.PathObjectHierarchy;
import qupath.imagej.tools.IJTools;
import qupath.lib.roi.ROIs;
import qupath.lib.objects.PathObjectFilter;
import qupath.lib.io.PointIO;
import qupath.lib.images.writers.ImageWriterTools;
import qupath.lib.images.servers.ImageServer;
import qupath.opencv.ml.pixel.PixelClassifierTools;
import qupath.lib.io.GsonTools;
import qupath.lib.common.ColorTools;
import qupath.lib.objects.classes.PathClassTools;
import qupath.lib.roi.GeometryTools;
import qupath.lib.projects.ProjectIO;
import qupath.lib.objects.PathObjectPredicates;
import qupath.lib.analysis.heatmaps.DensityMaps;
import qupath.lib.projects.Projects;
import qupath.lib.analysis.DelaunayTools;
import qupath.lib.images.ImageData;
import qupath.lib.scripting.QP;
import qupath.lib.objects.PathObjectTools;
import qupath.lib.images.servers.PixelType;
import qupath.lib.images.writers.TileExporter;
import qupath.lib.images.servers.ColorTransforms;
import qupath.opencv.tools.OpenCVTools;
import qupath.opencv.tools.GroovyCV;
import qupath.opencv.ops.ImageOps;
import qupath.lib.awt.common.AffineTransforms;
import qupath.lib.objects.CellTools;
import qupath.lib.io.PathIO;
import qupath.lib.objects.PathObject;
import qupath.lib.regions.Padding;
import qupath.opencv.dnn.DnnTools;
import qupath.lib.images.servers.ServerTools;
import java.awt.image.BufferedImage;
import qupath.lib.awt.common.BufferedImageTools;
import qupath.lib.analysis.DistanceTools;
import qupath.lib.io.UriUpdater;
import qupath.opencv.ml.BioimageIoTools;
import qupath.lib.images.servers.LabeledImageServer;
import qupath.lib.regions.RegionRequest;
import qupath.lib.objects.PathObjects;
import qupath.opencv.dnn.DnnModelParams;
import qupath.lib.objects.classes.PathClassFactory;
import qupath.lib.roi.RoiTools;
import qupath.opencv.tools.NumpyTools;
import qupath.lib.gui.QuPathGUI;
import qupath.lib.gui.dialogs.Dialogs;
import qupath.lib.gui.tools.GuiTools;
import qupath.lib.gui.charts.Charts;
import qupath.lib.gui.tools.MenuTools;
import qupath.lib.gui.tools.PaneTools;
import qupath.lib.gui.prefs.PathPrefs;
import qupath.lib.gui.logging.LogManager;
import javafx.application.Platform;
import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.InputStream;

import java.io.*
import org.util.*
import org.json.*

import org.apache.commons.io.IOUtils;
import java.nio.charset.StandardCharsets;

// PATHS
def ImageFileName = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())
pathOutput = buildFilePath(PROJECT_BASE_DIR, ImageFileName, 'Annotations')
mkdirs(pathOutput)
println("Output will be saved to $pathOutput")

// -------------------------------------------------------------
// Main OUTPUT file:
MATX_PATH = buildFilePath(pathOutput, 'intersections_matx.txt')

//optional csv of "Intersection of OBJ with SPOT: X.XX" format to assist review of matx
CSV_PATH = buildFilePath(pathOutput, 'intersections.csv')

// JSON for config file to read measurement settings from that determine spots drawn
filePath = "/Users/sydneyzink/Desktop/yale/projects/repos/dbittools/Analysis/config.json"


File file = new File(filePath);
InputStream in = new FileInputStream(file)
String contents = IOUtils.toString(in, StandardCharsets.UTF_8);
JSONTokener tokener = new JSONTokener(contents);
JSONObject root = new JSONObject(tokener);


// THESE VARIABLES SHOULD BE READ FROM CONFIG FILE
// # spots desired per dimension in DBiT area
int nChannels_horiz = root.nChannels_horiz as int //50
int nChannels_vert = root.nChannels_vert as int //50
// channel measurements in micrometers
double channel_size_horiz = root.channel_size_horiz as double //25
double channel_size_vert = root.channel_size_vert as double //25


// CLASS LABELS
//DBiT area should be a square (shift + draw rectangle obj) of any rotation
def DBiT_AREA_LABEL = getPathClass("DBiT-seq area")
//a spot is a smaller square within the DBiT SPOTS_PER_SIDE x SPOTS_PER_SIDE grid
def SPOT_LABEL = getPathClass('Spot')
//an intersection is an area of overlap between an annotation and a spot
def INTERSECTION_LABEL = getPathClass("Intersection")


//make sure to populate existing/default classes 1st so they aren't deleted
def newClasses = []
def pathClasses = getQuPath().getAvailablePathClasses()
newClasses.addAll(pathClasses)

if (pathClasses.contains(getPathClass('Spot')) == false) {
    println(pathClasses)
    newClasses.addAll([SPOT_LABEL, INTERSECTION_LABEL])
    newClasses = newClasses as Set
    //overwrite the available classes with your new list
    pathClasses.setAll(newClasses)
}

//now we're ready to start the intersection identification process

nonSpots = getAnnotationObjects().findAll { (it.getPathClass() != DBiT_AREA_LABEL) && (it.getPathClass() != SPOT_LABEL)}

//characterize the annotations throughout the img that the DBiT area might overlap
listOfClasses = [] 
nonSpots.each{ listOfClasses << it.getPathClass()}
nonspotClasses = listOfClasses as Set
print(nonspotClasses) //which types of annotations do we have geometrically?
for (classname in nonspotClasses) {
    classObjs = nonSpots.findAll { it.getPathClass() == classname}
    def idx = 0
    for (classobj in classObjs) {
        /*
        give each annotation within a geometric type an ascending ID
        e.g. "Rectangle 0," "Rectangle 1"
        so that we can tell which annotation is overlapping with a spot
        (otherwise each annotation falls under the same "null" identity)
        NOTE the DBiT seq area square is NOT counted among Rectangles
        */
        classobj.setName("${classname}_${idx}")
        idx = idx + 1
    }
}


/*
we need to store annotation objects' column IDs for the matrix's reference
so we create an Annotation Name: Column # dictionary, printed to console
*/
nonspot_inds = [:]
objColIdx = 0
for (obj in nonSpots) {
    nonspot_inds[obj.getName()] = objColIdx
    objColIdx = objColIdx + 1
}
println("The matrix columns correspond to the annotation objects as follows:")
println(nonspot_inds)
nonspot_names = nonspot_inds.keySet()

//next step is to fill a 2D arr w/ SPOTS_PER_SIDE x # objects' overlap proportions
def retarr = new Float [nChannels_horiz * nChannels_vert][nonSpots.size()]
def matx_rownames = []
//(optional) create csv of positive intersection values per spot & annotation object
def csv_header = "Intersection ID,Overlap Proportion"

def i_chan = nChannels_vert
def j_chan = nChannels_horiz

flip_flag = false
if (flip_flag) {
    i_chan = nChannels_horiz
    j_chan = nChannels_vert
}


def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
new File(CSV_PATH).withWriter { fw ->
   fw.writeLine(csv_header)
   for (int i = 1; i <= i_chan; i++) {
       for (int j = 1; j <= j_chan; j++) {
           
           ann01 = getAnnotationObjects().find { it.getName() == "${i} x ${j}"}
           ann01_roi = ann01.getROI()
           def plane = ann01_roi.getImagePlane()
           ann01_geo = ann01.getROI().getGeometry()
           def ann01Area = ann01_roi.getArea()
           
           
           nonSpots.each() {
               
               //first object is the spot in question
               //second object is the (non-spot, non-DBiT rectangle) annotation in question
               ann02 = it
                
               if (ann01 != ann02) {
                   matx_rownames << "${i}x${j}"
                   ann02_roi = ann02.getROI()
                   ann02_geo = ann02.getROI().getGeometry()
                   
                   intersect_geo = ann01_geo.intersection(ann02_geo)
                   intersect_roi = GeometryTools.geometryToROI(intersect_geo, plane)
                    
                   def intersectArea = intersect_roi.getArea()
                    
                   def ann02Name = ann02.getName()
                    
                   //store the proportion of overlap to matrix regardless of value
                   def rownum = ((i-1)*j_chan) + (j-1)
                   def colnum = nonspot_inds[ann02Name]
                   //overlap is simply intersection area over spot area
                   def cellval = intersectArea/ann01Area
                   retarr[rownum][colnum] = cellval
                    
                   //but only create intersection annotation objects if overlap exists
                   if (cellval > 0){
                       intersect_annotation = PathObjects.createAnnotationObject(intersect_roi)
                       intersect_annotation.setName("Intersection of ${ann02Name} with ${i} x ${j}")
                       intersect_annotation.setPathClass(INTERSECTION_LABEL)
                   
                       hierarchy.addPathObject(intersect_annotation) 
                        
                       //next 2 lines help to review intersection proportions per object & spot
                       //the csv otherwise exists for that reason as complement to the matrix txt file
                       //println(intersect_annotation.getName())
                       //println(intersectArea/ann01Area)
                       
                       fw.writeLine(intersect_annotation.getName() + "," + cellval + ",")
                   }    
               }
           }
       }
   }
}

def output_matx = new File(MATX_PATH)
def matx_header = nonspot_names.join(",")
output_matx.write("Row_ID," + matx_header + "\n") //use write erase existing contents before appending
def matx_rowname_idx = 0
for (int i = 0; i < retarr.length; i++) {
    for (int j = 0; j < retarr[i].length; j++) {
        if ((j == retarr[i].length-1) && (j == 0)) {
            output_matx.append(matx_rownames[matx_rowname_idx]  + "," + Float.toString(retarr[i][j]) + "\n");
        }
        else {
            if (j == retarr[i].length-1) {
                output_matx.append(Float.toString(retarr[i][j]) + "\n");
            }
            else {
                if (j == 0) {
                    output_matx.append(matx_rownames[matx_rowname_idx]  + "," + Float.toString(retarr[i][j]) + ",");
                }
                else {
                    output_matx.append(Float.toString(retarr[i][j]) + ",");
                }
            }
        }
        matx_rowname_idx = matx_rowname_idx + 1
  
    }
}

println("Wrote $MATX_PATH")

println("done")