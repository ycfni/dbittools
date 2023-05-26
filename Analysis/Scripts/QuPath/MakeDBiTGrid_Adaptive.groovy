// IMPORTS
import org.locationtech.jts.geom.Geometry
import qupath.lib.common.GeneralTools
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
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

/*
import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.JSONParser;
*/
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

import org.json.*;

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




filePath = "/Users/sydneyzink/Desktop/yale/projects/repos/dbittools/Analysis/config.json"
File file = new File(filePath);
InputStream in = new FileInputStream(file)
String contents = IOUtils.toString(in, StandardCharsets.UTF_8);
//URI uri = new URI("file://home/mmd47/Example.json");
JSONTokener tokener = new JSONTokener(contents);
JSONObject root = new JSONObject(tokener);
//println(root.channel_size_vert)

/*
// THESE VARIABLES SHOULD BE READ FROM CONFIG FILE
// # spots desired per dimension in DBiT area square
int nChannels_horiz = 10
int nChannels_vert = 8
int channel_size_horiz = 17
int channel_size_vert = 5
*/

int nChannels_horiz = root.nChannels_horiz as int //50
int nChannels_vert = root.nChannels_vert as int //50
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

//delete spots and intersections drawn in previous run, if applicable
def annots_to_remove = getAnnotationObjects().findAll { ((it.getPathClass() == SPOT_LABEL) || (it.getPathClass() == INTERSECTION_LABEL))}
removeObjects(annots_to_remove, true)
fireHierarchyUpdate()

//now we're ready to start the intersection identification process



// ------------------------------------------------------------------------------------------------
// Get information about the geometry of the user specified annotation 
// This is expected to be a rectangle (possibly rotated) of class matching DBiT_AREA_LABEL
//

def img_annotations = getAnnotationObjects()
//script selects the DBiT seq area to split *for* you, but has to be marked as such
def selected = img_annotations.findAll { it.getPathClass() == DBiT_AREA_LABEL }[0]
if (selected == null) {
    println 'No object selected! Did you assign the DBiT-seq area class to your annotation square?'
    return(0)
}

// Get points, removing duplicates
def roi = selected.getROI()


// Attempt conversion to calibrated units
def cal = getCurrentServer().getPixelCalibration()
if (cal.pixelWidth != cal.pixelHeight) {
    println "Pixel width != pixel height ($cal.pixelWidth vs. $cal.pixelHeight)"
    println "Distance measurements will be calibrated using the average of these"
}

// NOTE: points is in pixels
def points = new LinkedHashSet<>(roi.getAllPoints()) as List
if (points.size() != 4) {
    println 'I need a ROI with exactly 4 points'
    return
}

public static int indexOfSmallest(ArrayList array){
    if (array.size() == 0)
        return -1;
    int index = 0;
    int min = array[index];
    for (int i = 1; i < array.size(); i++){
        if (array[i] <= min){
        min = array[i];
        index = i;
        }
    }
    return index;
}

dist_orig = []
for(var i = 0; i < points.size(); i++){
   a = Math.pow(points.x[i],2) + Math.pow(points.y[i],2)
   dist_orig.add(Math.sqrt(a))
}

//print(points)
//print(dist_orig)

pointClosestToOrigin = indexOfSmallest(dist_orig)


//print("Point Closest to Origin")
//print(points[pointClosestToOrigin].getX() * cal.getAveragedPixelSize())
//print(points[pointClosestToOrigin].getY() * cal.getAveragedPixelSize())
//print("Next Point")
//print(points[(pointClosestToOrigin+1)%4].getX() * cal.getAveragedPixelSize())
//print(points[(pointClosestToOrigin+1)%4].getY() * cal.getAveragedPixelSize())


s1y = points[(pointClosestToOrigin+1)%4].getY() - points[pointClosestToOrigin].getY() 
s1x = points[(pointClosestToOrigin+1)%4].getX() - points[pointClosestToOrigin].getX()
s2y = points[(pointClosestToOrigin+1)%4].getY() - points[(pointClosestToOrigin+2)%4].getY()
s2x = points[(pointClosestToOrigin+1)%4].getX() - points[(pointClosestToOrigin+2)%4].getX()
angle_s1 = Math.atan(s1y/s1x)
angle_s2 = -Math.atan(s2y/s2x)


shape_height = points[pointClosestToOrigin].distance(points[(pointClosestToOrigin+1)%4])
shape_width =  points[(pointClosestToOrigin+1)%4].distance(points[(pointClosestToOrigin+2)%4])
if (angle_s1 == 0.0) {
    angle = angle_s1
    shape_width = points[pointClosestToOrigin].distance(points[(pointClosestToOrigin+1)%4])
    shape_height =  points[(pointClosestToOrigin+1)%4].distance(points[(pointClosestToOrigin+2)%4])
} else if (angle_s1 < angle_s2) {
    angle = angle_s1 + Math.PI/2
} else {    
    angle = -angle_s2 
} 
shape_width  = shape_width * cal.getAveragedPixelSize()
shape_height = shape_height * cal.getAveragedPixelSize()


// ------------------------------------------------------------------------------------------------
// Use the genometry information from the DBiT_AREA_LABEL annotation along with chip configuration 
// information in the config file to populate the rectangular spots
//

public static double[] linspace(double min, double max, int points) {  
    double[] d = new double[points];  
    for (int i = 0; i < points; i++){      
        def temp =  min + i * (max - min) / (points - 1);   
        d[i] = temp 
    }  
    return d;  
} 

//From config
def Nboxx = nChannels_horiz
def Nboxy = nChannels_vert

def dx = channel_size_horiz 
def dy = channel_size_vert 

// From drawing
def Lx = shape_width 
//println("Lx")
//println(Lx)
def Ly = shape_height 
//println("Ly")
//println(Ly)


def Ox = points[pointClosestToOrigin].x * cal.getAveragedPixelSize()
def Oy = points[pointClosestToOrigin].y * cal.getAveragedPixelSize()
//print("Oy")
//print(Oy)

//---------------------------------------------------------------------------------------
// Compute vertices of boxes
def dx2 = (Lx - (Nboxx * dx)) / (Nboxx - 1)
def dy2 = (Ly - (Nboxy * dy)) / (Nboxy - 1)

def original_x = linspace(Ox, Ox + Lx - dx, Nboxx)
//print("original_x")
//print(original_x)
def Px1 = new double[Nboxy * original_x.size()]
for (int i = 0; i < Nboxy; i++) {
    for (int j = 0; j < original_x.size(); j++) {
        Px1[(i*original_x.size()) + j] = original_x[j]
    }
}
//println("Px1")
//println(Px1)
def Px2 = new double[Px1.size()]
for (i = 0; i < Px1.size(); i++) {
    Px2[i] = Px1[i] + dx
}
//println("Px2")
//print(Px2)


def original_y = linspace(Oy, Oy + Ly - dy, Nboxy)
def repeated = new double[original_y.size() * Nboxx]
//println("original y - ")
//print(original_y)
def Py1 = new double[Nboxx * original_y.size()]
for (int i = 0; i < original_y.size(); i++) {
    for (int j = 0; j < Nboxx; j++) {
        Py1[i*Nboxx + j] = original_y[i]
    }
}
//println("Py1")
//println(Py1)

def Py2 = new double[Py1.size()]
for (i = 0; i < Py1.size(); i++) {
    Py2[i] = Py1[i] + dy
}
//println("py2")
//print(Py2)

def P1 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [Px1[it], Py1[it]] }.toArray()
def P2 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [Px2[it], Py1[it]] }.toArray()
def P3 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [Px1[it], Py2[it]] }.toArray()
def P4 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [Px2[it], Py2[it]] }.toArray()

P1no = P1.stream().map(i -> [i[0] - Ox, i[1] - Oy]).collect(Collectors.toList());
P2no = P2.stream().map(i -> [i[0] - Ox, i[1] - Oy]).collect(Collectors.toList());
P3no = P3.stream().map(i -> [i[0] - Ox, i[1] - Oy]).collect(Collectors.toList());
P4no = P4.stream().map(i -> [i[0] - Ox, i[1] - Oy]).collect(Collectors.toList());


theta = angle

def c = Math.cos(theta)
def s = Math.sin(theta)

def PPx1 = (P1no.collect { (it[0] * c - it[1] * s) + Ox } ).toArray()
def PPy1 = (P1no.collect { (it[0] * s + it[1] * c) + Oy } ).toArray()
def PPx2 = (P2no.collect { (it[0] * c - it[1] * s) + Ox } ).toArray()
def PPy2 = (P2no.collect { (it[0] * s + it[1] * c) + Oy } ).toArray()
def PPx3 = (P3no.collect { (it[0] * c - it[1] * s) + Ox } ).toArray()
def PPy3 = (P3no.collect { (it[0] * s + it[1] * c) + Oy } ).toArray()
def PPx4 = (P4no.collect { (it[0] * c - it[1] * s) + Ox } ).toArray()
def PPy4 = (P4no.collect { (it[0] * s + it[1] * c) + Oy } ).toArray()

def PP1 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [PPx1[it], PPy1[it]] }.toArray()
def PP2 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [PPx2[it], PPy2[it]] }.toArray()
def PP3 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [PPx3[it], PPy3[it]] }.toArray()
def PP4 = IntStream.range(0, Nboxx * Nboxy).mapToObj { [PPx4[it], PPy4[it]] }.toArray()

def patches = []
def ind_horiz = 1
def ind_vert = 1
for (a in [PP1, PP3, PP4, PP2].transpose().collect { [it[0], it[1], it[2], it[3]] }) {

    w1 = a[0][0] / cal.getAveragedPixelSize()
    w2 = a[1][0] / cal.getAveragedPixelSize()
    w3 = a[2][0] / cal.getAveragedPixelSize()
    w4 = a[3][0] / cal.getAveragedPixelSize()
    w5 = a[0][1] / cal.getAveragedPixelSize()
    w6 = a[1][1] / cal.getAveragedPixelSize()
    w7 = a[2][1] / cal.getAveragedPixelSize()
    w8 = a[3][1] / cal.getAveragedPixelSize()

    def spotpoly = ROIs.createPolygonROI(
        [w1, w2, w3, w4] as double[],
        [w5, w6, w7, w8] as double[],
        roi.getImagePlane()

    )
    def spotannot = PathObjects.createAnnotationObject(spotpoly)
    spotannot.setPathClass(SPOT_LABEL)
    //name the spot w/ AxB convention for matrix utility & orientation clarity
   
    spotannot.setName("${ind_horiz} x ${ind_vert}")
    
    patches << spotannot

    ind_vert += 1
    if (ind_vert % nChannels_horiz == 1) {
      ind_vert = 1
      ind_horiz += 1
    }
    
    
}

addObjects(patches)

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

def i_chan = nChannels_horiz
def j_chan = nChannels_vert

flip_flag = true
if (flip_flag) {
    i_chan = nChannels_vert
    j_chan = nChannels_horiz
}


def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
new File(CSV_PATH).withWriter { fw ->
   fw.writeLine(csv_header)
   for (int i = 1; i <= i_chan; i++) {
       for (int j = 1; j <= j_chan; j++) {
           
           //matx_rownames << "${i-1}x${j-1}"
           //print("matx rownames")
           //print(matx_rownames)
           ann01 = getAnnotationObjects().find { it.getName() == "${i} x ${j}"}
           ann01_roi = ann01.getROI()
           def plane = ann01_roi.getImagePlane()
           ann01_geo = ann01.getROI().getGeometry()
           def ann01Area = ann01_roi.getArea()
           
           
           nonSpots.each() {
               
             //matx_rownames << "${i-1}x${j-1}"
             //print("matx rownames")
             //print(matx_rownames)
               //first object is the spot in question
               //second object is the (non-spot, non-DBiT rectangle) annotation in question
               ann02 = it
                
               if (ann01 != ann02) {
                    matx_rownames << "${i-1}x${j-1}"
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
//googled around, is there really no "else if" to use??
//and using a switch runs into issues with multiple-cond checks (multiple satisfied)\
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