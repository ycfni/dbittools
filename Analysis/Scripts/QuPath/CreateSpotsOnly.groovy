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

//delete spots and intersections drawn in previous run, if applicable
def annots_to_remove = getAnnotationObjects().findAll { ((it.getPathClass() == SPOT_LABEL) || (it.getPathClass() == INTERSECTION_LABEL))}
removeObjects(annots_to_remove, true)
fireHierarchyUpdate()

//now we're ready to start the intersection identification process



// ------------------------------------------------------------------------------------------------
// Get information about the geometry of the user specified annotation 
// This is expected to be a rectangle (possibly rotated) of class matching DBiT_AREA_LABEL

println("Analyzing user-drawn/labeled DBiT-seq area for subsequent spot placement...")

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

pointClosestToOrigin = indexOfSmallest(dist_orig)

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

println("Starting spot creation...")

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
def Ly = shape_height 


def Ox = points[pointClosestToOrigin].x * cal.getAveragedPixelSize()
def Oy = points[pointClosestToOrigin].y * cal.getAveragedPixelSize()

//---------------------------------------------------------------------------------------
// Compute vertices of boxes
def dx2 = (Lx - (Nboxx * dx)) / (Nboxx - 1)
def dy2 = (Ly - (Nboxy * dy)) / (Nboxy - 1)

def original_x = linspace(Ox, Ox + Lx - dx, Nboxx)
def Px1 = new double[Nboxy * original_x.size()]
for (int i = 0; i < Nboxy; i++) {
    for (int j = 0; j < original_x.size(); j++) {
        Px1[(i*original_x.size()) + j] = original_x[j]
    }
}

def Px2 = new double[Px1.size()]
for (i = 0; i < Px1.size(); i++) {
    Px2[i] = Px1[i] + dx
}


def original_y = linspace(Oy, Oy + Ly - dy, Nboxy)
def Py1 = new double[Nboxx * original_y.size()]
for (int i = 0; i < original_y.size(); i++) {
    for (int j = 0; j < Nboxx; j++) {
        Py1[i*Nboxx + j] = original_y[i]
    }
}

def Py2 = new double[Py1.size()]
for (i = 0; i < Py1.size(); i++) {
    Py2[i] = Py1[i] + dy
}

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
println("done")