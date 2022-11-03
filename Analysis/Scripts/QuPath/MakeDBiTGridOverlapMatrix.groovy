/**
 * modified script from Pete Bankhead
 * to split a rectangle into equal parts along its longest dimension.
 * Written for https://forum.image.sc/t/how-to-divide-annotation-roi-into-equal-dimension-bins/51563/8
 *
 * @author of modications - Sydney Zink
 */
 
 
 /**
 * Usage:  
 * If they don't exist already, create the following annotation classes:
 *    -  'DBiT-seq area'
 *    -  'spot'
 *    -  'strip'
 *
 * Using the QuPath GUI, draw a sqaure annotation (hold SHIFT with rectangle tool).
 * Assign class 'DBiT-seq area' to the square annotation
 * Create other annotations with arbitrary class names, to segment image as you like
 * Run this script
 *
 */
 
// IMPORTS
import org.locationtech.jts.geom.Geometry
import qupath.lib.common.GeneralTools
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import static qupath.lib.gui.scripting.QPEx.*

// PATHS
//optional csv of "Intersection of OBJ with SPOT: X.XX" format to assist review of matx
CSV_PATH = '/Users/mmd47/Dropbox/Workspace/dbittools/Analysis/Scripts/QuPath/tmp.csv'
//final comma-separated matrix in .txt format for downstream ingestion as obs columns
MATX_PATH = '/Users/mmd47/Dropbox/Workspace/dbittools/Analysis/Scripts/QuPath/Annotations.txt'

// # spots desired per dimension in DBiT area square
int SPOTS_PER_SIDE = 50

// CLASS LABELS
//DBiT area should be a square (shift + draw rectangle obj) of any rotation
def DBiT_AREA_LABEL = getPathClass("DBiT-seq area")
//a spot is a smaller square within the DBiT SPOTS_PER_SIDE x SPOTS_PER_SIDE grid
def SPOT_LABEL = getPathClass('spot')
//a strip is DBiT area sliced into SPOTS_PER_SIDE-many rectangles as an intermediate step
//before slicing each of the SPOTS_PER_SIDE-many strips into SPOTS_PER_SIDE-many spots
def STRIP_LABEL = getPathClass('strip')

def img_annotations = getAnnotationObjects()
//script selects the DBiT seq area to split *for* you, but has to be marked as such
def selected = img_annotations.findAll { it.getPathClass() == DBiT_AREA_LABEL }[0]

if (selected == null) {
    println 'No object selected! Did you assign the DBiT-seq area class to your annotation square?'
}

// Get points, removing duplicates
def roi = selected.getROI()
def points = new LinkedHashSet<>(roi.getAllPoints()) as List
if (points.size() != 4) {
    println 'I need a ROI with exactly 4 points'
    return
}

// Get the side lengths
double d1 = points[1].distance(points[0])
double d2 = points[2].distance(points[1])
double d3 = points[3].distance(points[2])
double d4 = points[0].distance(points[3])

double eps = 0.01
if (Math.abs(d1 - d3) > eps || Math.abs(d4 - d2) > eps) {
    println 'Points do not appear to form a rectangle!'
    println '1st pass should be 1 square; 2nd is many rectangles (strips).'
    return
}


// Determination of where to begin slicing
int ind = 0
if (d1 < d2) {
    points.add(0, points.remove(3))
}
double x = points[ind].x
double y = points[ind].y
double dx = (points[ind+1].x - x) / SPOTS_PER_SIDE
double dy = (points[ind+1].y - y) / SPOTS_PER_SIDE
double dx2 = (points[ind+2].x - points[ind+1].x)
double dy2 = (points[ind+2].y - points[ind+1].y)

// Add annotations
def annotations = []
for (int i = 0; i < SPOTS_PER_SIDE; i++) {
    double originX = x + dx*i
    double originY = y + dy*i
    def polygon = ROIs.createPolygonROI(
        [originX, originX+dx, originX+dx+dx2, originX+dx2] as double[],
        [originY, originY+dy, originY+dy+dy2, originY+dy2] as double[],
        roi.getImagePlane()
    )
    def newAnnotation = PathObjects.createAnnotationObject(polygon)
    //our first pass slices the DBiT seq square into SPOTS_PER_SIDE-many strips
    newAnnotation.setPathClass(STRIP_LABEL)
    newAnnotation.setName("Strip ${i+1}")
    annotations << newAnnotation
}

//we need to label and add the strip annotations to the existing annotations
//so we can find them when slicing each strip into spots in the second pass
addObjects(annotations)

def stripCount = 0 //for iterating over each slice out of SPOTS_PER_SIDE many
def smallSquareAnnotations = []

/*
we don't want to keep the strips from the first annotations pass
we just want to split each of those strips (a Tile object) into its own strips
which will give small square spaces (spots), since our starting shape is square
*/

//don't run this second pass over the original DBiT area 
//(nor over irrelevant annotations outside of that area)
def secondary = getAnnotationObjects().findAll { it.getPathClass() == STRIP_LABEL }

for (annotation in secondary){

    def stripROI = annotation.getROI()
    def stripPoints = new LinkedHashSet<>(stripROI.getAllPoints()) as List
    if (stripPoints.size() != 4) {
        println 'I need a ROI with exactly 4 points'
        return
    }
    
    // Get the side lengths
    double stripD1 = stripPoints[1].distance(stripPoints[0])
    double stripD2 = stripPoints[2].distance(stripPoints[1])
    double stripD3 = stripPoints[3].distance(stripPoints[2])
    double stripD4 = stripPoints[0].distance(stripPoints[3])
    
    if (Math.abs(d1 - d3) > eps || Math.abs(d4 - d2) > eps) {
        println 'Points do not appear to form a rectangle!'
        return
    }
    
    // Get starting point based on longest side
    int stripInd = 0
    if (stripD1 < stripD2) {
        //we do want these to be longer strips, subsections of original square
        stripPoints.add(0, stripPoints.remove(3))
    }
    double strip_x = stripPoints[stripInd].x
    double strip_y = stripPoints[stripInd].y
    double strip_dx = (stripPoints[stripInd+1].x - strip_x) / SPOTS_PER_SIDE
    double strip_dy = (stripPoints[stripInd+1].y - strip_y) / SPOTS_PER_SIDE
    double strip_dx2 = (stripPoints[stripInd+2].x - stripPoints[stripInd+1].x)
    double strip_dy2 = (stripPoints[stripInd+2].y - stripPoints[stripInd+1].y)
    
    // Add annotations
    for (int i = 0; i < SPOTS_PER_SIDE; i++) {
        double strip_originX = strip_x + strip_dx*i
        double strip_originY = strip_y + strip_dy*i
        def strip_polygon = ROIs.createPolygonROI(
            [strip_originX, 
            strip_originX+strip_dx, 
            strip_originX+strip_dx+strip_dx2, 
            strip_originX+strip_dx2] as double[],
            [strip_originY, 
            strip_originY+strip_dy, 
            strip_originY+strip_dy+strip_dy2, 
            strip_originY+strip_dy2] as double[],
            stripROI.getImagePlane()
        )
        def strip_newAnnotation = PathObjects.createAnnotationObject(strip_polygon)
        strip_newAnnotation.setPathClass(SPOT_LABEL)
        //name the spot w/ AxB convention for matrix utility & orientation clarity
        strip_newAnnotation.setName("${stripCount} x ${i}")
        smallSquareAnnotations << strip_newAnnotation
    }
    stripCount = stripCount + 1

}

//we no longer need the strip intermediate annotations
removeObjects(annotations, true)
//though we do need to add the spot annotations
addObjects(smallSquareAnnotations)
//and for getting overlaps, we finally filter for non-spot, non-DBiT annotations
nonSpots = getAnnotationObjects().findAll { (it.getPathClass() != DBiT_AREA_LABEL) && (it.getPathClass() != SPOT_LABEL)}


//characterize the annotations throughout the img that the DBiT area might overlap
listOfClasses = [] 
nonSpots.each{ listOfClasses << it.getROI().getRoiName()}
nonspotGeometries = listOfClasses as Set
print(nonspotGeometries) //which types of annotations do we have geometrically?
for (shape in nonspotGeometries) {
    shapeObjs = nonSpots.findAll { it.getROI().getRoiName() == shape}
    def idx = 0
    for (obj in shapeObjs) {
        /*
        give each annotation within a geometric type an ascending ID
        e.g. "Rectangle 0," "Rectangle 1"
        so that we can tell which annotation is overlapping with a spot
        (otherwise each annotation falls under the same "null" identity)
        NOTE the DBiT seq area square is NOT counted among Rectangles
        */
        obj.setName("${shape} ${idx}")
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

//next step is to fill a 2D arr w/ SPOTS_PER_SIDE x # objects' overlap proportions
def retarr = new Float [SPOTS_PER_SIDE**2][nonSpots.size()]
//(optional) create csv of positive intersection values per spot & annotation object
def header = "Intersection ID,Overlap Proportion"
new File(CSV_PATH).withWriter { fw ->
   fw.writeLine(header)
   for (int i = 0; i < SPOTS_PER_SIDE; i++) {
       for (int j = 0; j < SPOTS_PER_SIDE; j++) {
           nonSpots.each() {
               //first object is the spot in question
               ann01 = getAnnotationObjects().find { it.getName() == "${i} x ${j}"}
               //second object is the (non-spot, non-DBiT rectangle) annotation in question
               ann02 = it
                
               if (ann01 != ann02) {
                   ann01_roi = ann01.getROI()
                   ann02_roi = ann02.getROI()
                   ann01_geo = ann01.getROI().getGeometry()
                   ann02_geo = ann02.getROI().getGeometry()
                    
                   def imageData = getCurrentImageData()
                   def hierarchy = imageData.getHierarchy()
                   def server = imageData.getServer()
                    
                   def plane = ann01_roi.getImagePlane()
                    
                   intersect_geo = ann01_geo.intersection(ann02_geo)
                   intersect_roi = GeometryTools.geometryToROI(intersect_geo, plane)
                    
                   def ann01Area = ann01_roi.getArea()
                   def intersectArea = intersect_roi.getArea()
                    
                   def ann02Name = ann02.getName()
                    
                   //store the proportion of overlap to matrix regardless of value
                   def rownum = (i*SPOTS_PER_SIDE) + j
                   def colnum = nonspot_inds[ann02Name]
                   //overlap is simply intersection area over spot area
                   def cellval = intersectArea/ann01Area
                   retarr[rownum][colnum] = cellval
                    
                   //but only create intersection annotation objects if overlap exists
                   if (cellval > 0){
                       intersect_annotation = PathObjects.createAnnotationObject(intersect_roi)
                       intersect_annotation.setName("Intersection of ${ann02Name} with ${i} x ${j}")
                   
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
output_matx.write("") //erase existing contents before writing
for (int i = 0; i < retarr.length; i++) {
    for (int j = 0; j < retarr[i].length; j++) {
        if(j == retarr[i].length-1) {
            output_matx.append(Float.toString(retarr[i][j]) + "\n");
        } 
        else {
            output_matx.append(Float.toString(retarr[i][j]) + ",");
        }
    }
}

println("done")
