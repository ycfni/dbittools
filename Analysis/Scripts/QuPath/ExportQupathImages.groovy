/* 
a script to export 3 TIF image categories: (1) the DBiT-seq region, (2) the full annotated image (3) the 2500 spots (tiles)
this script is intended to follow use of the MakeDBiTGridOverlapMatrix.groovy script generating said spots (among other entities) 
*/

// SETUP
def imageData = getCurrentImageData()
def imageName = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())

// Define output path (relative to project)
def outputDir = buildFilePath(PROJECT_BASE_DIR, 'export')
mkdirs(outputDir)
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def path = buildFilePath(outputDir, name + "-labels.tif")

// Define how much to downsample during export (may be required for large images)
double downsample = 8

/*

PART ONE:
Save the DBiT-seq area as a TIF file cropped square region of the full image
(not necessarily smaller; may encapsulate the full image, if the DBiT-seq area annotation is drawn as such)

make a directory called export
store as a TIF image the smallest bounding box that fits the DBiT-seq region
i.e., crop to the DBiT-seq square itself (a user-defined square cropped from the full image)
*/

def img_annotations = getAnnotationObjects()
def DBiT_AREA_LABEL = getPathClass("DBiT-seq area")
def dbit_region = img_annotations.findAll { it.getPathClass() == DBiT_AREA_LABEL }[0]

dbit_region.eachWithIndex{anno,x->
    //Name of the file and the path to where it goes in the Project
    //fileName = outputDir+"//"+imageName+"-"+anno.getPathClass()+"-"+x+".tif"
    fileName = outputDir+"//"+imageName+"_"+anno.getPathClass()+".tif"
    //For each annotation, we get its outline
    def roi = anno.getROI()
    //For each outline, we request the pixels within the bounding box of the annotation
    def requestROI = RegionRequest.createInstance(getCurrentServer().getPath(), 1, roi)
    //The 1 in the function above is the downsample, increase it for smaller images
    writeImageRegion(getCurrentServer(), requestROI, fileName)
}

/*

PART TWO: 
Save the 2500 spots as individual TIF image files

make a subdir of the export dir called spots
store each of the 2500 spot tiles as a small TIF image
*/

def SPOT_LABEL = getPathClass('Spot')
def spots = img_annotations.findAll { it.getPathClass() == SPOT_LABEL }

def spotsOutputDir = buildFilePath(outputDir, 'spots')
mkdirs(spotsOutputDir)

spots.eachWithIndex{anno,x->
    //Name of the file and the path to where it goes in the Project
    fileName = spotsOutputDir+"//"+imageName+"_"+anno.getName()+".tif"
    //For each annotation, we get its outline
    def roi = anno.getROI()
    //For each outline, we request the pixels within the bounding box of the annotation
    def requestROI = RegionRequest.createInstance(getCurrentServer().getPath(), 1, roi)
    //The 1 in the function above is the downsample, increase it for smaller images
    writeImageRegion(getCurrentServer(), requestROI, fileName)
}

/*

PART THREE:
Save the full image with all annotations marked

Ensure that QuPath is toggled to show/hide respective classes you want seen in the final image,
as well as toggled accordingly to show/not show names of said entities.
Write the full image, displaying objects according to how they are currently shown in the viewer
i.e., store the original image with all the annotations, spot squares, etc. drawn on it
if you are showing names/labels in the QuPath viewer, they will be written to the img, else not
*/

def viewer = getCurrentViewer()
path = buildFilePath(outputDir, name + "_viewerpane.tif")
writeRenderedImage(viewer, path)

println("done")