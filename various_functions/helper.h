/* Collection of useful functions */

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "../../common/barrel_field.h"

#ifndef HELPER
#define HELPER

class helper
{
	public: 
		// Statistics functions
		static double computeMean(std::list<double> list); 
		static double computeSTD(std::list<double> list); 
		static void computeMeanSTD(std::list< double> list, double& mean, double& std); 
		static void computeMeanSTDMinMax(std::list< double> list, double& mean, double& std, double& minVal, double& maxVal);
		static double computeCorr(std::list< double > list1 , std::list< double > list2); 
		static void computeCumSum(double cumsum[], double inputarray[], int length); 
		static double * randn(double mean, double SD, int numSamples); // Draws one sample from random normal Distribution with mean and SD
		static double * randn(double mean, double SD, int numSamples, gsl_rng * r); // Draws one sample from random normal Distribution with mean and SD
		static double norm(double vector[3]); // Normalize 3D vector
		static double * cross(double u[3], double v[3]); 
		
		// ImageDataVolume Functions
		static ImageDataPointerType correctVolume(ImageDataPointerType volume, double VoxelSize);  // So that extend runs from Minus over 0 to Plus
		static ImageDataPointerType loadVolume(const char* inputFilename);
		static ImageDataPointerType loadVolumeN(const char* inputFilename); 
		static void calculateExtent(double bounds[6], int extent[6],double VoxelSize);
		static ImageDataPointerType createImageVolume(int dims[6], double VoxelSize);
		static ImageDataPointerType createImageVolume(double bounds[6], double VoxelSize);
		static ImageDataPointerType createImageVolume(int dims[6], double VoxelSize, int NumScalarComponents); 
		static ImageDataPointerType createImageVolumeNeuroNet(double bounds[6], double VoxelSize);
		static void storeImageVolume(const char * outputFilename, ImageDataPointerType volume); 
		static bool isValidPt(int x, int y, int z, ImageDataPointerType volume); 
		static ImageDataPointerType drawFromHist(ImageDataPointerType volume, double VoxelSize, int binSz);
		static ImageDataPointerType dimPlusOneVolumne(ImageDataPointerType inputVolume);
		
		// BoundingBox
		static void getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph); 
		static void getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph); 
		static void getBoundingBox(double xBox[2], double yBox[2], double zBox[2], ImageDataPointerType volume, double VoxelSize); 
		static void getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph, int neuriteID);
		static void getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph, int neuriteID);
		
		// Miscellaneous
		static void getProfileAsTxt(ImageDataPointerType volume, const char * path, double VoxelSize);
		static void getProfileAsTxtImprov(ImageDataPointerType volume, ImageDataPointerType volumeOld, const char * path, double VoxelSize); 
		static std::string getRootFromPath(const char* inputFilenameList); 
		static std::string getFilenameFromPath(const char* inputPathFilename);
		static void printArray(double array[], int len); 
		static void convertImageData2txt(ImageDataPointerType volume, const char * outputFilename);
		static void mkDir(const char * path); 
		static std::string exec(char* cmd); 
		static std::vector<std::string> getReturnOfLsCmd(std::string path);
		static std::vector<std::string> returnNames(const char* inputFilenameList, const char* root); 
		static std::vector<unsigned int> returnCellIDs(const char * inputFilenameList);

		static bool isInhibitory(int celltype);
		static void extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath);
		static void extractMorphologies(const char * spatialGraphSetFilename, const char * CellIDsFilename, const char * outputpath);
		static void extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath, double bounds[6]);
		static TransformPointerType amiraToVTKTransform(double* amiraTransform);
		static void writeListToFile(std::list<double>,const char * filename);
		static double computeEuclideanDistance(double pt1[3], double pt2[3]);
		
		// Point Operations
		static bool isPtWithinBounds(double pt[3],double bounds[6]);
		static std::list<double> getIntersectionValue(double bounds[6], double pt1[3], double pt2[3]);
		static double computeIntersectionValue(double Dst1, double Dst2);
		static void getIntersectionPoint(double t, double pt1[3], double pt2[3], double intersectPt[3]);
		static bool isPtOnEdge(double pt[3], double bounds[6]);

		// AmiraSpatialGraph
		static AmiraSpatialGraph * getSpatialGraph(const char * filename); 
		static void zScaling(AmiraSpatialGraph * spatialGraph); 
		static void zScaling(AmiraSpatialGraph * spatialGraph, double thickness);
		static void centerOfSpatialGraph(int neuriteID, double centerPt[3], AmiraSpatialGraph * spatialGraph); 
		static void align(AmiraSpatialGraph * spatialGraph, int neuriteID); 
		static double computeSurfaceArea(AmiraSpatialGraph * spatialGraph, int neuriteID);
		static double computeSurfaceArea(AmiraSpatialGraph * SpatialGraph, int neuriteID, bool likeAmira);
		static std::list< int > getBranchPointIDs(AmiraSpatialGraph * spatialGraph, int neuriteID); 
		static std::map< int, int > getBranchOrder(AmiraSpatialGraph * spatialGraph, int neuriteID); 
		static void flipYSpatialGraph(AmiraSpatialGraph * spatialGraph); 
		static void switchYZ(AmiraSpatialGraph * spatialGraph, int d);
		static void scaleDendrites(AmiraSpatialGraph * spatialGraph, double scaleLat, double scaleVert, bool registrated);
		static double lengthCell(AmiraSpatialGraph * SpatialGraph, int neuriteID);
		static void checkCell(AmiraSpatialGraph * SpatialGraph, int neuriteID);
		static void changeSpatialGraphLabel(AmiraSpatialGraph * spatialGraph,int oldLabel, int newLabel);
};
#endif