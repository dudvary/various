/****************************************************************************/
/*                                                                          */
/* Program:   downSampleSections                                            */
/*                                                                          */
/* File:      downSampleSections.cpp                                        */
/*                                                                          */
/* Purpose:   downsamples images and generated amira						*/
/* 			  denisty file based on pixel values		                    */
/*                                                                          */
/* Author:    Daniel Udvary                                                 */
/****************************************************************************/


#include <execinfo.h>
#include <signal.h>

#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include <set>
#include <utility>

#include "itkAddImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkAffineTransform.h"

// IO
CalcImage2DType::Pointer readImage(const char * filename);
void writeImage(const char * filename, CalcImage2DType::Pointer image);
void writeImage8Bit(const char * filename, Image2DType::Pointer image);
void storeImageVolume(const char * outputFilename, ImageDataPointerType volume);

// Main functions
//void processDownSampling(const char * filepath);
void processDownSampling(const char * filepath, const char * outputpath, int opt);
void processAmiraDensity(const char * filepath);
void getProjectionImages(const char * filepath, const char * outputpath);
void getAmiraDensity(const char * filepath, const char * filepathSpatialGraphs,
						const char * filepathAmiraDensity);
void getAmiraDensity2(const char * filepath, const char * filepathTextFile,
						const char * filepathAmiraDensity);

// Projection Image
CalcImage2DType::Pointer addImages(CalcImage2DType::Pointer image1, CalcImage2DType::Pointer image2);
CalcImage2DType::Pointer maxImages(CalcImage2DType::Pointer image1, CalcImage2DType::Pointer image2);
CalcImage2DType::Pointer multiplyImageWithScalar(CalcImage2DType::Pointer image, double x);
CalcImage2DType::Pointer averageImageWithinSection(const char * filepath, int& numberOfZPlanes);
CalcImage2DType::Pointer maxImageWithinSection(const char * filepath, int& numberOfZPlanes);

void projectionImageWithinSection(const char * filepath, int& numberOfZPlanes,
		CalcImage2DType::Pointer& maxImg, CalcImage2DType::Pointer& avgImg);

// Downsampling
CalcImage2DType::Pointer resizeImage(CalcImage2DType::Pointer image, double x);
CalcImage2DType::Pointer downsampleImage(CalcImage2DType::Pointer image, double x, int opt);
ImageDataPointerType modifyVolumeWithImage(ImageDataPointerType volume,
											CalcImage2DType::Pointer image,
											int currentZidx);

void printUsage();
std::string getFilenameFromPath(const char* inputPathFilename);
TransformPointerType getTransformationMatrixOfSpatialGraph(const char * inputFilename);
void transformPt(TransformPointerType transform, double x, double y,
					double& xT, double& yT, double& zT);
bool isImgNonZero(CalcImage2DType::Pointer inputImg);
float maxValueFromImage(CalcImage2DType::Pointer inputImg);

// System
std::string exec(char * cmd);
std::vector<std::string> getReturnOfLsCmd(std::string path);
void mkDir(const char * path);

// Global variables
void initializeGlobalVariables();
void setImageResXY(double x);
void setSamplingXY(double x);
void handler(int sig);
std::string subFolderNameDwnSmplng;
double samplingXY;
double samplingZ;
double imageResXY;
double imageResZ;
double scalingXY; // scalingFactor -> 1 pixel

typedef std::map< int, std::map< MatrixIndexType, float> > DensityMap;
ImageDataPointerType createImageVolume(int dim[6], std::list<int> depthList,
										DensityMap amiraDensity);
std::list<int> updateDepthList(std::list<int> depthList);

/* Functionality is there,
 * awaiting further instructions for implementation
 * NOTE: ZSampling needs to be taken care of!
 */
int main(int argc , char * argv[])
{
	signal(SIGSEGV, handler);

	initializeGlobalVariables();

	if(argc == 3)
	{
		const char * pathToInputFolder = argv[1];
		const char * outputprefix = argv[2];

		getProjectionImages(pathToInputFolder,outputprefix);
	}
	else if (argc==4)
	{
		const char * filepathImages = argv[1];
		const char * filepathTextFile = argv[2];
		const char * filepathAmiraDensity = argv[3];
		setImageResXY(50.0);
		getAmiraDensity2(filepathImages,filepathTextFile,filepathAmiraDensity);
	}
	else
	{
		printUsage();
	}
}

void printUsage()
{
	std::cout << "USAGE: ./DownSampleSections [path/to/inputfolder] [outputprefix] (writes projection images)" << std::endl;
	std::cout << "USAGE: ./DownSampleSections [path/to/inputfolder] [path/to/csvfile.csv] [outpath/to/Density.am] (writes projection images)" << std::endl;
}

// Extract Transformation Matrix from csv file
// csv delimited by \t [tab]
// SectionID[S0XX]\tT(1)\tT(2)\tT(3)\t...\tT(15)\tT(16)
// last line should be empty
std::map<std::string, TransformPointerType> getTransformationMatrixFromText(const char * filepathTextFile)
{
	std::map<std::string,TransformPointerType> sectionTransformMap;
	std::ifstream inputStream(filepathTextFile);

	if (!inputStream.fail())
	{
		std::string currentLine;

		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				TransformPointerType transform = TransformPointerType::New();
				HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();

				size_t nameStart = 0;
				size_t nameEnd = currentLine.find_first_of("\t", nameStart);
				std::string sectionID(currentLine, nameStart, (nameEnd-nameStart));
				std::string::size_type loc1 = nameEnd+1;
				char * matChar = new char[currentLine.size()-loc1];
				currentLine.copy(matChar, currentLine.size()-loc1, loc1);
				double tmpMat[16];
				char ** endptr = new char*;
				tmpMat[0] = strtod(matChar, endptr);
				for(int ii = 1; ii < 16; ++ii)
					tmpMat[ii] = strtod(*endptr, endptr);

				for(int ii = 0; ii < 4; ++ii)
				{
					for(int jj = 0; jj < 4; ++jj)
					{
						mat->SetElement(jj, ii, tmpMat[ii*4+jj]);
					}
				}
				transform->SetMatrix(mat);
				transform->Update();
				sectionTransformMap[sectionID] = transform;

				//std::cout << sectionID << " " << tmpMat[0] << " " << tmpMat[1] << " ... ";
				//std::cout << tmpMat[14] << " " << tmpMat[15] << std::endl;
			}
		}
	}
	else
	{
		std::cout << "ERROR! Could not open file " << filepathTextFile << std::endl;
	}

	inputStream.close();
	return sectionTransformMap;
}

void getAmiraDensity2(const char * filepathImages, const char * filepathTextFile,
						const char * filepathAmiraDensity)
{
	// Get all the projection views
	std::string pathQuery = std::string(filepathImages) + "*.tif";
	std::vector<std::string> filenamesImages = getReturnOfLsCmd(pathQuery);
	int numberOfImages = filenamesImages.size();

	// Get Transformation Matrices
	std::map< std::string, TransformPointerType > sectionTransformMap =
			getTransformationMatrixFromText(filepathTextFile);
	int numberTransformationMatrices = sectionTransformMap.size();

	std::cout << "-------------------------------------" << std::endl;
	std::cout << "** COMPUTE AMIRA IMAGE DATA VOLUME **" << std::endl;
	std::cout << "  > found " << numberOfImages << " .tif projection images in " << filepathImages << std::endl;
	std::cout << "  > found " << numberTransformationMatrices << " transformation matrices in " << filepathTextFile << std::endl;

	if (numberTransformationMatrices != numberOfImages)
	{
		std::cout << "Number of Images does not match number of SpatialGraphs found!" << std::endl;
		return;
	}
	std::cout << "  Reading in " << numberOfImages << " projection images and their transformation matrices..." << std::endl;
	std::cout << "-------------------------------------" << std::endl;

	// Declare / Initialize variables
	std::vector<std::string>::iterator itFileImgs = filenamesImages.begin();
	int opt = 1; // 1: Average, 2: Max
	int bounds[6] = {1E8,-1E8,1E8,-1E8,1E8,-1E8};
	DensityMap amiraDensity; //std::map< int, std::map< MatrixIndexType, float> >

	int minDepth = 1E8;
	int maxDepth = -1E8;

	// Go through projection images
	for (; itFileImgs != filenamesImages.end(); ++itFileImgs)
	{
		std::cout << "  Reading file " << getFilenameFromPath((*itFileImgs).c_str()) << std::endl;

		// Declare variables
		std::map< MatrixIndexType, float> PixelPos2PixelValue;
		std::map< MatrixIndexType, int> PixelPos2Counter;
		double xT, yT; // Transform pixel position
		double depth = 0.0;

		// Check whether SectionID of Image and SpatialGraph match
		// Here image is S0XX while SpatialGraph is SXX
		std::string secName = getFilenameFromPath((*itFileImgs).c_str()).substr(0,4);

		// Get Transformation matrix
		if (sectionTransformMap.count(secName)==0)
		{
			std::cout << "WARNING! SectionID of Image " << secName;
			std::cout << " could not be found in secionTransformMap" << std::endl;
			continue;
		}

		TransformPointerType transform = sectionTransformMap[secName];
		transformPt(transform,0.0,0.0,xT,yT,depth);

		// Read in image and compute image size
		CalcImage2DType::Pointer tmpImg = readImage((*itFileImgs).c_str());
		CalcImage2DType::SizeType inputsize = tmpImg->GetLargestPossibleRegion().GetSize();

		float summedPixelValue = 0;

		bool WarningMessagePrinted = false;

		// Iterate over image
		for(int ii = 0; ii < inputsize[0]; ++ii)
		{
			for(int jj = 0; jj < inputsize[1]; ++jj)
			{
				CalcImage2DType::IndexType inputIndex;
				inputIndex[0] = ii;
				inputIndex[1] = jj;
				float pixelValue = tmpImg->GetPixel(inputIndex);

				summedPixelValue += pixelValue;

				// Transform pixel position
				double zT;
				transformPt(transform,double(ii)*samplingXY,double(jj)*samplingXY,xT,yT,zT);

				if (int(depth) != int(zT))
				{
					if (!WarningMessagePrinted)
					{
						std::cout << "  WARNING! Depth (z-value) not consistent! " << depth << " " << zT << std::endl;
						WarningMessagePrinted = true;
					}

				}
				// Index of pixel in downsampled image
				int outputIndex[3];
				outputIndex[0] = int(xT/imageResXY);
				outputIndex[1] = int(yT/imageResXY);
				outputIndex[2] = int(zT/imageResZ);
				MatrixIndexType newIndex(outputIndex[0],outputIndex[1]);

				// Check whether index already exists in map
				// If not, add it and update boundingbox
				if (PixelPos2PixelValue.count(newIndex)==0)
				{
					PixelPos2PixelValue[newIndex] = pixelValue;
					PixelPos2Counter[newIndex] = 1;

					// update bounding box
					bounds[0] = std::min(bounds[0],outputIndex[0]);
					bounds[1] = std::max(bounds[1],outputIndex[0]);
					bounds[2] = std::min(bounds[2],outputIndex[1]);
					bounds[3] = std::max(bounds[3],outputIndex[1]);
					bounds[4] = std::min(bounds[4],outputIndex[2]);
					bounds[5] = std::max(bounds[5],outputIndex[2]);

				} // If exists, add pixelValues up
				else
				{
					if (opt==1) // Average
					{
						PixelPos2PixelValue[newIndex] += pixelValue;
						PixelPos2Counter[newIndex] += 1;
					}
					else if (opt==2) // Max
					{
						PixelPos2PixelValue[newIndex] = std::max(pixelValue,PixelPos2PixelValue[newIndex]);
					}
					else
					{
						std::cout << "  ERROR! Only valid options are opt = 1 (avg) or opt = 2 (max)" << std::endl;
					}
				}
			}
		}

		// Delete image from Memory
		tmpImg = NULL;

		// Normalize pixel values if option for averaging is on (opt==1)
		if (opt==1)
		{
			for (std::map<MatrixIndexType,float>::iterator it=PixelPos2PixelValue.begin();
					it!=PixelPos2PixelValue.end(); ++it)
			{
				std::map<MatrixIndexType,int>::iterator itCount = PixelPos2Counter.find(it->first);

				if (itCount != PixelPos2Counter.end())
				{
					it->second = (it->second)/float(itCount->second);
				}
				else
				{
					std::cout << "ERROR! PixelPos not found! [" << it->first.first << " ";
					std::cout << it->first.second << "]" << std::endl;
				}
			}
		}

		// Add map to amiraDensity Map
		amiraDensity[int(depth)] = PixelPos2PixelValue;

		minDepth = std::min(minDepth,(int)depth);
		maxDepth = std::max(maxDepth,(int)depth);

		std::cout << "    Image Size = [" << inputsize[0] << " " << inputsize[1] << "]"<< std::endl;
		std::cout << "    Current BB = [" << bounds[0] << " " << bounds[1] << ", ";
		std::cout << bounds[2] << " " << bounds[3] << ", ";
		std::cout << bounds[4] << " " << bounds[5] << "]" << std::endl;
		std::cout << "    Added z plane at depth = " << depth << " to AmiraDensityMap" << std::endl;
		std::cout << "    AveragePixelValue = " << summedPixelValue/((float)(inputsize[0] * inputsize[1])) << std::endl;
	}

	// Finished reading in images
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "  Finished reading in all projection images" << std::endl;
	std::cout << "  > Final Dim = [" << bounds[0] << " " << bounds[1] << ", ";
	std::cout << bounds[2] << " " << bounds[3] << ", ";
	std::cout << bounds[4] << " " << bounds[5] << "]" << std::endl;
	std::cout << "  > Final Extent = [" << bounds[0]*imageResXY << " " << bounds[1]*imageResXY << ", ";
	std::cout << bounds[2]*imageResXY << " " << bounds[3]*imageResXY << ", ";
	std::cout << bounds[4]*imageResZ << " " << bounds[5]*imageResZ << "]" << std::endl;
	std::cout << "  > Range of in z (Depth) = [" << minDepth << " " << maxDepth << "]" << std::endl;
	std::cout << "  > Resolution/Spacing = [" << imageResXY << " " << imageResXY;
	std::cout << " " << imageResZ << "]" << std::endl;

	// Get depth list
	std::list<int> depthList;
	int z = minDepth;
	while (z<=maxDepth)
	{
		depthList.push_back(z);
		z += imageResZ;
	}

	// Add a zero top and bottom
	int dim[6] = {bounds[0]-1, bounds[1]+1, bounds[2]-1, bounds[3]+1, bounds[4]-1, bounds[5]+2};
	std::cout << "  > AmiraDensityMap has " << amiraDensity.size() << " z planes" << std::endl;
	std::cout << "  > Final Dim = [" << dim[0] << " " << dim[1] << ", ";
	std::cout << dim[2] << " " << dim[3] << ", ";
	std::cout << dim[4] << " " << dim[5] << "]" << std::endl;
	std::cout << "  > Final Extent = [" << dim[0]*imageResXY << " " << dim[1]*imageResXY << ", ";
	std::cout << dim[2]*imageResXY << " " << dim[3]*imageResXY << ", ";
	std::cout << dim[4]*imageResZ << " " << dim[5]*imageResZ << "]" << std::endl;
	std::cout << "  Create Amira Image Data Volume" << std::endl;

	// Create Amira Image Data Volume from amiraDensity
	ImageDataPointerType volume = createImageVolume(dim, depthList, amiraDensity);
	storeImageVolume(filepathAmiraDensity, volume);

	std::cout << "  Wrote Amira Image Data Volume to " << filepathAmiraDensity << std::endl;
	std::cout << "-------------------------------------" << std::endl;
}


void getAmiraDensity(const char * filepathImages, const char * filepathSpatialGraphs,
						const char * filepathAmiraDensity)
{
	// Get all the projection views
	std::string pathQuery = std::string(filepathImages) + "*.tif";
	std::vector<std::string> filenamesImages = getReturnOfLsCmd(pathQuery);
	int numberOfImages = filenamesImages.size();

	// Get all the SpatialGraphs for Transformation Matrix
	pathQuery = std::string(filepathSpatialGraphs) + "*.am";
	std::vector<std::string> filenamesSGs = getReturnOfLsCmd(pathQuery);
	int numberOfSGs = filenamesSGs.size();

	std::cout << "-------------------------------------" << std::endl;
	std::cout << "** COMPUTE AMIRA IMAGE DATA VOLUME **" << std::endl;
	std::cout << "  > found " << numberOfImages << " .tif projection images in " << filepathImages << std::endl;
	std::cout << "  > found " << numberOfSGs << " .am spatial graph files in " << filepathSpatialGraphs << std::endl;

	if (numberOfSGs != numberOfImages)
	{
		std::cout << "Number of Images does not match number of SpatialGraphs found!" << std::endl;
		return;
	}
	std::cout << "  Reading in " << numberOfSGs << " projection images and their transformation matrices..." << std::endl;
	std::cout << "-------------------------------------" << std::endl;

	// Declare / Initialize variables
	std::vector<std::string>::iterator itFileImgs = filenamesImages.begin();
	std::vector<std::string>::iterator itFileSGs = filenamesSGs.begin();
	int opt = 1; // 1: Average, 2: Max
	int bounds[6] = {1E8,-1E8,1E8,-1E8,1E8,-1E8};
	DensityMap amiraDensity; //  std::map< int, std::map< MatrixIndexType, float> >
	std::list<int> depthList;

	// Go through projection images
	for (; itFileImgs != filenamesImages.end(); ++itFileImgs, ++itFileSGs)
	{
		std::cout << "  Reading file " << getFilenameFromPath((*itFileImgs).c_str()) << std::endl;

		// Check whether SectionID of Image and SpatialGraph match
		// Here image is S0XX while SpatialGraph is SXX
		std::string f1 = getFilenameFromPath((*itFileImgs).c_str()).substr(2,2);
		std::string f2 = getFilenameFromPath((*itFileSGs).c_str()).substr(1,2);
		if (f1.compare(f2) != 0)
		{
			std::cout << "WARNING! SectionID of Image " << f1;
			std::cout << " does not match SectionID of SpatialGraph " << f2 << std::endl;
			std::cout << (*itFileImgs) << " vs " << (*itFileSGs) << std::endl;
			continue;
		}

		// Declare variables
		std::map< MatrixIndexType, float> PixelPos2PixelValue;
		std::map< MatrixIndexType, int> PixelPos2Counter;
		double xT, yT; // Transform pixel position
		double depth = 0.0;

		// Get Transformation matrix
		TransformPointerType transform = getTransformationMatrixOfSpatialGraph((*itFileSGs).c_str());

		// Read in image and compute image size
		CalcImage2DType::Pointer tmpImg = readImage((*itFileImgs).c_str());
		CalcImage2DType::SizeType inputsize = tmpImg->GetLargestPossibleRegion().GetSize();

		transformPt(transform,0.0,0.0,xT,yT,depth);
		depthList.push_back(int(depth));

		// Iterate over image
		for(int ii = 0; ii < inputsize[0]; ++ii)
		{
			for(int jj = 0; jj < inputsize[1]; ++jj)
			{
				CalcImage2DType::IndexType inputIndex;
				inputIndex[0] = ii;
				inputIndex[1] = jj;
				float pixelValue = tmpImg->GetPixel(inputIndex);

				// Transform pixel position
				double zT;
				transformPt(transform,double(ii)*samplingXY,double(jj)*samplingXY,xT,yT,zT);

				if (int(depth) != int(zT))
					std::cout << "  WARNING! Depth (z-value) not consistent! " << depth << " " << zT << std::endl;

				// Index of pixel in downsampled image
				int outputIndex[2];
				outputIndex[0] = int(xT/imageResXY);
				outputIndex[1] = int(yT/imageResXY);
				MatrixIndexType newIndex(outputIndex[0],outputIndex[1]);

				// Check whether index already exists in map
				// If not, add it and update boundingbox
				if (PixelPos2PixelValue.count(newIndex)==0)
				{
					PixelPos2PixelValue[newIndex] = pixelValue;
					PixelPos2Counter[newIndex] = 1;

					// update bounding box
					bounds[0] = std::min(bounds[0],outputIndex[0]);
					bounds[1] = std::max(bounds[1],outputIndex[0]);
					bounds[2] = std::min(bounds[2],outputIndex[1]);
					bounds[3] = std::max(bounds[3],outputIndex[1]);
					bounds[4] = std::min(bounds[4],int(depth));
					bounds[5] = std::max(bounds[5],int(depth));

				} // If exists, add pixelValues up
				else
				{
					if (opt==1) // Average
					{
						PixelPos2PixelValue[newIndex] += pixelValue;
						PixelPos2Counter[newIndex] += 1;
					}
					else if (opt==2) // Max
					{
						PixelPos2PixelValue[newIndex] = std::max(pixelValue,PixelPos2PixelValue[newIndex]);
					}
					else
					{
						std::cout << "  ERROR! Only valid options are opt = 1 (avg) or opt = 2 (max)" << std::endl;
					}
				}
			}
		}

		// Delete image from Memory
		tmpImg = NULL;

		// Normalize pixel values if option for averaging is on (opt==1)
		if (opt==1)
		{
			for (std::map<MatrixIndexType,float>::iterator it=PixelPos2PixelValue.begin();
					it!=PixelPos2PixelValue.end(); ++it)
			{
				std::map<MatrixIndexType,int>::iterator itCount = PixelPos2Counter.find(it->first);

				if (itCount != PixelPos2Counter.end())
				{
					it->second = (it->second)/float(itCount->second);
				}
				else
				{
					std::cout << "ERROR! PixelPos not found! [" << it->first.first << " ";
					std::cout << it->first.second << "]" << std::endl;
				}
			}
		}

		// Add map to amiraDensity Map
		amiraDensity[int(depth)] = PixelPos2PixelValue;

		std::cout << "    Current BB = [" << bounds[0] << " " << bounds[1] << ", ";
		std::cout << bounds[2] << " " << bounds[3] << ", ";
		std::cout << bounds[4] << " " << bounds[5] << "]" << std::endl;
		std::cout << "    Added z plane at depth = " << depth << " to AmiraDensityMap" << std::endl;
	}

	// Finished reading in images
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "  Finished reading in all projection images" << std::endl;
	std::cout << "  > Final BB = [" << bounds[0] << " " << bounds[1] << ", ";
	std::cout << bounds[2] << " " << bounds[3] << ", ";
	std::cout << bounds[4] << " " << bounds[5] << "]" << std::endl;
	std::cout << "  > Resolution/Spacing = [" << imageResXY << " " << imageResXY;
	std::cout << " " << imageResZ << "]" << std::endl;

	// Compute min and max difference between sections (sanity check)
	// Add sections in case distance between sections is twice the imageResolution in z
	depthList = updateDepthList(depthList);

	// Add a zero top and bottom
	int dim[6] = {bounds[0]-1, bounds[1]+1, bounds[2]-1, bounds[3]+1, -depthList.size()-1, 0};
	std::cout << "  > AmiraDensityMap has " << amiraDensity.size() << " z planes" << std::endl;
	std::cout << "  > Final Dim = [" << dim[0] << " " << dim[1] << ", ";
	std::cout << dim[2] << " " << dim[3] << ", ";
	std::cout << dim[4] << " " << dim[5] << "]" << std::endl;
	std::cout << "  Create Amira Image Data Volume" << std::endl;

	// Create Amira Image Data Volume from amiraDensity
	ImageDataPointerType volume = createImageVolume(dim, depthList, amiraDensity);
	storeImageVolume(filepathAmiraDensity, volume);

	std::cout << "  Wrote Amira Image Data Volume to " << filepathAmiraDensity << std::endl;
	std::cout << "-------------------------------------" << std::endl;
}

std::list<int> updateDepthList(std::list<int> depthList)
{
	depthList.sort();
	bool updateDepthList = true;
	int depthDist[2] = {1E8,-1E8};

	while (updateDepthList)
	{
		std::list<int> addDepthValues;
		depthDist[0] = 1E8;
		depthDist[1] = -1E8;
		std::list<int>::iterator itDepth = depthList.begin();
		int prevDepth = (*itDepth);
		itDepth++;

		for (; itDepth != depthList.end(); ++itDepth)
		{
			int currDepth = (*itDepth);
			int diff = currDepth-prevDepth;
			depthDist[0] = std::min(depthDist[0],diff);
			depthDist[1] = std::max(depthDist[1],diff);

			if (diff>2*imageResZ)
			{
				int newDepth = prevDepth + int(diff/2);

				std::cout << "  > Add section between depth " << prevDepth;
				std::cout << " and " << currDepth << " (diff = "  << diff;
				std::cout << ") at " << newDepth << std::endl;
				addDepthValues.push_back(newDepth);
			}

			prevDepth = currDepth;
		}

		if (addDepthValues.size()==0)
			updateDepthList = false;

		addDepthValues.sort();
		depthList.merge(addDepthValues);
		depthList.sort();
	}

	std::cout << "  > Range of differences in z between sections = [" << depthDist[0];
	std::cout << "  " << depthDist[1] << "]" << std::endl;
	std::cout << "  > Resolution in z set at " << imageResZ << std::endl;
	return depthList;
}

ImageDataPointerType createImageVolume(int dim[6], std::list<int> depthList, DensityMap amiraDensity)
{
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(imageResXY, imageResXY, imageResZ);
	volume->SetExtent(dim);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Iterate over image volume, sorted by depth
	std::list<int>::iterator itDepth = depthList.begin();

	for (int z = dim[4]; z <= dim[5]; ++z)
	{
		std::map<MatrixIndexType,float> PixelPos2PixelValue;

		// If not top or bottom z plane and depth values still exist
		if (z>dim[4] && z<dim[5] && itDepth != depthList.end())
		{
			if (amiraDensity.count(*itDepth)==0 && itDepth != depthList.begin())
			{
				std::list<int>::iterator itDepthPrev = itDepth;
				--itDepthPrev;
				std::list<int>::iterator itDepthNext = itDepth;
				++itDepthNext;

				std::cout << "  > No section found in depth " << (*itDepth);
				std::cout << " for " << z << " in amiraDensityMap." << std::endl;

				if (itDepthNext != depthList.end() &&
						amiraDensity.count(*itDepthPrev)>0 &&
						amiraDensity.count(*itDepthNext)>0)
				{
					std::cout << "Take average between previous section " << (*itDepthPrev);
					std::cout << " and subsequent section " << (*itDepthNext) << std::endl;

					PixelPos2PixelValue = amiraDensity[(*itDepthPrev)];
					std::map<MatrixIndexType,float> PixelPos2PixelValueNext = amiraDensity[(*itDepthNext)];

					// Get average by dividing PixelPos2PixelValue values by 2
					for (std::map<MatrixIndexType,float>::iterator it = PixelPos2PixelValue.begin();
							it != PixelPos2PixelValue.end(); ++it)
					{
						it->second = it->second/2;
					}

					// Get average dividing PixelPos2PixelValueNext values and adding them to PixelPos2PixelValue values
					for (std::map<MatrixIndexType,float>::iterator itNext = PixelPos2PixelValueNext.begin();
							itNext != PixelPos2PixelValueNext.end(); ++itNext)
					{
						float newValue = (itNext->second)/2;
						std::map<MatrixIndexType,float>::iterator it = PixelPos2PixelValue.find(itNext->first);
						if (it != PixelPos2PixelValue.end())
						{ // If exists, udpate value
							it->second += newValue;
						}
						else
						{ // If not exists, insert new pair
							PixelPos2PixelValue[itNext->first] = newValue;
						}
					}
				}
				else
				{
					std::cout << "No average between " << (*itDepthPrev) << " and " << (*itDepthNext);
					std::cout << " previousSection: " << amiraDensity.count(*itDepthPrev);
					std::cout << " subsequentSection: " << amiraDensity.count(*itDepthNext);
					std::cout << " because either previous or subsequent section does not exist!"<< std::endl;
				}
			}
			else
			{
				PixelPos2PixelValue = amiraDensity[(*itDepth)];
				//std::cout << "  > Section found at depth " << (*itDepth) << " and added." << std::endl;
			}
			itDepth++;
		}

		for(int x = dim[0]; x <= dim[1]; ++x)
		{
		   for(int y = dim[2]; y <= dim[3]; ++y)
		   {
			   double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));

				MatrixIndexType tmp(x,y);

				std::map<MatrixIndexType,float>::iterator it = PixelPos2PixelValue.find(tmp);
				if (it != PixelPos2PixelValue.end())
				{
					*px = it->second;
				}
				else
				{
					*px = 0.0;
				}
		   }
		}
	}

	// Shift Volume up
	double origin[3];
	volume->GetOrigin(origin);
	origin[2] -= imageResZ/2;
	//std::cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
	volume->SetOrigin(origin);

	return volume;
}

// Apply transformation to point ii and jj
// Return transformed ii (=xT), jj (=yT) and zT (assuming that z
void transformPt(TransformPointerType transform, double x, double y,
					double& xT, double& yT, double& zT)
{
	HomogeneousMatrixPointerType mat = transform->GetMatrix();
	//oldCoords[4], newCoords[4];
	double oldCoords[4] = {x, y, 0.0, 1.0};
	double newCoords[4] = {0.0, 0.0, 0.0, 1.0};

	for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			newCoords[ii] += (mat->GetElement(ii, jj))*oldCoords[jj];

	xT = newCoords[0];
	yT = newCoords[1];
	zT = newCoords[2];
}

TransformPointerType getTransformationMatrixOfSpatialGraph(const char * inputFilename)
{
	std::string numbers = "0123456789";
	std::string signs = "+-";
	std::ifstream inputStream(inputFilename);
	bool transformFound = 0;
	TransformPointerType transform = TransformPointerType::New();
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();

	if (!inputStream.fail())
	{
		std::string currentLine;

		while(!std::getline(inputStream, currentLine).eof() && !transformFound)
		{
			if(currentLine.size())
			{
				if(currentLine.find("TransformationMatrix ", 0) != std::string::npos)
				{
					std::string::size_type loc1, loc2;
					loc1 = currentLine.find_first_of(numbers, 0);
					loc2 = currentLine.find_first_of(signs, 0);
					if(loc2 != std::string::npos)
						if(loc2 < loc1)
							loc1 = loc2;
					char * matChar = new char[currentLine.size()-loc1];
					currentLine.copy(matChar, currentLine.size()-loc1, loc1);
					double tmpMat[16];
					char ** endptr = new char*;
					tmpMat[0] = strtod(matChar, endptr);
					for(int ii = 1; ii < 16; ++ii)
						tmpMat[ii] = strtod(*endptr, endptr);

					for(int ii = 0; ii < 4; ++ii)
						for(int jj = 0; jj < 4; ++jj)
							mat->SetElement(jj, ii, tmpMat[ii*4+jj]);

					transform->SetMatrix(mat);
					transform->Update();
					transformFound = 1;
				}

				if (currentLine.find("# Data section follows", 0) != std::string::npos)
				{
					break;
				}
			}
		}
	}
	else
	{
		std::cout << "ERROR! Could not open file " << inputFilename << std::endl;
	}

	inputStream.close();

	if (!transformFound)
	{
		std::cout << "No transformation matrix found in " << inputFilename << std::endl;
		std::cout << "Set to Identity matrix" << std::endl;

		for(int ii = 0; ii < 4; ++ii)
		{
			for(int jj = 0; jj < 4; ++jj)
			{
				double tmp = (ii == jj) ? 1 : 0;
				mat->SetElement(ii, jj, tmp);
			}
		}
		transform->SetMatrix(mat);
		transform->Update();
	}

	return transform;
}

void getProjectionImages(const char * filepath, const char * outputprefix)
{
	std::string filepathStr(filepath);
	int numberOfZPlanes;

	// Get Maximum and Average Projection Images
	// Save them to outputpath folder
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "Reading in .tif files in folder " << filepath << std::endl;

	CalcImage2DType::Pointer maxImg;
	CalcImage2DType::Pointer avgImg;
	projectionImageWithinSection(filepath, numberOfZPlanes, maxImg, avgImg);

	std::cout << "  Found " << numberOfZPlanes << " z planes (#images)" << std::endl;

	std::string outputFilenameMax = std::string(outputprefix) + "projectionImgMax.tif";
	writeImage(outputFilenameMax.c_str(), maxImg);
	std::cout << "Wrote " << outputFilenameMax << std::endl;

	std::string outputFilenameAvg = std::string(outputprefix) + "projectionImgAvg.tif";
	writeImage(outputFilenameAvg.c_str(), avgImg);
	std::cout << "Wrote " << outputFilenameAvg << std::endl;
	std::cout << "-------------------------------------" << std::endl;
}

void projectionImageWithinSection(const char * filepath, int& numberOfZPlanes,
		CalcImage2DType::Pointer& maxImg, CalcImage2DType::Pointer& avgImg)
{
	std::string mainpath(filepath);
	std::string pathQuery = mainpath + "*.tif";
	std::vector<std::string> filenames = getReturnOfLsCmd(pathQuery);

	numberOfZPlanes = filenames.size();

	std::vector<std::string>::iterator iteratorFilename = filenames.begin();
	maxImg = readImage((*iteratorFilename).c_str());
	avgImg = readImage((*iteratorFilename).c_str());
	iteratorFilename++;

	// Add up images
	for (; iteratorFilename != filenames.end(); ++iteratorFilename)
	{
		CalcImage2DType::Pointer tmpImg = readImage((*iteratorFilename).c_str());
		maxImg = maxImages(maxImg,tmpImg);
		avgImg = addImages(avgImg,tmpImg);
	}

	// Normalize pixel value by number of images to get average image (projection view)
	double norm = 1.0/numberOfZPlanes;
	avgImg = multiplyImageWithScalar(avgImg,norm);
}

std::string getFilenameFromPath(const char* inputPathFilename)
{
	// Save truncated Spatial Graph
	std::string tmpStr(inputPathFilename);
	unsigned found = tmpStr.find_last_of("/");
	std::string Filename = tmpStr.substr(found+1);
	return Filename;
}


/* Main function to create amira densities based on stack of images
 * - loads in images in folder one by one
 * - add pixel values of each image in the respective z plane to a
 * 		image data volume
 * - save image data volume to hard disk */
void processAmiraDensity(const char * filepath)
{
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "Create Amira Density of .tif files in folder";
	std::cout << filepath << std::endl;

	// Get all tif images from folder
	std::string filepathStr(filepath);
	std::string pathQuery = filepathStr + "*.tif";
	std::vector<std::string> filenames = getReturnOfLsCmd(pathQuery);
	int numOfImgs = filenames.size();
	std::cout << "  Found " << numOfImgs << " .tif images" << std::endl;

	// Load first image
	std::vector<std::string>::iterator iteratorFilename = filenames.begin();
	std::cout << "Reading " << getFilenameFromPath((*iteratorFilename).c_str()) << std::endl;
	CalcImage2DType::Pointer img1 = readImage((*iteratorFilename).c_str());
	iteratorFilename++;

	// Get dimensions for ImageDataVolume
	CalcImage2DType::SizeType size = img1->GetLargestPossibleRegion().GetSize();
	int dims[6];
	dims[0] = 0;
	dims[1] = size[0]-1;
	dims[2] = 0;
	dims[3] = size[1]-1;
	dims[4] = 0;
	dims[5] = numOfImgs-1;

	// Create ImageDataVolume
	std::cout << "  Create Amira Image Data Volume of dimensions [" << size[0];
	std::cout << " " << size[1] << " " << numOfImgs << "]" << std::endl;

	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(imageResXY, imageResXY, imageResZ);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Do first image
	int currentZidx = 0;
	volume = modifyVolumeWithImage(volume,img1,currentZidx);

	std::cout << "-------------------------------------" << std::endl;
	std::cout << "  Reading in remaining " << numOfImgs-1 << " images ..." << std::endl;

	// Load remaining images from folder
	for (; iteratorFilename != filenames.end(); ++iteratorFilename, ++currentZidx)
	{
		std::cout << "Reading " << getFilenameFromPath((*iteratorFilename).c_str()) << std::endl;
		img1 = readImage((*iteratorFilename).c_str());
		volume = modifyVolumeWithImage(volume,img1,currentZidx);
	}
	std::cout << "  Reading finished" << std::endl;
	std::cout << "-------------------------------------" << std::endl;

	// Save final density as amira file
	std::string outputfilename = filepathStr + "density.am";
	storeImageVolume(outputfilename.c_str(),volume);
	std::cout << "Wrote Amira Image Data Volume " << outputfilename << std::endl;
	std::cout << "-------------------------------------" << std::endl;
}

// Add pixel values from image to image data volume at respective z-level
// Assumes that volume and image are of same size
ImageDataPointerType modifyVolumeWithImage(ImageDataPointerType volume,
											CalcImage2DType::Pointer image,
											int currentZidx)
{
	CalcImage2DType::SizeType size = image->GetLargestPossibleRegion().GetSize();

	// Iterate over image volume
	for(int x = 0; x <= size[0]-1; ++x)
	{
	   for(int y = 0; y <= size[1]-1; ++y)
	   {
			double * px = static_cast< double * >(volume->GetScalarPointer(x, y, currentZidx));

			CalcImage2DType::IndexType inputIndex;
			inputIndex[0] = x;
			inputIndex[1] = y;
			float pixelValue = image->GetPixel(inputIndex);
			*px = (double)pixelValue;
	   }
	}
	return volume;
}

/* Save image data volume to hard disk given outputfilename */
void storeImageVolume(const char * outputFilename, ImageDataPointerType volume)
{
	std::string fname(outputFilename);
	if (fname.compare(fname.size()-3,3,".am") == 0)
		fname = fname.substr(0,fname.size()-3);

	Reader * Writer = new Reader(fname.c_str(), fname.c_str());
	Writer->writeScalarField(volume);
	delete Writer;
}

void initializeGlobalVariables()
{
	// Variables
	// True VALUES for microscope
	samplingXY = 0.092; // um
	samplingZ = 0.5; // um
	imageResXY = 50.0;
	imageResZ = 50.0;
	scalingXY = round(imageResXY / samplingXY); // scalingFactor -> 1 pixel
	subFolderNameDwnSmplng = "projectionImage";
}

void setImageResXY(double x)
{
	imageResXY = x;
	scalingXY = round(imageResXY / samplingXY);
}

void setSamplingXY(double x)
{
	samplingXY = x;
	scalingXY = round(imageResXY / samplingXY);
}

void processDownSampling(const char * filepath, const char * outputpath, int opt)
{
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "Image Settings: " << std::endl;
	std::cout << "  XY-Sampling: " << samplingXY << std::endl;
	std::cout << "  Z-Sampling: " << samplingZ << std::endl;
	std::cout << "  XY-Resolution of final image: " << imageResXY << std::endl;
	std::cout << "  XY-DownScalingFactor: " << scalingXY << std::endl;

	std::cout << "-------------------------------------" << std::endl;
	std::cout << "Resizing and downsampling of projection images" << std::endl;

	// Read in all images in folder and downsample them
	std::string pathQuery = std::string(filepath) + "*.tif";
	std::vector<std::string> filenames = getReturnOfLsCmd(pathQuery);
	int numberOfProjectionImages = filenames.size();

	std::cout << "Found " << numberOfProjectionImages << " Projection Images in ";
	std::cout << filepath << std::endl;
	std::cout << "-------------------------------------" << std::endl;

	// Convert ImageRes into String
	std::stringstream num2str;
	num2str << imageResXY;

	// Downsample individual images
	for (std::vector<std::string>::iterator iteratorFilename = filenames.begin();
			iteratorFilename != filenames.end(); ++iteratorFilename)
	{
		CalcImage2DType::Pointer projectionImg = readImage((*iteratorFilename).c_str());

		CalcImage2DType::Pointer resizedImg = resizeImage(projectionImg, scalingXY);
		projectionImg = NULL;

		CalcImage2DType::Pointer downsampledImg = downsampleImage(resizedImg, scalingXY, opt);
		resizedImg = NULL;

		std::string base = std::string(outputpath) + getFilenameFromPath((*iteratorFilename).c_str());
		std::string outputFilename = base.substr(0,base.size()-4) +
												"_Downsample_x" + num2str.str() + ".tif";

		writeImage(outputFilename.c_str(),downsampledImg);
		std::cout << "  wrote " << outputFilename << std::endl;
		std::cout << "-------------------------------------" << std::endl;
	}

	std::cout << "Finished resizing and downsampling" << std::endl;
}

bool isImgNonZero(CalcImage2DType::Pointer inputImg)
{
	Calc2DIteratorType it(inputImg, inputImg->GetLargestPossibleRegion());
	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		float pixelValue = it.Get();
		if (pixelValue>0)
			return true;
		++it;
	}
	return false;
}

float maxValueFromImage(CalcImage2DType::Pointer inputImg)
{
	float maxValue = 0;
	Calc2DIteratorType it(inputImg, inputImg->GetLargestPossibleRegion());
	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		float pixelValue = it.Get();
		maxValue = std::max(pixelValue,maxValue);
		++it;
	}
	return maxValue;
}


/* Downsample image uniformly in xy given Scaling factor x
 * Takes average across x pixels in both dimensions and writes them in
 * new image. */
CalcImage2DType::Pointer downsampleImage(CalcImage2DType::Pointer inputImg, double x, int opt)
{
	CalcImage2DType::SizeType inputsize = inputImg->GetLargestPossibleRegion().GetSize();
	CalcImage2DType::SizeType outputsize;
	for (int i = 0; i<2; i++)
	{
		outputsize[i] = int(std::floor(inputsize[i] / x)); // New image size
	}

	std::cout << "  downsampleImage: Resize image from " << inputsize << " to " << outputsize;
	std::cout << "  (XYScalingFactor = " << x << ")" << std::endl;

	// Create new image in downsampled size
	CalcImage2DType::RegionType outputRegion;
	outputRegion.SetSize(outputsize);
	CalcImage2DType::Pointer outputImg = CalcImage2DType::New();
	outputImg->SetRegions(outputRegion);
	outputImg->Allocate();
	outputImg->FillBuffer(0); // Image filled with zeros

	// Iterate over image
	Calc2DIteratorType it(inputImg, inputImg->GetLargestPossibleRegion());
	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		CalcImage2DType::IndexType inputIndex = it.GetIndex();

		// Index of pixel in downsampled image
		CalcImage2DType::IndexType outputIndex;

		for (int i = 0; i<2; i++)
		{
			outputIndex[i] = int(std::floor(inputIndex[i]/x));

			if (outputIndex[i]<0 || outputIndex[i]>outputsize[i])
			{
				std::cout << "ERROR! " << i << " = " << outputIndex[i];
				std::cout << " " << outputsize[i] << std::endl;
				return outputImg;
			}
		}

		float pixelValue = it.Get();

		// Set new Pixel value in downsampled image (already normalized!)
		float newPixelValue = -1;

		if (opt==1) // Average
		{
			newPixelValue = pixelValue/(pow(x,2.0)) + outputImg->GetPixel(outputIndex);
		}
		else if (opt==2) // Maximum
		{
			newPixelValue = std::max(pixelValue,outputImg->GetPixel(outputIndex));
		}
		else
		{
			std::cout << "ERROR! Only valid options are opt = 1 (avg) or opt = 2 (max)" << std::endl;
			return outputImg;
		}

		outputImg->SetPixel(outputIndex,newPixelValue);
		++it;
	}

	outputImg->Update();

	return outputImg;
}

/* Resize Image to appropiate dimensions given down scaling factor
 * adds if necessary pixels (value of white (255)) in x and/or y.
 * 		new_size = ceil(size / scalingFacor) * scalingFactor */
CalcImage2DType::Pointer resizeImage(CalcImage2DType::Pointer image, double x)
{
	CalcImage2DType::SizeType inputsize = image->GetLargestPossibleRegion().GetSize();

	// Calculates new image size
	CalcImage2DType::SizeType outputsize;
	CalcImage2DType::IndexType destinationIndex;
	for (int i = 0; i<2; i++)
	{
		int r = (int)(std::ceil(inputsize[i] / x));
		outputsize[i] = r * x; // New image size
		// calculate starting index for merging of images
		destinationIndex[i] = (outputsize[i]-inputsize[i])/2;
	}

	// Paste smaller image in new image
	// Check if resizing necessary
	if (destinationIndex[0]>0 || destinationIndex[1]>0)
	{
		// Create new (bigger) image
		CalcImage2DType::RegionType scaledRegion;
		scaledRegion.SetSize(outputsize);
		CalcImage2DType::Pointer imageNew = CalcImage2DType::New();
		imageNew->SetRegions(scaledRegion);
		imageNew->Allocate();
		imageNew->FillBuffer(0);

		// Paste smaller image into bigger one
		typedef itk::PasteImageFilter <CalcImage2DType, CalcImage2DType> PasteImageFilterType;
		PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
		pasteFilter->SetSourceImage(image);
		pasteFilter->SetDestinationImage(imageNew);
		pasteFilter->SetSourceRegion(image->GetLargestPossibleRegion());
		pasteFilter->SetDestinationIndex(destinationIndex);
		pasteFilter->Update();

		std::cout << "  resizeImage: Resize image from " << inputsize << " to " << outputsize << std::endl;

		return pasteFilter->GetOutput();
	}
	else // if no resizing necessary, return original image
	{
		return image;
	}
}

/* Go over all tif images in folder and return average image
 * also returns the number of Z planes (number over images read in and taken average of) */
CalcImage2DType::Pointer averageImageWithinSection(const char * filepath, int& numberOfZPlanes)
{
	std::string mainpath(filepath);
	std::string pathQuery = mainpath + "*.tif";
	std::vector<std::string> filenames = getReturnOfLsCmd(pathQuery);

	numberOfZPlanes = filenames.size();

	std::vector<std::string>::iterator iteratorFilename = filenames.begin();
	CalcImage2DType::Pointer avgImage = readImage((*iteratorFilename).c_str());
	iteratorFilename++;

	// Add up images
	for ( ; iteratorFilename != filenames.end(); ++iteratorFilename)
	{
		CalcImage2DType::Pointer tmpImg = readImage((*iteratorFilename).c_str());
		avgImage = addImages(avgImage,tmpImg);
	}

	// Normalize pixel value by number of images to get average image (projection view)
	double norm = 1.0/numberOfZPlanes;
	return multiplyImageWithScalar(avgImage,norm);
}

/* Go over all tif images in folder and return projection image with max values
 * also returns the number of Z planes */
CalcImage2DType::Pointer maxImageWithinSection(const char * filepath, int& numberOfZPlanes)
{
	std::string mainpath(filepath);
	std::string pathQuery = mainpath + "*.tif";
	std::vector<std::string> filenames = getReturnOfLsCmd(pathQuery);

	numberOfZPlanes = filenames.size();

	std::vector<std::string>::iterator iteratorFilename = filenames.begin();
	CalcImage2DType::Pointer maxImage = readImage((*iteratorFilename).c_str());
	iteratorFilename++;

	// Add up images
	for ( ; iteratorFilename != filenames.end(); ++iteratorFilename)
	{
		CalcImage2DType::Pointer tmpImg = readImage((*iteratorFilename).c_str());
		maxImage = maxImages(maxImage,tmpImg);
	}

	return maxImage;
}

/* Reads in image as double value */
CalcImage2DType::Pointer readImage(const char * filename)
{
	typedef itk::ImageFileReader< CalcImage2DType > CalcImage2DReaderType;

	CalcImage2DReaderType::Pointer imgReader = CalcImage2DReaderType::New();
	imgReader->SetFileName(filename);
	imgReader->Update();
	return imgReader->GetOutput();
}

/* Multiply pixels in image with double value */
CalcImage2DType::Pointer multiplyImageWithScalar(CalcImage2DType::Pointer image, double x)
{
	Calc2DIteratorType imageIterator(image,image->GetLargestPossibleRegion());
	imageIterator.GoToBegin();
	while(!imageIterator.IsAtEnd())
	{
		imageIterator.Set(imageIterator.Get() * x);
		++imageIterator;
	}
	return image;
}

/* Adds two CalcImage2DImages */
CalcImage2DType::Pointer addImages(CalcImage2DType::Pointer image1, CalcImage2DType::Pointer image2)
{
	typedef itk::AddImageFilter <CalcImage2DType, CalcImage2DType >  AddImageFilterType;
	AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

	addFilter->SetInput1(image1);
	addFilter->SetInput2(image2);
	addFilter->Update();
	return addFilter->GetOutput();
}

/* Maximum between two images */
CalcImage2DType::Pointer maxImages(CalcImage2DType::Pointer image1, CalcImage2DType::Pointer image2)
{
	typedef itk::MaximumImageFilter <CalcImage2DType, CalcImage2DType >  MaximumImageFilterType;
	MaximumImageFilterType::Pointer maxFilter = MaximumImageFilterType::New();
	maxFilter->SetInput1(image1);
	maxFilter->SetInput2(image2);
	maxFilter->Update();
	return maxFilter->GetOutput();
}

/* Write CalcImage2D as image (8bit image file). First converts it, then writes it */
void writeImage(const char * filename, CalcImage2DType::Pointer image)
{
	typedef itk::CastImageFilter< CalcImage2DType, Image2DType > CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(image);
	writeImage8Bit(filename,castFilter->GetOutput());
}

/* Write Image file (8bit)*/
void writeImage8Bit(const char * filename, Image2DType::Pointer image)
{
	Image2DWriterType::Pointer imgWriter = Image2DWriterType::New();
	imgWriter->SetFileName(filename);
	imgWriter->SetInput(image);
	imgWriter->Update();
}

/* executes command line */
std::string exec(char * cmd)
{
	FILE* pipe = popen(cmd, "r");

	if (!pipe)
		return "ERROR";

	char buffer[128];
	std::string result = "";

	while(!feof(pipe))
	{
		if(fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}

	pclose(pipe);
	return result;
}

/* Returns the output of the ls command to terminal
 * for example the list of all files that match criteria */
std::vector<std::string> getReturnOfLsCmd(std::string path)
{
	std::string cmd =  "ls " + path;
	char * cstr = &cmd[0u];
	std::string txt = exec(cstr);

	std::vector<std::string> filenames;

    std::istringstream iss(txt);
    do
    {
        std::string subs;
        iss >> subs;

        if (subs.size()>0)
        	filenames.push_back(subs);
    } while (iss);
    return filenames;
}

/* Creates Directory / Folder in case it does not exist already */
void mkDir(const char * path)
{
	std::string fname(path);
	std::string mkdir_cmd = "mkdir " + fname;

	if (access(fname.c_str(), F_OK) == -1)
		system(mkdir_cmd.c_str());
}

void handler(int sig)
{

	void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
