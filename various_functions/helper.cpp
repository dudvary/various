#include "helper.h"

/* Load ImageDataVolume */
ImageDataPointerType helper::loadVolume(const char* inputFilename)
{
	std::string fname(inputFilename);
	if (fname.compare(fname.size()-3,3,".am") != 0)
		fname += ".am";
    
	if((access(fname.c_str(), F_OK) != -1))
	{
		Reader * Reader1 = new Reader(fname.c_str(),fname.c_str()); 
		return Reader1->readScalarField(); 
	}
	else
	{
		std::cout << "Error! Could not load " << fname.c_str() << "!" << std::endl;
	}
	return 0; 
} 

/* Loads ImageDataVolume Histograms (more than one scalar value at each position) */
ImageDataPointerType helper::loadVolumeN(const char* inputFilename)
{
	std::string fname(inputFilename);
	if (fname.compare(fname.size()-3,3,".am") != 0)
		fname += ".am";
  
	if((access(fname.c_str(), F_OK) != -1))
	{
		Reader * Reader1 = new Reader(fname.c_str(),fname.c_str()); 
		return Reader1->readNVectorField(); 
	}
	else
	{
		std::cout << "Error! Could not load " << fname << "!" << std::endl;
	}
	return 0; 
} 

/* Computes Mean of list */
double helper::computeMean(std::list< double > list)
{
	std::list<double>::iterator it = list.begin(); 
	
	double mean = 0.0;
	
	if (list.size()>0)
	{
		for (; it!=list.end(); it++)
		{
			mean += (*it);
		}
		mean /= list.size();
	}
	else
	{
		std::cout << "Error! List has no elements! Mean cannot be computed!" << std::endl; 
	}
	return mean; 
}

/* Computes SD of list */
double helper::computeSTD(std::list< double > list)
{
	std::list<double>::iterator it = list.begin(); 
	
	double mean = computeMean(list); 
	double std = 0.0; 
	
	if (list.size()>0)
	{
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);
		}
		std = sqrt(std/list.size());
	}
	else
	{
		std::cout << "Error! List has no elements! STD cannot be computed!" << std::endl; 
	}
	return std; 
}

/* Computes Mean and SD of list */
void helper::computeMeanSTD(std::list< double> list, double& mean, double& std)
{
	if (list.size()>0)
	{
		mean = computeMean(list); 
		std = 0.0; 
		std::list<double>::iterator it = list.begin(); 
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);
		}
		std = sqrt(std/list.size()); 
	}
	else
	{
		std::cout << "Error! List has no elements! Mean and STD cannot be computed!" << std::endl; 
	}
}

void helper::computeMeanSTDMinMax(std::list< double> list, double& mean, double& std, double& minVal, double& maxVal)
{
	if (list.size()>0)
	{
		mean = computeMean(list);
		std = 0.0;
		minVal = std::numeric_limits<double>::infinity();
		maxVal = -std::numeric_limits<double>::infinity();

		std::list<double>::iterator it = list.begin();
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);

			if (minVal>(*it))
				minVal = (*it);

			if (maxVal<(*it))
				maxVal = (*it);
		}
		std = sqrt(std/list.size());
	}
	else
	{
		std::cout << "Error! List has no elements! Mean and STD cannot be computed!" << std::endl;
	}
}


/* Computes Pearson's Correlation between list1 and list2 */
double helper::computeCorr(std::list< double > list1, std::list< double > list2)
{
	double corrcoef = 0.0; 
	
	if (list1.size()==list2.size() && list1.size()>0)
	{
		double m1, m2, std1, std2;  
		computeMeanSTD(list1,m1,std1); 
		computeMeanSTD(list2,m2,std2); 
		
		std::list<double>::iterator it1 = list1.begin(); 
		std::list<double>::iterator it2 = list2.begin(); 
			
		for (; it1!=list1.end();it1++,it2++)
		{
			corrcoef += ((*it1)-m1)*((*it2)-m2); 
		}
		
		corrcoef /= list1.size(); 
		if (std1<1E-9 || std2<1E-9)
		{
			//std::cout << " Warning! STD is 0.0! Corrcoeff is set to 0.0!" << std::endl;
			corrcoef = 0.0;  
		}
		else
		{
			corrcoef /= (std1*std2);  
		}
	}
	else
	{
		std::cout << " Error! Correlation of two lists with different number of elements cannot be computed!" << std::endl; 
		std::cout << " 		list1.size() = " << list1.size() << " list2.size() = " << list2.size() << std::endl;
	}
	
	return corrcoef; 
}

/* computes Cumulative Sum of inputarray of length */
void helper::computeCumSum(double cumsum[], double inputarray[], int length)
{ 
	for(int i = 0; i<length; i++)
	{
		cumsum[i] = inputarray[i];
		if (i>0)
			cumsum[i] += cumsum[i-1]; 
	}
}

/* Stores ImageDataVolume as .txt file */
void helper::convertImageData2txt(ImageDataPointerType volume, const char* outputFilename)
{
	std::ofstream outStream(outputFilename);
	
	if(!outStream.fail())
	{	  
		int extent[6];
		double CoordOrigin[3]; 
		volume->GetOrigin(CoordOrigin);
		volume->GetExtent(extent);
		
		outStream << extent[0] << "," << extent[1] << "," << extent[2] << "," << extent[3] << "," << extent[4] << "," << extent[5] << std::endl;
		outStream << CoordOrigin[0] << "," << CoordOrigin[1] << "," << CoordOrigin[2] << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{  
			for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
			{
				for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
				{ 
					double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
					outStream << (*px);
					
					if (z==volume->GetExtent()[5])
					{
						outStream << "; ";
					}
					else
					{
						outStream << " "; 
					}
				}
			}
			outStream << "" << std::endl; 
		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  outputFilename << std::endl; 
	}
}

/* Corrects ImageDataVolume by setting Origin and extent correctly so that values run from - to + */
ImageDataPointerType helper::correctVolume(ImageDataPointerType volume, double VoxelSize)
{
	double OriginCoord[3];
	int extent[6];
	volume->GetOrigin(OriginCoord);
	volume->GetExtent(extent);
	
	extent[0] += OriginCoord[0]/VoxelSize;
	extent[1] += OriginCoord[0]/VoxelSize;
	extent[2] += OriginCoord[1]/VoxelSize;
	extent[3] += OriginCoord[1]/VoxelSize;
	extent[4] += OriginCoord[2]/VoxelSize;
	extent[5] += OriginCoord[2]/VoxelSize;

	volume->SetOrigin(0,0,0);
	volume->SetExtent(extent);
	
	return volume; 
}

/* Create Image Data Volume (VTK), Input: BoundingBox Coordinates, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(double bounds[6], double VoxelSize)
{
	double spacing[3];
	int dims[6];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;
	
	calculateExtent(bounds, dims, VoxelSize);
	
	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	{
	   for(int y = dims[2]; y <= dims[3]; ++y)
	   {
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	   }
	}
	volume->Update();
	return volume;
}

/* Create Image Data Volume (VTK), Input: BoundingBox Coordinates, Output: ImageDataVolume
 * NeuroNet: Voxels run from [-50 to 0; 0 to 50], origin is [-50], extends from [0 1] */
ImageDataPointerType helper::createImageVolumeNeuroNet(double bounds[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	int extent[6];
	extent[0] = 0;
	extent[1] = (bounds[1] - bounds[0])/spacing[0];
	extent[2] = 0;
	extent[3] = (bounds[3] - bounds[2])/spacing[1];
	extent[4] = 0;
	extent[5] = (bounds[5] - bounds[4])/spacing[2];

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetOrigin(bounds[0]+VoxelSize/2,bounds[2]+VoxelSize/2,bounds[4]+VoxelSize/2);
	volume->SetExtent(extent);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Initialize each Scalar with 0.0
	for(int z = extent[4]; z <= extent[5]; ++z)
	{
	   for(int y = extent[2]; y <= extent[3]; ++y)
	   {
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	   }
	}
	volume->Update();
	return volume;
}

ImageDataPointerType helper::dimPlusOneVolumne(ImageDataPointerType inputVolume)
{
	// Get Spacing
	double Spacing[3];
	inputVolume->GetSpacing(Spacing[0],Spacing[1],Spacing[2]);

	// Correct Spacing
	double maxSpacing = 0;
	if (Spacing[0]>0 && std::isnan(Spacing[1]) && std::isnan(Spacing[2]))
		maxSpacing = Spacing[0];
	else if (std::isnan(Spacing[0]) && Spacing[1]>0 && std::isnan(Spacing[2]))
		maxSpacing = Spacing[1];
	else if (std::isnan(Spacing[0]) && std::isnan(Spacing[1]) && Spacing[2]>0)
		maxSpacing = Spacing[2];
	else
	{
		std::cout << "ERROR in helper::dimPlusOneVolumne!" << std::endl;
		std::cout << "  Spacing = [" << Spacing[0] << "," << Spacing[1] << "," << Spacing[2] << "]" << std::endl;
		return inputVolume;
	}

	if (maxSpacing==0)
	{
		std::cout << "ERROR in helper::dimPlusOneVolumne! MaxSpacing is zero!" << std::endl;
		std::cout << "  Spacing = [" << Spacing[0] << "," << Spacing[1] << "," << Spacing[2] << "]" << std::endl;
		return inputVolume;
	}

	Spacing[0] = maxSpacing;
	Spacing[1] = maxSpacing;
	Spacing[2] = maxSpacing;

	// Correct Extent by adding one
	int extent[6];
	inputVolume->GetExtent(extent[0],extent[1],extent[2],extent[3],extent[4],extent[5]);
	extent[1] = extent[1] + 1;
	extent[3] = extent[3] + 1;
	extent[5] = extent[5] + 1;

	// Extract Origin
	double origin[3];
	inputVolume->GetOrigin(origin[0],origin[1],origin[2]);

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetExtent(extent);
	volume->SetSpacing(Spacing[0],Spacing[1],Spacing[2]);
	volume->SetOrigin(origin[0],origin[1],origin[2]);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Copy values into new density
	for(int z = extent[4]; z <= extent[5]; ++z)
	{
	   for(int y = extent[2]; y <= extent[3]; ++y)
	   {
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));

				// Copy values over from inputVolume
				if (z<extent[5] && y<extent[3] && x<extent[1])
				{
					double * px2 = static_cast< double * >(inputVolume->GetScalarPointer(x, y, z));
					*px = (*px2);
				}
				else // set new values in added dimensions to zero
				{
					*px = 0.0;
				}
			}
	   }
	}
	volume->Update();
	return volume;
}

/* Create Image Data Volume (VTK), Input: Dimensions, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(int dims[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	volume->Update();
	return volume;  
}

/* Create Image Data Volume (VTK), Input: Dimensions, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(int dims[6], double VoxelSize, int NumScalarComponents)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(NumScalarComponents);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				
				for(int i = 0; i<NumScalarComponents; i++)
					px[i] = 0;
			}
	volume->Update();
	return volume;  
}

/* Draw values from Histogram ImageDataVolume (if binsz set to 0, assume integer values) 
 * returns drawn values as ImageDataVolume
 */
ImageDataPointerType helper::drawFromHist(ImageDataPointerType volume, double VoxelSize, int binSz)
{
	int dim[6]; 
	double OriginCoord[3];
	volume->GetExtent(dim); 
	volume->GetOrigin(OriginCoord);
	// Create new ImageDataVolume
	ImageDataPointerType volumeReturn = createImageVolume(dim,VoxelSize); 
	volumeReturn->SetOrigin(OriginCoord[0],OriginCoord[1],OriginCoord[2]);
	int NumScalarComponents = volume->GetNumberOfScalarComponents(); 
	
	srand(time(NULL));
	int maxValue; 
	int randVal;
	
// 	volume->Print(std::cout); 
// 	volumeReturn->Print(std::cout); 
	
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{    
				double *pxHist = static_cast< double * >(volume->GetScalarPointer(x,y,z)); 	// Histogram
				double *px = static_cast< double * >(volumeReturn->GetScalarPointer(x,y,z)); 
				
				// Compute cumulative sum for histogram
				double cumsum[NumScalarComponents]; 
				computeCumSum(cumsum, pxHist, NumScalarComponents); 
				
				// Draw random number between 1 and maxValue
				maxValue = cumsum[NumScalarComponents-1]; 
				int i = 0; 
				
				if (maxValue>0) // otherwise length is set to zero
				{
					randVal = rand() % maxValue + 1;
					
					while(randVal>cumsum[i] && i<NumScalarComponents)	// find valid index
						i++; 
				}
				
				if (binSz==0)
					*px = double(i); 
				else
				{
					// Get Max value in this voxel (0 25 75 ... )
					*px = double(i*binSz)-binSz/2;
					if ((*px)<0)
						(*px) = 0; 
				}
			}
		}
	}
	
	volumeReturn->Update();
	return volumeReturn; 
}

/* Computes Extend of given Bounds and VoxelSize (Number of Voxels in each Dimensions +/-) */
void helper::calculateExtent(double bounds[6], int extent[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;
	
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = (bounds[0] - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = (bounds[1] + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = (bounds[2] - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = (bounds[3] + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = (bounds[4] - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = (bounds[5] + spacing[2])/spacing[2]/* + 0.5*/;
}

/* Compute BoundingBox along all dimensions (calls getBoundingBox); 
 * min in [0], max in [1] 
 */
void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph)
{
	getBoundingBox(X_COORD,xBox,SpatialGraph);
	getBoundingBox(Y_COORD,yBox,SpatialGraph);
	getBoundingBox(Z_COORD,zBox,SpatialGraph);
}

void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	getBoundingBox(X_COORD,xBox,SpatialGraph,neuriteID);
	getBoundingBox(Y_COORD,yBox,SpatialGraph,neuriteID);
	getBoundingBox(Z_COORD,zBox,SpatialGraph,neuriteID);
}


/* Get BoundingBox along one dimension.int Coordinate between 0 and 2 (X_COORD, Y_COORD, Z_COORD), double Coord returns min and max point */
void helper::getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph)
{
	double maxC = -1E9;
	double minC = 1E9;
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label >= Neuron && (*edgeIt)->label <= Soma)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double tmpC = (*edgeListIt)[Coordinate];
				if(tmpC > maxC)
					maxC = tmpC;
				if(tmpC < minC)
					minC = tmpC;
			}
		}
	}
	
	Coord[0] = minC;
	Coord[1] = maxC;   
}

/* Get BoundingBox along one dimension.int Coordinate between 0 and 2 (X_COORD, Y_COORD, Z_COORD), double Coord returns min and max point */
void helper::getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	double maxC = -1E9;
	double minC = 1E9;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == neuriteID)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double tmpC = (*edgeListIt)[Coordinate];
				if(tmpC > maxC)
					maxC = tmpC;
				if(tmpC < minC)
					minC = tmpC;
			}
		}
	}

	Coord[0] = minC;
	Coord[1] = maxC;
}

/* Stores ImageDataVolume under given name */
void helper::storeImageVolume(const char * outputFilename, ImageDataPointerType volume)
{
	std::string fname(outputFilename);
	if (fname.compare(fname.size()-3,3,".am") == 0)
		fname = fname.substr(0,fname.size()-3);
	
	Reader * Writer = new Reader(fname.c_str(), fname.c_str());
	
	int numScalar = volume->GetNumberOfScalarComponents(); 
	
	if (numScalar==1)
		Writer->writeScalarField(volume);
	else
		Writer->writeNVectorField(volume);
	
	delete Writer; 
}

/* Checks whether given point coordinates x y z are within ImageDataVolume */
bool helper::isValidPt(int x, int y, int z, ImageDataPointerType volume)
{
	int extent[6];
	volume->GetExtent(extent);	
	return (extent[0]<=x && extent[1]>=x && extent[2]<=y && extent[3]>=y && extent[4]<=z && extent[5]>=z);
}

/* Gets ImageDataPointerType as Input and creates .txt containing the normalized profile along x and z axis */
void helper::getProfileAsTxt(ImageDataPointerType volume, const char * path, double VoxelSize)
{   
	int dim[3]; 
	volume->GetDimensions(dim);

	// Get 2D-Array and Initialize with Zeros
	int x_off = volume->GetExtent()[0]; 
	double x_profile[dim[0]][2];
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		x_profile[x-x_off][1] = 0.0; 
		x_profile[x-x_off][0] = x*VoxelSize; 
	}
	
	int y_off = volume->GetExtent()[2]; 
	double y_profile[dim[1]][2];
	for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
	{
		y_profile[y-y_off][1] = 0.0;
		y_profile[y-y_off][0] = y*VoxelSize;
	}
	
	int z_off = volume->GetExtent()[4]; 
	double z_profile[dim[2]][2];
	for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
	{
		z_profile[z-z_off][1] = 0.0;
		z_profile[z-z_off][0] = z*VoxelSize;
	}
	
	// Sum up values in volume
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{ 
				double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
				x_profile[x-x_off][1] += (*px); 
				z_profile[z-z_off][1] += (*px); 
				y_profile[y-y_off][1] += (*px); 
			}
		}
	}
	
	// Normalize 
	double x_norm =  1E9/(VoxelSize*VoxelSize*VoxelSize*dim[1]*dim[2]); 	// (50*50*50*y_dim*z_dim)
	double z_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[1]); 	// (50*50*50*x_dim*y_dim)
	double y_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[2]); 	// (50*50*50*x_dim*z_dim)
	
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		x_profile[x-x_off][1] *= x_norm; 
	}	
	for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
	{
		y_profile[y-y_off][1] *= y_norm;
	}
	for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
	{
		z_profile[z-z_off][1] *= z_norm;
	}
	
	std::ofstream outStream(path);
	
	if(!outStream.fail())
	{	  
		outStream << "x_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{
			outStream << x_profile[x-x_off][0] << " " << x_profile[x-x_off][1] << std::endl; 
		}
		
		outStream << "\n" << std::endl; 
		outStream << "y_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
		
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			outStream << y_profile[y-y_off][0] << " " << y_profile[y-y_off][1] << std::endl;
 		}

 		outStream << "\n" << std::endl; 
		outStream << "z_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
		
		for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
		{
			outStream << z_profile[z-z_off][0] << " " << z_profile[z-z_off][1] << std::endl;
 		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  path << std::endl; 
	}
}

/* Computes directory one level up the given directory or file */
std::string helper::getRootFromPath(const char* inputFilenameList)
{
	std::string ofName(inputFilenameList); 
	unsigned found = ofName.find_last_of("/"); 
	std::string rootDir = ofName.substr(0,found);
	rootDir += "/";
	return rootDir; 
}

std::string helper::getFilenameFromPath(const char* inputPathFilename)
{
	// Save truncated Spatial Graph
	std::string tmpStr(inputPathFilename);
	unsigned found = tmpStr.find_last_of("/");
	std::string Filename = tmpStr.substr(found+1);
	return Filename;
}

/* Displays array of length len in console / shell */
void helper::printArray(double array[], int len)
{
	for (int i = 0; i<len; i++)
		std::cout << array[i] << " "; 
	std::cout << "" << std::endl; 
}

/* Compute BoundingBox of ImageDataVolume/Density for given VoxelSize (50) */
void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2], ImageDataPointerType volume, double VoxelSize)
{
	double box[6];
	volume->GetBounds(box); 
	
	xBox[0] = box[0]-VoxelSize/2; 
	xBox[1] = box[1]+VoxelSize/2; 
	yBox[0] = box[2]-VoxelSize/2; 
	yBox[1] = box[3]+VoxelSize/2; 
	zBox[0] = box[4]-VoxelSize/2; 
	zBox[1] = box[5]+VoxelSize/2; 
}

/* Computes length of spatialGraph neuriteID */
double helper::lengthCell(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{	
	// length of edges within clippedBox
	double len = 0; 
// 	double p = 0; 
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		// Check only certain label
		if((*edgeIt)->label == neuriteID)
		{   
			std::list< double * >::iterator edgeListIt;			
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;
			++edgeListIt;

			// Compute Distance between first Vertex and first Point (they should be identical!)
			int tmpVertexID = (*edgeIt)->edgeConnectivity[0];
			Vertex * tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			len += sqrt(pow((tmpVertex->coordinates[X_COORD]-previousPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-previousPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-previousPt[Z_COORD]),2.0));
			
			// Go through all points of each Edge starting with the second edge
			for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{  
				double * currentPt = *edgeListIt;
				len += sqrt(pow((currentPt[X_COORD]-previousPt[X_COORD]),2.0)+pow((currentPt[Y_COORD]-previousPt[Y_COORD]),2.0)+pow((currentPt[Z_COORD]-previousPt[Z_COORD]),2.0));  
				previousPt = currentPt;
			}
// 			p = p + (*edgeIt)->numEdgePoints;
// 			std::cout << (*edgeIt)->numEdgePoints << std::endl; 
		}
	}
	
// 	std::cout << " >> Edges " << e << " EdgePts " << p << std::endl; 
	return len;  
}

/* Function checks whether Vertex and Edge Points match! (Essential for calculation of Length) */
void helper::checkCell(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	double errFront = 0;
	double errBack = 0;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		// Check only certain label
		if((*edgeIt)->label == neuriteID)
		{
			std::list< double * >::iterator edgeListIt;

			// Compute Distance between first Vertex and first Point (they should be identical!)
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * firstPt = *edgeListIt;
			int tmpVertexID = (*edgeIt)->edgeConnectivity[0];
			Vertex * tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			errFront += sqrt(pow((tmpVertex->coordinates[X_COORD]-firstPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-firstPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-firstPt[Z_COORD]),2.0));

			// Compute Distance between first Vertex and last Point (they should be identical!)
			edgeListIt = (*edgeIt)->edgePointCoordinates.end();
			edgeListIt--;
			double * lastPt = *edgeListIt;
			tmpVertexID = (*edgeIt)->edgeConnectivity[1];
			tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			errBack += sqrt(pow((tmpVertex->coordinates[X_COORD]-lastPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-lastPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-lastPt[Z_COORD]),2.0));
		}
	}

	if (errFront > 1)
	{
		std::cout << "-- ISSUE FOUND -- " << std::endl;
		std::cout << "WARNING! The first Edge Point and the first Vertex Point do not match! Introduced length error is " << errFront << "!" << std::endl;
		std::cout << " helper::lengthCell _does_ takes care of that issue when calculating the length!" << std::endl;
	}

	if (errBack > 1)
	{
		std::cout << "-- ISSUE FOUND -- " << std::endl;
		std::cout << "WARNING! The last Edge Point and the last Vertex Point do not match! Introduced length error is " << errBack << "!" << std::endl;
		std::cout << " WARNING! helper::lengthCell _DOES_NOT_ take care of that issue when calculating the length!" << std::endl;
	}

	//std::cout << "Error (front) : " << errFront << " Error (back) : " << errBack << std::endl;
}


/* Returns spatialGraph from given directory */
AmiraSpatialGraph * helper::getSpatialGraph(const char * filename)
{
  	Reader * fileReader = new Reader(filename, filename);
	
	std::string fname(filename); 
	if((access(fname.c_str(), F_OK) != -1))
	{
		if (fname.compare(fname.size()-3,3,".am") == 0)
		{
			fileReader->readSpatialGraphFile(0);
			return fileReader->getSpatialGraph(); 
		}
		else if (fname.compare(fname.size()-4,4,".hoc") == 0)
		{
			fileReader->readHocFile();
			return fileReader->getSpatialGraph(); 
		}
		else
		{
			std::cout << "Inputfile is neither .hoc nor .am! Cannot read file" << std::endl; 
			return 0; 
		}
	}
	else
	{
		std::cout << "Inputfile not found: " << fname << std::endl; 
		return 0; 
	}
}

/* Gets two ImageDataPointerType as Input and creates .txt containing the normalized profile along x and z axis */
void helper::getProfileAsTxtImprov(ImageDataPointerType volume, ImageDataPointerType volumeOld, const char * path, double VoxelSize)
{   
	int dim[3]; 
	volume->GetDimensions(dim);
	
	int x_off = volume->GetExtent()[0]; 
	double profile[dim[0]][5];
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		profile[x-x_off][1] = 0.0; 
		profile[x-x_off][2] = 0.0; 
		profile[x-x_off][3] = 0.0; 
		profile[x-x_off][4] = 0.0; 
		profile[x-x_off][0] = x*VoxelSize; 
	}

	// Sum up values in volume
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{ 
				double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
				profile[x-x_off][1] += (*px); 
				profile[z-x_off][2] += (*px); 
				
				if (helper::isValidPt(x,y,z,volumeOld))
				{
					double *pxOld = static_cast< double * >(volumeOld->GetScalarPointer(x,y,z));
					profile[x-x_off][3] += (*pxOld); 
					profile[z-x_off+8][4] += (*pxOld);
				}
			}
		}
	}
	
	// Normalize 
	double x_norm =  1E9/(VoxelSize*VoxelSize*VoxelSize*dim[1]*dim[2]); // (50*50*50*y_dim*z_dim)	
	double z_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[1]); // (50*50*50*y_dim*x_dim)
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		profile[x-x_off][1] *= x_norm; 
		profile[x-x_off][2] *= z_norm;
		profile[x-x_off][3] *= x_norm; 
		profile[x-x_off][4] *= z_norm;
	}
	
	//std::ofstream outStream(name.c_str());
	std::ofstream outStream(path);
	
	if(!outStream.fail())
	{	  
		outStream << "range x z x_old z_old  Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{
			outStream << profile[x-x_off][0] << " " << profile[x-x_off][1] << " " << profile[x-x_off][2] << " " << profile[x-x_off][3] << " " << profile[x-x_off][4] << std::endl; 
		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  path << std::endl; 
	}
}

/* Scales spatialGraph along cut/shrunk z-dimension. Scales it back to 300um by scaling the BB of the spatialGraph */
void helper::zScaling(AmiraSpatialGraph * spatialGraph)
{
	std::cout << "WARNING! Using default scaling of 300um! L2/3 INs have 350um!" << std::endl;
	zScaling(spatialGraph,300);
};

/* Scales spatialGraph along cut/shrunk z-dimension. Scales it back to original slicing thickness by scaling the BB of the spatialGraph */
void helper::zScaling(AmiraSpatialGraph* spatialGraph, double thickness)
{
	// Get Z coordinates for scaling	
	double r[2]; 
	getBoundingBox(Z_COORD,r,spatialGraph);
	double minZ = r[0];
	double maxZ = r[1];
	
	TransformPointerType scal = TransformPointerType::New();
	
	// Scale z-axis coordinates to Slicing thickness in micrometer
	double zfactor = thickness/(maxZ-minZ);
	scal->Scale(1,1,zfactor);
	
	std::cout << "    z-scaling factor = " << zfactor << " (Slicing thickness = " << thickness << " um)" << std::endl;
	
	spatialGraph->setTransformation(scal);
	spatialGraph->applyTransformation();
}

/* Compute the center of a label in the spatialGraph (most likely center of soma) 
 * can then be aligned with helper::align
 */
void helper::centerOfSpatialGraph(int neuriteID, double centerPt[3], AmiraSpatialGraph * spatialGraph)
{
	PolyDataPointerType structure = PolyDataPointerType::New();
	if(!spatialGraph->extractLandmark(neuriteID, structure))
	{
		std::cout << "Error! Could not find structure with ID " << neuriteID <<" in SpatialGraph!" << std::endl;
		return;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
	structure->GetCell(0)->GetParametricCenter(pCoords);
	structure->GetCell(0)->EvaluateLocation(subID, pCoords, centerPt, weights);
	delete [] weights;     
}

/* Aligns spatialGraph to neuriteID (Soma) */
void helper::align(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	double centerPt[3]; 
	centerOfSpatialGraph(neuriteID, centerPt, spatialGraph); 
		
	// Subtract centerPt from spatialGraph
	for(int i=0; i<3; i++)
	{
	    centerPt[i] *= -1;
	}
  
	TransformPointerType sub = TransformPointerType::New();
	sub->Translate(centerPt);
	spatialGraph->setTransformation(sub);
	spatialGraph->applyTransformation();
}

/* Computes L2 Euclidean Norm */
double helper::norm(double vector[3])
{
	return (sqrt(pow(vector[0],2.0) + pow(vector[1],2.0) + pow(vector[2],2.0))); 
}

/* Draws #numSamples Gaussian Random Numbers given Mean and SD,
 * WARNING! Called in a loop will produce same output because random number generator is set within function
 */
double * helper::randn(double mean, double SD, int numSamples)
{
	gsl_rng_env_setup();	
	//const gsl_rng_type * T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set (r,time(NULL));
	
	double * randnSample = new double[numSamples];
	
	for (int i = 0; i < numSamples; i++) 
	{
		randnSample[i] = mean + gsl_ran_gaussian_ratio_method(r, SD);
	}
	
	gsl_rng_free(r);
	return randnSample; 
}

/* Draws #numSamples Gaussian Random Numbers given Mean and SD, and random number generator  */
double * helper::randn(double mean, double SD, int numSamples, gsl_rng * r)
{
	double * randnSample = new double[numSamples];
	
	for (int i = 0; i < numSamples; i++) 
	{
		randnSample[i] = mean + gsl_ran_gaussian_ratio_method(r, SD);
	}
	return randnSample; 
}

/* Computes Cross Product between two 3D vectors */
double * helper::cross(double u[3], double v[3])
{
	double * crossProduct = new double[3];
	
	crossProduct[0] = u[1]*v[2] - u[2]*v[1]; 
	crossProduct[1] = u[2]*v[0] - u[0]*v[2]; 
	crossProduct[2] = u[0]*v[1] - u[1]*v[0]; 
	
	return crossProduct; 
}

/* Returns list of IDs to Branch Points in given spatialGraph for given neurite (axon, dendrite, soma)
 * WARNING! Only works for .hoc files. Requires fatherID!
 */
 std::list< int > helper::getBranchPointIDs(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	std::map< int, int > nrOfChildBranches;
	std::list< int > branchPtIDs;
	std::vector< Edge * >::const_iterator edgeIt;
		
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == neuriteID)
		{  
			int fatherID = (*edgeIt)->fatherID;
			if(nrOfChildBranches.find(fatherID) != nrOfChildBranches.end()) // Check whether fatherID exists, if so, increase count by one
			{
				nrOfChildBranches[fatherID]++;
			}
			else // if not, insert fatherID with corresponding count of 1
			{
				nrOfChildBranches.insert(std::pair< int, int >(fatherID, 1));
			}
		}
	}
	
	std::map< int, int >::const_iterator nrBranchIt;
	for(nrBranchIt = nrOfChildBranches.begin(); nrBranchIt != nrOfChildBranches.end(); ++nrBranchIt)
	{
		if(nrBranchIt->second > 1) // If count of fatherID is larger than 1, then this is a branching point, put fatherID in list
		{
			branchPtIDs.push_back(nrBranchIt->first);
		}
	}
	
	return branchPtIDs;
}

bool helper::isPtOnEdge(double pt[3], double bounds[6])
{
	bool onEdge = (pt[X_COORD]==bounds[0] || pt[X_COORD]==bounds[1]
					   || pt[Y_COORD]==bounds[2] || pt[Y_COORD]==bounds[3]
					   || pt[Z_COORD]==bounds[4] || pt[Z_COORD]==bounds[5]);
    return onEdge;
}

/* Return whether 3D point is within 3D bounding box
 * [xmin xmax ymin ymax zmin zmax]
 */
bool helper::isPtWithinBounds(double pt[3], double bounds[6])
{
	bool outsideOfBounds = (pt[X_COORD]<bounds[0] || pt[X_COORD]>bounds[1]
					   || pt[Y_COORD]<bounds[2] || pt[Y_COORD]>bounds[3]
					   || pt[Z_COORD]<bounds[4] || pt[Z_COORD]>bounds[5]);

    return !outsideOfBounds;
}

// Computes t values between two points and a box (standing_vector + t*direction_vector)
std::list<double> helper::getIntersectionValue(double bounds[6],double pt1[3], double pt2[3])
{
	std::list < double> t_values;
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<3; j++)
		{
			double tmp = computeIntersectionValue(pt1[j]-bounds[j*2+i],pt2[j]-bounds[j*2+i]);
			if (tmp != 0.0)
				t_values.push_back(tmp);
		}
	}
	return t_values;
}

// Returns t-value (pt1 + t * (pt2 - pt1))
double helper::computeIntersectionValue(double Dst1, double Dst2)
{
	// Check whether Line is intersecting with given plane
	if ((Dst1*Dst2) >= 0.0)
		return 0.0;
	if (Dst1 ==  Dst2)
		return 0.0;

	// Compute intersecting Point via standing vetor plus weighted direction vector
	double tmp = (-Dst1/(Dst2-Dst1 /*+ 1E-8 */));

	if (tmp>0.0 && tmp<1.0)
		return tmp;
	else
		return 0.0;
}

// Compute intersection point between two points pt1 and pt2 given the relative distance measurement t
void helper::getIntersectionPoint(double t, double pt1[3], double pt2[3], double intersectPt[3])
{
	for (int i = 0; i<3; i++)
	{
		intersectPt[i] = pt1[i] + (pt2[i]-pt1[i]) * t;
	}
}

/* Rotate SpatialGraph so that Y becomes Z (column axis) 
 * Input: 
 * 	-spatialGraph which should be rotated
 * 	- int d: which direction: 
 * 		 1 : y -> z (FINAL, z is column axis)
 * 		-1 : z -> y (ORIGINAL, y is column axis)
 */
void helper::switchYZ(AmiraSpatialGraph * spatialGraph, int d)
{
	if (!(spatialGraph->isLabelInSpatialGraph(Soma)))
	{
		std::cout << "ERROR! Soma was not found in given SpatialGraph! Cannot compute Soma Location for Rotation!" << std::endl;
		return; 
	}
	
	if (d != 1 && d != -1)
	{
		std::cout << "ERROR! d has to be either 1 (make z column axis) or -1 (make y column axis), not " << d << std::endl; 
		return; 
	}
  
	double somaPt[3]; 
	double shiftPt[3]; 
	centerOfSpatialGraph(Soma, somaPt, spatialGraph);
  
	// Set Soma Position to 0 0 0
	TransformPointerType sub1 = TransformPointerType::New();
	shiftPt[0] = -somaPt[0]; 
	shiftPt[1] = -somaPt[1]; 
	shiftPt[2] = -somaPt[2]; 
	sub1->Translate(shiftPt);
	spatialGraph->setTransformation(sub1);
	spatialGraph->applyTransformation();
	
	// Rotate in x-dimension 90°
	TransformPointerType rot = TransformPointerType::New();
	double rotVec[3] = {d,0,0}; 
	rot->RotateWXYZ(90, rotVec);
	spatialGraph->setTransformation(rot);
	spatialGraph->applyTransformation();
	
	// Set back to original Soma Position
	TransformPointerType sub2 = TransformPointerType::New();
	shiftPt[0] = somaPt[0]; 
	shiftPt[1] = somaPt[2]; 
	shiftPt[2] = somaPt[1]; 
	sub2->Translate(shiftPt);
	spatialGraph->setTransformation(sub2);
	spatialGraph->applyTransformation();
}

/* Creates Directory / Folder in case it does not exist already */
void helper::mkDir(const char * path) 
{
	std::string fname(path); 
	std::string mkdir = "mkdir " + fname; 
	
	if (access(fname.c_str(), F_OK) == -1)
		system(mkdir.c_str());
}

/* executes command line
 * Example:
 * 	char * cmd = ""/bin/bash /home/dudvary/genAxons.sh yeyey""; 
 * 	std::string out = helper:: exec(cmd); 
 * 	std::cout << out << std::endl;
*/
std::string helper::exec(char* cmd)
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
std::vector<std::string> helper::getReturnOfLsCmd(std::string path)
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

/* Computes BranchOrder for each EdgeID
 * WARNING! Only works for .hoc input! Checks fatherID!
 */
std::map< int, int > helper::getBranchOrder(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	// Branch Orderfor each Edge/Segment
	std::map< int, int > edgeBO; 
	std::list< int > BranchPointIDs = getBranchPointIDs(spatialGraph, neuriteID); 
	
	int edgeID = 0; 
	for(std::vector< Edge * >::iterator edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt, ++edgeID)
	{ 
		if((*edgeIt)->label == neuriteID)
		{
			// Get Branch Order
			int branchOrder = 0; 
			int fatherID = (*edgeIt)->fatherID;
			while (fatherID != -1)
			{
				branchOrder++; 
				Edge * newEdge = (*(spatialGraph->edgesPointer()))[fatherID];
				fatherID = newEdge->fatherID; 
			}
			edgeBO[edgeID] = branchOrder; 
		}
	}
	return edgeBO; 
}

void helper::flipYSpatialGraph(AmiraSpatialGraph * spatialGraph)
{
  	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label >= Neuron && (*edgeIt)->label <= Soma)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				(*edgeListIt)[Y_COORD] = -(*edgeListIt)[Y_COORD];
			}
		}
	}
}

std::vector<unsigned int> helper::returnCellIDs(const char * inputFilenameList)
{
  	std::vector<unsigned int> CellIDs;
	std::string currentLine;
	std::ifstream inputStream(inputFilenameList);
	if(!inputStream.fail())
	{
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine);
			if (currentLine.empty())
				break;
			unsigned int CellID;
			std::istringstream(currentLine) >> CellID;
			CellIDs.push_back(CellID);
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Could not find " << inputFilenameList << std::endl;
	}

	return CellIDs;
}

std::vector<std::string> helper::returnNames(const char* inputFilenameList, const char* root)
{
  	std::vector<std::string> fnameList;
	std::string currentLine; 
	std::ifstream inputStream(inputFilenameList); 
	if(!inputStream.fail())
	{
		std::string rootStr(root); 
	  
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine); 
			if (currentLine.empty())
				break; 
			currentLine = rootStr + currentLine; 
			fnameList.push_back(currentLine);
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Could not find " << inputFilenameList << std::endl;
	}
	
	return fnameList; 
}

double helper::computeEuclideanDistance(double pt1[3], double pt2[3])
{
	double dist = sqrt( pow( ( pt1[X_COORD] - pt2[X_COORD] ), 2.0 ) +
						pow( ( pt1[Y_COORD] - pt2[Y_COORD] ), 2.0 ) +
						pow( ( pt1[Z_COORD] - pt2[Z_COORD] ), 2.0 ));
	return dist;
}

// Computes the Surface Area, (specifically for Dendrites)
// Connects primary dendrites with Soma.
// Returns Lateral Surface Area using A = PI * (r1 + r2) * sqrt(height^2 + (r1 - r2)^2)
// 	https://www.quora.com/What-is-the-area-of-a-cylinder-of-different-radius-at-top-and-bottom-respectively
// 			misses Top/Bottom Surface area at sides PI * (r1^2 + r2^2)
//			Total Surface Area = LateralSurfaceArea + Top/BottomSurfaceArea
// with height being the length between two points and r1 and r2 being the radius of the two points
// NOTE: In .hoc files the radius denotes the diameter, thus in amira files, the radius also denotes the diameter!
// If computing the "true" surface area the radius needs to be divided by two.
// If computing the surface area like the NetworkAnalyzer in Amira, the diameter is set to be the radius.
double helper::computeSurfaceArea(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	return helper::computeSurfaceArea(SpatialGraph, neuriteID, false);
}
// If likeAmira true, computes "surface area" as NetworkAnalyzer in Amira does it, by setting the radius = diameter
double helper::computeSurfaceArea(AmiraSpatialGraph * SpatialGraph, int neuriteID, bool likeAmira)
{
	double surfaceArea = 0.0;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{

		if((*edgeIt)->label == neuriteID)
		{
			// Point
			std::list< double * >::iterator edgeListIt;
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;
			++edgeListIt;

			// Radius
			std::list< double >::iterator edgeRadiusListIt;
			edgeRadiusListIt = (*edgeIt)->radiusList.begin();
			double previousRad = *edgeRadiusListIt;
			++edgeRadiusListIt;

			if (likeAmira)
				previousRad = previousRad/2.0;

			// Go through all points of each Edge starting with the second edge
			for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt, ++edgeRadiusListIt)
			{
				double * currentPt = *edgeListIt;
				double currentRad = *edgeRadiusListIt;

				if (likeAmira)
					currentRad = currentRad/2.0;

				double heigth = sqrt(pow((currentPt[X_COORD]-previousPt[X_COORD]),2.0) +
									pow((currentPt[Y_COORD]-previousPt[Y_COORD]),2.0) +
									pow((currentPt[Z_COORD]-previousPt[Z_COORD]),2.0));
				double lateralArea = PI * (currentRad+previousRad) * sqrt(pow(heigth,2.0) +
									pow(currentRad-previousRad,2.0));
				surfaceArea += lateralArea;

				previousPt = currentPt;
				previousRad = currentRad;
			}
		}
	}

	return surfaceArea;
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellType: CellType ID that should be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * stores morphologies in outputpath
 * stores soma.csv in outputpath (contains soma position of each morphology)
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath)
{
	double bounds[6];

	bounds[0] = -(std::numeric_limits<double>::infinity());
	bounds[1] = std::numeric_limits<double>::infinity();
	bounds[2] = -(std::numeric_limits<double>::infinity());
	bounds[3] = std::numeric_limits<double>::infinity();
	bounds[4] = -(std::numeric_limits<double>::infinity());
	bounds[5] = std::numeric_limits<double>::infinity();

	extractMorphologies(spatialGraphSetFilename, CellType, outputpath, bounds);
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellIDsFilename: path/to/CellIDs.txt file containing all CellIDs to be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * stores morphologies in outputpath
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, const char * CellIDsFilename, const char * outputpath)
{
	mkDir(outputpath);
	std::string outputpathstr(outputpath);
	std::ifstream inputStream(CellIDsFilename);

	if(!inputStream.fail())
	{
		// Read in SpatialGraphSet
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		Reader::readSpatialGraphSetFile(spatialGraphSetFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);

		// Read in CellIDs
		std::string currentLine;

		while(!inputStream.eof()) // go through all CELLIDs
		{
			getline(inputStream,currentLine);
			if (currentLine.empty())
				break;

			int CellID = atoi(currentLine.c_str());

			// Get Filename
			std::string originalName = originalGraphFiles[originalGraphIndices[CellID]];
			originalName = originalName.substr(1, originalName.size()-2);
			std::string loadName = getRootFromPath(spatialGraphSetFilename) + originalName;

			// Outputfilename
			std::ostringstream CellIDStr;   // stream used for the conversion
			CellIDStr << CellID;      // insert the textual representation of 'Number' in the characters in the stream
			std::string saveName = outputpathstr + CellIDStr.str();

			// Load Spatial Graph
			Reader * spatialGraphReader = new Reader(loadName.c_str(), saveName.c_str());
			spatialGraphReader->readSpatialGraphFile(0);
			AmiraSpatialGraph * sg = spatialGraphReader->getSpatialGraph();

			// Apply Transformation
			sg->setTransformation(amiraToVTKTransform(spatialGraphTransforms[CellID]));
			sg->applyTransformation();

			// Store SpatialGraph
			spatialGraphReader->writeSpatialGraphFile();

			delete sg;
			delete spatialGraphReader;
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Error! Reading file with CellIDs failed! Path: " <<  CellIDsFilename << std::endl;
	}
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellType: CellType ID that should be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * - bounds[6]: [xmin xmax ymin ymax zmin zmax]: bounding box of to be extracted morphologies.
 * 	 morphologies outside this box will not be extracted/stored
 * stores morphologies in outputpath
 * stores soma.csv in outputpath (contains soma position of each morphology)
 * stores cells.txt in outputpath (contains list of all morphology names)
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath, double bounds[6])
{
	mkDir(outputpath);
	std::string outputpathstr(outputpath);
	std::string outputFilenameSoma = outputpathstr + "soma.csv";
	std::ofstream outStreamSoma(outputFilenameSoma.c_str());

	std::string outputFilenameCells = outputpathstr + "cells.txt";
	std::ofstream outStreamCells(outputFilenameCells.c_str());

	if(!outStreamSoma.fail() && !outStreamCells.fail())
	{
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		Reader::readSpatialGraphSetFile(spatialGraphSetFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);
		int dendCounter = 0;

		outStreamSoma << "OriginalName,StoredName,SomaX,SomaY,SomaZ,CellID," << std::endl;

		for(int i = 0; i < cellTypeIDs.size(); ++i)
		{
			if (cellTypeIDs[i]==CellType)
			{
				// Get Filename
				std::string originalName = originalGraphFiles[originalGraphIndices[i]];
				originalName = originalName.substr(1, originalName.size()-2);

				std::ostringstream dendCounterStr;   // stream used for the conversion
				dendCounterStr << dendCounter;      // insert the textual representation of 'Number' in the characters in the stream

				std::string saveName = outputpathstr + "Dend_" + dendCounterStr.str();
				std::string loadName = getRootFromPath(spatialGraphSetFilename) + originalName;

				// Load Spatial Graph
				Reader * spatialGraphReader = new Reader(loadName.c_str(), saveName.c_str());
				spatialGraphReader->readSpatialGraphFile(0);
				AmiraSpatialGraph * sg = spatialGraphReader->getSpatialGraph();

				// Apply Transformation
				sg->setTransformation(amiraToVTKTransform(spatialGraphTransforms[i]));
				sg->applyTransformation();

				// Get Soma Position
				double centerPt[3];
				centerOfSpatialGraph(Soma, centerPt, sg);

				// Check whether Soma is within bounding box
				if (bounds[0]<=centerPt[0] && bounds[1]>=centerPt[0] && bounds[2]<=centerPt[1] && bounds[3]>=centerPt[1] && bounds[4]<=centerPt[2] && bounds[5]>=centerPt[2])
				{
					// Write out name and soma position
					outStreamSoma << originalName << "," << saveName << "," << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << "," << i << "," << std::endl;
					outStreamCells << "Dend_" << dendCounterStr.str() << ".am" << std::endl;

					// Store SpatialGraph
					spatialGraphReader->writeSpatialGraphFile();
					dendCounter++;
				}

				delete sg;
				delete spatialGraphReader;
			}
		}
		outStreamSoma.close();
		outStreamCells.close();
	}
	else
	{
		std::cout << "Error! Writing soma.csv failed! Path: " <<  outputFilenameSoma << std::endl;
		std::cout << "Error! Writing cells.txt failed! Path: " <<  outputFilenameCells << std::endl;
	}
}

// Change Label in SpatialGraph (Edge and Vertex)
// for example change BasalDendrite to Dendrite
void helper::changeSpatialGraphLabel(AmiraSpatialGraph * spatialGraph,int oldLabel, int newLabel)
{
	// Change EdgeLabels
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == oldLabel)
		{
			(*edgeIt)->label = newLabel;
		}
	}

	// Change VertexLabels
	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = spatialGraph->verticesBegin(); vertexIt != spatialGraph->verticesEnd(); ++vertexIt)
	{
		if((*vertexIt)->label == oldLabel)
		{
			(*vertexIt)->label = newLabel;
		}
	}
}

TransformPointerType helper::amiraToVTKTransform(double* amiraTransform)
{
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
		{
			mat->SetElement(j, i, amiraTransform[i*4+j]);
		}

	TransformPointerType vtkTransform = TransformPointerType::New();
	vtkTransform->SetMatrix(mat);
	return vtkTransform;
}

void helper::writeListToFile(std::list<double> valList,const char * filename)
{
	std::ofstream TXTWriter;
	TXTWriter.open(filename);
	if(!TXTWriter.fail())
	{
		for (std::list<double>::iterator it=valList.begin(); it != valList.end(); ++it)
			TXTWriter << *it << std::endl;
	}
	else
	{
		std::cout << "ERROR! Writing " << filename << " failed! " << std::endl;
	}
	TXTWriter.close();
}