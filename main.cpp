#include "mainCLP.h"

#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkImageIOBase.h>
#include <itkVariableLengthVector.h>
#include "itkArray2D.h"
#include "itkImageDuplicator.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include "itkTimeProbe.h"
#include "itkVariableSizeMatrix.h"
#include "itkAffineTransform.h"
#include <itksys/SystemTools.hxx>
#include <itkImage.h>
#include <itkImageIOBase.h> 
#include <itkVectorImage.h> 
#include <itkVariableLengthVector.h>
#include <itkNrrdImageIO.h>
#include <itkMetaDataObject.h> 
#include <itkVector.h>
#include <itkMatrix.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkMath.h>
#include <vtkTensorGlyph.h>
#include <vtkProperty.h>
#include <vtkDataSet.h>
#include <vtkPolynomialSolversUnivariate.h>



#include <sstream>
#include <stdio.h> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES
#include "math.h"
#include <vcl_cmath.h>
#include <vxl_config.h>
#include <vnl/vnl_config.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

#define VOLUME_DIMENSION 3
#define ImageDimension 3

using namespace std; 

struct FOMeta{
    double origin[3];
    double spacing[3];
    unsigned int size[3];
};


int readFOMeta(const char* foImgFileName, FOMeta* m);

template <class ImagePointer>
		int getFOImageHandler(ImagePointer &foImage, const char* foImgFileName);

double estimateHinderedDiffusion(std::vector< double > EigenValue, itk::Vector<double, 3> bkdirection,double WidthPulseGradient,double DiffTime);

int generationdwi(string dwiImgFilename,string T2ImgFilename,string OutFilename,string EigenFilename,string foImgFilename,string EigenHinderedImgFilename,double Timetoecho,double DiffTime,double WidthPulseGradient,double MagnitudeG,float fH,float fR,double noiseSigma);


int main(int argc, char *argv[])
{
	PARSE_ARGS;
	
	
	generationdwi(dwiImgFilename,T2ImgFilename,OutFilename,EigenFilename,foImgFilename,EigenHinderedImgFilename,Timetoecho,DiffTime, WidthPulseGradient,MagnitudeG,fH, fR,noiseSigma);
	
	return 0;
}


int generationdwi(string dwiImgFilename,string T2ImgFilename,string OutFilename,string EigenFilename,string foImgFilename,string EigenHinderedImgFilename,double Timetoecho,double DiffTime,double WidthPulseGradient,double MagnitudeG,float fH,float fR,double noiseSigma)
{
	//We define dwi parameters
	typedef double RealType;
	float Radius=0.008;
	double gyroRad=267.5;
	double gyro=42.576;
	double pi=M_PI;   
	double t=Timetoecho/2;
	double DPa;
	double DPe;
	std::string Inputname = dwiImgFilename;
	std::string T2Name = T2ImgFilename;
	const char * EigenFileName=EigenFilename.c_str();
	const char *EigenHinderedimgFileName=EigenHinderedImgFilename.c_str();
	const char *foimgFileName=foImgFilename.c_str();	
	
	typedef double PixelType;
	
	
	//Read dwi original and store gradients:
		
	typedef itk::Image< PixelType , 3 > ImageType ;
	typedef itk::ImageFileReader< ImageType > FileReaderType ; 
	typedef itk::VectorImage< PixelType , 3 > VectorImageType ; 
	ImageType::Pointer image ;
	std::vector< ImageType::Pointer > vectorOfImage ;
	itk::MetaDataDictionary dico ;
	itk::VectorImage< PixelType, 3 >::Pointer dwi_template ;
	dwi_template = itk::VectorImage< PixelType , 3 >::New() ; 
	
	itk::ImageFileReader< VectorImageType >::Pointer reader ;
	reader = itk::ImageFileReader< VectorImageType >::New() ;
	reader->SetFileName( Inputname) ;
	reader->Update() ;
	dwi_template = reader->GetOutput();
	//Save metadata dictionary
	dico = reader->GetOutput()->GetMetaDataDictionary() ;
	//on recupere le vecteur directions de :TransformGradients(dico)

	//define variables to add Rician Noise
	//RealType noiseSigma = argv[12];
	/*typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
	typename RandomizerType::Pointer randomizer = RandomizerType::New();
	randomizer->Initialize();*/
	
	//we take the gradient directions from the original dwi thanks to its metadatadictionnary then we store it in a vector of vector
	std::vector<itk::Vector<double, 3> > directions;
	itk::Vector<double, 3> direction;
	
	
	double b_value;
	itk::Vector<double> b_values;
	int i=0;
	typedef itk::MetaDataObject< std::string > MetaDataStringType ;
	itk::MetaDataDictionary::ConstIterator itr = dico.Begin() ;
	itk::MetaDataDictionary::ConstIterator end = dico.End() ; 
	while( itr != end )
	{
		itk::MetaDataObjectBase::Pointer entry = itr->second ;
		MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType* >( entry.GetPointer() ) ;
		if( entryvalue )
		{
			//get the gradient directions
			int pos = itr->first.find( "DWMRI_gradient" ) ;
			int pos2 = itr->first.find( "DWMRI_b-value" ) ;
			if( pos2 != -1 )//we find the b-value from original dwi metadictionnary
			{
				std::string tagvalue = entryvalue->GetMetaDataObjectValue() ;
				std::istringstream iss( tagvalue ) ;
				iss >> b_value;
				b_values[i]=b_value;
				++i;
				
			}  
			else if( pos != -1 )//we find the gradient directions from original dwi metadictionnary
			{
				std::string tagvalue = entryvalue->GetMetaDataObjectValue() ;
				itk::Vector< double , 3 > vec ;
				std::istringstream iss( tagvalue ) ;
				iss >> vec[ 0 ] >> vec[ 1 ] >> vec[ 2 ] ;//we copy the metavalue in an itk::vector
				direction[0]=vec[ 0 ];direction[1]=vec[ 1 ];direction[2]=vec[ 2 ];
				if( iss.fail() )
				{
					iss.str( tagvalue ) ;
					iss.clear() ;
					std::string trash ;
					iss >> vec[ 0 ] >> trash >> vec[ 1 ] >> trash >> vec[ 2 ] ;//in case the separator between the values is something else than spaces
					direction[0]=vec[ 0 ];direction[1]=vec[ 1 ];direction[2]=vec[ 2 ];
					if( iss.fail() )//problem reading the gradient values
					{
						std::cerr << "Error reading a DWMRI gradient value" << std::endl ;
					}
				} 
				directions.push_back(direction);
			}
		} 
		++itr ;
	} 
   
	/*estimate max b_value*/
	double maxBValue=0;
	int numBValue=b_values.GetNumberOfComponents();
	for (int i=0; i < numBValue; ++i)
	{
		if (b_values[i] >= maxBValue)
		{
			maxBValue = b_values[i];
		}
	} 
	
	/*estimate max gradient norm */
	int numGradients=directions.size();
	double normGradMax=0;
	for (int i=0; i <numGradients;i++)
	{
		itk::Vector<double, 3> bi = directions[i];
		if (bi.GetNorm()>normGradMax)
		{
			normGradMax=bi.GetNorm();
		}
	}
	
	/*estimate b of each gradient and then estimate q */
	std::vector<itk::Vector<double, 4> > directions2;
	
	std::cout<<"MagnitudeG is "<<MagnitudeG<<std::endl;
	for (int i=0; i <numGradients;i++)
	{
		itk::Vector<double, 3> bi = directions[i];
		double b=(pow((bi.GetNorm()),2)/pow(normGradMax,2))*maxBValue;
		float DiffTime=0.032020;
		
		/*estimation of the gradient normed with its magnitude*/
		itk::Vector<double, 4> direction2;
		for (int j=0;j<3;j++)
		{
			if(bi.GetNorm()!=0)
			{ 
			direction2[j]=(bi[j]/(bi.GetNorm()))*MagnitudeG*WidthPulseGradient*gyro;//estimate q
			}
			else
			{
			direction2[j]=0;
			}
		}
		direction2[3]=DiffTime;
		directions2.push_back(direction2);/*we store in directions2 vectors containing the three components of the normed gradient and at the end the DiffTime*/
	} 
  
	//Read the baseline
	typedef itk::ImageFileReader< ImageType  > ImageReaderType;
	ImageReaderType::Pointer  imageReader = ImageReaderType::New();
	ImageType::Pointer b0_img = ImageType::New();
	imageReader->SetFileName(T2Name);
	try{
		imageReader->Update();
		b0_img = imageReader->GetOutput();
	}
	catch (itk::ExceptionObject &ex){
	std::cout << ex << std::endl;
	return EXIT_FAILURE;
	}
	itk::ImageRegionIterator<ImageType> b0_img_it (b0_img, b0_img->GetLargestPossibleRegion());
	typedef itk::ImageRegionIterator< VectorImageType > IteratorType;
	IteratorType dwi_template_it( dwi_template, dwi_template->GetLargestPossibleRegion().GetSize() );
	
	/*Read FOImage, file containing fiber orientations for each voxel*/
	const unsigned int FODimension = 3;
	typedef std::vector< double > FOPixelType;
	typedef itk::Image< FOPixelType, FODimension > FOImageType;
	typedef itk::ImageRegionConstIterator< FOImageType > FOConstIteratorType;
	typedef itk::Vector<double, 3> VectorType;
	FOImageType::Pointer foImage = FOImageType::New();
	getFOImageHandler(foImage, foimgFileName); /*function to read txt file containing fiber orientation vectors*/
	FOConstIteratorType fo_it( foImage, foImage->GetLargestPossibleRegion());
	FOImageType::PixelType foValue;
	
	//signal estimation for each gradient direction and walking through each voxel of the new dwi
	
	itk::VariableLengthVector<PixelType> dwi_val = dwi_template_it.Get();
	//for each voxel
	//walking through every location on template_dwi and b0 image
	/*file in which we will store hindered parameters*/
	const unsigned int EigenHinderedDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenHinderedDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	EigenImageType::Pointer EigenHinderedImage = EigenImageType::New();
	getFOImageHandler(EigenHinderedImage, EigenHinderedimgFileName); /*function to read txt file containing hindered diffusion coefficients*/	
	EigenConstIteratorType eigen_hindered_it(EigenHinderedImage, EigenHinderedImage->GetLargestPossibleRegion());
	EigenImageType::PixelType EigenValue;
	
	eigen_hindered_it.GoToBegin();
	b0_img_it.GoToBegin();
	dwi_template_it.GoToBegin();
	fo_it.GoToBegin();
	while(!dwi_template_it.IsAtEnd() && !b0_img_it.IsAtEnd() && !fo_it.IsAtEnd() && !eigen_hindered_it.IsAtEnd()){
		std::cout<<"reading hindered diffusion of this voxel"<<std::endl;
		EigenValue=eigen_hindered_it.Get(); // a list that has 12 components of three eigenvalues and three eigenvectors
		
		std::cout << "Reading Fiber Orientations per Voxel Image"<< std::endl;
		foValue = fo_it.Get();
		unsigned int count = foValue.size(); //size of fiber orientation file
		VectorType RPD; //restricted principal direction
		VectorType FPD; //fiber principal direction which is normalized RPD
		if(count!=0){
			RPD[0] = foValue[i];//we take the x component of current fiber orientation vector
			RPD[1] = foValue[i+1];//we take the y component of current fiber orientation vector
			RPD[2] = foValue[i+2];//we take the z component of current fiber orientation vector
			double normPD =RPD.GetNorm();//fiber orientation's norm
			FPD = RPD/normPD;
		}
					
		/*we reset the signal to zero as we change of voxel*/	
		RealType signal=0;
		for( unsigned int grad_no = 0; grad_no < directions2.size(); grad_no++ )//For each gradient direction
		{	
			
			itk::Vector<double, 4> bktemp= directions2[grad_no];
			double DiffTime = bktemp[3];
			VectorType bkdirection;/*we store the current gradient direction in bkdirection*/
			for(int i=0;i<3;i++)
			{
				bkdirection[i]=bktemp[i];
			}
			//Initialization of signals
			//std::cout << "  Applying direction " << grad_no << " of " <<directions2.size()-1 << "): [" << bkdirection << "]" <<std::endl;
			
			RealType SignalHindered=0;
			RealType SignalRestrictedPa=0;
			RealType SignalRestrictedPe =0;
			RealType SignalRestricted=0;
					
			std::cout<<"new direction"<<std::endl;

			//we walk through the image containing fiber orientation vectors for each voxel, 

			DPa = 0.001;
			DPe = 0.00001;
			//there is fiber in this voxel
			//in the case of no fiber or 0,0,0 grad direction, hindered diffusion is zero
			SignalHindered = estimateHinderedDiffusion(EigenValue,bkdirection,WidthPulseGradient,DiffTime);
			if (SignalHindered!=1){
				std::cout<<"hindered diffusion is "<<SignalHindered<<std::endl;
			}
			int fiber_count = 0;
			if(count!=0 && bkdirection[0]!=0 && bkdirection[1]!=0 && bkdirection[2]!=0){
				SignalRestricted = 0;
				//go through all the fiber orientations
				for(unsigned int i = 0; i < count; i += 3){ /*for each RPD(=vector=fiber direction)*/
					VectorType FPD; //fiber principal direction which is normalized RPD
					RPD[0] = foValue[i];//we take the x component of current fiber orientation vector
					RPD[1] = foValue[i+1];//we take the y component of current fiber orientation vector
					RPD[2] = foValue[i+2];//we take the z component of current fiber orientation vector
					double normPD =RPD.GetNorm();//fiber orientation's norm
					FPD = RPD/normPD;
									

					//Projection of gradient direction on the vector representing fiber direction
					double dotproduct_bk_FPD = bkdirection*FPD;//dot product of gradient direction and the vector representing fiber 
					// q parallel and perpendicular estimation //
					double normd=bkdirection.GetNorm();//gradient direction's norm
					double qpa = fabs(dotproduct_bk_FPD);//projection result
					double qpe = sqrt((pow (normd,2))-(pow (qpa,2)));
					/*Estimation paralell component of restricted diffusion coefficient(DPa) and perpendicular component of restricted diffusion coefficient (DPe)*/
					

					/*estimation of the parallel restricted signal*/
					SignalRestrictedPa = vcl_exp((-4) * pow(pi,2) * (pow(qpa,2)) * (DiffTime - (WidthPulseGradient/3)) * DPa);
	
					/*estimation of the perpendicular restricted signal*/
					SignalRestrictedPe = vcl_exp(-(4 * pow(pi,2)*pow(Radius,4)*(pow(qpe,2))/DPe*t)*(7/96)*(2-(99/112)*(pow(Radius,2)/DPe*t)));
	
					/*estimation of the total restricted signal*/
					SignalRestricted += SignalRestrictedPa * SignalRestrictedPe;
					fiber_count++;	
				}
					
				//Final estimation of signal for one voxel//
				
				/*fR is re-normalized based on how many of restricted diffusion components are presented*/
				signal = (fR * SignalRestricted)/fiber_count+ fH*SignalHindered; //temporially get rid of hindered component
			}
			//there is no fiber in this voxel
			else{
				signal = SignalHindered; //1 for bkdirection = 0,0,0
			}
			
		
			dwi_val[grad_no]=signal * b0_img_it.Get();
			
			
			//std::cout<<"signal is dwi_val[0]--"<<dwi_val[0]<<"--end--"<< dwi_val[grad_no]<<std::endl;
			
			
			++fo_it;
			++dwi_template_it;
			++b0_img_it;
			++eigen_hindered_it;
			//////////////ADD RICIAN NOISE/////////////////////
			/* RealType realNoise = 0.0;
			RealType imagNoise = 0.0;
			if( noiseSigma > 0.0 )
			{
			realNoise = randomizer->GetNormalVariate( 0.0,
			vnl_math_sqr( noiseSigma ) );
			imagNoise = randomizer->GetNormalVariate( 0.0,
			vnl_math_sqr( noiseSigma ) );
			}
			RealType realSignal = signal + realNoise;
			RealType imagSignal = imagNoise;

			vcl_complex<RealType> noisySignal( realSignal, imagSignal );

			RealType finalSignal = vcl_sqrt( vcl_norm( noisySignal ) );*/
				
				
		}
	}
	dwi_template_it.Set(dwi_val);	
	/*we write the new dwi*/
	typedef itk::ImageFileWriter<VectorImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( dwi_template );
	writer->UseCompressionOn();
	writer->SetFileName(OutFilename);
	writer->Update();

	return 0;
}


/*specific function to take information from txt files created in the first pipeline : eigenImage, EigenHinderedImage, FiberOrientationImage*/
template <class ImagePointer>
int getFOImageHandler(ImagePointer &foImage, const char* foImgFileName){
    //Read the Fiber Orientation Image
    const unsigned int FODimension = 3;
    typedef std::vector< double > FOPixelType;
    typedef itk::Image< FOPixelType, FODimension > FOImageType;

    typedef itk::ImageFileWriter< FOImageType > FOWriterType;
    typedef itk::ImageRegionConstIterator< FOImageType > FOConstIteratorType;

    FOImageType::PointType foOrigin;
    FOImageType::SpacingType foSpacing;
    FOImageType::SizeType foSize;

    FOMeta *meta = new FOMeta();
    readFOMeta(foImgFileName, meta);

    for(int i = 0; i < VOLUME_DIMENSION; i++){
        foOrigin[i] = meta->origin[i];
        foSpacing[i] = meta->spacing[i];
        foSize[i] = meta->size[i];
    }
    
   // FOImageType::Pointer foImage = FOImageType::New();
    FOImageType::RegionType foRegion;
    foRegion.SetSize(foSize);
    foImage->SetSpacing(foSpacing);
    foImage->SetOrigin(foOrigin);
    foImage->SetRegions(foRegion);
    foImage->Allocate();

    FOImageType::IndexType pixelIndex;
    FOImageType::PixelType pixelValue;

    pixelIndex[0] = 1; pixelIndex[1] = 1; pixelIndex[2] = 1;
    pixelValue = foImage->GetPixel(pixelIndex);
    pixelValue.push_back(0.234);
    foImage->SetPixel(pixelIndex, pixelValue);

    ifstream inFile(foImgFileName);
    std::string line;
    bool headerFinish = false;
    bool keyFinish = false;
    unsigned int i, count;
    int j, k;
    char value[100];
	
    while(getline(inFile, line)){
        if(line.find("end header") != std::string::npos){
            headerFinish = true; 
        }
        if(!headerFinish){
            continue;
        }
        
        //Process voxels
        if(line.find("voxel") != std::string::npos){
            keyFinish = false;
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   pixelIndex[k++] = (unsigned int)atoi(value);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            pixelIndex[k] = (unsigned int)atoi(value);

            pixelValue = foImage->GetPixel(pixelIndex);

            //Compute size 
            getline(inFile, line);
            if(line.find("size") != std::string::npos){
                keyFinish = false;
                j = 0;
                for(i = 0; i < line.length(); i++){
                   if(line[i] == '='){
                       keyFinish = true;
                       continue;
                   }
                   if(!keyFinish){
                       continue;
                   }
                   value[j++] = line[i];
                }
                value[j] = '\0';
                count = (unsigned int)atoi(value);
            }
            while(count > 0){
                getline(inFile, line);
                pixelValue.push_back(strtod(line.c_str(), NULL));
                count--;
            }
            foImage->SetPixel(pixelIndex, pixelValue);
        }
    }
    return 0;
}

/*specific function to read information from txt files created in the first pipeline : eigenImage, EigenHinderedImage, FiberOrientationImage*/
int readFOMeta(const char* foImgFileName, FOMeta* m){
    ifstream inFile(foImgFileName);
    std::string line;
    bool keyFinish = false;
    unsigned int i;
    int j, k;
    char value[100];

    while(getline(inFile, line)){
        if(line.find("size") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->size[k++] = (unsigned int)atoi(value);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->size[k] = (unsigned int)atoi(value);
        }

        if(line.find("origin") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i<line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->origin[k++] = strtod(value, NULL);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->origin[k] = strtod(value, NULL);
        }
        if(line.find("spacing") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->spacing[k++] = strtod(value, NULL);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->spacing[k] = strtod(value, NULL);
        }
        if(line.find("end header") != std::string::npos){
            break;
        }
    }
    return 0;
}

double estimateHinderedDiffusion(std::vector< double > EigenValue, itk::Vector<double, 3> bkdirection,double WidthPulseGradient,double DiffTime)
{
	const unsigned int EigenDimension = 3;	
	double pi=M_PI;
	
	
	vnl_matrix<double> eigenVectormatrix(3,3);
	vnl_matrix<double> eigenValuematrix(3,3);
	vnl_matrix<double> eigenVectortransmatrix(3,3);
	
	vnl_vector<double> q(3);
	q = bkdirection.GetVnlVector();

	double HinderedSignal = 0; //initialize hindered signal
	//initialize
	fill(eigenVectormatrix.begin(), eigenVectormatrix.end(), 0.0);
	fill(eigenValuematrix.begin(), eigenValuematrix.end(), 0.0);
	fill(eigenVectortransmatrix.begin(), eigenVectortransmatrix.end(), 0.0);
		 
	//temporary solution is to compute the tensor based on the eigenvalue and eigenvactor
	
	unsigned int count= EigenValue.size();
	
	if(count>0 && bkdirection[0]!=0 && bkdirection[1]!=0 && bkdirection[2]!=0){
		for(unsigned int i = 0; i < 3; i++){
			/*we take the eigenvectors and store it in v1, v2, v3*/
			eigenValuematrix(i,i) = EigenValue[i];
			std::cout<<"eigenvalue("<<i<<","<<i<<") is "<<eigenValuematrix(i,i)<<std::endl;
			for(unsigned int j = 0; j < 3; j++){
				eigenVectormatrix(i,j) = EigenValue[3*i+j+3];
				eigenVectortransmatrix(j,i) = EigenValue[3*i+j+3];
			}
		}
		vnl_matrix<double> tensor(3,3);
		tensor = eigenVectormatrix*eigenValuematrix* eigenVectortransmatrix;
		/*we estimate hindered diffusion signal*/
		vnl_vector<double> temp_qD(3);
		temp_qD.post_multiply(tensor);
		double temp_qDq = dot_product(temp_qD,q);
		HinderedSignal = vcl_exp(( -4 )* pow(pi,2) * (DiffTime - (WidthPulseGradient/3)) * temp_qDq);
	}
	/*average tensor is zeros, attenuation is 1*/
	else{
		HinderedSignal = 1;
	}
	return HinderedSignal;
}






