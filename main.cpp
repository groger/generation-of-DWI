#include "mainCLP.h"
//testing
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
template <class ImagePointer>
int estimateHinderedDiffusion(ImagePointer &EigenHinderedImage,const char* EigenHinderedImgFileName, std::vector<itk::Vector<double, 2> > &HinderedD,itk::Vector<double, 3> FPD);
int generationdwi(string dwiImgFilename,string T2ImgFilename,string OutFilename,string EigenFilename,string foImgFilename,string EigenHinderedImgFilename,float Timetoecho,float DiffTime,double WidthPulseGradient,double MagnitudeG,float fH,float fR,double noiseSigma);


int main(int argc, char *argv[])
{
PARSE_ARGS;


generationdwi(dwiImgFilename,T2ImgFilename,OutFilename,EigenFilename,foImgFilename,EigenHinderedImgFilename,Timetoecho,DiffTime, WidthPulseGradient,MagnitudeG,fH, fR, noiseSigma);

return 0;
}


int generationdwi(string dwiImgFilename,string T2ImgFilename,string OutFilename,string EigenFilename,string foImgFilename,string EigenHinderedImgFilename,float Timetoecho,float DiffTime,double WidthPulseGradient,double MagnitudeG,float fH,float fR,double noiseSigma)
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
double lambdaPe=1;
double lambdaPa=1;
std::string Inputname = dwiImgFilename;
std::string T2Name = T2ImgFilename;
const char * EigenFileName=EigenFilename.c_str();
const char *EigenHinderedImgFileName=EigenHinderedImgFilename.c_str();
const char *foImgFileName=foImgFilename.c_str();	 
typedef float      PixelType;


 //Read dwi original and store gradients:
	
	typedef itk::Image< PixelType , 3 > ImageType ;
 	typedef itk::ImageFileReader< ImageType > FileReaderType ; 
	typedef itk::VectorImage< PixelType , 3 > VectorImageType ; 
	ImageType::Pointer image ;
	ImageType::IndexType Index;
 	std::vector< ImageType::Pointer > vectorOfImage ;
 	itk::MetaDataDictionary dico ;
	itk::VectorImage< PixelType, 3 >::Pointer olddwi ;
 	olddwi = itk::VectorImage< PixelType , 3 >::New() ; 
	
		itk::ImageFileReader< VectorImageType >::Pointer reader ;
		reader = itk::ImageFileReader< VectorImageType >::New() ;
		reader->SetFileName( Inputname) ;
		reader->Update() ;
		olddwi = reader->GetOutput();
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
ImageType::Pointer img = ImageType::New();
imageReader->SetFileName(  T2Name );
  try{
        imageReader->Update();
        img = imageReader->GetOutput();
}
  catch (itk::ExceptionObject &ex){
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }
itk::ImageRegionIterator<ImageType> img_it (img, img->GetLargestPossibleRegion());

//write the baseline 
typedef itk::Image< double, ImageDimension >  ScalarImageType;
ScalarImageType::Pointer b0 = ScalarImageType::New();
  b0->CopyInformation(img);
  b0->SetRegions(b0->GetLargestPossibleRegion());
  b0->Allocate();
  b0->Update();
  std::cout << b0 << std::endl;
  itk::ImageRegionIterator<ScalarImageType> itB0(b0, b0->GetLargestPossibleRegion());



     //signal estimation for each gradient direction and walking through each voxel of the new dwi
     for( unsigned int d = 0; d < directions2.size(); d++ )//For each gradient direction
      {
	itk::Vector<double, 4> bktemp= directions2[d];
	typedef itk::Vector<double, 3> VectorType;
	VectorType bkdirection;/*we store the current gradient direction in bkdirection*/
	for(int i=0;i<3;i++)
	{
	  bkdirection[i]=bktemp[i];
	}
	
	std::cout << "  Applying direction " << d << " of " <<directions2.size()-1 << "): [" << bkdirection << "]" <<std::endl;
      
	itk::Vector<double, 3> bk0;
	bk0[0]=0;bk0[1]=0;bk0[2]=0;

        /*file in which we will store hindered parameters*/
        std::cout << "Reading EigenValue Hindered part per Voxel Image"<< std::endl;
	const unsigned int EigenHinderedDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenHinderedDimension > EigenImageType;
	
	EigenImageType::Pointer EigenHinderedImage = EigenImageType::New();
	getFOImageHandler(EigenHinderedImage, EigenHinderedImgFileName); /*function to read txt file containing hindered diffusion coefficients*/
	std::vector<itk::Vector<double, 2> > HinderedD;//vector of vector to store hindered parameters
	itk::Vector<double, 2> f;f[0]=0;f[1]=0;
	std::fill( HinderedD.begin(), HinderedD.end(), f );
	

	/*Read FOImage, file containing fiber orientations for each voxel*/
	std::cout << "Reading Fiber Orientations per Voxel Image"<< std::endl;
	const unsigned int FODimension = 3;
	typedef std::vector< double > FOPixelType;
	typedef itk::Image< FOPixelType, FODimension > FOImageType;
	typedef itk::ImageRegionConstIterator< FOImageType > FOConstIteratorType;
	
	FOImageType::Pointer foImage = FOImageType::New();
	getFOImageHandler(foImage, foImgFileName); /*function to read txt file containing fiber orientation vectors*/
	FOConstIteratorType fo_it( foImage, foImage->GetLargestPossibleRegion());
	FOImageType::PixelType foValue;
	FOImageType::IndexType foIndex;
	fo_it.GoToBegin();
	
	itk::Vector<double, 3> RPD; //restricted principal direction
	unsigned int count;

	/*Read EigenValue Image,file containing fiber's eigenvalues for each voxel : restricted diffusion coefficients */
	std::cout << "Reading EigenValue per Voxel Image"<< std::endl;
	const unsigned int EigenDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	
	EigenImageType::Pointer EigenImage = EigenImageType::New();
	getFOImageHandler(EigenImage, EigenFileName);/*function to read txt file containing restricted diffusion coefficients*/ 
	EigenConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
	EigenImageType::PixelType EigenValue;
	EigenImageType::IndexType eigenIndex;
	Eigen_it.GoToBegin();
	typedef itk::Vector<double, 3> VectorType;
	VectorType eigen;
	
	//define iterator on the original dwi
	typedef itk::ImageRegionIterator< VectorImageType > IteratorType;
	IteratorType olddwi_it( olddwi, olddwi->GetLargestPossibleRegion().GetSize() );
	olddwi_it.GoToBegin();
	
	//Initialization of signals
	RealType SignalHindered=0;
	RealType SignalRestrictedPa=0;
	RealType SignalRestrictedPe =0;
	RealType SignalRestricted=0;
	RealType signal=0;
	double DiffTime = bktemp[3];
	
	int j=0;
	int compteur=0;

	std::cout<<"new direction"<<std::endl;

	//Write Baseline on gradient 0
	if(bkdirection==bk0){//baseline
		img_it.GoToBegin();
		while(!img_it.IsAtEnd() && !olddwi_it.IsAtEnd()){
		  
			itk::VariableLengthVector<PixelType> val=olddwi_it.Get();
			val[0]=img_it.Get();
			olddwi_it.Set(val);
			++img_it;
			++olddwi_it;
		}
	}




	else{//write in all other gradient directions than 0
			
		
/*for each voxel*/while(!fo_it.IsAtEnd() && !Eigen_it.IsAtEnd() && !olddwi_it.IsAtEnd()){//we walk through the image containing fiber orientation vectors for each voxel, the image containing restricted diffusion coefficients and the image containing hindered diffusion, the original dwi and we store the signal in the dwi(nrrd file) 
			
			Index=fo_it.GetIndex();
			foValue = fo_it.Get();
			count = foValue.size();
			EigenValue=Eigen_it.Get();
			int compt=0;

			if(count!=0){
		       	for(unsigned int i = 0; i < count; i += 3){ /*for each RPD(=vector=fiber direction)*/
				VectorType FPD; //fiber principal direction which is normalized RPD
				RPD[0] = foValue[i];//we take the x component of current fiber orientation vector
				RPD[1] = foValue[i+1];//we take the y component of current fiber orientation vector
				RPD[2] = foValue[i+2];//we take the z component of current fiber orientation vector
				double normPD =RPD.GetNorm();//fiber orientation's norm
				
				for (unsigned int temp_k = 0; temp_k<3;temp_k++){
					FPD[temp_k] = RPD[temp_k]/normPD;
					}				

				//Projection of gradient direction on the vector representing fiber direction
				double dotproduct_bk_FPD = bkdirection*FPD;//dot product of gradient direction and the vector representing fiber direction
				//v is the vector of a point from a fiber in a voxel
				estimateHinderedDiffusion(EigenHinderedImage,EigenHinderedImgFileName,HinderedD,FPD);
				itk::Vector<double, 2> lambda = HinderedD[j];
				lambdaPa=lambda[0];/*Hindered Diffusion parallel component*/
				lambdaPe=lambda[1];/*Hindered Diffusion perpendicular component*/
				// q parallel and perpendicular estimation //
				double normd=bkdirection.GetNorm();//gradient direction's norm
				double qpa = fabs(dotproduct_bk_FPD);//projection result
				double qpe = sqrt((pow (normd,2))-(pow (qpa,2)));//estimation of the vector perpendicular to projection of gradient direction
				
				/*std::cout << "qpa: "<<qpa<<" qpe: "<<qpe<<"  "<<std::endl;*/
				//Take restricted parameters from EigenFile
				eigen[0] =EigenValue[i];
				eigen[1] =EigenValue[i+1];
				eigen[2] =EigenValue[i+2];

				/*Estimation paralell component of restricted diffusion coefficient(DPa) and perpendicular component of restricted diffusion coefficient (DPe)*/
				//DPa=eigen[0];
				DPa = 0.001;
				DPe = 0.00001;
				//DPe=((eigen[1]+eigen[2])/2);

				/*prints for test*/
				if(compt==0 ){
					std::cout << "DPa: "<<DPa<<" DPe:   "<<DPe<<"  "<<std::endl;
					std::cout <<"Difftime : "<<DiffTime<<std::endl;
					std::cout<<"Width pulse gradient : "<<WidthPulseGradient<<std::endl;
					std::cout<<"t : "<<t<<std::endl;
					std::cout<<"Radius : "<<Radius<<std::endl;
					std::cout << "lambdaPa "<<lambdaPa<< std::endl;
					std::cout << "lambdaPe "<<lambdaPe<< std::endl;
					}
				
				


				/*estimation of the parallel restricted signal*/
				SignalRestrictedPa = vcl_exp((-4) * pow(pi,2) * (pow(qpa,2)) * (DiffTime - (WidthPulseGradient/3)) * DPa);

				/*estimation of the perpendicular restricted signal*/
				SignalRestrictedPe = vcl_exp(-(4 * pow(pi,2)*pow(Radius,4)*(pow(qpe,2))/DPe*t)*(7/96)*(2-(99/112)*(pow(Radius,2)/DPe*t)));

				/*estimation of the total restricted signal*/
				SignalRestricted += SignalRestrictedPa * SignalRestrictedPe;
				
				/*Estimation of hindered signal*/
				SignalHindered += vcl_exp(( -4 )* pow(pi,2) * (DiffTime - (WidthPulseGradient/3)) * ((pow(qpa,2)) * lambdaPa + (pow(qpe,2)) * lambdaPe));
	 
				++compt;
				}
				
				//Final estimation of signal for one voxel//
				if (compt == 0){
					signal = SignalHindered/compt;
				}
				else{	
					std::cout << "SignalR: "<< SignalRestricted<<std::endl;
					std::cout << "there are "<< compt<< " fiber bundles in this voxel"<<std::endl;
					std::cout << "SignalH: "<< SignalHindered<<std::endl;
					/*fR is re-normalized based on how many of restricted diffusion components are presented*/
					signal = (fH * SignalHindered + fR * SignalRestricted)/compt; 
				}
			 
		        }
			
			else{
				signal =0;
			}//empty voxel so signal = 0
			//on remplit la composante numero "componentToInsert"(=numero du gradient) du vecteur contenu dans le pixel numero Index
			int componentToInsert=d;
			itk::VariableLengthVector<PixelType> val=olddwi_it.Get();
			val[componentToInsert]=signal;

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
				
				  
			//We set the value of final signal  
			olddwi_it.Set(val);
				
			++fo_it;
			++olddwi_it;
			++Eigen_it;
			++j;
			if(signal!= 0)
			{ ++compteur;}
			
			/*we reset the signal to zero as we change of voxel*/	
			SignalHindered=0;
			SignalRestrictedPa=0;
			SignalRestrictedPe =0;
			SignalRestricted=0;
			signal=0;
		}
	}
				
	}
/*we write the new dwi*/
typedef itk::ImageFileWriter<VectorImageType> WriterType;
WriterType::Pointer writer = WriterType::New();
writer->SetInput( olddwi );
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

template <class ImagePointer>
int estimateHinderedDiffusion(ImagePointer &EigenHinderedImage,const char* EigenHinderedImgFileName, std::vector<itk::Vector<double, 2> > &HinderedD,itk::Vector<double, 3> FPD){

	HinderedD.clear();
	const unsigned int EigenDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	EigenImageType::IndexType eigenIndex;
	EigenConstIteratorType eigen_hindered_it( EigenHinderedImage, EigenHinderedImage->GetLargestPossibleRegion());
	EigenImageType::PixelType EigenValue;
	eigen_hindered_it.GoToBegin();
	itk::Vector<double, 2> hinderedD; 
	hinderedD[0]=0;hinderedD[1]=0;/*initialization*/
	itk::Vector<double, 2> temp;
	itk::Vector<double, 3> v1;
	itk::Vector<double, 3> v2;
	itk::Vector<double, 3> v3;
	int count2=0;
	int compteur=0;
	while(!eigen_hindered_it.IsAtEnd()){/*walk through txt file containing eigenvalues and eigenvectors of average tensor*/
		EigenValue=eigen_hindered_it.Get();
		unsigned int count= EigenValue.size();
		if(count>0 && FPD[0]!=0 && FPD[1]!=0 && FPD[2]!=0){
			for(unsigned int i = 0; i < 3; i++){
			  /*we take the eigenvectors and store it in v1, v2, v3*/
			  v1[i]=EigenValue[i+3];
			  v2[i]=EigenValue[i+6];
			  v3[i]=EigenValue[i+9]; 
			} 
		
			++compteur;
			/*we estimate hindered diffusion coefficient : the parallel one :temp[0] and the perpendicular one : temp[1]*/
			//HINDERED DIFFUSION coefficient projected along the gradient directions
			//
			temp[0]=EigenValue[0]*(fabs(FPD*v1))+EigenValue[1]*(fabs(FPD*v2))+EigenValue[2]*(fabs(FPD*v3));
			temp[1]=sqrt(((pow(EigenValue[0],2)+pow(EigenValue[1],2)+pow(EigenValue[2],2))-pow(temp[0],2))/2);
		}
		else{hinderedD[0]=0;hinderedD[1]=0;temp[0]=0;temp[1]=0;}
		HinderedD.push_back(temp);/*we pushback hindered diffusion coefficients in HinderedD*/
		++eigen_hindered_it;++count2;
	}		
	return 0;
}






