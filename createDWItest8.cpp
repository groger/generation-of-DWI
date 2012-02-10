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
//#include "dll.h"
#include <vxl_config.h>
#include <vnl/vnl_config.h>
//#include <vnl_numeric_traits.h>
//#include <vnl_complex.h>

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
int estimateHinderedDiffusion(ImagePointer &EigenHinderedImage,const char* EigenHinderedImgFileName, std::vector<itk::Vector<double, 2> > &HinderedD,itk::Vector<double, 3> bkdirection);
//Fiber Orientation Image meta info struct 



int main(int argc, char *argv[])
{

//We define dwi parameters
typedef double RealType;
//const char *T2ImgFileName= "/home/gwendo/taf/FICHIERS_CHARMED/T2_Gwen.nrrd";
const char *T2ImgFileName= "/home/gwendo/taf/FICHIERS_CHARMED/GwenDWI_QCed_B0.nrrd";
const char *foImgFileName = "/home/gwendo/taf/FICHIERS_CHARMED/fiberOrientationsPerVoxelTEST26janv.txt";
const char *dwiImgFileName = "/home/gwendo/taf/FICHIERS_CHARMED/GwenDWI_QCed.nhdr";
const char *EigenFileName = "/home/gwendo/taf/FICHIERS_CHARMED/EigenValuePerVoxelTEST26janv.txt";
const char *EigenHinderedImgFileName="/home/gwendo/taf/FICHIERS_CHARMED/EigenValuePerVoxelHinderedPartABSNegSingVal31janv.txt";
//float TE=0.2;
float TE=0.12;
//float DiffTime = 0.15;
float DiffTime =0.00476;
//double WidthPulseGradient=0.04;
//double WidthPulseGradient=0.05528;
double MagnitudeG=26.98;
float fH=0.7;
float fR=0.3;
float Radius= 0.002;
double gyro=42.576;
double pi=M_PI;
//double multGrad=(gyro*WidthPulseGradient)/(2*pi);

double t=TE/2;
double DPa;
double DPe;
double lambdaPe=1;
double lambdaPa=1;

std::string Inputname = dwiImgFileName;
std::string T2Name = T2ImgFileName;
std::string OutFileName = "/home/gwendo/taf/FICHIERS_CHARMED/testNewCorrespondanceVoxel08febrABSandNewEstimationOfq.nrrd";	 
typedef float      PixelType;
//typedef itk::Vector< double, 36 >    PixelType;

 //Read dwi original and store gradients:
	
	typedef itk::Image< PixelType , 3 > ImageType ;
 	typedef itk::ImageFileReader< ImageType > FileReaderType ; 
	typedef itk::VectorImage< PixelType , 3 > VectorImageType ; 
	/*typename */ImageType::Pointer image ;
	ImageType::IndexType Index;
 	std::vector< /*typename*/ ImageType::Pointer > vectorOfImage ;
 	itk::MetaDataDictionary dico ;
	itk::VectorImage< PixelType, 3 >::Pointer olddwi ;
 	olddwi = itk::VectorImage< PixelType , 3 >::New() ; 
	
		/*typename*/ itk::ImageFileReader< VectorImageType >::Pointer reader ;
		reader = itk::ImageFileReader< VectorImageType >::New() ;
		reader->SetFileName( Inputname) ;
		reader->Update() ;
		olddwi = reader->GetOutput();
		//Save metadata dictionary
		dico = reader->GetOutput()->GetMetaDataDictionary() ;
		//on recupere le vecteur directions de :TransformGradients(dico)

  //define variables to add Rician Noise
  RealType noiseSigma = 0;
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
  typename RandomizerType::Pointer randomizer = RandomizerType::New();
  randomizer->Initialize();
  
  
//int TransformGradients( itk::MetaDataDictionary &dico){//on recupere les gradients dans directions
  //on recupere les gradients du dwi original en passant son metadatadictionnary puis on le stocke dans un vecteur de vecteur
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
		else if( pos != -1 )
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
 std::cout<<"max b value: "<<maxBValue<<std::endl;
 
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
std::cout<<"norme max des gradients"<<normGradMax<<std::endl;

/*estimate b of each gradient and the WidthPulseGradient correspondant and store it */
std::vector<itk::Vector<double, 4> > directions2;
  
for (int i=0; i <numGradients;i++)
{
  itk::Vector<double, 3> bi = directions[i];
  double b=(pow((bi.GetNorm()),2)/pow(normGradMax,2))*maxBValue;
  std::cout<<pow((bi.GetNorm()),2)/pow(normGradMax,2)<<std::endl;
  std::cout<<"b with g normed and multiplied by maxbvalue "<<b<<std::endl;
  double   r1, r2, r3;/*roots*/
  r1=0;r2=0;r3=0;/*initialize roots at zero*/
  int    num_roots; /* Number Of Unique Roots To Equation */
  /*declaration of the coefficient of cubic equation to solve to find WidthPulseGradient*/
 /* double c1=(1/3); 
  double c2=(-1)*DiffTime;
  double c3=0;
  double c4=b/(pow(gyro,2)*pow(MagnitudeG,2));*/
   double c1=1; 
  double c2=(-1);
  double c3=(-4);
  double c4=4;
  //vtkPolynomialSolversUnivariate::SolveCubic( c1, c2, c3, c4, &r1, &r2, &r3, &num_roots );/*solve the equation to estimate WidthPulseGradient from b value*/
  int result=vtkMath::SolveCubic( c1, c2, c3, c4, &r1, &r2, &r3, &num_roots );
  double WidthPulseGradient=0;
std::cout<<"result width pulse gradient : "<<r1<<" "<<r2<<" "<<r3<<" "<<num_roots<<" "<<result<<std::endl;
  /*WidthPulseGradient is the positive solution*/
  if (r1>0){ WidthPulseGradient=r1;}
  if (r2>0 && r2>r1){ WidthPulseGradient=r2;}
  if (r3>0 && r3>r2){ WidthPulseGradient=r3;}
  //std::cout<<"width pulse gradient "<<WidthPulseGradient<<std::endl;
  /*estimation of the gradient normed with its magnitude*/
   itk::Vector<double, 4> direction2;
   for (int j=0;j<3;j++)
   {
     if(bi.GetNorm()!=0)
     { //std::cout<<"g normed "<<bi[j]/(bi.GetNorm())<<std::endl;
       //std::cout<<"to multiplied: "<<MagnitudeG*WidthPulseGradient*gyro<<std::endl;
	direction2[j]=(bi[j]/(bi.GetNorm()))*MagnitudeG*WidthPulseGradient*gyro;
	//std::cout<<"g normed multiplied by ..."<<direction2[j]<<std::endl;
     }
     else
     {
       direction2[j]=0;
     }
   }
   direction2[3]=WidthPulseGradient;
   directions2.push_back(direction2);/*we store in e vector vectors containing the three components of the normed gradient and at the end the WidthPulseGradient*/
} 
   
//estimate WidthPulseGradient from b value
//std::cout<<"b value  "<<b_values[0]<<std::endl;
//double   r1, r2, r3;/*roots*/
//r1=0;r2=0;r3=0;/*initialize roots at zero*/
//int    num_roots; /* Number Of Unique Roots To Equation */
//double c1=(1/3); double c2=(-1)*DiffTime;double c3=0;double c4=b_values[0]/(pow(gyro,2)*pow(MagnitudeG,2));
//vtkPolynomialSolversUnivariate::SolveCubic( c1, c2, c3, c4, &r1, &r2, &r3, &num_roots );/*solve the equation to estimate WidthPulseGradient from b value*/
//double WidthPulseGradient=0;

/*WidthPulseGradient is the positive solution*/
//if (r1>0){ WidthPulseGradient=r1;}
//if (r2>0){ WidthPulseGradient=r2;}
//if (r3>0 && num_roots==3){ WidthPulseGradient=r3;}
//double multGrad=(gyro*WidthPulseGradient); 
  
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


 
	


     //signal estimation for each voxel of the new dwi
     for( unsigned int d = 0; d < directions2.size(); d++ )//For each gradient direction
      {
	itk::Vector<double, 4> bktemp= directions2[d];
	itk::Vector<double, 3> bkdirection;
	for(int i=0;i<3;i++)
	{
	  bkdirection[i]=bktemp[i];
	}
	
	std::cout << "  Applying direction " << d << " of " <<directions2.size()-1 << "): [" << bkdirection << "]" <<std::endl;
	itk::Vector<double, 3> bk0;
	bk0[0]=0;bk0[1]=0;bk0[2]=0;

        //file in which we will store hindered parameters
        std::cout << "Reading EigenValue Hindered part per Voxel Image"<< std::endl;
	const unsigned int EigenHinderedDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenHinderedDimension > EigenImageType;
	//typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	
	EigenImageType::Pointer EigenHinderedImage = EigenImageType::New();
	getFOImageHandler(EigenHinderedImage, EigenHinderedImgFileName); 
	//EigenConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
	//FOImageType::IndexType foIndex;
	//EigenImageType::PixelType EigenValue;
	/*Eigen_it.GoToBegin();*/
	/*typedef itk::Vector<double, 3> VectorType;
	VectorType eigen;*/
	std::vector<itk::Vector<double, 2> > HinderedD;//vector of vector to store hindered parameters
	itk::Vector<double, 2> f;f[0]=0;f[1]=0;
	std::fill( HinderedD.begin(), HinderedD.end(), f );
	estimateHinderedDiffusion(EigenHinderedImage,EigenHinderedImgFileName,HinderedD,bkdirection);
	
	
	
	
	//Read FOImage, file containing fiber orientations for each voxel
	std::cout << "Reading Fiber Orientations per Voxel Image"<< std::endl;
	const unsigned int FODimension = 3;
	typedef std::vector< double > FOPixelType;
	typedef itk::Image< FOPixelType, FODimension > FOImageType;
	typedef itk::ImageRegionConstIterator< FOImageType > FOConstIteratorType;
	
	FOImageType::Pointer foImage = FOImageType::New();
	getFOImageHandler(foImage, foImgFileName); 
	FOConstIteratorType fo_it( foImage, foImage->GetLargestPossibleRegion());
	//FOImageType::IndexType foIndex;
	FOImageType::PixelType foValue;
	FOImageType::IndexType foIndex;
	fo_it.GoToBegin();
	typedef itk::Vector<double, 3> VectorType;
	VectorType v;
	unsigned int count;

	//Read EigenValue Image,file containing fiber's eigenvalues for each voxel 
	std::cout << "Reading EigenValue per Voxel Image"<< std::endl;
	const unsigned int EigenDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	
	EigenImageType::Pointer EigenImage = EigenImageType::New();
	getFOImageHandler(EigenImage, EigenFileName); 
	EigenConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
	//FOImageType::IndexType foIndex;
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
	double WidthPulseGradient = bktemp[3];
	//std::cout<<"verif WidthPulseGradient "<<WidthPulseGradient<<std::endl;
	int j=0;
	int compteur=0;

	std::cout<<"new direction"<<std::endl;
	
	/*while(!fo_it.IsAtEnd() && !Eigen_it.IsAtEnd() && !olddwi_it.IsAtEnd()){
			Index=fo_it.GetIndex();std::cout << "foindex "<<Index<<std::endl;
			Index=Eigen_it.GetIndex();std::cout << "eigenindex "<<Index<<std::endl;
			Index=olddwi_it.GetIndex();std::cout << "dwiindex "<<Index<<std::endl;
				++fo_it;
				++olddwi_it;
				++Eigen_it;
}*/

	//Write Baseline on gradient 0
	if(bkdirection==bk0){//baseline
		//itB0.GoToBegin(); std::cout << "hello2 "<< std::endl;
		img_it.GoToBegin();
		//while(!itB0.IsAtEnd() && !olddwi_it.IsAtEnd()){//std::cout << "hello22 "<< std::endl;
		while(!img_it.IsAtEnd() && !olddwi_it.IsAtEnd()){
		  
			itk::VariableLengthVector<PixelType> val=olddwi_it.Get();
			val[0]=img_it.Get();
			//if(val[0]!=0)
			  //std::cout << "B0 val: " << val[0] << std::endl;
			olddwi_it.Set(val);
			++img_it;
			++olddwi_it;
		}
	}




	else{//write in all other gradient directions than 0
			
			std::cout << "hello "<< std::endl;
/*for each voxel*/while(!fo_it.IsAtEnd() && !Eigen_it.IsAtEnd() && !olddwi_it.IsAtEnd()){//on parcourt  l'image contenant les vecteurs directions des fibres par voxel et l'image contenant les eigenvalues et on inscrit le signal dans le dwi par voxel
			//std::cout << "new voxel "<<endl;
			Index=fo_it.GetIndex();
			/*std::cout << "foindex "<<Index<<std::endl;
			if(Index[0]==1 && Index[1]==1 && Index[2]==1)
			{Index=Eigen_it.GetIndex();std::cout << "eigenindex "<<Index<<std::endl;
			Index=olddwi_it.GetIndex();std::cout << "dwiindex "<<Index<<std::endl;
			 std::cout << "size: "<<count<< std::endl;} */
			
			foValue = fo_it.Get();
			count = foValue.size();
			//if (count!=0){std::cout << "size: "<<count<< std::endl;}
			itk::Vector<double, 2> lambda = HinderedD[j];
			lambdaPa=lambda[0];//if(lambdaPa!=0){std::cout << "lambdaPa "<<lambdaPa<< std::endl;}//Hindered Diffusion parallel component
			lambdaPe=lambda[1];//if(lambdaPe!=0){std::cout << "lambdaPe "<<lambdaPe<< std::endl;}//Hindered Diffusion perpendicular component
			EigenValue=Eigen_it.Get();
			//unsigned int count2= EigenValue.size();
			//if(count2!=0){std::cout << "hello22 "<< std::endl;std::cout << "size2: "<<count2<< std::endl;}
			//unsigned int count2= EigenValue.size();
			//if(count2!=0){std::cout << "hello22 "<< std::endl;std::cout << "size2: "<<count2<< std::endl;}
			//Index=olddwi_it.GetIndex();//on recupere l'index du voxel en cours
			int compt=0;
			//if(lambdaPa!=0 || (lambdaPe!=0)){std::cout << "lambdaPa: "<<lambdaPa<< "lambdaPe: "<<lambdaPe<< std::endl;}
		       if(count!=0){	
		       for(unsigned int i = 0; i < count; i += 3){ /*for each v(=vector=fiber direction)*/
				
				//if(foValue[i] != 0 && foValue[i+1] != 0 && foValue[i+2] != 0){
					v[0] = foValue[i];//on recupere la valeur de x de chaque vecteur
					//std::cout << "hello1 "<<v[0]<< std::endl;
					v[1] = foValue[i+1];//on recupere la valeur de y de chaque vecteur
					//std::cout << "hello2 "<<v[1]<< std::endl;
					v[2] = foValue[i+2];//on recupere la valeur de z de chaque vecteur
					//std::cout << "hello22 "<<v[2]<< std::endl;

					//std::cout<<"verif gradient direction "<<bkdirection[0]<<" "<<bkdirection[1]<<" "<<bkdirection[2]<<std::endl;
					//Projection of gradient direction on the vector representing fiber direction
					double dotproduct = bkdirection*v;//dot product of gradient direction and the vector representing fiber direction//v is the vector of a point from a fiber in a voxel
					double normd=bkdirection.GetNorm();//gradient direction's norm
					//std::cout << "norm grad verif: "<<normd<<std::endl;
					//double norm2=sqrt( (pow(bk[0],2))+(pow(bk[1],2))+(pow(bk[2],2)));
					//std::cout << "q value 2: "<<norm2<<std::endl;
					double normv =v.GetNorm();//fiber orientation's norm
		
					// q parallel and perpendicular estimation //
					//VectorType v1;VectorType v2;VectorType v3;
					// v1[0]=-0.998972 ;v1[1]=0.0339967 ;v1[2]=-0.0299989;v2[0]=-0.998922;v2[1]=0.0355988;v2[2]=-0.0298004;v3[0]=-0.998864;v3[1]=0.0377998;v3[2]=-0.0290012;
					double qpa= fabs(dotproduct)/normv;//projection result
					//std::cout<<" tempqpa "<<tempqpa<<std::endl;
					/*if( d==1 && (v==v1 || v==v2 || v==v3) )
					{
					  std::cout<<"v "<<v[0]<<v[1]<<v[2]<<std::endl;
					  std::cout<<" tempqpa "<<tempqpa<<std::endl;
					}*/
					//VectorType tempqpe= directions[d]-tempqpa;
					double qpe=sqrt((pow (normd,2))-(pow (qpa,2)));//estimation of the vector perpendicular to projection of gradient direction
					//std::cout<<" tempqpa "<<tempqpa<<std::endl;
					
					//if (compteur==1 && i==0){
					
					
					//Take restricted parameters from EigenFile
					eigen[0] =EigenValue[i];
					eigen[1] =EigenValue[i+1];
					eigen[2] =EigenValue[i+2];
					/*if(eigen[0]==0.234){Index=Eigen_it.GetIndex();std::cout << "eigenindex "<<Index<<std::endl;
					  std::cout<<"eigenvalues : "<<eigen[0]<<"	"<<eigen[1]<<"	"<<eigen[2]<<std::endl;}*/
					//Estimation Dparallel and Dperpendicular//
					DPa=eigen[0];
					DPe=(eigen[1]+eigen[2])/2;
					if(compt==0 ){std::cout << "DPa: "<<DPa<<" DPe:   "<<DPe<<"  "<<std::endl;
					  std::cout <<"Difftime : "<<DiffTime<<std::endl;
					std::cout<<"Width pulse gradient : "<<WidthPulseGradient<<std::endl;
					std::cout<<"t : "<<t<<std::endl;
					std::cout<<"Radius : "<<Radius<<std::endl;
					std::cout << "lambdaPa "<<lambdaPa<< std::endl;
					std::cout << "lambdaPe "<<lambdaPe<< std::endl;
					std::cout << "qpa: "<<qpa<<" qpe: "<<qpe<<"  "<<std::endl;}
					//b=(-4)*pi*pi*(pow(normd,2))*(DiffTime - (WidthPulseGradient/3))*multGrad*multGrad;
					//b=(-4)*pi*pi*(pow(qpa,2)+pow(qpe,2));
					//std::cout<<"test b value"<<b<<std::endl;
		
					//Estimation of signal hindered and restricted//
					SignalHindered += vcl_exp(( -4 )* pow(pi,2) * (DiffTime - (WidthPulseGradient/3)) * ((pow(qpa,2)) * lambdaPa + (pow(qpe,2)) * lambdaPe));
					if(compt==0){
					  std::cout <<"For Hindered signal : "<<"  pow(pi,2): "<<pow(pi,2)<<" (DiffTime - (WidthPulseGradient/3) : "<<(DiffTime - (WidthPulseGradient/3))<<"	(pow(qpa,2)) : "<<(pow(qpa,2))<<"	(pow(qpe,2)) : "<<(pow(qpe,2))<<std::endl;
					  std::cout<<"	(pow(qpa,2)) * lambdaPa : "<<((pow(qpa,2)) * lambdaPa)<<"	(pow(qpe,2)) * lambdaPe) : "<<((pow(qpe,2)) * lambdaPe)<<std::endl;
					  std::cout<<"everything multiplied : "<<(( -4 )* pow(pi,2) * (DiffTime - (WidthPulseGradient/3)) * ((pow(qpa,2)) * lambdaPa + (pow(qpe,2)) * lambdaPe))<<std::endl;
					  std::cout << "SignalH: "<< SignalHindered<<std::endl;}
					SignalRestrictedPa = vcl_exp((-4) * pow(pi,2) * (pow(qpa,2)) * (DiffTime - (WidthPulseGradient/3)) * DPa);
					SignalRestrictedPe = vcl_exp(-(4 * pow(pi,2)*pow(Radius,4)*(pow(qpe,2))/DPe*t)*(7/96)*(2-(99/112)*(pow(Radius,2)/DPe*t)));
					SignalRestricted += SignalRestrictedPa * SignalRestrictedPe;
					if(compt==0){std::cout << "SignalR: "<< SignalRestricted<<std::endl;
					  std::cout<<"                        "<<std::endl;
					}
					//} 
			++compt;
				}
				
				//Final estimation of signal for one voxel//
				signal = fH *SignalHindered + fR * SignalRestricted;//std::cout << "Signal: "<< signal<<std::endl;
			 
		      }
			
			else{signal =0;/*std::cout << "Signal: "<< signal<<std::endl;*/}//empty voxel so signal = 0
				//std::cout << "Signal: "<< signal<<std::endl;
				
				//on remplit la composante numero "componentToInsert"(=numero du gradient) du vecteur contenu dans le pixel numero Index
				int componentToInsert=d;
				itk::VariableLengthVector<PixelType> val=olddwi_it.Get();
				val[componentToInsert]=signal;
				//if(signal != 0)
				 // std::cout << "Signal: " << signal << std::endl;
				
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
				
				SignalHindered=0;
				SignalRestrictedPa=0;
				SignalRestrictedPe =0;
				SignalRestricted=0;
				signal=0;
		}
	}
				
	}
typedef itk::ImageFileWriter<VectorImageType> WriterType;
WriterType::Pointer writer = WriterType::New();
writer->SetInput( olddwi );
writer->UseCompressionOn();
writer->SetFileName(OutFileName);
writer->Update();


return 0;
}





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
int estimateHinderedDiffusion(ImagePointer &EigenHinderedImage,const char* EigenHinderedImgFileName, std::vector<itk::Vector<double, 2> > &HinderedD,itk::Vector<double, 3> bkdirection){
	HinderedD.clear();
	std::cout<< "je suis dedans "<< std::endl;
	std::cout << "Reading EigenValue per Voxel Image"<< std::endl;
	const unsigned int EigenDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	EigenImageType::IndexType eigenIndex;
	EigenConstIteratorType eigen_hindered_it( EigenHinderedImage, EigenHinderedImage->GetLargestPossibleRegion());
	EigenImageType::PixelType EigenValue;
	eigen_hindered_it.GoToBegin();
	itk::Vector<double, 2> hinderedD;
	itk::Vector<double, 2> temp;
	itk::Vector<double, 3> v1;
	itk::Vector<double, 3> v2;
	itk::Vector<double, 3> v3;
	int count2=0;
	int compteur=0;
	double lambdaPe=1;
  	double lambdaPa=1;
	while(!eigen_hindered_it.IsAtEnd()){//std::cout<< "je suis dans eigen "<< std::endl;
		EigenValue=eigen_hindered_it.Get();
		unsigned int count= EigenValue.size();hinderedD[0]=0;hinderedD[1]=0;
		//if(count!=0){std::cout<<"count:"<<count<< std::endl;}
		if(count>0 && bkdirection[0]!=0 && bkdirection[1]!=0 && bkdirection[2]!=0 ){//std::cout<< "count ok "<< std::endl;
			for(unsigned int i = 0; i < 3; i++){
			  
			  v1[i]=EigenValue[i+3];
			  v2[i]=EigenValue[i+6];
			  v3[i]=EigenValue[i+9]; 
			} 
		/*if (EigenValue[0]==0.234){//std::cout<<"direction gradient : "<<bkdirection[0]<<"	"<<bkdirection[1]<<"	"<<bkdirection[2]<<std::endl;
		  eigenIndex=eigen_hindered_it.GetIndex();std::cout << "eigenHinderedindex "<<eigenIndex<<std::endl;
		  std::cout<<"EigenValue[0]: "<<EigenValue[0]<<"EigenValue[1]: "<<EigenValue[1]<<"EigenValue[2]: "<<EigenValue[2]<<std::endl;}*/
		
		++compteur;
		//std::cout<<"lambda1 "<<EigenValue[0]<<"lambda2 "<<EigenValue[1]<<"lambda3 "<<EigenValue[2]<<std::endl;	
		//std::cout<<"bk direction"<<bkdirection[0]<<" "<<bkdirection[1]<<" "<<bkdirection[2]<<std::endl;
		//std::cout<<"v1 "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<std::endl;
		//std::cout<<"v2 "<<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<std::endl;
		//std::cout<<"v3 "<<v3[0]<<" "<<v3[1]<<" "<<v3[2]<<std::endl;
		temp[0]=EigenValue[0]*(fabs(bkdirection*v1))/(bkdirection.GetNorm())+EigenValue[1]*(fabs(bkdirection*v2))/(bkdirection.GetNorm())+EigenValue[2]*(fabs(bkdirection*v3))/(bkdirection.GetNorm());
		//std::cout<<"first term "<<pow(temp[0],2)<<std::endl;
		//std::cout<<"norm grad "<<bkdirection.GetNorm()<<std::endl;
		//std::cout<<"second term "<<(bkdirection.GetNorm())*(pow(EigenValue[0],2)+pow(EigenValue[1],2)+pow(EigenValue[2],2))<<std::endl;
		//temp[1]=sqrt((((bkdirection.GetNorm())*(pow(EigenValue[0],2)+pow(EigenValue[1],2)+pow(EigenValue[2],2)))-pow(temp[0],2))/2);*/
		//std::cout<<"premier calcul: "<<temp[0]<<std::endl;
		/*if(pow(temp[0],2)>(pow(EigenValue[0],2)+pow(EigenValue[1],2)+pow(EigenValue[2],2))){std::cout<<"c'est la merde "<<"EigenValue[0]: "<<EigenValue[0]<<"EigenValue[1]: "<<EigenValue[1]<<"EigenValue[2]: "<<EigenValue[2]<<std::endl;
		  std::cout<<"EigenValue[3]: "<<EigenValue[3]<<"EigenValue[4]: "<<EigenValue[4]<<"EigenValue[5]: "<<EigenValue[5]<<std::endl;
		std::cout<<"EigenValue[6]: "<<EigenValue[6]<<"EigenValue[7]: "<<EigenValue[7]<<"EigenValue[8]: "<<EigenValue[8]<<std::endl;
		}*/
		temp[1]=sqrt(((pow(EigenValue[0],2)+pow(EigenValue[1],2)+pow(EigenValue[2],2))-pow(temp[0],2))/2);
		double v11[3];double v21[3];double v31[3];double direction1[3];double cross1[3];double cross2[3];double cross3[3];
		for(int i=0;i<3;i++){
		  v11[i]= v1[i];
		  v21[i]= v2[i];
		  v31[i]= v3[i];
		  direction1[i]= bkdirection[i];}
		vtkMath::Cross( v11,direction1,cross1);
		vtkMath::Cross( v21,direction1,cross2);
		vtkMath::Cross( v31,direction1,cross3);
		double cross1norm=sqrt(pow(cross1[0],2)+pow(cross1[1],2)+pow(cross1[2],2));double cross2norm=sqrt(pow(cross2[0],2)+pow(cross2[1],2)+pow(cross2[2],2));double cross3norm=sqrt(pow(cross3[0],2)+pow(cross3[1],2)+pow(cross3[2],2));
		//std::cout<<"norm "<<cross1norm<<" "<<cross2norm<<" "<<cross3norm<<std::endl;
		//double temp1=EigenValue[0]*cross1norm+EigenValue[1]*cross2norm+EigenValue[2]*cross3norm;
		//std::cout<<"premier calcul: "<<temp[1]<<"deuxieme calcul (cross product): "<<temp1<<std::endl;
		//std::cout<<"deuxieme calcul: "<<temp[1]<<std::endl;}
		}
		else{hinderedD[0]=0;hinderedD[1]=0;temp[0]=0;temp[1]=0;}
		//if(temp[0]!=0 && temp[1]!=0){std::cout<< "ok2 "<<temp[0]<<temp[1]<<std::endl;}
		HinderedD.push_back(temp);
		itk::Vector<double, 2> lambda = HinderedD[count2];
		lambdaPa=lambda[0];//if(lambdaPa!=0){std::cout << "hello21 "<<lambdaPa<< std::endl;}
		lambdaPe=lambda[1];//if(lambdaPe!=0){std::cout << "hello31 "<<lambdaPe<< std::endl;}
		++eigen_hindered_it;++count2;
	}		
return 0;
}

/*std::cout << "Reading EigenValue per Voxel Image"<< std::endl;
	const unsigned int EigenDimension = 3;
	typedef std::vector< double > EigenPixelType;
	typedef itk::Image< EigenPixelType, EigenDimension > EigenImageType;
	typedef itk::ImageRegionConstIterator< EigenImageType > EigenConstIteratorType;
	
	EigenImageType::Pointer EigenImage = EigenImageType::New();
	getFOImageHandler(EigenImage, EigenFileName); 
	EigenConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
	//FOImageType::IndexType foIndex;
	EigenImageType::PixelType EigenValue;
	Eigen_it.GoToBegin();
	typedef itk::Vector<double, 3> VectorType;
	VectorType eigen;*/





