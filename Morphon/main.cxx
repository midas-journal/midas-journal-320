
/*
This file creates a morphon registration pipeline.
*/

#include "itkMorphonPipe.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"

#include "complex"

const unsigned int Dimension = 2;

int main( int argc, char ** argv )
{
    
    std::cout<<" Welcom in the Morphon pipeline \n"<<std::endl;

    if( argc < 5 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Minimum usage: " << argv[0]<<"\n";
        std::cerr << " fixedImageFile movingImageFile outputImageFileName [PARAM]"<<"\n";
        std::cerr << "[PARAM] = [level] ; [level] | [empty]" << "\n";
        std::cerr << "[level] = [iteration number] [variance]" <<"\n";
        std::cerr << "[empty] = NULL"<<"\n";
        return 1;
    }    /* Constants */

    int nlevel; /* the level number */
    if(argc>4)
    {
    	nlevel=(int)floor((argc-4)/2);
    }// nbr levels depends on number of pairs of parameters (nbr it - variances) put after the image inputs arguments.
    else
    {
        nlevel=1;
    }

    unsigned int noit[nlevel]; /*number of iterations on each level */
    float nofdev[nlevel]; /* value of deviation on each level*/
    if(argc>4)
    {
    	for(int i=0;i<nlevel;i++)
    	{
    	  noit[i]=atoi(argv[2*i+4]);
          nofdev[i]=atof(argv[2*i+5]);
          std::cout << "Level " << nlevel-i << " with " << noit[i] << " iterations and " << nofdev[i] << " as variance" << std::endl;
        }
    }
    else
    {
        noit[0]=1;
        nofdev[0]=1.0;
    }

    /*The files' names*/
    std::string imageFixedName = argv[1];
    std::string imageMovingName = argv[2];
    std::string imageOutputName = argv[3];

    size_t ext_pos = imageOutputName.find_last_of( '.' );
    std::string ext( imageOutputName, ext_pos );
    
    typedef unsigned char OutputPixelType;
 
    
    typedef double                                           PixelType;
    typedef float                                                     InternalPixelType;
    typedef itk::Image < InternalPixelType, Dimension >    InternalImageType;
        InternalImageType::Pointer CertaintyImage;
        InternalImageType::Pointer InternalMovingImage;
    typedef itk::Vector< float, Dimension >                     VectorPixelType;
    typedef itk::Image < VectorPixelType, Dimension >     DeformationFieldType;
    typedef itk::Image< OutputPixelType, Dimension >      OutputImageType;    
    typedef itk::Image < PixelType, Dimension >               FixedImageType;
    typedef itk::Image < PixelType, Dimension >               MovingImageType;
    typedef itk::ImageFileReader < FixedImageType >        FixedImageReaderType;
        FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
        fixedImageReader->SetFileName( imageFixedName );
        fixedImageReader->Update();
        FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    typedef itk::ImageFileReader < MovingImageType >      MovingImageReaderType;
        MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
        movingImageReader->SetFileName( imageMovingName );
        movingImageReader->Update();
        MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
    typedef itk::CastImageFilter < MovingImageType, InternalImageType >         MovingCasterType;
        MovingCasterType::Pointer movingImageCaster = MovingCasterType::New();
        movingImageCaster->SetInput( movingImageReader->GetOutput() );
        movingImageCaster->Update();
    typedef itk::CastImageFilter < FixedImageType, InternalImageType >          FixedCasterType;
        FixedCasterType::Pointer fixedImageCaster = FixedCasterType::New();
        fixedImageCaster->SetInput( fixedImageReader->GetOutput() ) ;
        fixedImageCaster->Update();

    /* The morphons pipeline */
    typedef itk::MorphonPipe< InternalImageType, InternalImageType,OutputImageType> MorphonType;
    

    MorphonType::Pointer morphon = MorphonType::New();
        morphon->SetMovingImage( movingImageCaster->GetOutput() );
        InternalImageType::IndexType index;
        morphon->SetFixedImage( fixedImageCaster->GetOutput() );
        morphon->SetLevelNumber ( nlevel );
        morphon->SetPadding( true );
        morphon->SetIterationsNumber( noit );
        morphon->SetStandardDeviation ( nofdev );
        morphon->SetOutputName( imageOutputName );
        try
        {
            morphon->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception catched in the main in the morphons Update()" << std::endl;
            std::cerr << excep << std::endl;
        }

        OutputImageType::Pointer out_morphon = morphon->GetOutput();
        
}//end main
