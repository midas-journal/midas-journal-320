/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonPipe.txx,v $
  Language:  C++
  Date:      $Date: 13th February 2009$
  Version:   $Revision: $

  Copyright (c) 2009 Plumat Jerome, Benoit Macq
  Laboratoire de Télécommunications et Télédétection 
  Ecole Polytechnique de Louvain 
  Université Catholique de Louvain

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkMorphonPipe_txx
#define __itkMorphonPipe_txx

#include "itkMorphonPipe.h"

#include "itkCommand.h"
#include "itkImageToImageFilter.h"

#include "itkMorphonRegistrationFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkVectorResampleImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "itkImageFileWriter.h"

namespace itk
{
  template<class TFixedImage, class TMovingImage, class TOutputImage>
  void
  MorphonPipe<TFixedImage, TMovingImage, TOutputImage>
  ::GenerateData()
  {
    std::cout<<" Start the pyramide \n"<<std::endl;
        
    int resol = 0, iterLevel;
        
    fixedXImageSmooth = DiscretGaussianType::New();
    fixedYImageSmooth = DiscretGaussianType::New();
    if( ImageDimension==3 ) { fixedZImageSmooth = DiscretGaussianType::New(); }
        
    movingXImageSmooth = DiscretGaussianType::New();
    movingYImageSmooth = DiscretGaussianType::New();
    if( ImageDimension==3 ) { movingZImageSmooth = DiscretGaussianType::New(); }
        
    certaintyXImageSmooth = DiscretGaussianType::New();
    certaintyYImageSmooth = DiscretGaussianType::New();	
    if( ImageDimension==3 ) { certaintyZImageSmooth = DiscretGaussianType::New(); }
        
    /** The deformations field */
    deformationField = DeformationFieldType::New();
        
    /** PYRAMIDE */
    for(iterLevel=m_nlevel; iterLevel>0; iterLevel--)
    {
      std::cout << "\nLevel in the Pyramid: " << iterLevel << std::endl;

      /** Compute the reduce factor in this level */
      /** used to compute the dimensions over the pyramide */
      unsigned int factor[ImageDimension];
      for(unsigned j=0;j<ImageDimension;j++)
      {
        factor[j] = static_cast<int>(pow(2,static_cast<float>(iterLevel-1)/resampling));
      }

      /** Smoothing */
      float variance[ImageDimension]; /* variance for the smoothing*/
      std::cout<<"\t Setting variances for smoothing ";
      for (unsigned int i=0; i<ImageDimension; i++)
      {
        variance[i] = vnl_math_sqr(0.5*factor[i] ) ;
        std::cout<<variance[i]<<" ";
      }
      std::cout<<std::endl;
	     
     /** Smoothing */
     fixedXImageSmooth->SetSigma( variance[0] );
     fixedXImageSmooth->SetDirection( 0 );
     fixedXImageSmooth->SetNormalizeAcrossScale( false );
     fixedYImageSmooth->SetSigma( variance[1] );
     fixedYImageSmooth->SetDirection( 1 );
     fixedYImageSmooth->SetNormalizeAcrossScale( false );
     if( ImageDimension==3 ) 
     {
       fixedZImageSmooth->SetSigma( variance[2] );
       fixedZImageSmooth->SetDirection( 2 );
       fixedZImageSmooth->SetNormalizeAcrossScale( false );
     }
          
     fixedXImageSmooth->SetInput( this->GetInput(1) );
     fixedYImageSmooth->SetInput( fixedXImageSmooth->GetOutput() );
     if( ImageDimension==3 )
     {
       fixedZImageSmooth->SetInput( fixedYImageSmooth->GetOutput() );
     }

     movingXImageSmooth->SetSigma( variance[0] );
     movingXImageSmooth->SetDirection( 0 );
     movingXImageSmooth->SetNormalizeAcrossScale( false );
     movingYImageSmooth->SetSigma( variance[1] );
     movingYImageSmooth->SetDirection( 1 );
     movingYImageSmooth->SetNormalizeAcrossScale( false );
     if( ImageDimension==3 )
     {
       movingZImageSmooth->SetSigma( variance[2] );
       movingZImageSmooth->SetDirection( 2 );
       movingZImageSmooth->SetNormalizeAcrossScale( false );
     }
          
    movingXImageSmooth->SetInput( this->GetInput(0) );
    movingYImageSmooth->SetInput( movingXImageSmooth->GetOutput() );
    if( ImageDimension==3 )
    {
      movingZImageSmooth->SetInput( movingYImageSmooth->GetOutput() ); 
    }


     /** Get Spacing */
     typename FixedImageType::SpacingType spacing = this->GetInput(1)->GetSpacing();
     std::cout<< "\t Final fixed image pacing in all dimensions are "<<spacing<<std::endl;
     SizeType size = this->GetInput(1)->GetLargestPossibleRegion().GetSize();


     /** Shrinking */
     idTransform = TransformType::New();
     idTransform->SetIdentity();
     movingResampler = ResampleFilterType::New();
     fixedResampler = ResampleFilterType::New();
     interpolator = ImageInterpolatorType::New();
     typename FixedImageType::SpacingType spacing_shrink;
     SizeType size_shrink;

     for (unsigned int i=0; i<ImageDimension; i++)
     {
       spacing_shrink[i] = spacing[i] * pow(2,static_cast<float>(iterLevel-1)/resampling);
       size_shrink[i] = (long unsigned int) floor(size[i] / pow(2,static_cast<float>(iterLevel-1)/resampling));
     }

     std::cout << "\t Spacing after resampling" << spacing_shrink << "  -  Size " << size_shrink << std::endl;

     fixedResampler->SetOutputSpacing(spacing_shrink);
     fixedResampler->SetSize(size_shrink);
     if( ImageDimension == 3 )
     { 
       fixedResampler->SetOutputOrigin( fixedZImageSmooth->GetOutput()->GetOrigin() );
       fixedResampler->SetInput(fixedZImageSmooth->GetOutput());
     }
     else
     {
       fixedResampler->SetOutputOrigin( fixedYImageSmooth->GetOutput()->GetOrigin() );
       fixedResampler->SetInput(fixedYImageSmooth->GetOutput());
     }
     fixedResampler->SetTransform(idTransform);
     try
     {
       fixedResampler->Update();
     }
     catch( itk::ExceptionObject & excep )
     {
       std::cerr << "Exception catched in the fixed image resampler!" << std::endl;
       std::cerr << excep << std::endl;
     }

     movingResampler->SetOutputSpacing(spacing_shrink);
     movingResampler->SetSize(size_shrink);
     if( ImageDimension == 3 )
     { 
       movingResampler->SetOutputOrigin( movingZImageSmooth->GetOutput()->GetOrigin() );
       movingResampler->SetInput(movingZImageSmooth->GetOutput());
     }
     else
     {
       movingResampler->SetOutputOrigin( movingYImageSmooth->GetOutput()->GetOrigin() );
       movingResampler->SetInput(movingYImageSmooth->GetOutput());
     }
     movingResampler->SetInterpolator( interpolator );            
     movingResampler->SetTransform(idTransform);
     try
     {
       movingResampler->Update();
     }
     catch( itk::ExceptionObject & excep )
     {
       std::cerr << "Exception catched in the moving image resampler!" << std::endl;
       std::cerr << excep << std::endl;
     }

     std::cout << "Resampling ok" << std::endl;

     /** 
      * Inintiate the certainty image with zeros in the first level
      */
     typename MovingImageType::SizeType Imsize = movingResampler->GetOutput()->GetRequestedRegion().GetSize();
     int numPt = 1;
     for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
     typename MovingImageType:: IndexType zerosIndex;
     for( unsigned int i=0; i<ImageDimension; i++) { zerosIndex[i]=0; }
           
     /** Certainly Image*/
     if( iterLevel == this->m_nlevel )
     {
       CertaintyImage = MovingImageType::New();
       CertaintyImage->SetRegions( movingResampler->GetOutput()->GetRequestedRegion() );
       CertaintyImage->SetSpacing( movingResampler->GetOutput()->GetSpacing());
       CertaintyImage->SetOrigin( movingResampler->GetOutput()->GetOrigin());
       CertaintyImage->Allocate();

       /** Initialize the certainly image */
       typename MovingImageType::PixelType *cert = &( CertaintyImage->GetPixel( zerosIndex ) );
       for( unsigned int i=0; i<numPt; i++)
	{
	  *cert = 0.0;
	   cert += 1;
	}//end for
      }//end if

      /** Morphon Pipe */
      observer = CommandIterationUpdate::New();
      filter = RegistrationFilterType::New();
        filter->SetFilters(m_PhaseFilterOperator);/* pass the pointer to have an access to the filters */
        filter->AddObserver( itk::IterationEvent(), observer );
	filter->SetFixedImage( fixedResampler->GetOutput() );
        filter->SetMovingImage( movingResampler->GetOutput() );
        filter->SetPadding( this->GetPadding() );
        if( iterLevel == m_nlevel )
        {
          /* for the first level, initiate the certainty image*/
          filter->SetCertaintyImage( CertaintyImage );
        }
        filter->SetNumberOfIterations( m_noit[resol] );
        filter->SetStandardDeviations( m_nofdev[resol] );

      /** Make a fixed moving image copy, used to deforme */
      typename MovingImageType ::Pointer rescaledPrototype = MovingImageType::New();
        rescaledPrototype->SetRegions( movingResampler->GetOutput()->GetRequestedRegion() );
        rescaledPrototype->SetSpacing( movingResampler->GetOutput()->GetSpacing());
        rescaledPrototype->SetOrigin( movingResampler->GetOutput()->GetOrigin());
        rescaledPrototype->Allocate();

     /** Initialize the rescaled prototype as the moving image type. */
     typename MovingImageType::PixelType *rescaledProtPixel = &( rescaledPrototype->GetPixel( zerosIndex ) );
     typename MovingImageType::PixelType *movingPixel = &( movingResampler->GetOutput()->GetPixel( zerosIndex ) );
     for( unsigned int i=0; i<numPt; i++)
     {
       *rescaledProtPixel = *movingPixel;
	rescaledProtPixel += 1;
	movingPixel += 1;
     }//end for
     filter->SetRescaledPrototype( rescaledPrototype );

     /** The deformations and the certainties are computed , take it in consideration */
     if( iterLevel != m_nlevel )
     {
       typename DeformationFieldType::SpacingType spacing_defor = deformationField->GetSpacing();
       std::cout<< " Deformation field's spacing befor reinterpolation : " << spacing_defor << std::endl;

       SizeType size_resol = fixedResampler->GetOutput()->GetLargestPossibleRegion().GetSize();

       /** Take the compute deformation field and transform it : smoothing and resize*/
       /** Resize */
       typename ResampleFilterDeformationType::Pointer resampler = ResampleFilterDeformationType::New();

       typename TransformType::Pointer transform = TransformType::New();
         transform->SetIdentity();

       typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

       for (unsigned int i=0; i<ImageDimension; i++)
       {
         spacing_defor[i] = spacing_defor[i] / pow(2,1.0/resampling);
       }
       std::cout << "Spacing deformation reinterpolator in all directions are: " << spacing_defor << std::endl;

       resampler->SetOutputSpacing(spacing_defor);
       resampler->SetOutputOrigin( fixedResampler->GetOutput()->GetOrigin() );
       resampler->SetInterpolator( interpolator );
       resampler->SetSize(size_resol);
       resampler->SetInput(deformationField);
       resampler->SetTransform(transform);
       try
       {
         resampler->Update();
       }
       catch( itk::ExceptionObject & excep )
       {
         std::cerr << "Exception catched in the deformation field resampler!" << std::endl;
         std::cerr << excep << std::endl;
       }

       /** Set the new initial deformation field */
       deformationField = resampler->GetOutput();
       filter->SetInitialDeformationField(resampler->GetOutput());

       /** Take the certainly image and transform it : smoothing and resize*/
       /** Smoothing */
       certaintyXImageSmooth->SetSigma( variance[0] );
       certaintyXImageSmooth->SetDirection( 0 );
       certaintyXImageSmooth->SetNormalizeAcrossScale( false );
       certaintyYImageSmooth->SetSigma( variance[1] );
       certaintyYImageSmooth->SetDirection( 1 );
       certaintyYImageSmooth->SetNormalizeAcrossScale( false );
       if( ImageDimension == 3 )
       { 
         certaintyZImageSmooth->SetSigma( variance[2] );
         certaintyYImageSmooth->SetDirection( 2 );
         certaintyYImageSmooth->SetNormalizeAcrossScale( false );  
       }
       certaintyXImageSmooth->SetInput( CertaintyImage );
       certaintyYImageSmooth->SetInput( certaintyXImageSmooth->GetOutput() );
       if( ImageDimension == 3 )
       {
         certaintyZImageSmooth->SetInput( certaintyYImageSmooth->GetOutput() ); 
       }
                
       /** Resize the certainty image */
       certaintyResampler = ResampleFilterType::New();
       certaintyInterpolator = ImageInterpolatorType::New();
       certaintyResampler->SetSize( size_resol );
       certaintyResampler->SetTransform( transform );
       if( ImageDimension == 3 )
       {
         certaintyResampler->SetInput( certaintyZImageSmooth->GetOutput() );
       }
       else 
       {
         certaintyResampler->SetInput( certaintyYImageSmooth->GetOutput() );
       }
       certaintyResampler->SetOutputSpacing( spacing_defor );
       certaintyResampler->SetOutputOrigin( fixedResampler->GetOutput()->GetOrigin() );
       certaintyResampler->SetInterpolator( certaintyInterpolator );
       try
       {
         certaintyResampler->Update();
       }
       catch( itk::ExceptionObject & excep )
       {
         std::cerr << "Exception catched in the certainty resampler!" << std::endl;
         std::cerr << excep << std::endl;
       }
       
       /** Set the new initial certainty image */
       filter->SetCertaintyImage( certaintyResampler->GetOutput() );
       CertaintyImage = certaintyResampler->GetOutput();

       std::cout << std::endl << "********************************************************* " << std::endl;
       std::cout << std::endl;

     }//end if

     /** Compute the morphons' algorithm */
     try
     {
       filter->Update();
     }
     catch( itk::ExceptionObject & excep )
     {
       std::cerr << "Exception catched in the morphons update!" << std::endl;
       std::cerr << excep << std::endl;
     }

     deformationField = filter->GetOutput();
     CertaintyImage = filter->GetCertaintyImage();
     MovingImage = const_cast< MovingImageType *>( filter->GetMovingImage() );

     if (iterLevel == 1)
     {
     /**
      *  WRITING RESULTS
      */

       typedef itk::WarpImageFilter<MovingImageType,MovingImageType,DeformationFieldType  >     WarperType;
       typedef itk::LinearInterpolateImageFunction<MovingImageType,double > InterpolatorType;
           typename WarperType::Pointer warper = WarperType::New();
                      
       /** Warpe the original moving image */
       warper->SetInput( movingResampler->GetOutput() );//OutMovingMirrorPad->GetOutput() );
       warper->SetInterpolator( interpolator );
       warper->SetOutputSpacing( movingResampler->GetOutput()->GetSpacing() );
       warper->SetOutputOrigin( movingResampler->GetOutput()->GetOrigin() );
       warper->SetDeformationField( filter->GetDeformationField() );

       typename WriterType::Pointer      writer =  WriterType::New();
       typename CastFilterType::Pointer  caster =  CastFilterType::New();

       std::cout<< "Saving output - image - File out: " << m_OutputName << std::endl;
       /** Write the output image and the deformation field */
       writer->SetFileName( m_OutputName );
       caster->SetInput( filter->GetMovingImage() );
       writer->SetInput( caster->GetOutput()   );
       try
       {
         writer->Update();
       }
       catch( itk::ExceptionObject & excep )
       {
          std::cerr << "Exception catched in the morphons update!" << std::endl;
          std::cerr << excep << std::endl;
       }

       std::string fieldOutputName = "DeformationField.mhd";
       std::cout<< "Saving output - deformation field - File out: " << fieldOutputName << std::endl;

       typedef itk::ImageFileWriter< DeformationFieldType > FieldWriterType;
       typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
            fieldWriter->SetFileName( fieldOutputName );
            fieldWriter->SetInput( filter->GetOutput() );
            fieldWriter->Update();

       this->GraftOutput( caster->GetOutput()  );

     }//end if
     
     /** level is finish, compute the nex one*/
     resol = resol +1;
   }//end pyramide
 }//end GenerateData

}//end namespace


#endif

