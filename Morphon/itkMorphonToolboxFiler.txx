/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonToolboxFiler.txx,v $
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

#ifndef _itkMorphonToolboxFiler_txx
#define _itkMorphonToolboxFiler_txx

#include "itkMorphonToolboxFiler.h"

#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkMorphonToolbooxFunction.h"

#include "itkExceptionObject.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkDataObject.h"

#include "itkGaussianOperator.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkNeighborhoodOperatorImageFilter.h"

#include "vnl/vnl_math.h"


#include "itkDiscreteGaussianImageFilter.h"

namespace itk {
  /**
   * void SetFixedImage( const FixedImageType * ptr )
   * Set the fixed image
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage, TMovingImage, TDeformationField>
  ::SetFixedImage( const FixedImageType * ptr )
    {
        this->ProcessObject::SetNthInput( 1, const_cast< FixedImageType * >( ptr ) );
    }
    
  /**
   * FixedImageType * GetFixedImage(void) const
   * Return the fixed image
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  const typename MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::FixedImageType *
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GetFixedImage(void) const
  {
      return dynamic_cast< const FixedImageType * > ( this->ProcessObject::GetInput( 1 ) );
  }
  
  /**
   * void SetMovingImage(const MovingImageType * ptr )
   * Set the Moving image
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetMovingImage( const MovingImageType * ptr )
  {
      this->ProcessObject::SetNthInput( 2, const_cast< MovingImageType * >( ptr ) );
  }
  
  /**
   * void SetCertaintyImage(const MovingImageType * ptr )
   * Set the Certainty matrix
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetCertaintyImage( MovingImageType * ptr )
  {
      m_CertaintyImage = ptr;
  }
  
  /**
   * void SetRescaledPrototype(const MovingImageType * ptr )
   * Set the Rescaked Prototype
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetRescaledPrototype( MovingImageType * ptr )
  {
      m_RescaledPrototype = ptr;
  }
  
  /**
   * MovingImageType * GetMovingImage(void) const
   * Return the Moving image
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  const typename MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::MovingImageType *
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GetMovingImage(void) const
  {
      return dynamic_cast< const MovingImageType * >     ( this->ProcessObject::GetInput( 2 ) );
  }
  
  /**
   * MovingImageType * GetCertaintyImage(void) const
   * Return the Certainty Matrix
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  typename MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::MovingImageType *
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GetCertaintyImage(void) const
  {
      return m_CertaintyImage;
  }
  
  /**
   * MovingImageType * GetRescaledPrototype(void) const
   * Return the Rescaled Prototype
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  typename MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::MovingImageType*
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GetRescaledPrototype(void) const
  {
      return m_RescaledPrototype ;
  }
  
  /**
   * void SetInitialDeformationField( DeformationFieldType * ptr )
   * Set the initial Deformation Field
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetInitialDeformationField( DeformationFieldType * ptr )
  { this->SetInput( ptr );}

  /**
   * DeformationFieldType * GetDeformationField()
   * Return the Deformation Field
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  typename MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::DeformationFieldType *
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GetDeformationField()
  { return this->GetOutput(); }
  
  /**
   * void SetStandardDeviations( double value )
   * Set the standard deviations
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetStandardDeviations( double value )
  {
      unsigned int j;
      for( j = 0; j < ImageDimension; j++ )
        {
        if( value != m_StandardDeviations[j] )
          {
          break;
          }
        }
      if( j < ImageDimension )
        {
        this->Modified();
        for( j = 0; j < ImageDimension; j++ )
          {
          m_StandardDeviations[j] = value;
          }
        }
  }
  
  /**
   * void SetUpdateFieldStandardDeviations( double value )
   * Use to update the deformation field's standard deviation
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SetUpdateFieldStandardDeviations( double value )
  {
      unsigned int j;
      for( j = 0; j < ImageDimension; j++ )
        {
        if( value != m_UpdateFieldStandardDeviations[j] )
          {
          break;
          }
        }
      if( j < ImageDimension )
        {
        this->Modified();
        for( j = 0; j < ImageDimension; j++ )
          {
          m_UpdateFieldStandardDeviations[j] = value;
          }
        }
  }
  
  /**
   * MorphonToolboxFiler()
   * Constructor
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::MorphonToolboxFiler()
  {
      this->SetNumberOfRequiredInputs(2);
      this->SetNumberOfIterations(10);

      unsigned int j;
      for( j = 0; j < ImageDimension; j++ )
        {
        m_StandardDeviations[j] = 1.0;
        m_UpdateFieldStandardDeviations[j] = 1.0;
        }

      m_TempField = DeformationFieldType::New();
      m_MaximumError = 0.1;
      m_MaximumKernelWidth = 90;
      m_StopRegistrationFlag = false;

      m_SmoothDeformationField = true;
      m_SmoothUpdateField = false;

  }
  
  
  /**
   * void PrintSelf(std::ostream& os, Indent indent) const
   * Use to print local labrary's information
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
      Superclass::PrintSelf(os, indent);
      os << indent << "Smooth deformation field: "
         << (m_SmoothDeformationField ? "on" : "off") << std::endl;
      os << indent << "Standard deviations: [";
      unsigned int j;
      for( j = 0; j < ImageDimension - 1; j++ )
        {
        os << m_StandardDeviations[j] << ", ";
        }
      os << m_StandardDeviations[j] << "]" << std::endl;
      os << indent << "Smooth update field: "
         << (m_SmoothUpdateField ? "on" : "off") << std::endl;
      os << indent << "Update field standard deviations: [";
      for( j = 0; j < ImageDimension - 1; j++ )
        {
        os << m_UpdateFieldStandardDeviations[j] << ", ";
        }
      os << m_UpdateFieldStandardDeviations[j] << "]" << std::endl;
      os << indent << "StopRegistrationFlag: ";
      os << m_StopRegistrationFlag << std::endl;
      os << indent << "MaximumError: ";
      os << m_MaximumError << std::endl;
      os << indent << "MaximumKernelWidth: ";
      os << m_MaximumKernelWidth << std::endl;
  }
  
  /**
   * bool Halt()
   * Supplies the halting criteria for this class of filters.  The
   * algorithm will stop after a user-specified number of iterations.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  bool
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::Halt()
    {

    if ( m_StopRegistrationFlag )
      {
      return true;
      }

    return this->Superclass::Halt();
    }
  
  /**
   * void CopyInputToOutput()
   * A simple method to copy the data from the input to the output.
   * If the input does not exist, a zero field is written to the output.
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::CopyInputToOutput()
  {
      typename Superclass::InputImageType::ConstPointer  inputPtr  = this->GetInput();

  if( inputPtr )
    {
    this->Superclass::CopyInputToOutput();
    }
  else
    {
    typename Superclass::PixelType zeros;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      zeros[j] = 0;
      }

    typename OutputImageType::Pointer output = this->GetOutput();

    ImageRegionIterator<OutputImageType> out(output, output->GetRequestedRegion());

    while( ! out.IsAtEnd() )
      {
      out.Value() =  zeros;
      ++out;
      }
    }
  }
  
  /** 
   * void InitializeIteration()
   * Initialize the state of filter and equation before each iteration.
   * Progress feeback is implemented as part of this method.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::InitializeIteration()
  {
      MovingImageConstPointer movingPtr = this->GetMovingImage();
      FixedImageConstPointer fixedPtr = this->GetFixedImage();
      MovingImagePointer certaintyPtr = this->GetCertaintyImage();
      if( !movingPtr || !fixedPtr || !certaintyPtr )
        {
        itkExceptionMacro( << "Fixed and/or moving and/or certainly image not set" );
        }

      // update variables in the equation object
      PDEDeformableRegistrationFunctionType *f =
        dynamic_cast<PDEDeformableRegistrationFunctionType *>
        (this->GetDifferenceFunction().GetPointer());
      if ( !f )
        {
        itkExceptionMacro(<<"FiniteDifferenceFunction not of type PDEDeformableRegistrationFilterFunction");
        }

      f->SetFixedImage( fixedPtr );
      f->SetMovingImage( movingPtr );
      f->SetCertaintyImage( certaintyPtr );
      f->SetRescaledPrototype( this->GetRescaledPrototype() );
      this->Superclass::InitializeIteration();
        /** Now, set the deformed moving image as the new moving and  certainty image.
            *   This is necessary to take previous iterations into consideration.
            */
      this->SetMovingImage( f->GetMovingImage() );
      this->SetCertaintyImage( f->GetCertaintyImage() );
  }
  
  /** 
   * void SmoothDeformationField()
   * Utility to smooth the deformation field (represented in the Output)
   * using a Guassian operator. The amount of smoothing can be specified
   * by setting the StandardDeviations.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SmoothDeformationField()
  {

    DeformationFieldPointer field = this->GetOutput();
    typename DeformationFieldType::SpacingType spacing;
    spacing = this->GetDeformationField()->GetSpacing();
    typedef itk::RecursiveGaussianImageFilter < DeformationFieldType, DeformationFieldType > RecursiveGaussianType;
    typename RecursiveGaussianType::Pointer SmoothX = RecursiveGaussianType::New();
      SmoothX->SetSigma( m_StandardDeviations[0]*spacing[0] );
      SmoothX->SetDirection( 0 );
      SmoothX->SetNormalizeAcrossScale( false );
      SmoothX->SetInput( this->GetDeformationField() );
    typename RecursiveGaussianType::Pointer SmoothY = RecursiveGaussianType::New();
      SmoothY->SetSigma( m_StandardDeviations[1]*spacing[1] );
      SmoothY->SetDirection( 1 );
      SmoothY->SetNormalizeAcrossScale( false );
      SmoothY->SetInput( SmoothX->GetOutput() );
   typename RecursiveGaussianType::Pointer SmoothZ;
    if( ImageDimension == 3 ){
      SmoothZ = RecursiveGaussianType::New();
      SmoothZ->SetSigma( m_StandardDeviations[2]*spacing[2] );
      SmoothZ->SetDirection( 2 );
      SmoothZ->SetNormalizeAcrossScale( false );
      SmoothZ->SetInput( SmoothY->GetOutput() );
      try{
        SmoothZ->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr<< "Exception catched in smoother Update in SmoothDeformationField() function"<<std::endl;
        std::cerr<< excep<<std::endl;
      }
      field->SetPixelContainer( SmoothZ->GetOutput()->GetPixelContainer() );
      field->SetRequestedRegion( SmoothZ->GetOutput()->GetRequestedRegion() );
      field->SetBufferedRegion( SmoothZ->GetOutput()->GetBufferedRegion() );
      field->SetLargestPossibleRegion( SmoothZ->GetOutput()->GetLargestPossibleRegion() );
      field->CopyInformation( SmoothZ->GetOutput() );
    }
    else{
      SmoothY->Update();     
     field->SetPixelContainer( SmoothY->GetOutput()->GetPixelContainer() );
     field->SetRequestedRegion( SmoothY->GetOutput()->GetRequestedRegion() );
     field->SetBufferedRegion( SmoothY->GetOutput()->GetBufferedRegion() );
     field->SetLargestPossibleRegion( SmoothY->GetOutput()->GetLargestPossibleRegion() );
     field->CopyInformation( SmoothY->GetOutput() );
     
    }
  }
  
  /**
   *    void SmoothCertaintyImage()
   *     Utility to smooth the certainty image
   *    using a Guassian operator. The amount of smoothing can be specified
   *    by setting the StandardDeviations.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::SmoothCertaintyImage()
  {
        /** Smooth the certainty image */
    typename MovingImageType::SpacingType spacing;
    spacing = this->GetCertaintyImage()->GetSpacing();
    typedef itk::RecursiveGaussianImageFilter < TMovingImage, TMovingImage > RecursiveGaussianType;
    typename RecursiveGaussianType::Pointer SmoothX = RecursiveGaussianType::New();
      SmoothX->SetSigma( m_StandardDeviations[0]*spacing[0] );
      SmoothX->SetDirection( 0 );
      SmoothX->SetNormalizeAcrossScale( false );
      SmoothX->SetInput( this->GetCertaintyImage() );
    typename RecursiveGaussianType::Pointer SmoothY = RecursiveGaussianType::New();
      SmoothY->SetSigma( m_StandardDeviations[1]*spacing[1] );
      SmoothY->SetDirection( 1 );
      SmoothY->SetNormalizeAcrossScale( false );
      SmoothY->SetInput( SmoothX->GetOutput() );
   typename RecursiveGaussianType::Pointer SmoothZ;
    if( ImageDimension == 3 ){
      SmoothZ = RecursiveGaussianType::New();
      SmoothZ->SetSigma( m_StandardDeviations[2]*spacing[2] );
      SmoothZ->SetDirection( 2 );
      SmoothZ->SetNormalizeAcrossScale( false );
      SmoothZ->SetInput( SmoothY->GetOutput() );
      try{
        SmoothZ->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr<< "Exception catched in smoother Update in SmoothCertaintyImage() function"<<std::endl;
        std::cerr<< excep<<std::endl;
      }
      m_CertaintyImage->SetPixelContainer( SmoothZ->GetOutput()->GetPixelContainer() );
      m_CertaintyImage->SetRequestedRegion( SmoothZ->GetOutput()->GetRequestedRegion() );
      m_CertaintyImage->SetBufferedRegion( SmoothZ->GetOutput()->GetBufferedRegion() );
      m_CertaintyImage->SetLargestPossibleRegion( SmoothZ->GetOutput()->GetLargestPossibleRegion() );
      m_CertaintyImage->CopyInformation( SmoothZ->GetOutput() );
    }
    else{
      SmoothY->Update();     
     m_CertaintyImage->SetPixelContainer( SmoothY->GetOutput()->GetPixelContainer() );
     m_CertaintyImage->SetRequestedRegion( SmoothY->GetOutput()->GetRequestedRegion() );
     m_CertaintyImage->SetBufferedRegion( SmoothY->GetOutput()->GetBufferedRegion() );
     m_CertaintyImage->SetLargestPossibleRegion( SmoothY->GetOutput()->GetLargestPossibleRegion() );
     m_CertaintyImage->CopyInformation( SmoothY->GetOutput() );
     
    }
        
    
  }
  
  /**
   *   void LocalSmooth( itk::Image< float, ImageDimension > * ptr )
   *   This function is mostly used to smooth an unknown prior image
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::LocalSmooth( itk::Image< float, ImageDimension > * ptr )
  {
   
    typename MovingImageType::SpacingType spacing;
    spacing = ptr->GetSpacing();
    typedef itk::RecursiveGaussianImageFilter < TMovingImage, TMovingImage > RecursiveGaussianType;
    typename RecursiveGaussianType::Pointer SmoothX = RecursiveGaussianType::New();
      SmoothX->SetSigma( m_StandardDeviations[0]*spacing[0] );
      SmoothX->SetDirection( 0 );
      SmoothX->SetNormalizeAcrossScale( false );
      SmoothX->SetInput( ptr );
    typename RecursiveGaussianType::Pointer SmoothY = RecursiveGaussianType::New();
      SmoothY->SetSigma( m_StandardDeviations[1]*spacing[1] );
      SmoothY->SetDirection( 1 );
      SmoothY->SetNormalizeAcrossScale( false );
      SmoothY->SetInput( SmoothX->GetOutput() );
   typename RecursiveGaussianType::Pointer SmoothZ;
    if( ImageDimension == 3 ){
      SmoothZ = RecursiveGaussianType::New();
      SmoothZ->SetSigma( m_StandardDeviations[2]*spacing[2] );
      SmoothZ->SetDirection( 2 );
      SmoothZ->SetNormalizeAcrossScale( false );
      SmoothZ->SetInput( SmoothY->GetOutput() );
      try{
        SmoothZ->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr<< "Exception catched in smoother Update in LocalSmooth() function"<<std::endl;
        std::cerr<< excep<<std::endl;
      }
      ptr->SetPixelContainer( SmoothZ->GetOutput()->GetPixelContainer() );
      ptr->SetRequestedRegion( SmoothZ->GetOutput()->GetRequestedRegion() );
      ptr->SetBufferedRegion( SmoothZ->GetOutput()->GetBufferedRegion() );
      ptr->SetLargestPossibleRegion( SmoothZ->GetOutput()->GetLargestPossibleRegion() );
      ptr->CopyInformation( SmoothZ->GetOutput() );
    }
    else{
      SmoothY->Update();     
     ptr->SetPixelContainer( SmoothY->GetOutput()->GetPixelContainer() );
     ptr->SetRequestedRegion( SmoothY->GetOutput()->GetRequestedRegion() );
     ptr->SetBufferedRegion( SmoothY->GetOutput()->GetBufferedRegion() );
     ptr->SetLargestPossibleRegion( SmoothY->GetOutput()->GetLargestPossibleRegion() );
     ptr->CopyInformation( SmoothY->GetOutput() );
     
    }
  }
  
  /**
   * virtual void GenerateOutputInformation();
   * By default the output deformation field has the same Spacing, Origin
   * and LargestPossibleRegion as the input/initial deformation field.  If
   * the initial deformation field is not set, the output information is
   * copied from the fixed image.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GenerateOutputInformation()
  {
       typename DataObject::Pointer output;

      if( this->GetInput(0) )
        {
        // Initial deformation field is set.
        // Copy information from initial field.
        this->Superclass::GenerateOutputInformation();

        }
      else if( this->GetFixedImage() )
        {
        // Initial deforamtion field is not set.
        // Copy information from the fixed image.
        for (unsigned int idx = 0; idx <
               this->GetNumberOfOutputs(); ++idx )
          {
          output = this->GetOutput(idx);
          if (output)
            {
            output->CopyInformation(this->GetFixedImage());
            }
          }

        }
  }
  
  /**
   * virtual void GenerateInputRequestedRegion()
   * It is difficult to compute in advance the input moving image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole moving image.
   *
   * For the fixed image and deformation field, the input requested region
   * set to be the same as that of the output requested region.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void 
  MorphonToolboxFiler<TFixedImage,TMovingImage,TDeformationField>
  ::GenerateInputRequestedRegion()
  {
      // call the superclass's implementation
      Superclass::GenerateInputRequestedRegion();

      // request the largest possible region for the moving image
      MovingImagePointer movingPtr = const_cast< MovingImageType * >( this->GetMovingImage() );
      if( movingPtr )
        {
        movingPtr->SetRequestedRegionToLargestPossibleRegion();
        }

      // just propagate up the output requested region for
      // the fixed image and initial deformation field.
      DeformationFieldPointer inputPtr =
        const_cast< DeformationFieldType * >( this->GetInput() );
      DeformationFieldPointer outputPtr = this->GetOutput();
      FixedImagePointer fixedPtr =
        const_cast< FixedImageType *>( this->GetFixedImage() );

      if( inputPtr )
        {
        inputPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
        }

      if( fixedPtr )
        {
        fixedPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
        }
  }
  
 }//end namespace
 
 #endif
 
 
