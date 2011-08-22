/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonRegistrationFilter.txx,v $
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

#ifndef _itkMorphonRegistrationFilter_txx
#define _itkMorphonRegistrationFilter_txx

#include "itkMorphonRegistrationFilter.h"
#include "itkMorphonToolboxFiler.h"
#include "itkComputingMorphonDeformationField.h"

namespace itk {

  /**
   * double GetMetric()
   * This function is used to have an access to a metric.
   * Actually, the used metric is the mean square difference in intensity
   * between the fixed image and transforming moving image computed over
   * the the overlapping region between the two images. This is done for the
   * current iteration at a current level.
   * But, in this case, this is not the best metric we can use
   * because of the multi-modality.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  double
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::GetMetric() const
  {
    MorphonsRegistrationFunctionType *mrfp =
      dynamic_cast<MorphonsRegistrationFunctionType *>
     (this->GetDifferenceFunction().GetPointer());
      
    if( !mrfp )
    {
     itkExceptionMacro( <<
       "Could not cast difference function to SymmetricForcesDemonsRegistrationFunction" );
    }
    return mrfp->GetMetric();
  }
    
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "Intensity difference threshold: " << std::endl;
  }
     	
  /**
   * void InitializeIteration()
   * 
   * Initialize the iteration. The first step is to deform the moving image.
   *
   * If it's the first iteration, we have to convolute the moving AND the fixed 
   * image, if it's the second (or higher) than, we just have to convolute the 
   * moving image (because it was deformed with the computed deformation field).
   *
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::InitializeIteration()
  {
    typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;
      // update variables in the equation object
    MovingImageConstPointer movingPtr = this->GetMovingImage();

    MorphonsRegistrationFunctionType *f =
      dynamic_cast<MorphonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
    if( IterationNum==1) 
    {
      tempCert->SetRegions( this->GetCertaintyImage()->GetRequestedRegion() );
      tempCert->SetSpacing( this->GetCertaintyImage()->GetSpacing());
      tempCert->SetOrigin( this->GetCertaintyImage()->GetOrigin());
      try
      {
        tempCert->Allocate();
      }
      catch( itk::ExceptionObject & excep )
      {
         std::cerr << "Exception catched in the memory allocation" << std::endl;
         std::cerr << excep << std::endl;
      }
    }

    if ( !f )
    {
      itkExceptionMacro(<<"FiniteDifferenceFunction not of type MorphonsRegistrationFunctionType");
    }
    
    f->SetFilters(this->GetFilters());
    f->SetDeformationField( this->GetDeformationField() );

    f->SetCertaintyImage( this->GetCertaintyImage() );
    f->SetMovingImage( movingPtr );
    f->setIterationNum( IterationNum );
    /* Set the quadrature coefficients */
    f->SetPadding( this->GetPadding() );
    f->SetQuadCoef( m_FixedConv, m_MovingConv, m_qqin );

    IterationNum = IterationNum +1;

    // call the superclass  implementation
    Superclass::InitializeIteration();
    /* Take back the deformed moving image. */
    this->SetMovingImage( f->GetMovingImage() );
    /* Take back the certainty image. */
    this->SetCertaintyImage( f->GetCertaintyImage() );
    /* take back the quadrature coefficients. */
    f->GetQuadCoef( m_FixedConv, m_MovingConv, m_qqin );
  }
  
  /**
   * void ApplyUpdate(TimeStepType dt)
   * This function is used to call the superClass' applyUpdate (as it's usually
   * done) and to make a normalized convolution over the deformation field and 
   * the certainty image (as explain in Normalized Convolution: {A} Technique 
   * for Filtering Incomplete and Uncertain Data, H. Knutsson and C-F. Westin, 
   * 1993)
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::ApplyUpdate(TimeStepType dt)
  {
    this->Superclass::ApplyUpdate(dt);
    MorphonsRegistrationFunctionType *mrfp =
      dynamic_cast<MorphonsRegistrationFunctionType *>
        (this->GetDifferenceFunction().GetPointer());

    if( !mrfp )
    {
      itkExceptionMacro( <<
         "Could not cast difference function to MorphonsRegistrationFunction" );
    }

    /**
     *  Morphon method doesn't only smooth the certainty image and the deformation field
     */
        
     typename MovingImageType:: IndexType index;
     for( unsigned int i=0; i<ImageDimension; i++) { index[i]=0; }
     
     /** copy the certainty matrix in a temp matrix */
     typename MovingImageType::SizeType Imsize = this->GetCertaintyImage()->GetRequestedRegion().GetSize();
     int numPt = 1;
     for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i];}
     
     typename MovingImageType::PixelType* certPoint = &( this->GetCertaintyImage()->GetPixel( index ) );
     typename MovingImageType::PixelType* tempCerPt = &( tempCert->GetPixel( index ) );
     for(unsigned int i=0; i<numPt; i++)
     {
        *tempCerPt = *certPoint;
        tempCerPt += 1;
        certPoint += 1;   
     }//end iterator*/
     
     /** normalized averaging smoothing */
     //C = liss_gauss(cert)
     this->LocalSmooth( tempCert.GetPointer() );

     // A = def*cert
     DefMult( this->GetDeformationField(), this->GetDeformationField(),this->GetCertaintyImage() );
     
     // B = liss_gauss(A)
     this->SmoothDeformationField();
     //this->LocalSmooth( &(this->GetDeformationField()) );
     
     // D = B/C;
     DefDiv ( this->GetDeformationField(), this->GetDeformationField(),tempCert );
     
     // E = cert * cert
     MatrixMult( this->GetCertaintyImage(), this->GetCertaintyImage(),this->GetCertaintyImage() );
     
     // F = liss_gauss(E)
     this->SmoothCertaintyImage();
     
     // G = F/C;
     MatrixDiv( this->GetCertaintyImage(), this->GetCertaintyImage(),tempCert );
     
  }//end ApplyUpdate
  
  
  /**
   * void DefMult(
   *    itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > result,
   *    const itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > a,
   *    const itk::SmartPointer < itk::Image< float, ImageDimension > > b
   *             )
   *
   * Used to compute the division between a deformation field and an certainty image.
   * result->GetPixel(index) = a->GetPixel(index)*b->GetPixel(b)
   * All matrix must have the same dimension.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  inline void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::DefMult
        (
            DeformationFieldPointer result,
            const DeformationFieldPointer a,
            const MovingImagePointer b
        )
  {
    typedef typename MovingImageType::PixelType	          PixelType;
    typedef itk::ConstNeighborhoodIterator< MovingImageType > NeighborhoodIteratorType;
    typedef itk::ImageRegionIterator< DeformationFieldType >       IteratorType;
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
       
    IteratorType iterator(a, a->GetRequestedRegion());
    typename MovingImageType::IndexType index;
    for( unsigned int i=0; i<ImageDimension; i++) { index[i] = 0; }

    /* Multiplication */
    typename MovingImageType::SizeType Imsize = b->GetRequestedRegion().GetSize();
    int numPt = 1;
    for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
    
    PixelType* localB = &( b->GetPixel( index ) );
    VectorType* defA = & ( a->GetPixel( index ) ), *resultP = & ( result->GetPixel( index ) );
    for( unsigned int k=0; k<numPt; k++)
    {
      for( unsigned int i=0; i<ImageDimension; i++)
      {
        if( *localB != 0 & *localB != 1 )
        {
          (*resultP)[i] = (*defA)[i] * (*localB);
        }
        if( *localB==0 )
        {
          (*resultP)[i] = 0;
        }
      }//end for
      localB += 1;
      defA += 1 ;
      resultP += 1;
    }//end for
  }//end DefMult
  
  /**
    *       void DefDiv(
    *       itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > result,
    *       itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > a,
    *       itk::SmartPointer < itk::Image< float, DIMENSION > > b,
    *        )
    *
    *       Used to compute the multiplication between the deformation field and the certainty image.
    *                   result->GetPixel(index) = a->GetPixel(index)/b->GetPixel(b)
    *       All matrix must have the same dimension.
    **/
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  inline void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::DefDiv
        (
            DeformationFieldPointer result,
            const DeformationFieldPointer a,
            const MovingImagePointer b
        )
  {
    typedef typename MovingImageType::PixelType	          PixelType;
    typedef itk::ConstNeighborhoodIterator< MovingImageType > NeighborhoodIteratorType;
    typedef itk::ImageRegionIterator< DeformationFieldType >       IteratorType;
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
        
    IteratorType iterator(a, a->GetRequestedRegion());
    typename MovingImageType:: IndexType index;
    for( unsigned int i=0; i<ImageDimension; i++) { index[i] = 0; }
    
    /* Division */
    typename MovingImageType::SizeType Imsize = b->GetRequestedRegion().GetSize();
    int numPt = 1;
    for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
    
    PixelType* localB = &( b->GetPixel( index ) );
    VectorType* defA = & ( a->GetPixel( index ) ), *resultP = & ( result->GetPixel( index ) );
    for( unsigned int k=0; k<numPt; k++)
    {
      for( unsigned int i=0; i<ImageDimension; i++)
      {
        if( *localB != 0 & *localB != 1 )
        {
          (*resultP)[i] = (*defA)[i] / (*localB);
        }
        if( *localB==0 )
        {
         (*resultP)[i] = 0;
        }
      }//end for
      localB += 1;
      defA += 1 ;
      resultP += 1;
    }//end for
  }//end DefDiv
    
  /**
   *       void MatrixMult(
   *        itk::SmartPointer < itk::Image< float, ImageDimension > > result,
   *        const itk::SmartPointer < itk::Image< float, ImageDimension > > a,
   *       const itk::SmartPointer < itk::Image< float, ImageDimension > > b
   *        )
   *
   *       Used to compute the matrix multiplication :
   *                   result->GetPixel(index) = a->GetPixel(index)*b->GetPixel(b)
   *       All matrix must have the same dimension.
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  inline void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::MatrixMult
        (
            MovingImagePointer result,
            const MovingImagePointer a,
            const MovingImagePointer b
        )
  {
    typename MovingImageType::SizeType Imsize = b->GetRequestedRegion().GetSize();
    int numPt = 1;
    for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
    
    typename MovingImageType:: IndexType index;
    for( unsigned int i=0; i<ImageDimension; i++) { index[i] = 0; }
    
    typename MovingImageType::PixelType *APixel, *BPixel, *ResultPixel;
    APixel = &( a->GetPixel( index ) );
    BPixel = &( b->GetPixel( index ) );
    ResultPixel = & ( result->GetPixel( index ) );
    for( unsigned int k=0; k<numPt; k++)
    {
      *ResultPixel = (*APixel) * (*BPixel);
      ResultPixel += 1;
      APixel += 1;
      BPixel += 1;
    }//end for
  }//end MatrixMult
    
  /**
   *       inline void MatrixDiv
   *       (
   *         MovingImagePointer result,
   *         const MovingImagePointer a,
   *        const MovingImagePointer b
   *       );
   *       Used to compute the a matrix division :
   *                   result->GetPixel(index) = a->GetPixel(index)/b->GetPixel(b)
   *       All matrix must have the same dimension.
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  inline void
  MorphonRegistrationFilter<TFixedImage, TMovingImage, TDeformationField>
  ::MatrixDiv
        (
           MovingImagePointer result,
           const MovingImagePointer a,
           const MovingImagePointer b
        )
  {
    typename MovingImageType::SizeType Imsize = b->GetRequestedRegion().GetSize();
    int numPt = 1;
    for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
    typename MovingImageType:: IndexType index;
    for( unsigned int i=0; i<ImageDimension; i++) { index[i] = 0; }
    
    typename MovingImageType::PixelType *APixel, *BPixel, *ResultPixel;
    APixel = &( a->GetPixel( index ) );
    BPixel = &( b->GetPixel( index ) );
    ResultPixel = & ( result->GetPixel( index ) );
    for( unsigned int k=0; k<numPt; k++)
    {
      if( *BPixel == 0 ) { *ResultPixel = 1; }
      if( *BPixel != 1 ) { *ResultPixel = (*APixel) / (*BPixel); }
      ResultPixel += 1;
      APixel += 1;
      BPixel += 1;
    }//end for
  }//end MatrixDiv

} // end namespace itk


#endif
