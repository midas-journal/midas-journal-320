/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputingMorphonDeformationField.txx,v $
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

#ifndef _ComputingMorphonDeformationField_txx
#define _ComputingMorphonDeformationField_txx

#include "itkComputingMorphonDeformationField.h"

#include "itkMorphonToolbooxFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

#include "itkNeighborhoodAlgorithm.h"

#include "itkBuildingMorphonFilters.h"

#include <string.h>

#include "itkWarpImageFilter.h"

#include <itkMirrorPadImageFilter.h>

namespace itk {

  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void*
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::GetGlobalDataPointer() const
  {
    GlobalDataStruct *global = new GlobalDataStruct();
    global->m_SumOfSquaredDifference  = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange      = 0;
    return global;
  }
  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::ReleaseGlobalDataPointer( void *gd ) const
  {
    GlobalDataStruct * globalData = (GlobalDataStruct *) gd;

   m_MetricCalculationLock.Lock();
   m_SumOfSquaredDifference  += globalData->m_SumOfSquaredDifference;
   m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
   m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
   if ( m_NumberOfPixelsProcessed )
   {
     m_Metric = m_SumOfSquaredDifference / static_cast<double>( m_NumberOfPixelsProcessed );
   }
   m_MetricCalculationLock.Unlock();
   delete globalData;
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
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::InitializeIteration()
  {
    if( !this->GetMovingImage() || !this->GetFixedImage() )
    {
      itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    }

    // cache fixed image information
    m_FixedImageSpacing    = this->GetFixedImage()->GetSpacing();
    m_FixedImageOrigin     = this->GetFixedImage()->GetOrigin();

    // initialize metric computation variables
    m_SumOfSquaredDifference  = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_SumOfSquaredChange      = 0.0;
      
    /* Get the images' size */
    m_imagesSize = this->GetMovingImage()->GetLargestPossibleRegion().GetSize();
    /* Get the pointer of the first pixel of the certainties image */
    m_FirstCertaintyPixel = &( this->GetCertaintyImage()->GetPixel( m_zerosIndex ) );
    /* Get the pointer of the first pixel of moving the image */
    m_FirstMovingPixel = &( this->GetMovingImage()->GetPixel( m_zerosIndex ) );
    /* Get the pointer of the first pixel of fixed the image */
    m_FirstFixedPixel = &( this->GetFixedImage()->GetPixel( m_zerosIndex ) );
    /* Get the images' spacing*/
    m_imagesSpacing = this->GetFixedImage()->GetSpacing();
      
    /**
     *      Compute the convolutions in all directions for the fixed and moving image
     */
    /* take the number of quadrature filters */
    int totNumFilter = m_Filters->GetNumMatrix();

    typedef itk::Image < float, ImageDimension >                                InternalImageType;
    typedef itk::SmartPointer<Image < float, ImageDimension > >                 InternalImagePointer;

    /**
     * Before compute the convolutions and if there is an deformation field, deforme the moving image.
     */
    if ( this->GetDeformationField() )
    {
      typedef itk::WarpImageFilter < MovingImageType, MovingImageType, DeformationFieldType> WarperType;
        typename WarperType::Pointer warper = WarperType::New();
      typedef itk::LinearInterpolateImageFunction< MovingImageType, double >        InterpolatorType;
        typename InterpolatorType::Pointer interpolatorSec = InterpolatorType::New();
	      
      warper->SetInput( this->GetRescaledPrototype() );//MovingMirrorPad->GetOutput() );
      warper->SetInterpolator( interpolatorSec );
      warper->SetOutputSpacing(this->GetMovingImage()->GetSpacing() );
      warper->SetOutputOrigin( this->GetMovingImage()->GetOrigin() );
      warper->SetDeformationField( this->GetDeformationField() );

      try
      {
        warper->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr << "Exception catched in moving image warper" << std::endl;
        std::cerr << excep << std::endl;
        throw(EXIT_FAILURE);
      }

      this->SetMovingImage( warper->GetOutput() );
    }
    else
    {
      std::cerr<<"Error in the itkSymmetricMorphonsRegistrationFunction"<<std::endl;
      std::cerr<<"We lost the deformation field..."<<std::endl;
      throw(EXIT_FAILURE);
    }
    if( IterationNum==1 )
    {
      /**
       *  Make convolution for fixed and moving image.
       *  Because it's the first iteration, it's necessary to 
       *  compute convolutions for the fixed and moving image.
       */
       for(int i=0; i<totNumFilter; i++)
       {
         /** i-th convolution with the fixed image. */
         m_FixedConv[i] = ComplexMatrixType::New();
         m_FixedConv[i]->SetRegions( this->GetFixedImage()->GetRequestedRegion() );
         m_FixedConv[i]->SetSpacing( this->GetFixedImage()->GetSpacing());
         m_FixedConv[i]->SetOrigin( this->GetFixedImage()->GetOrigin());
         m_FixedConv[i]->Allocate();
         Convolution(i, this->GetFixedImage(), m_FixedConv[i] );

         /** i-th convolution with the moving image. */
         m_MovingConv[i] = ComplexMatrixType::New();
         m_MovingConv[i]->SetRegions( this->GetMovingImage()->GetRequestedRegion() );
         m_MovingConv[i]->SetSpacing( this->GetMovingImage()->GetSpacing());
         m_MovingConv[i]->SetOrigin( this->GetMovingImage()->GetOrigin());
         m_MovingConv[i]->Allocate();
         Convolution(i, this->GetMovingImage(), m_MovingConv[i] );

         /** Compute the i-th complex mutliplication */
         m_qqin[i] = ComplexMatrixType::New();
         m_qqin[i]->SetRegions( this->GetFixedImage()->GetRequestedRegion() );
         m_qqin[i]->SetSpacing( this->GetFixedImage()->GetSpacing());
         m_qqin[i]->SetOrigin( this->GetFixedImage()->GetOrigin());
         m_qqin[i]->Allocate();
         ComplexMultConj(m_qqin[i], m_FixedConv[i], m_MovingConv[i] );
       }
    }
    else
    {
     /**
      * Because fixed image doesn't change only the moving 
      * image has to be convoluted
      */
      for(int i=0; i<totNumFilter; i++)
      {
        /** i-th real convolution with the moving image. */
        Convolution(i, this->GetMovingImage(), m_MovingConv[i] );
        /** i-th compute complex mutliplication */
        ComplexMultConj(m_qqin[i], m_FixedConv[i], m_MovingConv[i] );
      }//end for
    }//end if
  }//end InitializeIteration()

 /**
  * virtual PixelType  ComputeUpdate
  *        (
  *         const NeighborhoodType &it,
  *         void *gd,
  *         const FloatOffsetType &offset = FloatOffsetType(0.0)
  *        )
  * This method is called by a finite difference solver image filter at 
  * each pixel that does not lie on a data set boundary
  * The Morphon method compute update for the deformation field and for the
  * certainty image.
  */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  typename ComputingMorphonDeformationField<TFixedImage,TMovingImage,TDeformationField>
  ::PixelType
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::ComputeUpdate(const NeighborhoodType &it,
                                   void *gd,
                                   const FloatOffsetType &offset)
  {
    GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
    PixelType update;
    IndexType index = it.GetIndex();

    double fixedValue = (double) this->GetFixedImage()->GetPixel( index );
    double movingValue= (double) this->GetMovingImage()->GetPixel( index );
    double localDeff[ImageDimension];
    double oldCertitude=0.0, newCertitude=0.0;
    double m_vk, m_ck2, m_a11, m_a22, m_a12, m_a13, m_a23, m_a33, m_b1, m_b2, m_b3, m_Deta;
    if( ImageDimension == 3)
    {
      /* Initialisation */
      int totNumMatrix = m_Filters->GetNumMatrix();
      m_vk = 0.0; m_ck2 = 0.0; m_a11 = 0.0; m_a22 = 0.0; m_a12 = 0.0; m_a13 = 0.0; m_a23 = 0.0; m_a33 = 0.0; m_b1 =0.0; m_b2 = 0.0; m_b3 =0.0; m_Deta=0.0;
      const double ETA = -1;
      typename ComplexMatrixType::PixelType local;
	
      /* Compute the coefficients */ 
      for (unsigned int i=0; i<totNumMatrix; i++)
      {
        local = m_qqin[i]->GetPixel( index );
        m_vk = local.imag() / ( abs(local) + eps ) ; 
        m_ck2 = sqrt( abs(local) ) * pow( cos(m_vk/2), 4 );
        m_a11 = m_a11 + m_ck2 * m_directionnalCoef[i][0] * m_directionnalCoef[i][0];
        m_a12 = m_a12 + m_ck2 * m_directionnalCoef[i][0] * m_directionnalCoef[i][1];
        m_a22 = m_a22 + m_ck2 * m_directionnalCoef[i][1] * m_directionnalCoef[i][1];
        m_a13 = m_a13 + m_ck2 * m_directionnalCoef[i][0] * m_directionnalCoef[i][2];
        m_a23 = m_a23 + m_ck2 * m_directionnalCoef[i][1] * m_directionnalCoef[i][2];
        m_a33 = m_a33 + m_ck2 * m_directionnalCoef[i][2] * m_directionnalCoef[i][2];
        m_b1 = m_b1 + m_ck2 * m_directionnalCoef[i][0] * m_vk;
        m_b2 = m_b2 + m_ck2 * m_directionnalCoef[i][1] * m_vk;
        m_b3 = m_b3 + m_ck2 * m_directionnalCoef[i][2] * m_vk;
      }//end for
	    
      newCertitude =  m_a11 + m_a22 + m_a33;
      m_Deta = ( -pow(m_a13,2)*m_a22 + (2*m_a12*m_a13 *m_a23) - m_a11*pow(m_a23,2)
                 -pow(m_a12,2)*m_a33 + (m_a11*m_a22*m_a33) );
                             
      /** computes certainties */ 
      oldCertitude = *(m_FirstCertaintyPixel + ((index[1])*m_imagesSize[0]) + (index[2]*m_imagesSize[1]*m_imagesSize[0]) + index[0]);
      double tempFactor = newCertitude + oldCertitude + eps;
      if( m_Deta ==0 )
      {
        *(m_FirstCertaintyPixel+index[0]*m_imagesSize[0]+index[1]*m_imagesSize[1]+index[2]*m_imagesSize[2])=0;
      }
      else
      {
        *(m_FirstCertaintyPixel + ((index[1])*m_imagesSize[0]) + (index[2]*m_imagesSize[1]*m_imagesSize[0]) + index[0]) =
            static_cast<typename MovingImageType::PixelType> ( ( pow(oldCertitude,2) + pow(newCertitude, 2) ) / tempFactor );
      }    
      /** computes deformations */    
      if( m_Deta == 0 )
      {
        localDeff[0] = 0;
        localDeff[1] = 0;
        localDeff[2] = 0;
      }
      else
      {
        m_Deta = 1 / m_Deta;
        localDeff[0] = 
          (
           m_Deta *
           (m_a13*m_a23*m_b1 -
            m_a12*m_a33*m_b1 -
            pow(m_a13,2)*m_b2 +
            m_a11*m_a33*m_b2 +
            m_a12*m_a13*m_b3 -
            m_a11*m_a23*m_b3 )
            * ETA
           );
        localDeff[1] = 
          (
           m_Deta *
           (-pow(m_a23,2)*m_b1 +
           m_a22*m_a33*m_b1 +
           m_a13*m_a23*m_b2 -
           m_a12*m_a33*m_b2 -
           m_a13*m_a22*m_b3 +
           m_a12*m_a23*m_b3 )
           * ETA
         );
          localDeff[2] = 
           (
            m_Deta *
            (-m_a13*m_a22*m_b1 +
            m_a12*m_a23*m_b1 +
            m_a12*m_a13*m_b2 -
            m_a11*m_a23*m_b2 -
            pow(m_a12,2)*m_b3 +
            m_a11*m_a22*m_b3 )
            * ETA 
           );
      }
      update[0] = static_cast < typename MovingImageType::PixelType > 
                  ( ( newCertitude  * localDeff[0]  / tempFactor ) );
      update[1] = static_cast < typename MovingImageType::PixelType > 
                  ( ( newCertitude  * localDeff[1]  / tempFactor ) );
      update[2] = static_cast < typename MovingImageType::PixelType > 
                  ( ( newCertitude  * localDeff[2]  / tempFactor ) );
          
      update[0] = update[0] * m_imagesSpacing[0];
      update[1] = update[1] * m_imagesSpacing[1];
      update[2] = update[2] * m_imagesSpacing[2];

    }//end if
    else if(ImageDimension == 2)
    {           
      /* Initialisation */         
      const double ETA =2.0; /* ETA is an unchanged constant */
      int totNumMatrix = m_Filters->GetNumMatrix();
      m_a11 = 0; m_a22 = 0; m_a12 = 0; m_b1 =0; m_b2 = 0;
      std::complex<double> local;
      /* Compute the coefficients */  
      for (unsigned int i=0; i<totNumMatrix; i++)
      {
        local = m_qqin[i]->GetPixel( index );
        m_vk = local.imag() / ( abs(local) + eps ) ;
        m_ck2 = ( sqrt( abs(local) )
                  *pow( cos ( atan2( (local.imag()), local.real()) * 0.5) ,4) )
                  + eps  ;

        m_a11 = m_a11 + m_ck2 * m_directionnalCoef[i][0] * m_directionnalCoef[i][0];
        m_a12 = m_a12 + m_ck2 * m_directionnalCoef[i][1] * m_directionnalCoef[i][0];
        m_a22 = m_a22 + m_ck2 * m_directionnalCoef[i][1] * m_directionnalCoef[i][1];
        m_b1 = m_b1 + ( m_ck2 * m_directionnalCoef[i][0] * m_vk );
        m_b2 = m_b2 + ( m_ck2 * m_directionnalCoef[i][1] * m_vk );
      }//end for
      m_Deta = (m_a11 * m_a22 - pow(m_a12,2) ) ;
      if( m_Deta == 0 ) { m_Deta=1; }
      else { m_Deta = 1 / m_Deta; }
      newCertitude = m_a11 + m_a22;
      oldCertitude = *(m_FirstCertaintyPixel + (index[1]*m_imagesSize[0]) + index[0]);
      double tempFactor = newCertitude + oldCertitude + eps;

       
       *(m_FirstCertaintyPixel + (index[1]*m_imagesSize[0]) +index[0])=
                 ( pow(oldCertitude,2) + pow(newCertitude, 2) ) / tempFactor;
       

       localDeff[0] = m_Deta * (-m_a12 * m_b1 + m_a11 * m_b2 ) * ETA ;
       localDeff[1] = m_Deta * ( m_a22 * m_b1 - m_a12 * m_b2 ) * ETA ;

       update[0] = static_cast < typename MovingImageType::PixelType >
                  ( ( newCertitude  * localDeff[0]  / tempFactor ) );
       update[1] = static_cast < typename MovingImageType::PixelType >
                  ( ( newCertitude  * localDeff[1]  / tempFactor ) );
		
       /* Compute a deformation in the real image  */
       update[0] = update[0] * m_imagesSpacing[0];
       update[1] = update[1] * m_imagesSpacing[1];
            
    }//end if
    else
    {
      itkExceptionMacro(<< "Dimension must be 2 or 3. " );
      return EXIT_FAILURE;
    }

    globalData->m_SumOfSquaredDifference += vnl_math_sqr( fixedValue - movingValue );
    globalData->m_NumberOfPixelsProcessed += 1;
    return update;
  }//end ComputeUpdate
    
  /**
   * void SetQuadCoef
   *     (
   *       ComplexMatrixPointerType * FixedRealConv,
   *       ComplexMatrixPointerType * MovingRealConv,
   *       ComplexMatrixPointerType * qqinReal,
   *     )
   * Used to set the quadrature coefficients.
   * This fonction is usefull to take previous iteration coefficients because the
   * itkComputingMorphonDeformationField object is created in each iteration.
   */  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::SetQuadCoef(
           ComplexMatrixPointerType* FixedConv,
           ComplexMatrixPointerType* MovingConv,
           ComplexMatrixPointerType* qqin
              )
    {
      m_FixedConv     =   FixedConv;
      m_MovingConv    =   MovingConv;
      m_qqin          =   qqin;
    }//end SetQuadCoef
    
  /**
   * void GetQuadCoef
   *           (
   *            ComplexMatrixPointerType * FixedRealConv,
   *            ComplexMatrixPointerType * MovingRealConv,
   *            ComplexMatrixPointerType * qqinReal,
   *           )
   * Used to save the quadrature coefficients out this local object.
   * This fonction is usefull to save this iteration coefficients because the
   * itkComputingMorphonDeformationField object is created in each iteration.
   */ 
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::GetQuadCoef(
           ComplexMatrixPointerType * FixedConv,
           ComplexMatrixPointerType * MovingConv,
           ComplexMatrixPointerType * qqin
            )
    {
      FixedConv     = m_FixedConv;
      MovingConv    = m_MovingConv;
      qqin          = m_qqin;
    }//end GetQuadCoef

  /**
   * ComputingMorphonDeformationField()
   * constructor
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::ComputingMorphonDeformationField()
  {
    RadiusType r;
    unsigned int j;
    for( j = 0; j < ImageDimension; j++ ) { r[j] = 0; }
    this->SetRadius(r);

    m_TimeStep = 1.0;
    this->SetMovingImage( NULL );
    this->SetFixedImage(NULL);
    m_FixedImageSpacing.Fill( 1.0 );
    m_FixedImageOrigin.Fill( 0.0 );

    typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

    m_Metric = NumericTraits<double>::max();
    m_SumOfSquaredDifference = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_SumOfSquaredChange = 0.0;

    IterationNum = 1;
      
    m_Filters = NULL;
    for( unsigned int i=0; i<ImageDimension; i++) { m_zerosIndex[i]=0; }
      
  }
    
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "Metric: ";
    os << m_Metric << std::endl;
    os << indent << "SumOfSquaredDifference: ";
    os << m_SumOfSquaredDifference << std::endl;
    os << indent << "NumberOfPixelsProcessed: ";
    os << m_NumberOfPixelsProcessed << std::endl;
    os << indent << "SumOfSquaredChange: ";
    os << m_SumOfSquaredChange << std::endl;
    os << indent << " iteration's number : ";
    os << IterationNum << std::endl;
    os << indent << " number of quad. filters : ";
    os << m_Filters->GetNumMatrix() << std::endl;

  }
    
  /**
   * inline void ComplexMultConj(
   *            (
   *             ComplexMatrixPointerType result,
   *             ComplexMatrixPointerType a,
   *             ComplexMatrixPointerType b
   *            )
   *
   * Used to compute the complex multiplication.
   *             result = (a)*conj(b)
   * with conj(c + j d) = (c - j d)
   *
   * Warning : All matrix must have the same dimension.
   */
  template<class TFixedImage, class TMovingImage, class TDeformationField>        
  inline void 
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>  
  ::ComplexMultConj
             (
              ComplexMatrixPointerType result,
              ComplexMatrixPointerType a,
              ComplexMatrixPointerType b
             )
  {
   typename ComplexMatrixType::IndexType index;
   for( unsigned int i=0; i<ImageDimension; i++) { index[i]=0; }
    
   typename ComplexMatrixType::SizeType Imsize = a->GetRequestedRegion().GetSize();
   int numPt = 1;
   for( unsigned int i=0; i<ImageDimension; i++) { numPt = numPt * Imsize[i]; }
   
   typename ComplexMatrixType::PixelType* APixel = &( a->GetPixel( index ) );
   typename ComplexMatrixType::PixelType* BPixel = &( b->GetPixel( index ) );
   typename ComplexMatrixType::PixelType* RPixel = &( result->GetPixel( index ) );
   for(unsigned int i=0; i<numPt; i++)
   {
     *RPixel = (*APixel) * ( conj(*BPixel) );
     RPixel += 1;
     APixel += 1;
     BPixel += 1;
   }//end for
 }//end ComplexMultConj
 
 /**
   * void Convolution
   *     (
   *      int choice,
   *      const MovingImageType * InImage,
   *      ComplexMatrixPointerType OutImage
   *     )
   *
   * This function is used to compute the convolutions.
   * choice defined the filter number choice.
   * InImage and OutImage are respectivelly the input and output image
   *
   * Warning : All matrix must have the same dimension.
   */      
  template<class TFixedImage, class TMovingImage, class TDeformationField>        
  inline void 
  ComputingMorphonDeformationField<TFixedImage, TMovingImage, TDeformationField>  
  ::Convolution
  (
    int choice,
    const MovingImageType * InImage,
    ComplexMatrixPointerType OutImage
  )
  {
    typedef itk::ConstNeighborhoodIterator< MovingImageType >  NeighborhoodIteratorType;
    typedef itk::ImageRegionIterator<ComplexMatrixType> ComplexNeighborhoodOperatorType;
      ComplexNeighborhoodOperatorType out;
    typedef typename FilterType::ComplexOperator   ComplexOpType;
      itk::NeighborhoodInnerProduct< MovingImageType, ComplexType, ComplexType > innerProduct;
      
    /* Quadrature Filters */     
    ComplexOpType* Complexfilter = m_Filters->GetComplexOp(choice) ;
    NeighborhoodIteratorType InIT = NeighborhoodIteratorType( (*Complexfilter).GetRadius() , InImage, InImage->GetRequestedRegion() );
     
    std::complex<double> temp;        
    out = ComplexNeighborhoodOperatorType( OutImage, OutImage->GetRequestedRegion() );
    
    for(InIT.GoToBegin(), out.GoToBegin(); !InIT.IsAtEnd(); ++InIT, ++out )
    {
      temp = innerProduct ( InIT, (*Complexfilter) ) ;
      out.Set( conj(temp) );
    }//end for
  }//end Convolution

};//end class

#endif
