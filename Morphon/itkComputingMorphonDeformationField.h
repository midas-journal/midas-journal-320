/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputingMorphonDeformationField.h,v $
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

#ifndef _itkComputingMorphonDeformationField_h_
#define _itkComputingMorphonDeformationField_h_

#include "itkMorphonToolbooxFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodOperator.h"

#include "itkBuildingMorphonFilters.h"

#include <string.h>

#include "itkWarpImageFilter.h"


namespace itk {

/**
 * \class ComputingMorphonDeformationField
 *
 * This class encapsulate the Morphon method, the deformation field is here 
 * compute in the ComputeUpdate function. The deformation field is computed such
 * as described in the [5].
 * The convolution with all quadrature filters are performed in the InitializeIteration.
 *
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type, 
 * certainty image type and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa MorphonToolboxFunction
 * \ingroup FiniteDifferenceFunctions
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT ComputingMorphonDeformationField :
    public MorphonToolbooxFunction< TFixedImage,
                                              TMovingImage, TDeformationField>
{
public:
 /** Standard class typedefs. */
 typedef ComputingMorphonDeformationField    Self;
 typedef MorphonToolbooxFunction< TFixedImage, TMovingImage, TDeformationField >    Superclass;
 typedef SmartPointer<Self> Pointer;
 typedef SmartPointer<const Self> ConstPointer;

 /** Method for creation through the object factory. */
 itkNewMacro(Self);

 /** Run-time type information (and related methods). */
 itkTypeMacro( ComputingMorphonDeformationField, MorphonToolbooxFunction );

 /** MovingImage image type. */
 typedef typename Superclass::MovingImageType     MovingImageType;
 typedef typename Superclass::MovingImagePointer  MovingImagePointer;

 /** FixedImage image type. */
 typedef typename Superclass::FixedImageType     FixedImageType;
 typedef typename Superclass::FixedImagePointer  FixedImagePointer;
 typedef typename FixedImageType::IndexType      IndexType;
 typedef typename FixedImageType::SizeType       SizeType;
 typedef typename FixedImageType::SpacingType    SpacingType;

 /** Deformation field type. */
 typedef typename Superclass::DeformationFieldType          DeformationFieldType;
 typedef typename Superclass::DeformationFieldTypePointer   DeformationFieldTypePointer;

 /** Inherit some enums from the superclass. */
 itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

 /** Internal Image type */
 typedef itk::SmartPointer<Image < float, ImageDimension > > InternalImagePointer;

 /** Inherit some enums from the superclass. */
 typedef typename Superclass::PixelType          PixelType;
 typedef typename Superclass::RadiusType         RadiusType;
 typedef typename Superclass::NeighborhoodType   NeighborhoodType;
 typedef typename Superclass::FloatOffsetType    FloatOffsetType;
 typedef typename Superclass::TimeStepType       TimeStepType;

 /** Interpolator type. */
 typedef double CoordRepType;
 typedef InterpolateImageFunction<MovingImageType,CoordRepType>       InterpolatorType;
 typedef typename InterpolatorType::Pointer                           InterpolatorPointer;
 typedef typename InterpolatorType::PointType                         PointType;
 typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType> DefaultInterpolatorType;

 /** Covariant vector type. */
 typedef CovariantVector<double,itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

 /** Gradient calculator type. */
 typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
 typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;
 
 /** To import the filters */
 typedef itk::BuildingMorphonFilters< MovingImageType, FixedImageType > FilterType;
 typedef typename FilterType::Pointer 		 	   FilterPointerType;
  
 /** Complex type*/
 typedef typename FilterType::ComplexType		   ComplexType;
 typedef itk::Image< ComplexType, ImageDimension >        ComplexMatrixType;
 typedef typename FilterType::VectorCoeffType 		   VectorCoeffType;
 typedef typename ComplexMatrixType::Pointer              ComplexMatrixPointerType;

 /** This class uses a constant timestep of 1. */
 virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const
 { return m_TimeStep; }

 /** Return a pointer to a global data structure that is passed to
  * this object from the solver at each calculation.  */
 virtual void *GetGlobalDataPointer() const;

 /** Release memory for global data structure. */
 virtual void ReleaseGlobalDataPointer( void *gd ) const;

 /** Used to count the iteration's number. */
 int IterationNum;
    
 /** To set the filters */
 void SetFilters( FilterPointerType filters)
 {
   m_Filters = filters;
   /** initiate the directionnal coefficients */
   m_directionnalCoef = m_Filters->GetDirectionnalCoeff();
 }

 /** To get the filters */
 FilterType* GetFilters(){ return m_Filters; }

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
 virtual void InitializeIteration();

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
  inline virtual PixelType  ComputeUpdate(const NeighborhoodType &it,
                                   void *gd,
                                   const FloatOffsetType &offset = FloatOffsetType(0.0));

  /** 
   * Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images. 
   */
  virtual double GetMetric() const { return m_Metric; }

  /**
   * void setIterationNum(int iter)
   * Used to fixe the iteration's number.
   */
  void setIterationNum(int iter){ this->IterationNum = iter;}

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
  void SetQuadCoef(
           ComplexMatrixPointerType* FixedRealConv,
           ComplexMatrixPointerType* MovingRealConv,
           ComplexMatrixPointerType* qqinReal
           );

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
  void GetQuadCoef(
           ComplexMatrixPointerType* FixedRealConv,
           ComplexMatrixPointerType* MovingRealConv,
           ComplexMatrixPointerType* qqinReal
            );
    
 /**
  * void SetPadding( bool padding )
  *  enable (if padding==true) or disable the padding during the Morphon algorithm
  */
 void SetPadding( bool padding ) { m_padding = padding; }
    
 /**
  * bool GetPadding()
  *  return the padding value
  */
  bool GetPadding() { return m_padding; }

protected:
  /**
   * Default constructor
   */
  ComputingMorphonDeformationField();

  ~ComputingMorphonDeformationField() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** 
   * A global data type for this class of equation. Used to store
   * iterators for the fixed image. 
   */
  struct GlobalDataStruct
  {
    double          m_SumOfSquaredDifference;
    unsigned long   m_NumberOfPixelsProcessed;
    double          m_SumOfSquaredChange;
  };

private:

  ComputingMorphonDeformationField(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
 /** A pointer to have an access to the filters */
 FilterPointerType m_Filters;
 double** m_directionnalCoef;
 static const double eps = 2.2204e-016; /*a very small constant*/
 
 /** padding option*/
  bool m_padding;
 
 /*save global information to do not compute them at each time*/
 typename MovingImageType::IndexType m_zerosIndex;
 typename MovingImageType::SizeType  m_imagesSize;
 typename MovingImageType::PixelType *m_FirstCertaintyPixel;
 const typename MovingImageType::PixelType *m_FirstMovingPixel;
 const typename FixedImageType::PixelType *m_FirstFixedPixel;
 typename MovingImageType::SpacingType m_imagesSpacing;
 
 
 /** Cache fixed image information. */
 SpacingType                     m_FixedImageSpacing;
 PointType                       m_FixedImageOrigin;

 /** The global timestep. */
 TimeStepType                    m_TimeStep;

 /** Threshold below which two intensity value are assumed to match. */
 double                          m_IntensityDifferenceThreshold;

 /** 
   * The metric value is the mean square difference in intensity between
   * the fixed image and transforming moving image computed over the
   * the overlapping region between the two images.
   */
  mutable double                  m_Metric;
  mutable double                  m_SumOfSquaredDifference;
  mutable unsigned long           m_NumberOfPixelsProcessed;
  mutable double                  m_SumOfSquaredChange;

  /** Pointers to set the convolutions  Results from the previous iterations.*/
  ComplexMatrixPointerType* m_FixedConv;
  ComplexMatrixPointerType* m_MovingConv;
  ComplexMatrixPointerType* m_qqin;

  /** Mutex lock to protect modification to metric. */
  mutable SimpleFastMutexLock     m_MetricCalculationLock;

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
  inline void Convolution
             (
               int choice,
               const MovingImageType * InImage,
               ComplexMatrixPointerType OutImage
             );

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
          
  inline void ComplexMultConj
             (
              ComplexMatrixPointerType result,
              ComplexMatrixPointerType a,
              ComplexMatrixPointerType b
             );

};//end class


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputingMorphonDeformationField.txx"
#endif

#endif
