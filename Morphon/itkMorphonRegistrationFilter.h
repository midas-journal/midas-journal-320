/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonRegistrationFilter.h,v $
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

#ifndef _itkMorphonRegistrationFilter_h_
#define _itkMorphonRegistrationFilter_h_

#include "itkMorphonToolboxFiler.h"
#include "itkComputingMorphonDeformationField.h"

#include "itkBuildingMorphonFilters.h"

namespace itk {

/** \class itkMorphonRegistrationFilter
 * \brief Deformably register two images using the Morphon algorithm.
 *
 * SymmetricForcesDemonsRegistrationFilter implements the demons deformable algorithm that
 * register two images by computing the deformation field which will map a
 * moving image onto a fixed image.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation field type.
 *
 * The input fixed and moving images and are set via methods 
 * SetFixedImage, SetMovingImage and respectively.
 * An initial deformation field maybe set via SetInitialDeformationField or
 * SetInput. If no initial field is set, a zero field is used as the initial
 * condition.
 *
 * The output deformation field can be obtained via methods GetOutput
 * or GetDeformationField.
 *
 * This class make use of the finite difference solver hierarchy. Update
 * for each iteration is computed in itkComputingMorphonDeformationField.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa itkMorphonRegistrationFilter
 * \ingroup DeformableImageRegistration MultiThreaded
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT MorphonRegistrationFilter :
    public MorphonToolboxFiler< TFixedImage, TMovingImage,
                                            TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef MorphonRegistrationFilter    Self;
  typedef MorphonToolboxFiler<
    TFixedImage, TMovingImage,TDeformationField>    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( MorphonRegistrationFilter, MorphonToolboxFiler );

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType       FixedImageType;
  typedef typename Superclass::FixedImagePointer    FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;

  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer DeformationFieldPointer;

  /** FiniteDifferenceFunction type. */
  typedef typename Superclass::FiniteDifferenceFunctionType  FiniteDifferenceFunctionType;

  /** Take timestep type from the FiniteDifferenceFunction. */
  typedef typename FiniteDifferenceFunctionType::TimeStepType  TimeStepType;


  /** DemonsRegistrationFilterFunction type. */
  typedef ComputingMorphonDeformationField<FixedImageType,MovingImageType,
                                     DeformationFieldType>  MorphonsRegistrationFunctionType;

  /** Inherit some enums and typedefs from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);
  
  /** To import the filters */
  typedef itk::BuildingMorphonFilters< MovingImageType, FixedImageType > FilterType;
  typedef typename FilterType::Pointer 		 	   FilterPointerType;
  
  /** Complex type*/
  typedef typename FilterType::ComplexType		   ComplexType;
  typedef itk::Image< ComplexType, ImageDimension >        ComplexMatrixType;
  typedef typename ComplexMatrixType::Pointer              ComplexMatrixPointerType;

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
  virtual double GetMetric() const;
  
 /** To set the filters */
 void SetFilters( FilterPointerType filter){ m_Filters=filter;}
 /** To get the filters */
 FilterPointerType GetFilters(){ return m_Filters;}
 
 /*
  * void SetPadding( bool padding )
  *  enable (if padding==true) or disable the padding during the Morphon algorithm
  */
  void SetPadding( bool padding ) { m_padding = padding; }
    
 /*
  * bool GetPadding()
  *  return the padding value
  */
  bool GetPadding() { return m_padding; }

protected:
  /*
   * Default constructor
   */
  MorphonRegistrationFilter()
  {
    m_Filters = NULL;
    typename MorphonsRegistrationFunctionType::Pointer mrfp;
    mrfp = MorphonsRegistrationFunctionType::New();
    tempCert = MovingImageType::New();

    this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(mrfp.GetPointer() ) );

    IterationNum = 1;
    
   }

  ~MorphonRegistrationFilter() {}

  void PrintSelf(std::ostream& os, Indent indent) const;

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
   * void ApplyUpdate(TimeStepType dt)
   * This function is used to call the superClass' applyUpdate (as it's usually
   * done) and to make a normalized convolution over the deformation field and 
   * the certainty image (as explain in Normalized Convolution: {A} Technique 
   * for Filtering Incomplete and Uncertain Data, H. Knutsson and C-F. Westin, 
   * 1993)
   */
  virtual void ApplyUpdate(TimeStepType dt);

private:
  /** padding option*/
  bool m_padding;
  
  inline void MatrixDiv
        (
            MovingImagePointer result,
            const MovingImagePointer a,
            const MovingImagePointer b
        );
  /** Filters */
  FilterPointerType m_Filters;
  
  /** Used to count the iteration's number.*/
  int IterationNum;
  
  typename MovingImageType::Pointer tempCert;

  /** Results come from the convolutions */
  ComplexMatrixPointerType m_FixedConv[ImageDimension*2];
  ComplexMatrixPointerType m_MovingConv[ImageDimension*2];
  ComplexMatrixPointerType m_qqin[ImageDimension*2];
  
  MorphonRegistrationFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

   /**
    *       void DefMult(
    *        itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > result,
    *        const itk::SmartPointer < itk::Image< itk::Vector<float, ImageDimension>, ImageDimension > > a,
    *       const itk::SmartPointer < itk::Image< float, ImageDimension > > b
    *        )
    *
    *       Used to compute the division between a deformation field and an certainty image.
    *                   result->GetPixel(index) = a->GetPixel(index)*b->GetPixel(b)
    *       All matrix must have the same dimension.
    **/

    inline void DefMult
        (
            DeformationFieldPointer result,
            const DeformationFieldPointer a,
            const MovingImagePointer b
        );

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

    inline void DefDiv
        (
            DeformationFieldPointer result,
            const DeformationFieldPointer a,
            const MovingImagePointer b
        );

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
    **/

    inline void MatrixMult
        (
            MovingImagePointer result,
            const MovingImagePointer a,
            const MovingImagePointer b
        );

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
    **/

    inline void voidMatrixDiv
        (
           MovingImagePointer result,
           const MovingImagePointer a,
           const MovingImagePointer b
        );

};//end class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMorphonRegistrationFilter.txx"
#endif


#endif
