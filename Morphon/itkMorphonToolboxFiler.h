/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonToolboxFiler.h,v $
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

#ifndef _itkMorphonToolboxFiler_h_
#define _itkMorphonToolboxFiler_h_

#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkMorphonToolbooxFunction.h"

#include "itkExceptionObject.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkDataObject.h"

#include "itkGaussianOperator.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkNeighborhoodOperatorImageFilter.h"

#include "vnl/vnl_math.h"


#include "itkDiscreteGaussianImageFilter.h"

namespace itk {

/**
 * \class MorphonToolboxFiler
 * \brief The Morphon algorithm toolbox over fixed and moving images, certainty 
 * matrix and the deformation field.
 *
 * MorphonToolboxFiler is a base case for filter implementing.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation Field type.
 *
 * The input fixed and moving images are set via methods SetFixedImage
 * and SetMovingImage respectively. Because Morphon algorithm uses a certainty
 * matrix, this matrix need to be set with the SetCertaintyImage. 
 * An initial deformation field maybe set via
 * SetInitialDeformationField or SetInput. If no initial field is set,
 * a zero field is used as the initial condition.
 *
 * The output deformation field can be obtained via methods GetOutput
 * or GetDeformationField.
 *
 * To compute the normalized smmoothing easily, three functions could be :
 * SmoothDeformationField, SmoothCertaintyImage and LocalSmooth.
 *
 *
 * The algorithm is run for a user defined number of iterations.
 * Typically the algorithm requires period Gaussin smoothing of the
 * deformation field to enforce an elastic-like condition. The amount
 * of smoothing is governed by a set of user defined standard deviations
 * (one for each dimension).
 *
 * In terms of memory, this filter keeps two internal buffers: one for storing
 * the intermediate updates to the field and one for double-buffering when
 * smoothing the deformation field. Both buffers are the same type and size as the
 * output deformation field.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa MorphonToolboxFiler.
 * \ingroup DeformableImageRegistration
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT MorphonToolboxFiler :
    public DenseFiniteDifferenceImageFilter<TDeformationField,TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef MorphonToolboxFiler    Self;
  typedef DenseFiniteDifferenceImageFilter<TDeformationField,TDeformationField>    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( itkMorphonToolboxFiler,
                DenseFiniteDifferenceImageFilter );

  /** FixedImage image type. */
  typedef TFixedImage   FixedImageType;
  typedef typename FixedImageType::Pointer  FixedImagePointer;
  typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;

  /** MovingImage image type. */
  typedef TMovingImage    MovingImageType;
  typedef typename MovingImageType::Pointer  MovingImagePointer;
  typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;

  /** Deformation field type. */
  typedef TDeformationField    DeformationFieldType;
  typedef typename DeformationFieldType::Pointer  DeformationFieldPointer;

  /** Types inherithed from the superclass */
  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef std::complex<float>			  ComplexPixelType;

  /** FiniteDifferenceFunction type. */
  typedef typename Superclass::FiniteDifferenceFunctionType
  FiniteDifferenceFunctionType;

  /** PDEDeformableRegistrationFilterFunction type. */
  typedef MorphonToolbooxFunction<FixedImageType,MovingImageType,
                                            DeformationFieldType>  PDEDeformableRegistrationFunctionType;

  /** Inherit some enums and typedefs from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /**
   * void SetFixedImage( const FixedImageType * ptr )
   * Set the fixed image
   */
  void SetFixedImage( const FixedImageType * ptr );

  /**
   * FixedImageType * GetFixedImage(void) const
   * Return the fixed image
   */
  const FixedImageType * GetFixedImage(void) const;

  /**
   * void SetMovingImage(const MovingImageType * ptr )
   * Set the Moving image
   */
  void SetMovingImage( const MovingImageType * ptr );

  /**
   * void SetCertaintyImage(const MovingImageType * ptr )
   * Set the Certainty matrix
   */
  void SetCertaintyImage( MovingImageType * ptr );

  /**
   * void SetRescaledPrototype(const MovingImageType * ptr )
   * Set the Rescaked Prototype
   */
  void SetRescaledPrototype( MovingImageType * ptr );

  /**
   * MovingImageType * GetMovingImage(void) const
   * Return the Moving image
   */
  const MovingImageType * GetMovingImage(void) const;

  /**
   * MovingImageType * GetCertaintyImage(void) const
   * Return the Certainty Matrix
   */
  MovingImageType * GetCertaintyImage(void) const;

  /**
   * MovingImageType * GetRescaledPrototype(void) const
   * Return the Rescaled Prototype
   */
  MovingImageType * GetRescaledPrototype(void) const;

  /**
   * void SetInitialDeformationField( DeformationFieldType * ptr )
   * Set the initial Deformation Field
   */
  void SetInitialDeformationField( DeformationFieldType * ptr );

  /**
   * DeformationFieldType * GetDeformationField()
   * Return the Deformation Field
   */
  DeformationFieldType * GetDeformationField();

  /** Get the number of valid inputs.  For PDEDeformableRegistration,
   * this checks whether the fixed and moving images have been
   * set. While PDEDeformableRegistration can take a third input as an
   * initial deformation field, this input is not a required input.
   */
  virtual std::vector<SmartPointer<DataObject> >::size_type GetNumberOfValidRequiredInputs() const
  {
        typename std::vector<SmartPointer<DataObject> >::size_type num = 0;

          if (this->GetFixedImage())
            {
            num++;
            }

          if (this->GetMovingImage())
            {
            num++;
            }
          if (this->GetCertaintyImage())
            {
            num++;
            }
          return num;
  }

  /** Set/Get whether the deformation field is smoothed
   * (regularized). If SmoothDeformationField is on, then the
   * deformation field is smoothed with a Gaussian whose standard
   * deviations are specified with SetStandardDeviations() */
  itkSetMacro( SmoothDeformationField, bool );
  itkGetMacro( SmoothDeformationField, bool );
  itkBooleanMacro( SmoothDeformationField );

  /** Set the Gaussian smoothing standard deviations for the
   * deformation field. The values are set with respect to pixel
   * coordinates. */
  itkSetVectorMacro( StandardDeviations, double, ImageDimension );
  /**
   * void SetStandardDeviations( double value )
   * Set the standard deviations
   */
  virtual void SetStandardDeviations( double value );

  /** Get the Gaussian smoothing standard deviations use for smoothing
   * the deformation field. */
  const double * GetStandardDeviations(void)
    { return (double *) m_StandardDeviations; }

  /** Set/Get whether the update field is smoothed
   * (regularized). Smoothing the update field yields a solution
   * viscous in nature. If SmoothUpdateField is on, then the
   * update field is smoothed with a Gaussian whose standard
   * deviations are specified with SetUpdateFieldStandardDeviations() */
  itkSetMacro( SmoothUpdateField, bool );
  itkGetMacro( SmoothUpdateField, bool );
  itkBooleanMacro( SmoothUpdateField );

  /** Set the Gaussian smoothing standard deviations for the update
   * field. The values are set with respect to pixel coordinates. */
  itkSetVectorMacro( UpdateFieldStandardDeviations, double, ImageDimension );
  /**
   * void SetUpdateFieldStandardDeviations( double value )
   * Use to update the deformation field's standard deviation
   */
  virtual void SetUpdateFieldStandardDeviations( double value );

  /** Get the Gaussian smoothing standard deviations used for
   * smoothing the update field. */
  const double * GetUpdateFieldStandardDeviations(void)
    { return (double *) m_UpdateFieldStandardDeviations; }



  /** Stop the registration after the current iteration. */
  virtual void StopRegistration()
    { m_StopRegistrationFlag = true; }

  /** Set/Get the desired maximum error of the Guassian kernel approximate.
   * \sa GaussianOperator. */
  itkSetMacro( MaximumError, double );
  itkGetMacro( MaximumError, double );

  /** Set/Get the desired limits of the Gaussian kernel width.
   * \sa GaussianOperator. */
  itkSetMacro( MaximumKernelWidth, unsigned int );
  itkGetMacro( MaximumKernelWidth, unsigned int );


  /** Certainty image. */
  MovingImageType * m_CertaintyImage;

  /** Rescaled Prototype */
  MovingImageType * m_RescaledPrototype;


protected:
  /**
   * MorphonToolboxFiler()
   * Constructor
   */
  MorphonToolboxFiler();

  ~MorphonToolboxFiler() {}
  
  /**
   * void PrintSelf(std::ostream& os, Indent indent) const
   * Use to print local labrary's information
   */
  void PrintSelf(std::ostream& os, Indent indent) const;


  /**
   * bool Halt()
   * Supplies the halting criteria for this class of filters.  The
   * algorithm will stop after a user-specified number of iterations.
   */
  virtual bool Halt();

  /**
   * void CopyInputToOutput()
   * A simple method to copy the data from the input to the output.
   * If the input does not exist, a zero field is written to the output.
   */
  virtual void CopyInputToOutput();

  /** 
   * void InitializeIteration()
   * Initialize the state of filter and equation before each iteration.
   * Progress feeback is implemented as part of this method.
   */
  virtual void InitializeIteration();

  /** 
   * void SmoothDeformationField()
   * Utility to smooth the deformation field (represented in the Output)
   * using a Guassian operator. The amount of smoothing can be specified
   * by setting the StandardDeviations.
   */
  virtual void SmoothDeformationField();

   /**
   *    void SmoothCertaintyImage()
   *     Utility to smooth the certainty image
   *    using a Guassian operator. The amount of smoothing can be specified
   *    by setting the StandardDeviations.
   */
  void SmoothCertaintyImage();

  /**
   *   void LocalSmooth( itk::Image< float, ImageDimension > * ptr )
   *   This function is mostly used to smooth an unknown prior image
   */
  void LocalSmooth( itk::Image< float, ImageDimension > * ptr );

  /** This method is called after the solution has been generated. In this case,
   * the filter release the memory of the internal buffers. */
  virtual void PostProcessOutput()
  {
        this->Superclass::PostProcessOutput();
        m_TempField->Initialize();
  }

  /** 
   * virtual void Initialize()
   * This method is called before iterating the solution.
   */
  virtual void Initialize()
  {
      this->Superclass::Initialize();
      m_StopRegistrationFlag = false;
  }

  /**
   * virtual void GenerateOutputInformation();
   * By default the output deformation field has the same Spacing, Origin
   * and LargestPossibleRegion as the input/initial deformation field.  If
   * the initial deformation field is not set, the output information is
   * copied from the fixed image.
   */
  virtual void GenerateOutputInformation();

  /**
   * virtual void GenerateInputRequestedRegion()
   * It is difficult to compute in advance the input moving image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole moving image.
   *
   * For the fixed image and deformation field, the input requested region
   * set to be the same as that of the output requested region.
   */
  virtual void GenerateInputRequestedRegion();

private:
  MorphonToolboxFiler(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Standard deviation for Gaussian smoothing */
  double                   m_StandardDeviations[ImageDimension];
  double                   m_UpdateFieldStandardDeviations[ImageDimension];

  /** Modes to control smoothing of the update and deformation fields */
  bool m_SmoothDeformationField;
  bool m_SmoothUpdateField;


  /** Temporary deformation field use for smoothing the
   * the deformation field. */
  DeformationFieldPointer   m_TempField;

private:
  /** Maximum error for Gaussian operator approximation. */
  double                    m_MaximumError;

  /** Limits of Guassian kernel width. */
  unsigned int              m_MaximumKernelWidth;

  /** Flag to indicate user stop registration request. */
  bool                      m_StopRegistrationFlag;
  

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMorphonToolboxFiler.txx"
#endif


#endif

