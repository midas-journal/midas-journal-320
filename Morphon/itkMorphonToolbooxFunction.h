/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonToolbooxFunction.h,v $
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

#ifndef _itkMorphonToolbooxFunction_h_
#define _itkMorphonToolbooxFunction_h_

#include "itkFiniteDifferenceFunction.h"

namespace itk {

/** \class MorphonToolbooxFunction
 *
 * This is an abstract base class for all PDE functions which drives a
 * deformable registration algorithm. It is used by
 * PDEDeformationRegistrationFilter subclasses to compute the
 * output deformation field which will map a moving image onto
 * a fixed image.
 *
 * This class is used to stock a rescaled, but not defomed, moving image,
 * this particular moving image type is used to produce the deformation field.
 *
 * This class is templated over the fixed image type, moving image type, 
 * certainty image type and the deformation field type.
 *
 * \sa MorphonToolbooxFunction
 * \ingroup FiniteDifferenceFunctions
 */

template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT MorphonToolbooxFunction :
  public FiniteDifferenceFunction<TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef MorphonToolbooxFunction    Self;
  typedef FiniteDifferenceFunction<TDeformationField>    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MorphonToolbooxFunction,
    FiniteDifferenceFunction );

  /** MovingImage image type. */
  typedef TMovingImage   MovingImageType;
  typedef typename MovingImageType::ConstPointer  MovingImagePointer;

  /** FixedImage image type. */
  typedef TFixedImage    FixedImageType;
  typedef typename FixedImageType::ConstPointer  FixedImagePointer;

  /** Deformation field type. */
  typedef TDeformationField    DeformationFieldType;
  typedef typename DeformationFieldType::Pointer  DeformationFieldTypePointer;

  /** Set the moving image.  */
  void SetMovingImage( const MovingImageType * ptr )
    { m_MovingImage = ptr; }

    /** Set the certainty image */
    void SetCertaintyImage( MovingImageType * ptr )
    { m_CertaintyImage = ptr; }

    /** Set the certainty image. */
    void SetRescaledPrototype( MovingImageType * ptr )
    {
        m_RescaledPrototype = ptr;
    }

  /** Get the moving image. */
  const MovingImageType * GetMovingImage(void) const
    { return m_MovingImage; }

    /** Get the certainty image. */
  MovingImageType * GetCertaintyImage(void) const
      { return m_CertaintyImage; }

    /** Get the rescaled prototye. */
    MovingImageType * GetRescaledPrototype(void) const
    {
        return m_RescaledPrototype ;
    }

  /** Set the fixed image. */
  void SetFixedImage( const FixedImageType * ptr )
    { m_FixedImage = ptr; }

  /** Get the fixed image. */
  const FixedImageType * GetFixedImage(void) const
    { return m_FixedImage; }

  /** Set the defromation field. */
  void SetDeformationField(  DeformationFieldTypePointer ptr )
    { m_DeformationField = ptr; }

  /** Get the Deformation field. */
  DeformationFieldTypePointer GetDeformationField(void)
    { return m_DeformationField; }


  void SetEnergy( double e) { m_Energy=e;}
  double GetEnergy( ) const { return m_Energy;}
  void SetGradientStep( double e) { m_GradientStep = e;}
  double GetGradientStep( ) const { return m_GradientStep ;}
  void SetNormalizeGradient( bool e) { m_NormalizeGradient=e;}
  bool GetNormalizeGradient( ) const { return m_NormalizeGradient;}

  protected:
    /** Default Constructor */
  MorphonToolbooxFunction();
/** Default Destructor */
  ~MorphonToolbooxFunction() {}

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** The moving image. */
  MovingImagePointer                m_MovingImage;

  /** The fixed image. */
  FixedImagePointer                   m_FixedImage;

  /** The deformation field. */
  DeformationFieldTypePointer                   m_DeformationField;

  /** The certainly image. */
  MovingImageType *  m_CertaintyImage;

  /** Rescaled Prototype */
  MovingImageType * m_RescaledPrototype;

  mutable double                          m_Energy;
  bool                                          m_NormalizeGradient;
  mutable double                          m_GradientStep;

private:
  MorphonToolbooxFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMorphonToolbooxFunction.txx"
#endif


#endif
