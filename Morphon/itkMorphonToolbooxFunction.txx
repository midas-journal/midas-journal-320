/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonToolbooxFunction.txx,v $
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

#ifndef _itkMorphonToolbooxFunction_txx
#define _itkMorphonToolbooxFunction_txx

#include "itkFiniteDifferenceFunction.h"

namespace itk {

  template<class TFixedImage, class TMovingImage, class TDeformationField>
  MorphonToolbooxFunction<TFixedImage, TMovingImage, TDeformationField>
  ::MorphonToolbooxFunction()
    {
      m_MovingImage = NULL;
      m_FixedImage = NULL;
      m_CertaintyImage = NULL;
      m_RescaledPrototype = NULL;
      m_DeformationField = NULL;
      m_Energy = 0.0;
      m_NormalizeGradient = true;
      m_GradientStep = 1.0;
    }
  
  template<class TFixedImage, class TMovingImage, class TDeformationField>
  void
  MorphonToolbooxFunction<TFixedImage, TMovingImage, TDeformationField>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "MovingImage: ";
    os << m_MovingImage.GetPointer() << std::endl;
    os << indent << "FixedImage: ";
    os << m_FixedImage.GetPointer() << std::endl;
    os<<indent<<"CertaintyImage";
    os<< m_CertaintyImage <<std::endl;
    os<<m_RescaledPrototype<<std::endl;
  }


} // end namespace itk

#endif
