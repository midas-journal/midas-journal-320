/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadFilter.h,v $
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

#ifndef __itkQuadFilter_h
#define __itkQuadFilter_h

#include "itkComplexNeighborhoodOperator.h"
#include "itkNeighborhoodOperator.h"

namespace itk{

  template<class ComplexType, unsigned int ImageDimension=2, class TAllocator = NeighborhoodAllocator<ComplexType> >
  class ITK_EXPORT QuadFilter: 
  public ComplexNeighborhoodOperator<ComplexType, ImageDimension, TAllocator>
  {
    public :
    typedef QuadFilter Self;
    typedef ComplexNeighborhoodOperator<ComplexType, ImageDimension, TAllocator> SuperClass;
    typedef SmartPointer<Self> Pointer;
    typedef std::vector<ComplexType> VectorCoeffType;
        	
    itkTypeMacro(local_class, itk::NeighborhoodOperator);
    
    
    Self &operator=(const Self& other)
    {
      SuperClass::operator=(other);
      return *this;
    }
        	
    /** constructor */
    ~QuadFilter(){}
    QuadFilter( const Self& other) : ComplexNeighborhoodOperator<ComplexType, ImageDimension, TAllocator>(other){}
    
    QuadFilter() {m_FilterNumber = -1; m_coeff=NULL;}
    QuadFilter(int FilterNumber, VectorCoeffType* local_coeff)
    {
      if( FilterNumber>(ImageDimension*2) )
      {
        std::cerr<<"Error in local_function constructor : only"<<(ImageDimension*2)<<" filters.";
        std::cerr<<"Please specify one of them"<<std::endl;
        throw(-1);
      }
      m_FilterNumber = FilterNumber;
      m_coeff = (local_coeff);
    }
        	
    /**
     * void CreateDirectionnal(int radius)
     * Creates the operator with length only in the specified direction.
     */
    void CreateDirectionnal(int radius) { this->CreateToRadius(radius); }
        	
    void SetVectorCoeff( VectorCoeffType* local_coeff ) { m_coeff = local_coeff ; }
        	
    void SetNumFilter ( int FilterNumber ) { m_FilterNumber = FilterNumber; }
    
    int GetNumFilter() { return m_FilterNumber;}
        	
    void CreateDirectional()
    {
     this->CreateToRadius(4);
    }
        	
    /**
     * CoefficientVectorType GenerateCoefficients() 
     * Calculates the operator coefficients.
     */
    virtual VectorCoeffType GenerateCoefficients() 
    {
      if( m_FilterNumber==(-1) | m_coeff==NULL)
      {
        std::cerr<<"Error in GenerateCoefficients() : need to specified the number of the filter"<<std::endl;
        std::cerr<<"FilterNumber = "<<m_FilterNumber<<std::endl;
        std::cerr<<"m_coeff = "<<m_coeff<<std::endl;
        throw(-1);
      }
      return (*m_coeff);
    }
              
    /**
     *    void Fill(const CoefficientVector &coeff)
     *    Arranges coefficients spatially in the memory buffer.
     */
     virtual void Fill( const VectorCoeffType &coeff )
     {
     
       this->InitializeToZero();
	// Note that this code is only good for 2d and 3d operators.  Places the
	// coefficients in the exact center of the neighborhood
    	unsigned int i;
    	int x,y,z, pos;
    	unsigned int center = this->GetCenterNeighborhoodIndex();

    	if (ImageDimension == 3)
	{
	  i = 0;
	  for (z = -4; z <= 4; z++)
	  {
	    for (y = -4; y <= 4; y++ )
	    {
              for (x = -4; x <= 4; x++)
	      {
		pos = center + z * this->GetStride(2) + y * this->GetStride(1) +
		    x * this->GetStride(0);
		this->operator[](pos) = static_cast<ComplexType>(coeff[i]);
		i++;
	      }
	    }
           }
	 }
	 if (ImageDimension == 2)
	 {
	   i = 0;
	   for (y = -4; y <= 4; y++ )
	   {
	     for (x = -4; x <= 4; x++)
	     {
	       pos = center + y * this->GetStride(1) + x * this->GetStride(0);
	       this->operator[](pos) = static_cast<ComplexType>(coeff[i]);
	       i++;
	    }
	 }
       }
     }//end Fill

   protected :
   private :
   int m_FilterNumber;
   VectorCoeffType* m_coeff;

};//end of the local_class	
} // namespace itk	    
#endif
	   
