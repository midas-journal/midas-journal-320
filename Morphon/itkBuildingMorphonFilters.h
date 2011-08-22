/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBuildingMorphonFilters.h,v $
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

#ifndef __itkBuildingMorphonFilters_h
#define __itkBuildingMorphonFilters_h

#include "itkImage.h"
#include "itkExceptionObject.h"
#include "itkNeighborhoodOperator.h"
#include "itkComplexNeighborhoodOperator.h"
#include "itkImageToImageFilter.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"

#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include "itkQuadFilter.h"

#include "itkShrinkImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"

namespace itk {

template< class TMovingImage, class TFixedImage >
class ITK_EXPORT BuildingMorphonFilters: public ImageToImageFilter< TMovingImage, TFixedImage >
{
    
  public:
  /** Standard typedefs */
  static const unsigned int ImageDimension = TMovingImage::ImageDimension;
  typedef BuildingMorphonFilters Self;
  typedef ImageToImageFilter< TMovingImage, TFixedImage >  Superclass;
  typedef SmartPointer<Self> 	 Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  typedef std::complex<double> 	 ComplexType;
        
  typedef itk::Image< ComplexType, ImageDimension> ComplexImageType;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  itkTypeMacro(myFiler, ImageToImageFilter);
        
  typedef itk::NeighborhoodOperator<ComplexType, ImageDimension> 	 NeighborhoodOperatorType;

  typedef std::vector<ComplexType> 			 VectorCoeffType;
        
  typedef itk::QuadFilter<ComplexType, ImageDimension >	 ComplexOperator;
  typedef typename ComplexOperator::Pointer 		 ComplexOperatorPointerType;
        
  
  /**
   * Assignment operator
   */
  Self &operator=(const Self& other)
  {
    Superclass::operator=(other);
    return *this;
  }

  /**
   * ComplexOperator* GetComplexOp( int ID )
   * Get a specific filter defined by ID
   * ID must be inferior to ImageDimension*2, otherwise, a error is return
   */
  ComplexOperator* GetComplexOp( int ID )
  {
    if( ID>=(ImageDimension*2) )
    {
      std::cerr<<"Error in the GetComplexOp( int ID ) function"<<std::endl;
      std::cerr<<"The filter's ID must be in 0 and "<<(ImageDimension*2)<<" and you specify "<<ID<<std::endl;
      throw(-1);
    }
    return (ComplexOp+ID);
  }

  /**
   *       int GetNumMatrix()
   *       used to get the numer of matrices thanks to the image's dimension
   **/
  int GetNumMatrix(){ return m_NumMatrix; }

  void PrintSelf(std::ostream &os, Indent i) const;
  /*{
    os << i << "Quadratic Filters operator { this=" << this  << "}" << std::endl;
    Superclass::PrintSelf(os, i.GetNextIndent());   }*/
	
	
  /**
   *  getChoicefor( unsigned int i=0;
   *  used to get the filter choice
   */
  const int GetChoice() { return Choice;}

  /**
   *  void SetChoice(int choice)
   *  used to set the filter choice
   */
  void SetChoice(int choice){ this->Choice = choice; }


  /**
   * int returnDirectionnalCoeff(int l, int c)
   * return the specific coefficinet 
   * l correspond to the line and must be <(ImageDimension*2)
   * c correspond to the column and must be <(ImageDiension)
   */
  double ReturnDirectionnalCoeff(int l, int c)
  {
    if( m_DirectionnalCoeff == NULL )
    {
      std::cerr<<"Error in the returnDirectionnalCoeff()\n Some problems append during the construction : no coefficients matrix"<< std::endl;
      throw(EXIT_FAILURE);  
    }
    if( l>=(ImageDimension*2) | c>=ImageDimension )
    {
      std::cerr<<"Error in the returnDirectionnalCoeff(int l, int c), we need to pass correct parameters l and c :\n l must be<(ImageDimension*2)\n  c must be <(ImageDiension)"<< std::endl;
      throw(EXIT_FAILURE);
    }
    else
      return m_DirectionnalCoeff[l][c];
  }
         
  /**
   * int returnDirectionnalCoeff()
   * return the matrix of directionnal coefficinets (wich has a size [ImageDimensin*2][ImageDimension] 
   */
  double** GetDirectionnalCoeff()
  {
    if( m_DirectionnalCoeff == NULL )
    {
      std::cerr<<"Error in the returnDirectionnalCoeff()\nSome problems append during the construction : no coefficients matrix"<< std::endl;
      throw(EXIT_FAILURE);  
    }
      return m_DirectionnalCoeff;
  }
         
  /**
   * int returnFiltersDimension()
   * return the dimension of the filters
   */
  int ReturnFIltersDimension(){ return FiltersSize;}
         
  /**
   * int returnNumberOfFilters()
   * return the number of filters
   */
  int ReturnNumberOfFilters(){ return m_NumMatrix; }
        
  /**
   * VectorCoeffType* GetSpecificCoeffVector(int n)
   * use to return a specific coefficient of a vector
   */
  VectorCoeffType* GetSpecificCoeffVector(int n){ return (coeff+n);}
        
  /**
   * VectorCoeffType* ReturnFilterCoeff(int n)
   * use to return the coefficient of a specific filter
   */
  VectorCoeffType* ReturnFilterCoeff(int n) {return (coeff+n);}  
        
        
         
  /* Protected function */
  protected:
  
  /** Constructor*/
  BuildingMorphonFilters();
  /* Destructor */
  ~BuildingMorphonFilters() { };

  /**
   *    CoefficientVector GenerateCoefficients()
   *    Calculates operator coefficients by reading two files.
   *	With no changes, the files are "2D_QuadPhaseFilter.csv" and "3D_QuadPhaseFilter.csv"
   *	defined in the class constructor.
   */
   void GenerateCoefficients();
          
          
 /**
  * void generateDirectionnalCoeff();
  * create the coefficient regard to the direction, 
  * usefull in the deformations computation
  */
  void generateDirectionnalCoeff();
          
 
  /* Private functions */
  private:
    /* Choice is used to select the matrix */
    int Choice;
    /* The name of the file witch contains the coefficient*/
    std::string CoeffFileName;    /* the number of phase filters, it's a function of ImageDimension*/
    int m_NumMatrix, FiltersSize, m_radius; 
    /* The coefficients of the phase filters
     * it's a matrix wich contains all of them instancied one time
     */
   VectorCoeffType* coeff; 
   double** m_DirectionnalCoeff;
   ComplexOperator* ComplexOp;

};//end Class


}//end namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBuildingMorphonFilters.txx"
#endif

#endif //end librairy
