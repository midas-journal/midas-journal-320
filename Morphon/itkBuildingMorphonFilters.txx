/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBuildingMorphonFilters.txx,v $
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

#ifndef __itkBuildingMorphonFilters_txx
#define __itkBuildingMorphonFilters_txx

#include "itkBuildingMorphonFilters.h"

namespace itk {

  /**
   * Constructor
   */
  template< class TMovingImage, class TFixedImage >
  BuildingMorphonFilters< TMovingImage, TFixedImage >
  ::BuildingMorphonFilters()
  {
    m_NumMatrix = ImageDimension*2;
    FiltersSize = 9;
    this->coeff = new VectorCoeffType[ImageDimension*2];
    /* Reserve the place for the coeffecients of each filter */
    for( unsigned int i=0; i<(ImageDimension*2); i++ )
    {
      (*(coeff+i)).reserve( vcl_pow(static_cast<double>(FiltersSize),static_cast<double>(ImageDimension)) );
    }
    
    /* the files' names which contain the quadrature filters' coefficients */
    if(ImageDimension==2) {CoeffFileName = "2D_QuadPhaseFilter.csv";}
    else {CoeffFileName = "3D_QuadPhaseFilter.csv";}
    
    /* Read the files to generate the coefficients */
    this->GenerateCoefficients(); 
           
    m_DirectionnalCoeff = new double*[ImageDimension*2];
    for( unsigned int i=0; i<(ImageDimension*2); i++ )
    {  
      m_DirectionnalCoeff[i] = new double[ImageDimension]; 
    }
    
    /* Generate the directionnal coefficients */
    generateDirectionnalCoeff();
           
    /**
     * Creates the QuadFilter, inherites to Neighboorhod class 
     * We need one for each filter beacause we need to convolute each image with each filter
     */           
    m_radius = static_cast<int>(FiltersSize/2);
           
    ComplexOp = new ComplexOperator[ImageDimension*2];
    for( unsigned int i=0; i<(ImageDimension*2); i++ )
    {
      (ComplexOp[i]).SetNumFilter(i);
      ((ComplexOp)[i]).SetVectorCoeff( coeff+i );
      (ComplexOp[i]).SetDirection(0);
      (ComplexOp[i]).CreateDirectional();
     }
   }
     
  /**
   * void generateDirectionnalCoeff();
   * create the coefficient regard to the direction, 
   * usefull in the deformations computation
   */
  template< class TMovingImage, class TFixedImage >
  void
  BuildingMorphonFilters< TMovingImage, TFixedImage >
  ::generateDirectionnalCoeff()
  {
    if( ImageDimension == 3 )
    {
      m_DirectionnalCoeff[0][0] = 0.0;       m_DirectionnalCoeff[0][1] = 0.5257;   m_DirectionnalCoeff[0][2] = 0.8507;
      m_DirectionnalCoeff[1][0] = 0.0;       m_DirectionnalCoeff[1][1] = -0.5257;  m_DirectionnalCoeff[1][2] = 0.8507;
      m_DirectionnalCoeff[2][0] = 0.5257;    m_DirectionnalCoeff[2][1] = 0.8507;   m_DirectionnalCoeff[2][2] = 0.0;
      m_DirectionnalCoeff[3][0] = -0.5257;   m_DirectionnalCoeff[3][1] = 0.8507;   m_DirectionnalCoeff[3][2] = 0.0;
      m_DirectionnalCoeff[4][0] = 0.8507;    m_DirectionnalCoeff[4][1] = 0.0;      m_DirectionnalCoeff[4][2] = 0.5257;
      m_DirectionnalCoeff[5][0] = 0.8507;    m_DirectionnalCoeff[5][1] = 0.0;      m_DirectionnalCoeff[5][2] = -0.5257;
    }
    if(ImageDimension == 2 )
    {
      m_DirectionnalCoeff[0][0] = 0.0;       m_DirectionnalCoeff[0][1] = 1.0; 
      m_DirectionnalCoeff[1][0] = 0.70711 ;  m_DirectionnalCoeff[1][1] = 0.70711 ;
      m_DirectionnalCoeff[2][0] = 1.0 ;        m_DirectionnalCoeff[2][1] =  6.1232e-017;
      m_DirectionnalCoeff[3][0] = 0.70711 ;  m_DirectionnalCoeff[3][1] = -0.70711 ; 
    }
  }//end generateDirectionnalCoeff
  
  /**
   *    CoefficientVector GenerateCoefficients()
   *    Calculates operator coefficients by reading two files.
   *	With no changes, the files are "2D_QuadPhaseFilter.csv" and "3D_QuadPhaseFilter.csv"
   *	defined in the class constructor.
   */
  template< class TMovingImage, class TFixedImage >
  void
  BuildingMorphonFilters< TMovingImage, TFixedImage >
  ::GenerateCoefficients()
  {
    std::string filter,       //will contain a filter as a string
               local_coeff,   //will contain a coeff as a string
               local_real, local_complex;
    std::istringstream MyStr;
    int local_filters_size=0, local_matrixNum=0;
    std::ifstream Coefffile(CoeffFileName.c_str(), std::ifstream::out);            
    if ( !Coefffile ) // test if the file is correctly open
    {
      std::cerr<<"Error in the itkBuildingMorphonFilters class"<<std::endl;
      std::cerr<<"Error while trainnig to open the file which contains filters' coefficients"<<std::endl;
      throw(-1);
    }
    	    
    while ( std::getline( Coefffile, filter ) )
    {
      // each line in the csv file corresponds to ( FiltersSize*m_NumfiltersSize ) coefficients.
      // all of them are the FiltersSize first elements of each filters
      local_filters_size = 0;
      local_matrixNum = 0;
              
      std::istringstream MyStr(filter, std::istringstream::in);
      while ( std::getline( MyStr, local_coeff, ',' ) )
      { //warning : must have an exponential in the number so the more general is +/-a*e-X +/-b*e-Y
	        
        size_t found=0;
  	bool foundEndReal=false;
  	while (found!=std::string::npos)
  	{
  	  if( local_coeff[found] == '+' ){ foundEndReal=true; found+=1; }
  	  if( local_coeff[found] == '-'  & local_coeff[found-1] != 'e' & found>=1  ) foundEndReal=true;
  	  if( foundEndReal == false ) local_real.push_back( local_coeff[found] );
  	  if( foundEndReal == true &  local_coeff[found] != 'i' ) local_complex.push_back( local_coeff[found] );
  	  if( local_coeff[found] == 'i' ) break;
  	  found += 1;
  	}
		
        ComplexType a ( atof(local_real.c_str()) , atof(local_complex.c_str()) );
        if( local_filters_size >= ( vcl_pow(static_cast<double>(FiltersSize),static_cast<double>(ImageDimension-1)) ) )
        {
          local_matrixNum += 1;
          local_filters_size = 0;
          if ( local_matrixNum >= m_NumMatrix )
          std::cout<<"ERROR : find more filers than declared !! "<<"\nDeclare : "<< m_NumMatrix << " found "<< local_matrixNum <<std::endl;
        }
        local_filters_size += 1;
        (*(coeff+local_matrixNum)).push_back( a );
        local_real = "";
        local_complex = "";
      }
    }
    std::cout<<"file read without any error"<<std::endl;
    for( unsigned int i=0; i<(ImageDimension*2); i++)
    std::cout<<(i+1)<<"th filter loaded with "<<(*(coeff+i)).size()<<" elements\n"	;
    std::cout<<std::endl;
  }//end GenerateCoefficients
  
  /**
   *  Printself
   */
  template< class TMovingImage, class TFixedImage >
  void
  BuildingMorphonFilters< TMovingImage, TFixedImage >
  ::PrintSelf(std::ostream &os, Indent i) const
  {
    os << i << "Quadrature filters constructor { this=" << this  << "}" << std::endl;
    os << i << "Number of filters : " << m_NumMatrix << std::endl;
    os << i << "Sizes (in all dimensions) : " << FiltersSize << std::endl;
    os << i << "Radius : " << m_radius << std::endl;
    Superclass::PrintSelf(os, i.GetNextIndent());   }
  
} // end namespace itk

#endif

