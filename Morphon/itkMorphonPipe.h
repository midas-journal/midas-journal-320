/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMorphonPipe.h,v $
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

#ifndef __itkMorphonPipe_h
#define __itkMorphonPipe_h

#include "complex"

#include "itkCommand.h"
#include "itkImageToImageFilter.h"
#include "itkBuildingMorphonFilters.h"

#include "itkMorphonRegistrationFilter.h"
#include "itkCastImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkVectorResampleImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "itkWarpImageFilter.h"
#include "itkImageFileWriter.h"


namespace itk
{
/**
 * =========================================================
 *       This library create the morphon pipeline.
 * 
 *       The algorithm has many parameters: 
 *         - The inputs: fixed and moving images
 *         - The ouput's name: the name of the output image
 *         - A list of [number\_of\_iterations, variances] couple which defines
 *               the number of levels. At least, one couple is necessary to work
 *               correctly.
 *
 *       For each level the algorithm
 *           - makes a smothing of the certainty matrix and the deformation field.
 *           - rezises all the images, matrices and filed to the level's size (i.e. size = 2^[(level_num-1)/2] ).
 *           - makes a rescaled prototype aims to be deformed.
 *           - sets the variance and the number of iteration for this level.
 *           
 *      For each iterations the algorithm
 *			- convolues the fixed and prototype images with all the filtes.
 *			- makes a complex multiplication between all of them.
 *			- computes the updates : the displacements and the certainties.
 *			- smoothes the certainty matrix, the deformation field and the deformed prototype.
 *
 *       References :
 *       [1] Andreas Wrangsjö and Johanna Pettersson and Hans Knutsson,
 *           "Non-Rigid Registration Using Morphons",
 *       [2] Hans Knutsson and Mats Andersson,
 *           "Morphons: segmentation using elastic canvas and paint on priors",
 *
 *
 *       NB1 : the moving and fixed image must have the same dimensions and
 *		number of dimensions
 *
 *       NB2 : the functions
 *                   - SetMovingImage
 *                   - SetFixedImage
 *                   - SetLevelNumber
 *                   - SetIterationsNumber
 *                   - SetVariance
 *                   - SetOutputName
 *               have to be called BEFORE the the Update() function call
 *       
 *  \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 * =========================================================
 */
template <class TFixedImage, class TMovingImage, class TOutputImage>
class ITK_EXPORT MorphonPipe :
public ImageToImageFilter< TMovingImage, TOutputImage >
{
  public:
  /** Standard "Self" & Superclass typedef. */
  typedef MorphonPipe Self;
  typedef ImageToImageFilter< TMovingImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MorphonsFilter, ImageToImageFilter);

  /** Extract some information from the image types. */
  typedef typename TOutputImage::PixelType            OutputPixelType;
  typedef std::complex<float>			  ComplexPixelType;
  typedef typename TOutputImage::InternalPixelType    OutputInternalPixelType;
  typedef typename  TFixedImage::PixelType            FixedPixelType;
  typedef typename  TFixedImage::SizeType             SizeType;
  typedef typename  TMovingImage::PixelType           MovingPixellType;
  typedef typename  TMovingImage::InternalPixelType   InputInternalPixelType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. 
   */
   itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension);

  /** Image typedef support. */
  typedef TFixedImage    FixedImageType;
  typedef TMovingImage  MovingImageType;
  typedef TOutputImage  OutputImageType;


  /** Typedef **/
  typedef itk::Vector< float, ImageDimension >            VectorPixelType;
  typedef itk::Image < VectorPixelType, ImageDimension >  DeformationFieldType;

  typedef itk::RecursiveGaussianImageFilter < MovingImageType, MovingImageType >        DiscretGaussianType;
    
  typedef itk::VectorLinearInterpolateImageFunction <DeformationFieldType, double>     InterpolatorType;
  typedef itk::CastImageFilter<MovingImageType,OutputImageType >                       CastFilterType;

  typedef itk::ResampleImageFilter<MovingImageType, MovingImageType >                 ResampleFilterType;
  typedef itk::IdentityTransform< double, ImageDimension >                            TransformType;
  typedef itk::NearestNeighborInterpolateImageFunction<MovingImageType, double >      ImageInterpolatorType;
  typedef itk::VectorResampleImageFilter<DeformationFieldType, DeformationFieldType > ResampleFilterDeformationType;

  /** The writer */
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  
  /** Declare toolbok used to build quadrature filters*/  
  typedef itk::BuildingMorphonFilters< MovingImageType, FixedImageType > FilterType;
  /** The Morphon's Filter */
  typedef itk::MorphonRegistrationFilter <
                                             MovingImageType,
                                             MovingImageType,
                                             DeformationFieldType
                                           > RegistrationFilterType;

    
  /** iterators */  
  typedef itk::ConstNeighborhoodIterator< MovingImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator < MovingImageType >      IteratorType;

  /**
   *   void SetMovingImage ( MovingImageType* MovingImage )
   *   Set the moving image
   */
  void SetMovingImage ( MovingImageType* MovingImage ) { this->SetNthInput ( 0, MovingImage ); }
     

  /**
   *   void SetFixedImage ( FixedImageType* FixedImage )
   *   Set the fixed image
   */
      void SetFixedImage ( FixedImageType* FixedImage ) { this->SetInput ( 1, FixedImage ); }

  /**
   *   void SetLevelNumber( int nlevel)
   *   Set the number of level in the morphons algorithm
   */
  void SetLevelNumber( int nlevel){ this->m_nlevel = nlevel; }

  /**
   *   void SetIterationsNumber( int *noit )
   *   Set the number of iteration for each level.
   *   Warning : The dimension of the tab has to be equal to 'm_level'
   */
  void SetIterationsNumber( unsigned int *noit ){this->m_noit = noit;}

  /**
   *   void StandardDeviations(int *nofdev)
   *   Set the value of the standard deviation for each level in pixel space.
   *   Warning : The dimension of the tab has to be equal to 'm_level'
   */
  void SetStandardDeviation( float *nofdev){ this->m_nofdev = nofdev;}

  /**
   *   void SetOutputName( std::string & OutputName )
   *   Set the output name of the file to write the output image.
   */
  void SetOutputName( std::string & OutputName ){ m_OutputName = OutputName;}
    
  /*
   * void SetPadding( bool padding )
   *  enable (if padding==true) or disable the padding during the Morphon algorithm.
   *  Default value is false.
   *  Warning : this function is yet used for the moment.
   */
  void SetPadding( bool padding ) { m_padding = padding; }
    
  /*
   * bool GetPadding()
   *  return the padding value
   */
  bool GetPadding() { return m_padding; }
    

  protected :
  ~MorphonPipe() { }
    
  /* Default constructor */
  MorphonPipe()
  {
    /* the default padding option is false */
    m_padding = false;
    /** initialise the filters */
    std::cout<<"prepare to initiate the quadrature filters"<<std::endl;
    m_PhaseFilterOperator = FilterType::New();
	
    std::cout<<"Initialisation done"<<std::endl;

  }
      
      
  void PrintSelf(std::ostream& os, Indent indent)
  {
    Superclass::PrintSelf(os,indent);
  }
      
  /* computes the morphon's algorithm */
  void GenerateData();

  private :
  int m_nlevel;
  unsigned int * m_noit;
  float *m_nofdev;
  bool m_padding;
  std::string m_OutputName;
  typename FilterType::Pointer m_PhaseFilterOperator; //the coefficients of the phase filters
  MorphonPipe(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
    
  /** 
   * Variables use in the main pipe 
   */
  // the resampling factor
  static const float resampling = 2.0;
  // Smooth operations
  typename DiscretGaussianType::Pointer fixedXImageSmooth;
  typename DiscretGaussianType::Pointer fixedYImageSmooth;
  typename DiscretGaussianType::Pointer fixedZImageSmooth;
  
  typename DiscretGaussianType::Pointer movingXImageSmooth;
  typename DiscretGaussianType::Pointer movingYImageSmooth;
  typename DiscretGaussianType::Pointer movingZImageSmooth;
    
  typename DiscretGaussianType::Pointer certaintyXImageSmooth;
  typename DiscretGaussianType::Pointer certaintyYImageSmooth;
  typename DiscretGaussianType::Pointer certaintyZImageSmooth;

  /** Images, matrices and deformations field */
  /** the moving image */
  typename MovingImageType::Pointer MovingImage;

  /** the certainty image*/
  typename MovingImageType::Pointer CertaintyImage;

  /** the deformation field */
  typename DeformationFieldType::Pointer deformationField;
    
  /** resamplers and transformations*/
  typename ResampleFilterType::Pointer fixedResampler;
  typename ResampleFilterType::Pointer movingResampler;
  typename TransformType::Pointer idTransform;
    
  typename ImageInterpolatorType::Pointer interpolator;
    
  typename ResampleFilterType::Pointer certaintyResampler;
  typename ImageInterpolatorType::Pointer certaintyInterpolator;
    
  /** Filters */
  typename RegistrationFilterType::Pointer filter;

  class CommandIterationUpdate : public itk::Command
  {
    public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef  itk::SmartPointer<CommandIterationUpdate>  Pointer;
    itkNewMacro( CommandIterationUpdate );
    protected:
    CommandIterationUpdate() {};
    typedef itk::Image< float, ImageDimension >              InternalImageType;
    typedef itk::Vector< float, ImageDimension >             VectorPixelType;
    typedef itk::Image<  VectorPixelType, ImageDimension >   DeformationFieldType;

    typedef itk::MorphonRegistrationFilter<
                                    InternalImageType,
                                    InternalImageType,
                                    DeformationFieldType>   RegistrationFilterType;
    public:
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      const RegistrationFilterType * filter =
         dynamic_cast< const RegistrationFilterType * >( object );
      if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
        return;
      }
      std::cout << filter->GetMetric() << std::endl;
    }
  };//end command class
      
 typename CommandIterationUpdate::Pointer observer;

};//end MorphonPipe class
} //end itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMorphonPipe.txx"
#endif

#endif
