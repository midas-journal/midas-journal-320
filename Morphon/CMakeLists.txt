cmake_minimum_required(VERSION 2.4)
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

# This project is designed to be built outside the Insight source tree.
PROJECT(itkMorophon)

# declare the needed package
SET(Required_Packages
ITK
)

# load the needed packages
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

# Set the libraries
#INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)
SET(Libraries
ITKIO
ITKBasicFilters
ITKNumerics
ITKCommon
ITKAlgorithms
ITKStatistics
)

SET(CurrentExe "itkMorphon")
ADD_EXECUTABLE(${CurrentExe} main.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

ADD_TEST(RegisterImage itkMorphon ${CMAKE_SOURCE_DIR}/hand_indata.png  ${CMAKE_SOURCE_DIR}/hand_prot.png ${CMAKE_SOURCE_DIR}/out.png 10 1.5 10 1.5 10 4.5 10 4.5 10 4.5 10 2 10 2 10 2 10 2 10 2 10 2 10 3.5 10 2.5)


