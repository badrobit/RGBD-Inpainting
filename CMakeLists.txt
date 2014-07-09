CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( RGBD-Inpainting ) 

#.....: Compiler Options :........................................................................:
#ADD_DEFINITIONS( -std=c++11 )
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )

#.....: ITK Library :.............................................................................:
FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

#.....: Libraries :...............................................................................:
ADD_LIBRARY( DepthDerivative src/DepthDerivative.cpp )
TARGET_LINK_LIBRARIES( DepthDerivative ${ITK_LIBRARIES} )

#.....: Demo Files :..............................................................................:
#ADD_EXECUTABLE( InpaintingDemo test/InpaintingDemo.cpp )
#TARGET_LINK_LIBRARIES( InpaintingDemo ${OpenCV_LIBS} )