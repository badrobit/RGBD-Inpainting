CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( RGBD-Inpainting ) 

#.....: Compiler Options :........................................................................:
ADD_DEFINITIONS( -std=c++11 )

#.....: OpenCV :..................................................................................:
FIND_PACKAGE( OpenCV REQUIRED )
INCLUDE_DIRECTORIES( ${OpenCV_INVCLUDE_DIRS} )

#.....: Libraries :...............................................................................:

#.....: Demo Files :..............................................................................:
ADD_EXECUTABLE( InpaintingDemo test/InpaintingDemo.cpp )
TARGET_LINK_LIBRARIES( InpaintingDemo ${OpenCV_LIBS} )