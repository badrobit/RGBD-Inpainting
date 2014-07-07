CMAKE_MINIMUM_REQUIRED( VERSION 2.8 FATAL_ERROR )
PROJECT( RGB-D Inpainting ) 

#.....: Compiler Options :........................................................................:
ADD_DEFINITIONS( -std=c++11 )

#.....: OpenCV :..................................................................................:
FIND_PACKAGE( OpenCV REQUIRED )

#.....: Libraries :...............................................................................:
