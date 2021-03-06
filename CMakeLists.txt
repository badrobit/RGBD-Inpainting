CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( RGBD-Inpainting )

#.....: CMAKE MODULES :...........................................................................:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake_modules")

#.....: Compiler Options :........................................................................:
ADD_DEFINITIONS( "-std=c++11 -fpermissive -fopenmp -w" )
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Modules )

#.....: Google Logging :..........................................................................:
find_package(GLog REQUIRED)
include_directories( ${GLOG_INCLUDE_DIRS} )
add_definitions(-DGOOGLE_GLOG_DLL_DECL=)
set( GOOGLE_LIBS ${GOOGLE_LIBS} libglog.so )

#.....: ITK Library :.............................................................................:
FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

#.....: Libraries :...............................................................................:
ADD_LIBRARY( Helpers Modules/Helpers/Helpers.hpp )
SET_TARGET_PROPERTIES( Helpers PROPERTIES LINKER_LANGUAGE CXX )
TARGET_LINK_LIBRARIES( Helpers ${ITK_LIBRARIES} ${GOOGLE_LIBS} )

ADD_LIBRARY( DepthMask Modules/DepthMask/DepthMask.cpp )
SET_TARGET_PROPERTIES( DepthMask PROPERTIES LINKER_LANGUAGE CXX )
TARGET_LINK_LIBRARIES( DepthMask ${ITK_LIBRARIES} Helpers )

ADD_LIBRARY( DepthDerivative Modules/DepthDerivative/DepthDerivative.cpp )
SET_TARGET_PROPERTIES( DepthDerivative PROPERTIES LINKER_LANGUAGE CXX )
TARGET_LINK_LIBRARIES( DepthDerivative Helpers ${ITK_LIBRARIES} ${GOOGLE_LIBS} )

#.....: Demo Files :..............................................................................:
ADD_EXECUTABLE( InpaintingDemo test/InpaintingDemo.cpp )
TARGET_LINK_LIBRARIES( InpaintingDemo DepthDerivative )