/*
 * DepthDerrivative.h
 *
 *  Created on: Jul 8, 2014
 *      Author: badrobit
 */
#ifndef DEPTHDERIVATIVE_H_
#define DEPTHDERIVATIVE_H_

// ITK Includes
#include <itkComposeImageFilter.h>
#include <itkDerivativeImageFilter.h>
#include <itkForwardDifferenceOperator.h>
#include <itkNeighborhoodOperatorImageFilter.h>

// Google Logging
#include <glog/logging.h>

// Project Helpers
#include <Helpers/Helpers.hpp>

using namespace google;

namespace RGBDInpainting
{

class DepthDerivative
{
public:
	DepthDerivative();

	virtual ~DepthDerivative();

	void SetInputDepthMap( ScalarImage* input_depth_map );

    void Update();

    GradientImagePointer GetOutput();

    void ForwardDifferenceDerivative( const ScalarImage::Pointer input_depth_map, GradientImage* const gradient_image );

private:


protected:
	ScalarImagePointer original_depth_map;
    GradientImagePointer m_gradient_image;
};

}
#endif /* DEPTHDERRIVATIVE_H_ */
