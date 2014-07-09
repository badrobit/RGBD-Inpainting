/*
 * DepthDerrivative.h
 *
 *  Created on: Jul 8, 2014
 *      Author: badrobit
 */

// ITK Includes
#include <itkImage.h>
#include <itkDerivativeImageFilter.h>
#include <itkForwardDifferenceOperator.h>
#include "itkNeighborhoodOperatorImageFilter.h"

#ifndef DEPTHDERIVATIVE_H_
#define DEPTHDERIVATIVE_H_

namespace RGBDInpainting
{

typedef itk::Image<unsigned char, 2> ScalarImage;
typedef ScalarImage::Pointer ScalarImagePointer;

typedef itk::Image<itk::CovariantVector<float, 2>, 2> GradientImage;
typedef GradientImage::Pointer GradientImagePointer;

class DepthDerivative
{
public:
	DepthDerivative();

	virtual ~DepthDerivative();

	void SetInputDepthMap( ScalarImagePointer input_depth_map );

private:
	void ForwardDifferenceDerivative( const GradientImagePointer gradient_image );

protected:
	ScalarImagePointer original_depth_map;
};

}
#endif /* DEPTHDERRIVATIVE_H_ */
