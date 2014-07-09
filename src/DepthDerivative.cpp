/*
 * DepthDerrivative.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: badrobit
 */

#include "DepthDerivative.h"

namespace RGBDInpainting
{

DepthDerivative::DepthDerivative()
{
	// TODO Auto-generated constructor stub

}

DepthDerivative::~DepthDerivative()
{
	// TODO Auto-generated destructor stub
}

void
DepthDerivative::SetInputDepthMap( ScalarImagePointer input_depth_map )
{
	this->original_depth_map = input_depth_map;
}

/**
 *
 * http://www.itk.org/Wiki/ITK/Examples/Operators/ForwardDifferenceOperator
 * http://www.itk.org/Wiki/ITK/Examples/Images/NeighborhoodOperatorImageFilter
 */
void DepthDerivative::ForwardDifferenceDerivative( const GradientImagePointer gradient_image )
{
	// Setup the Derivative Operator, and set its radius to 1. This will give us a kernal that is
	// 3px x 3px.
	typedef itk::ForwardDifferenceOperator<float, 2> DerivativeOperator;
	itk::Size<2> derivative_kernel_radius;
	derivative_kernel_radius.Fill( 1 );

	// Compute the derivative in the X direction.
	DerivativeOperator DerivativeX;
	DerivativeX.SetDirection( 0 );
	DerivativeX.CreateToRadius( derivative_kernel_radius );

	typedef itk::NeighborhoodOperatorImageFilter<ScalarImage, ScalarImage, float> NeighborhoodImageFilter;
	typename NeighborhoodImageFilter::Pointer DerivativeXFilter = NeighborhoodImageFilter::New();
	DerivativeXFilter->SetOperator( DerivativeX );
	DerivativeXFilter->SetInput( original_depth_map );
	DerivativeXFilter->Update();
}


} /* namespace RGBDInpainting */
