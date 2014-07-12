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
	InitGoogleLogging( "DepthDerivative" );
	LOG(INFO) << "Google logging initialized for DepthDerivative class";
	FLAGS_stderrthreshold = 2;

#ifndef NDEBUG
FLAGS_stderrthreshold = 0;
DLOG(INFO) << "Debugging Enabled, Logs will also output to Terminal";
#endif

	original_depth_map = ScalarImage::New();
	m_gradient_image = GradientImage::New();

	DLOG(INFO) << "DephDerivative Initialized";
}

DepthDerivative::~DepthDerivative()
{
}

void
DepthDerivative::SetInputDepthMap( ScalarImage* input_depth_map )
{
	this->original_depth_map = input_depth_map;
	DLOG(INFO) << "Depth Map Set.";
}

void DepthDerivative::Update()
{
	DLOG(INFO) << "Computing Depth Derivatives";
	this->ForwardDifferenceDerivative( original_depth_map.GetPointer(), m_gradient_image.GetPointer() );
}

GradientImagePointer DepthDerivative::GetOutput()
{
	return m_gradient_image;
}

/**
 * http://www.itk.org/Wiki/ITK/Examples/Operators/ForwardDifferenceOperator
 * http://www.itk.org/Wiki/ITK/Examples/Images/NeighborhoodOperatorImageFilter
 */
void
DepthDerivative::ForwardDifferenceDerivative( const ScalarImage::Pointer input_depth_map, GradientImage* const gradient_image )
{
	DLOG(INFO) << "Setting up Forward Difference Derivative";
	// Setup the Derivative Operator, and set its radius to 1. This will give us a kernel that is
	// 3px x 3px.
	typedef itk::ForwardDifferenceOperator<unsigned char, 2> DerivativeOperator;
	itk::Size<2> derivative_kernel_radius;
	derivative_kernel_radius.Fill( 1 );

	DLOG(INFO) << "Creating X Derivative Operator";
	// Compute the derivative in the X direction.
	DerivativeOperator DerivativeX;
	DerivativeX.SetDirection( 0 );
	DerivativeX.CreateToRadius( derivative_kernel_radius );

	DLOG(INFO) << "Computing X Derivative";
	typedef itk::NeighborhoodOperatorImageFilter<ScalarImage, ScalarImage> NeighborhoodImageFilter;
	NeighborhoodImageFilter::Pointer DerivativeXFilter = NeighborhoodImageFilter::New();
	DerivativeXFilter->SetOperator( DerivativeX );
	DerivativeXFilter->SetInput( input_depth_map );
	DLOG(INFO) << "Primary is set: " << (bool)DerivativeXFilter->HasInput( "Primary" );
	DerivativeXFilter->Update();

#ifndef NDEBUG
	DLOG(INFO) << "Saving X Derivative to File";
	SaveImage( DerivativeXFilter->GetOutput(), "xDerivative.png" );
#endif

// 	DLOG(INFO) << "computing Y derivative";
// 	// Compute the derivative in the Y direction.
// 	DerivativeOperator DerivativeY;
// 	DerivativeY.SetDirection( 1 ) ;
// 	DerivativeY.CreateToRadius( derivative_kernel_radius );

// 	NeighborhoodImageFilter::Pointer DerivativeYFilter = NeighborhoodImageFilter::New();
// 	DerivativeYFilter->SetOperator( DerivativeY );
// 	DerivativeYFilter->SetInput( this->original_depth_map );
// 	DerivativeYFilter->Update();

// #ifndef NDEBUG
// SaveImage( DerivativeYFilter->GetOutput(), "yDerivative.png" );
// #endif

// 	// Combine our derivative images together to form a gradient image
// 	typedef itk::ComposeImageFilter<ScalarImage, GradientImage> ImageCompositionFilter;
// 	typename ImageCompositionFilter::Pointer imageComposer = ImageCompositionFilter::New();
// 	imageComposer->SetInput( 0, DerivativeXFilter->GetOutput() );
// 	imageComposer->SetInput( 1, DerivativeYFilter->GetOutput() );
// 	imageComposer->Update();

// #ifndef NDEBUG
// SaveImage( imageComposer->GetOutput(), "gradientImage.png" );
// #endif

// 	DeepCopyImage( imageComposer->GetOutput(), gradient_image );
}



} /* namespace RGBDInpainting */
