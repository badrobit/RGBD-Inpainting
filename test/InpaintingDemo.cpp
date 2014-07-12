#include <DepthDerivative.h>
#include <Helpers.hpp>

using namespace RGBDInpainting;

int main( int argc, char** argv )
{
    DepthDerivative depthDerivative;

    //RGBImage::Pointer image;
    ScalarImage::Pointer image;
    ReadImage( "data/test.png", image.GetPointer() );

    GradientImagePointer gradient;
    gradient = GradientImage::New();

    //depthDerivative.SetInputDepthMap( image.GetPointer() );
    //depthDerivative.Update();
    depthDerivative.ForwardDifferenceDerivative( image.GetPointer(), gradient.GetPointer() );

}