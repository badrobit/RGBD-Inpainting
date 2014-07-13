#ifndef ITKImageTypes_H
#define ITKImageTypes_H

// ITK
#include <itkImage.h>
#include <itkRGBPixel.h>
#include <itkVectorImage.h>

namespace RGBDInpainting
{
  /** Scalar Image Types. */
  typedef itk::Image<unsigned char, 2> ScalarImage;
  typedef ScalarImage::Pointer ScalarImagePointer;

  /** RGB Image Types. */
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImage;
  typedef RGBImage::Pointer RGBImagePointer;

  /** Gradient Image Types */
  typedef itk::Image<itk::CovariantVector<float, 2>, 2> GradientImage;
  typedef GradientImage::Pointer GradientImagePointer;
}
#endif