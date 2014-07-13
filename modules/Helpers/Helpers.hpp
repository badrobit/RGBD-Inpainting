// SubModules
#include <Helpers/ITKImageTypes.h>

// Google Logging
#include <glog/logging.h>

// ITK Includes
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>


#ifndef HELPERS_HPP_
#define HELPERS_HPP_

using namespace google;

namespace RGBDInpainting
{

template< typename InputImageType, typename OutputImageType >
void DeepCopyImageRegions( const InputImageType* input_image, const itk::ImageRegion<2>& image_region, OutputImageType* output_image )
{
    DLOG(INFO) << "Starting Image Region Copy";
    // Create two iterator objects one for the source one for the target and set them up.
    itk::ImageRegionConstIterator<InputImageType> inputImageIterator( input_image, image_region );
    itk::ImageRegionIterator<OutputImageType> outputImageIterator( output_image, image_region );

    // Iterate through the images and copy them over.
    while( !inputImageIterator.IsAtEnd() )
    {
        outputImageIterator.Set( inputImageIterator.Get() );
        ++inputImageIterator;
        ++outputImageIterator;
    }
    DLOG(INFO) << "Image Region Copy Completed";
}

template< typename InputPixelType, typename OutputPixelType >
void DeepCopyImage( const itk::Image<InputPixelType, 2>* const input_image, itk::Image<OutputPixelType, 2>* const output_image )
{
    DLOG(INFO) << "Starting Image Deep Copy";
    // If the output image does not have enough memory or any memory allocated do so :)
    if( output_image->GetLargestPossibleRegion() != input_image->GetLargestPossibleRegion() )
    {
        output_image->SetRegions( input_image->GetLargestPossibleRegion() );
        output_image->Allocate();
    }

    DLOG(INFO) << "TEST";

    DeepCopyImageRegions( input_image, input_image->GetLargestPossibleRegion(), output_image );
    DLOG(INFO) << "Image Deep Copy Completed";
}

template< typename ImageType >
void ReadImage( const std::string& file_name, ImageType* output_image )
{
    DLOG(INFO) << "Starting Image File Read of: " << file_name;
    typedef itk::ImageFileReader<ImageType> ImageReader;
    typename ImageReader::Pointer imageReader = ImageReader::New();
    imageReader->SetFileName( file_name );
    try
    {
        imageReader->Update();
    }
    catch( ... )
    {
        LOG(FATAL) << "Image Reading Failed";
    }

    output_image = imageReader->GetOutput();
    DLOG(INFO) << "Reading in of Image: " << file_name << " is completed";
}

template< typename ImageType >
void SaveImage( const ImageType* const image_to_save, const std::string& file_name )
{
    try
    {
        typename itk::ImageFileWriter<ImageType>::Pointer ImageWriter = itk::ImageFileWriter<ImageType>::New();
        ImageWriter->SetFileName( file_name );
        ImageWriter->SetInput( image_to_save );
        ImageWriter->Update();
    }
    catch( itk::ExceptionObject& error )
    {
        LOG(ERROR) << "Image save has failed! Region was: " << image_to_save->GetLargestPossibleRegion() << "Original Error: " << error;
        //throw std::runtime_error( error );
    }
}


}

#endif /* HELPERS_HPP_ */