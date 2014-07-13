#ifndef HELPERS_HPP_
#define HELPERS_HPP_

// SubModules
#include <Helpers/ITKImageTypes.h>

// Google Logging
#include <glog/logging.h>

// ITK Includes
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

#include <itkPasteImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkRegionOfInterestImageFilter.h>

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

// ITKHELPERS

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


unsigned int GetNumberOfComponentsPerPixelInFile(const std::string& filename)
{
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO( filename.c_str(), itk::ImageIOFactory::ReadMode );
  return imageIO->GetNumberOfComponents();
}

std::string GetFileExtension(const std::string& fileName)
{
  return fileName.substr(fileName.find_last_of(".") + 1);
}

std::string GetPath(const std::string& fileName)
{
  return fileName.substr( 0, fileName.find_last_of('/') + 1);
}

itk::ImageRegion<2> GetRegionInRadiusAroundPixel(const itk::Index<2>& pixel, const unsigned int radius)
{
  // This function returns a Region with the specified 'radius' centered at 'pixel'. By the definition of the radius of a square patch, the output region is (radius*2 + 1)x(radius*2 + 1).
  // Note: This region is not necessarily entirely inside the image!

  // The "index" is the lower left corner, so we need to subtract the radius from the center to obtain it
  itk::Index<2> lowerLeft;
  lowerLeft[0] = pixel[0] - radius;
  lowerLeft[1] = pixel[1] - radius;

  itk::ImageRegion<2> region;
  region.SetIndex(lowerLeft);
  itk::Size<2> size;
  size[0] = radius*2 + 1;
  size[1] = radius*2 + 1;
  region.SetSize(size);

  return region;
}

template<typename TImage>
bool HasNeighborWithValueOtherThan(const itk::Index<2>& pixel, const TImage* const image,
                                   const typename TImage::PixelType& value)
{
  std::vector<itk::Offset<2> > offsets = Get8NeighborOffsets();

  for(auto offset : offsets)
  {
    if(image->GetPixel(pixel + offset) != value)
    {
      return true;
    }
  }

  return false;
}

template<typename TImage>
std::vector<itk::Index<2> > GetPixelsWithValueInRegion(const TImage* const image,
                                                       itk::ImageRegion<2> region,
                                                       const typename TImage::PixelType& value)
{
  region.Crop(image->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > pixelsWithValue;

  itk::ImageRegionConstIterator<TImage> regionIterator(image, region);
  while(!regionIterator.IsAtEnd())
  {
    if(regionIterator.Get() == value)
    {
      pixelsWithValue.push_back(regionIterator.GetIndex());
    }
    ++regionIterator;
  }

  return pixelsWithValue;
}

std::vector<itk::Offset<2> > IndicesToOffsets(const std::vector<itk::Index<2> >& indices, const itk::Index<2>& referenceIndex)
{
  std::vector<itk::Offset<2> > offsets;
  for(auto index : indices)
  {
    offsets.push_back(index - referenceIndex);
  }
  return offsets;
}

template<typename TImage>
bool AllPixelsEqualTo(const TImage* const image, const itk::ImageRegion<2>& region,
                      const typename TImage::PixelType& value)
{
  itk::ImageRegionConstIteratorWithIndex<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
  {
    if(imageIterator.Get() != value)
    {
      return false;
    }

    ++imageIterator;
  }

  // If none of the pixels were not equal to 'value', then all of them are equal to 'value'
  return true;
}

/** Get a list of the neighbors of a 'pixel' with a specified 'value'.*/
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsWithValue(const itk::Index<2>& pixel, const TImage* const image,
                                                       const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > neighborsWithValue;

  // Construct a 1x1 region (a single pixel)
  typename TImage::SizeType regionSize = {{1,1}};

  typename TImage::RegionType region(pixel, regionSize);

  // Construct a 1x1 radius (to make a 3x3 patch, or the 8-neighborhood)
  typename TImage::SizeType radius = {{1,1}};

  itk::ConstNeighborhoodIterator<TImage> neighborhoodIterator(radius, image, region);

  typename TImage::SizeType neighborhoodSize = neighborhoodIterator.GetSize();

  unsigned int numberOfPixelsInNeighborhood = neighborhoodSize[0] * neighborhoodSize[1];

  // Do not need to wrap this in a while(!iterator.IsAtEnd()) because we only want to visit the single pixel in 'region'
  for(unsigned int i = 0; i < numberOfPixelsInNeighborhood; i++)
  {
    if(i != neighborhoodIterator.GetCenterNeighborhoodIndex()) // Skip the center, as it is not a neighbor, but the query pixel itself
    {
      bool inBounds = false;
      if(neighborhoodIterator.GetPixel(i, inBounds) == value)
      {
        if(inBounds)
        {
          neighborsWithValue.push_back(neighborhoodIterator.GetIndex(i));
        }
      }
    }
  }

  return neighborsWithValue;
}

std::vector<itk::Offset<2> > Get8NeighborOffsets()
{
  std::vector<itk::Offset<2> > offsets;

  for(int i = -1; i <= 1; ++i)
  {
    for(int j = -1; j <= 1; ++j)
    {
      if(i == 0 && j == 0)
      {
        continue;
      }
      itk::Offset<2> offset;
      offset[0] = i;
      offset[1] = j;
      offsets.push_back(offset);
    }
  }
  return offsets;
}

std::vector<itk::Index<2> > Get8NeighborsInRegion(const itk::ImageRegion<2>& region, const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighborsInRegion;

  std::vector<itk::Offset<2> > neighborOffsets = Get8NeighborOffsets();
  for(auto offset : neighborOffsets)
  {
    itk::Index<2> index = pixel + offset;
    if(region.IsInside(index))
    {
      neighborsInRegion.push_back(index);
    }
  }
  return neighborsInRegion;
}

/** Get a list of the neighbors of a 'pixel' with a specified 'value' that are inside 'region'.
  * If 'region' is image->LargestPossibleRegion(), Get8NeighborsWithValue() will do this check automatically.
  */
// This technique is 2x slower than the one below (which uses a NeighborhoodIterator)
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsInRegionWithValue(const itk::Index<2>& pixel, const TImage* const image,
                                                           const itk::ImageRegion<2>& region,
                                                           const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > neighbors = Get8NeighborsInRegion(region, pixel);
  std::vector<itk::Index<2> > neighborsWithValue;
  for(auto neighbor : neighbors)
  {
    if(image->GetPixel(neighbor) == value)
    {
      neighborsWithValue.push_back(neighbor);
    }
  }
  return neighborsWithValue;
}

template <class TImage>
void CopyPatchIntoImage(const TImage* patch, TImage* const image, const itk::Index<2>& centerPixel)
{
  // This function copies 'patch' into 'image' centered at 'position'.

  // The PasteFilter expects the lower left corner of the destination position, but we have passed the center pixel.
  itk::Index<2> cornerPixel;
  cornerPixel[0] = centerPixel[0] - patch->GetLargestPossibleRegion().GetSize()[0]/2;
  cornerPixel[1] = centerPixel[1] - patch->GetLargestPossibleRegion().GetSize()[1]/2;

  typedef itk::PasteImageFilter <TImage, TImage> PasteImageFilterType;

  typename PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
  pasteFilter->SetInput(0, image);
  pasteFilter->SetInput(1, patch);
  pasteFilter->SetSourceRegion(patch->GetLargestPossibleRegion());
  pasteFilter->SetDestinationIndex(cornerPixel);
  pasteFilter->InPlaceOn();
  pasteFilter->Update();

  image->Graft(pasteFilter->GetOutput());

}

template <class TImage>
void CopyRegion(const TImage* sourceImage, TImage* targetImage,
               const itk::Index<2>& sourcePosition, const itk::Index<2>& targetPosition, const unsigned int radius)
{
  // Copy a patch of radius 'radius' centered at 'sourcePosition' from 'sourceImage' to 'targetImage' centered at 'targetPosition'
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(GetRegionInRadiusAroundPixel(sourcePosition, radius));
  extractFilter->SetInput(sourceImage);
  extractFilter->Update();

  CopyPatchIntoImage<TImage>(extractFilter->GetOutput(), targetImage, targetPosition);
}

template <class TImage>
void CopyRegion(const TImage* sourceImage, TImage* targetImage,
                const itk::ImageRegion<2>& sourceRegion, const itk::ImageRegion<2>& targetRegion)
{
  if(targetRegion.GetSize() != sourceRegion.GetSize())
  {
    std::cerr << "Can't copy regions that aren't the same size!" << std::endl
              << "Target region size " << targetRegion.GetSize() << std::endl
              << "Source region size " << sourceRegion.GetSize() << std::endl;
    return;
  }

  itk::ImageRegionConstIterator<TImage> sourceIterator(sourceImage, sourceRegion);
  itk::ImageRegionIterator<TImage> targetIterator(targetImage, targetRegion);

  while(!sourceIterator.IsAtEnd())
  {
    targetIterator.Set(sourceIterator.Get());

    ++sourceIterator;
    ++targetIterator;
  }
}

std::vector<itk::Index<2> > Get4NeighborIndicesInsideRegion(const itk::Index<2>& pixel,
                                                            const itk::ImageRegion<2>& region)
{
  std::vector<itk::Index<2> > indices;

  itk::Offset<2> offset;
  offset[0] = -1;
  offset[1] = 0;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 1;
  offset[1] = 0;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 0;
  offset[1] = -1;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 0;
  offset[1] = 1;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  return indices;
}

template<typename TImage>
std::vector<itk::Index<2> > Get4NeighborsWithValue(const TImage* const image,
                                                   const itk::Index<2>& pixel,
                                                   const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > neighborsWithValue;

  std::vector<itk::Index<2> > neighbors =
      Get4NeighborIndicesInsideRegion(pixel, image->GetLargestPossibleRegion());

  for (auto neighbor : neighbors)
  {
    if(image->GetPixel(neighbor) == value)
    {
      neighborsWithValue.push_back(neighbor);
    }
  }

  return neighborsWithValue;
}

template<typename TImage>
bool Has4NeighborsWithValue(const TImage* const image,
                            const itk::Index<2>& pixel,
                            const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > neighbors = Get4NeighborsWithValue(image, pixel, value);
  if(neighbors.size() > 0)
  {
    return true;
  }

  return false;
}

template<typename TInputImage, typename TOutputImage>
void DeepCopyInRegion(const TInputImage* input, const itk::ImageRegion<2>& region, TOutputImage* output)
{
  // This function assumes that the size of input and output are the same.

  itk::ImageRegionConstIterator<TInputImage> inputIterator(input, region);
  itk::ImageRegionIterator<TOutputImage> outputIterator(output, region);

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

/** Copy the input to the output*/
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::Image<TInputPixel, 2>* const input, itk::Image<TOutputPixel, 2>* const output)
{
  //std::cout << "DeepCopy()" << std::endl;
  if(output->GetLargestPossibleRegion() != input->GetLargestPossibleRegion())
    {
    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();
    }
  DeepCopyInRegion(input, input->GetLargestPossibleRegion(), output);
}

// This is a specialization that ensures that the number of pixels per component also matches.
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::VectorImage<TInputPixel, 2>* const input, itk::VectorImage<TOutputPixel, 2>* const output)
{
  //std::cout << "DeepCopy<FloatVectorImageType>()" << std::endl;
  bool changed = false;
  if(input->GetNumberOfComponentsPerPixel() != output->GetNumberOfComponentsPerPixel())
    {
    output->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
    //std::cout << "Set output NumberOfComponentsPerPixel to "
    //          << input->GetNumberOfComponentsPerPixel() << std::endl;
    changed = true;
    }

  if(input->GetLargestPossibleRegion() != output->GetLargestPossibleRegion())
    {
    output->SetRegions(input->GetLargestPossibleRegion());
    changed = true;
    }
  if(changed)
    {
    output->Allocate();
    }

  DeepCopyInRegion(input, input->GetLargestPossibleRegion(), output);
}

template<typename TImage>
void WriteImage(const TImage* const image, const std::string& filename)
{
  try
  {
    // This is a convenience function so that images can be written in 1 line instead of 4.
    typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();
    writer->SetFileName(filename);
    writer->SetInput(image);
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::stringstream ss;
    ss << "ITKHelpers::WriteImage failed! Region is: " << image->GetLargestPossibleRegion()
       << "Original ITK exception: " << err;
    throw std::runtime_error(ss.str());
  }
}

template<typename TImage>
void SetRegionToConstant(TImage* const image, const itk::ImageRegion<2>& region,
                         const typename TImage::PixelType& value)
{
  typename itk::ImageRegionIterator<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(value);

    ++imageIterator;
    }
}








}

#endif /* HELPERS_HPP_ */