/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkBinaryThinningIterationImageFilter_h
#define itkBinaryThinningIterationImageFilter_h

#include "itkNeighborhoodIterator.h"
#include "itkImageToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
/**
 *\class BinaryThinningIterationImageFilter
 *
 * \brief This filter computes Skeletonization iterations untill
 * we get one pixel wide skeleton or reach maxIterations 
 *
 * This class is a minor modification over itkBinaryThinningImageFilter
 * \sphinx
 * \sphinxexample{Filtering/BinaryMathematicalMorphology/ThinImage,Thin Image}
 * \endsphinx
 */

template <typename TInputImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT BinaryThinningIterationImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    //ITK_DISALLOW_COPY_AND_MOVE(BinaryThinningIterationImageFilter);

  /** Standard class type aliases. */
  using Self = BinaryThinningIterationImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BinaryThinningIterationImageFilter, ImageToImageFilter);

  /** Type for input image. */
  using InputImageType = TInputImage;

  /** Type for output image: Skelenton of the object.  */
  using OutputImageType = TOutputImage;

  /** Type for the region of the input image. */
  using RegionType = typename InputImageType::RegionType;

  /** Type for the index of the input image. */
  using IndexType = typename RegionType::IndexType;

  /** Type for the index of the input image. */
  using PixelType = typename InputImageType::PixelType;

  /** Type for the size of the input image. */
  using SizeType = typename RegionType::SizeType;

  /** Pointer Type for input image. */
  using InputImagePointer = typename InputImageType::ConstPointer;

  /** Pointer Type for the output image. */
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass type alias. */
  using OutputImagePixelType = typename Superclass::OutputImagePixelType;

  /** Neighborhood iterator type */
  using NeighborhoodIteratorType = NeighborhoodIterator<TInputImage>;

  /** Get Skelenton by thinning image. */
  OutputImageType *
  GetThinning();

  /** ImageDimension enumeration   */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  itkSetMacro(MaxIterations, unsigned);
  itkGetConstMacro(MaxIterations, unsigned);

  itkGetConstMacro(RemoveCount, double);

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(SameDimensionCheck, (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(InputAdditiveOperatorsCheck, (Concept::AdditiveOperators<PixelType>));
  itkConceptMacro(InputConvertibleToIntCheck, (Concept::Convertible<PixelType, int>));
  itkConceptMacro(IntConvertibleToInputCheck, (Concept::Convertible<int, PixelType>));
  itkConceptMacro(SameTypeCheck, (Concept::SameType<PixelType, typename TOutputImage::PixelType>));
  // End concept checking
#endif

protected:
  BinaryThinningIterationImageFilter();
  ~BinaryThinningIterationImageFilter() override = default;

  /** Compute thinning Image. */
  void
  GenerateData() override;

  /** Prepare data. */
  void
  PrepareData();

  /**  Compute thinning Image.
   * Returns number of pixels removed */
  double
  ComputeThinImage();

private:
  unsigned m_MaxIterations;
  double m_RemoveCount;
}; // end of BinaryThinningIterationImageFilter
   // class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryThinningIterationImageFilter.hxx"
#endif

#endif
