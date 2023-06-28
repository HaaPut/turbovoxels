//**********************************************************
// Copyright 2021 Tabish Syed
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//**********************************************************
//
// Created by tabish on 2021-05-24.
//

#ifndef itkBinaryPruningImageFilterMy_h
#define itkBinaryPruningImageFilterMy_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"

namespace itk
{
/**
 *\class BinaryPruningImageFilterMy
 *
 * \brief This filter removes "spurs" of less than a certain
 * length in the input image.
 *
 * This class is parameterized over the type of the input image
 * and the type of the output image.
 *
 * The input is assumed to be a binary image.
 *
 * This filter is a sequential pruning algorithm and known to be computational time
 * dependable of the image size.  The algorithm is the N-dimensional version
 * of that given for two dimensions in:
 *
 * Rafael C. Gonzales and Richard E. Woods.
 * Digital Image Processing.
 * Addison Wesley, 491-494, (1993).
 *
 * \sa MorphologyImageFilter
 * \sa BinaryErodeImageFilter
 * \sa BinaryDilateImageFilter
 * \sa BinaryThinningImageFilter
 * \ingroup ImageEnhancement MathematicalMorphologyImageFilters
 * \ingroup ITKBinaryMathematicalMorphology
 *
 * \sphinx
 * \sphinxexample{Filtering/BinaryMathematicalMorphology/PruneBinaryImage,Prune Binary Image}
 * \endsphinx
 */

template <typename TInputImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT BinaryPruningImageFilterMy : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    //ITK_DISALLOW_COPY_AND_MOVE(BinaryPruningImageFilterMy);

  /** Standard class type aliases. */
  using Self = BinaryPruningImageFilterMy;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BinaryPruningImageFilterMy, ImageToImageFilter);

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

  /** Neighborhood iterator type */
  using NeighborhoodIteratorType = NeighborhoodIterator<TInputImage>;

  /** Offset Type*/
  using OffsetType = typename NeighborhoodIteratorType::OffsetType;

  /** Get Skelenton by thinning image. */
  OutputImageType *
  GetPruning();


  /** Set/Get the iteration value */
  itkSetMacro(Iteration, unsigned int);
  itkGetConstMacro(Iteration, unsigned int);

  /** Set/Get Foreground value */
  itkSetMacro(ForegroundValue, PixelType);
  itkGetConstMacro(ForegroundValue, PixelType);

  /** Set/Get Background value */
  itkSetMacro(BackgroundValue, PixelType);
  itkGetConstMacro(BackgroundValue, PixelType);

  /** ImageDimension enumeration   */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(SameDimensionCheck, (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(SameTypeCheck, (Concept::SameType<PixelType, typename TOutputImage::PixelType>));
  itkConceptMacro(AdditiveOperatorsCheck, (Concept::AdditiveOperators<PixelType>));
  itkConceptMacro(IntConvertibleToPixelTypeCheck, (Concept::Convertible<int, PixelType>));
  itkConceptMacro(PixelLessThanIntCheck, (Concept::LessThanComparable<PixelType, int>));
  // End concept checking
#endif

protected:
  BinaryPruningImageFilterMy();
  ~BinaryPruningImageFilterMy() override = default;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Compute thinning Image. */
  void
  GenerateData() override;

  /** Prepare data. */
  void
  PrepareData();

  /**  Compute thinning Image. */
  void
  ComputePruneImage();

private:
  unsigned
      find(unsigned idx, std::vector<unsigned>& parent_vec);

  bool
      IsRemovable(NeighborhoodIteratorType& ot);


  unsigned int m_Iteration;
  PixelType m_ForegroundValue;
  PixelType m_BackgroundValue;
  OutputImagePointer m_NeighborhoodImage;
}; // end of BinaryThinningImageFilter class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryPruningImageFilterMy.hxx"

#endif

#endif
