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
// Created by tabish on 2021-04-13.
//

#ifndef itkBradleyAdaptiveThresholdImageCalculator_h
#define itkBradleyAdaptiveThresholdImageCalculator_h

#include <itkImageToImageFilter.h>
#include <itkIntTypes.h>
#include <itkMacro.h>
#include <itkNumericTraits.h>
#include "itkIntegralImageFilter.h"

namespace itk
{
/**
 *\class BradleyAdaptiveThresholdImageCalculator
 * \brief Computes Bradley Adaptive threshold for input image
 * Reference:
 *  Bradley, D., G. Roth, "Adapting Thresholding Using the Integral Image,"
 *  Journal of Graphics Tools. Vol. 12, No. 2, 2007, pp.13–21.
 * \ingroup IntensityImageFilters
 * \ingroup ITKThresholding
 */
template <typename TInputImage>
class ITK_TEMPLATE_EXPORT BradleyAdaptiveThresholdImageCalculator
    : public ImageToImageFilter<TInputImage,
                                Image<float, TInputImage::ImageDimension>> {
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(BradleyAdaptiveThresholdImageCalculator);

  /** Number of dimensions. */
  static constexpr unsigned int NDimensions = TInputImage::ImageDimension;

  /** Standard class type aliases. */
  using Self = BradleyAdaptiveThresholdImageCalculator;
  using OutputImageType = Image<float, TInputImage::ImageDimension>;

  /** Standard class type aliases. */
  using Superclass = ImageToImageFilter<TInputImage, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BradleyAdaptiveThresholdImageCalculator, ImageToImageFilter);

  /** Image size type alias. */
  using SizeType = Size<Self::NDimensions>;

  /**Typedef for input image*/
  using InputImageType = TInputImage;
  using InputImageRegionType = typename TInputImage::RegionType;

  /** Image pixel value type alias. */
  using InputPixelType = typename TInputImage::PixelType;

  /** Typedef for integral image. */
  using InternalPixelType = float;
  using InternalImageType = Image<InternalPixelType, NDimensions>;

  /**Typedef for output image*/

  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputPixelValueType = typename OutputImageType::PixelType;
  using ImageIndexType = typename OutputImageType::IndexType;

  /**Type of internal integral image filter*/
  using IntegralImageFilterType =
      IntegralImageFilter<InputImageType, InternalImageType>;

  /** Getters and setters for parameters*/
  itkGetConstMacro(Sensitivity, double);
  itkSetMacro(Sensitivity, double);

  itkGetConstMacro(WindowHalfWidth, SizeValueType);
  void SetWindowHalfWidth(SizeValueType w) {
    if (w < 1) {
      itkWarningMacro("Window Size should be > 0");
    } else {
      m_WindowHalfWidth = w;
    }
  }
  typename InternalImageType::Pointer GetIntegralImage(){
      return m_integralImage;
  }
protected:
  BradleyAdaptiveThresholdImageCalculator();
  ~BradleyAdaptiveThresholdImageCalculator() override = default;
  void PrintSelf(std::ostream &os, Indent indent) const override;
  void GenerateData() override;

  InternalPixelType integrate(ImageIndexType, ImageIndexType);

private:
  typename IntegralImageFilterType::Pointer m_integralImageFilter;
  SizeValueType m_WindowHalfWidth;
  double m_Sensitivity;
  typename InternalImageType::Pointer m_integralImage;
        };
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBradleyAdaptiveThresholdImageCalculator.hxx"
#endif

#endif
