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

#ifndef itkIntegralImageFilter_h
#define itkIntegralImageFilter_h

#include <itkImageToImageFilter.h>


namespace itk
{
/**
 *\class IntegralImageFilter
 * \brief Computes Integral Image
 *
 * \ingroup ImageFeatures
 */
    template <typename TInputImage, typename TOutputImage = TInputImage>
        class ITK_TEMPLATE_EXPORT IntegralImageFilter : public ImageToImageFilter< TInputImage, TOutputImage>
        {
        public:
        ITK_DISALLOW_COPY_AND_ASSIGN(IntegralImageFilter);

        /** Number of dimensions. */
        static constexpr unsigned int NDimensions = TInputImage::ImageDimension;

        /** Standard class type aliases. */
        using Self = IntegralImageFilter;

        /** Standard class type aliases. */
        using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(IntegralImageFilter, ImageToImageFilter);

        /** Image size type alias. */
        using SizeType = Size<Self::NDimensions>;

        /**Typedef for input image*/
        using InputImageType = TInputImage;
        using InputImageRegionType = typename TInputImage::RegionType;

        /** Image pixel value type alias. */
        using InputPixelType = typename TInputImage::PixelType;

        /** Typedef for output image. */
        using OutputImageType = TOutputImage;
        using OutputImageRegionType = typename TOutputImage::RegionType;
        using OutputPixelType = typename OutputImageType::PixelType;

        protected:
        IntegralImageFilter();// = default;
        ~IntegralImageFilter() override = default;
        void
        PrintSelf(std::ostream & os, Indent indent) const override;
        void
        GenerateData() override;
        };
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIntegralImageFilter.hxx"
#endif

#endif
