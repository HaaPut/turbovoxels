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
// Created by tabish on 2021-04-14.
//

#ifndef itkAdaptiveThresholdImageFilter_h
#define itkAdaptiveThresholdImageFilter_h

#include "itkBinaryGeneratorImageFilter.h"
#include "itkNumericTraits.h"
#include "itkMath.h"

namespace itk
{
    namespace Functor
    {
    /**
     * \class AdaptiveThreshold
     * \brief
     * \ingroup ITKImageThresholding
     */
        template <typename TInput, typename TThreshold, typename TOutput = TInput>
        class ThresholdAbove
        {
        public:
            ThresholdAbove(){
                    m_InsideValue = NumericTraits<TOutput>::OneValue();
                    m_OutsideValue = NumericTraits<TOutput>::ZeroValue();
            }
            ~ThresholdAbove() = default;
            bool
            operator!=(const ThresholdAbove &) const{
                return false;
            }
            bool
            operator==(const ThresholdAbove & other) const{
                return !(*this != other);
            }

            inline TOutput
            operator()(const TInput & A, const TThreshold & B) const {
                if (A >= B){
                    return m_InsideValue;
                }
                else{
                    return m_OutsideValue;
                }
            }

            /** Method to explicitly set the outside value of the mask */
            void
            SetOutsideValue(const TOutput & outsideValue){
                m_OutsideValue = outsideValue;
            }

            /** Method to get the outside value of the mask */
            const TOutput &
            GetOutsideValue() const{
                    return m_OutsideValue;
            }
            /** Method to explicitly set the inside value of the mask */
            void
            SetInsideValue(const TOutput & insideValue){
                m_InsideValue = insideValue;
            }

            /** Method to get the inside value of the mask */
            const TOutput &
            GetInsideValue() const{
                return m_InsideValue;
            }

        private:
            TOutput m_InsideValue;
            TOutput m_OutsideValue;
        };
    } // namespace Functor

    /** \class AdaptiveThresholdImageFilter
     * \brief Thresholds images pixel-wise using a threshold image or value.
     \code
     output_pixel = static_cast<OutputPixelType>( input1_pixel >= threshold_value )
     \endcode
     *
     *
     *
     * \ingroup ITKImageThresholding
     */
    template <typename TInputImage1, typename TInputImage2, typename TOutputImage>
    class AdaptiveThresholdImageFilter : public BinaryGeneratorImageFilter<TInputImage1, TInputImage2, TOutputImage>
    {
    public:
        //ITK_DISALLOW_COPY_AND_MOVE(AdaptiveThresholdImageFilter);

        /** Standard class type aliases. */
        using Self = AdaptiveThresholdImageFilter;
        using Superclass = BinaryGeneratorImageFilter<TInputImage1, TInputImage2, TOutputImage>;


        using FunctorType = Functor::
            ThresholdAbove<typename TInputImage1::PixelType, typename TInputImage2::PixelType, typename TOutputImage::PixelType>;

        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        //using RealType = typename FunctorType::RealType;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Runtime information support. */
        itkTypeMacro(AdaptiveThresholdImageFilter, BinaryGeneratorImageFilter);


        /** Method to explicitly set the outside value of output. Defaults to 0 */
        void
        SetOutsideValue(const typename TOutputImage::PixelType & outsideValue){
            if (Math::NotExactlyEquals(this->GetOutsideValue(), outsideValue)){
                this->Modified();
                this->GetFunctor().SetOutsideValue(outsideValue);
            }
        }
        const typename TOutputImage::PixelType &
        GetOutsideValue() const{
            return this->GetFunctor().GetOutsideValue();
        }


        /** Method to explicitly set the inside value of output. Defaults to 1 */
        void
        SetInsideValue(const typename TOutputImage::PixelType & insideValue){
            if (Math::NotExactlyEquals(this->GetInsideValue(), insideValue)){
                this->Modified();
                this->GetFunctor().SetInsideValue(insideValue);
            }
        }
        const typename TOutputImage::PixelType &
        GetInsideValue() const{
            return this->GetFunctor().GetInsideValue();
        }



#ifdef ITK_USE_CONCEPT_CHECKING
        // Begin concept checking
        itkConceptMacro(Input1HasNumericTraitsCheck, (Concept::HasNumericTraits<typename TInputImage1::PixelType>));
        // End concept checking
#endif

protected:
        AdaptiveThresholdImageFilter() = default;
        ~AdaptiveThresholdImageFilter() override = default;

        void
        BeforeThreadedGenerateData() override{
                this->SetFunctor(this->GetFunctor());
        }
private:
        itkGetConstReferenceMacro(Functor, FunctorType);
        FunctorType &
        GetFunctor(){
            return m_Functor;
        }

        FunctorType m_Functor;
    };
} // end namespace itk

#endif

//
// #ifndef itkAdaptiveThresholdImageFilter_h
// #define itkAdaptiveThresholdImageFilter_h

// #include <itkImageToImageFilter.h>
// #include <itkIntTypes.h>
// #include <itkNumericTraits.h>
// #include "itkIntegralImageFilter.h"

// namespace itk
// {
// /**
//  *\class AdaptiveThresholdImageFilter
//  * \brief Threshold input image using adaptive threshold
//  *
//  * \ingroup ImageFeatures
//  */
//     template <typename TInputImage, typename TOutputImage>
//         class ITK_TEMPLATE_EXPORT AdaptiveThresholdImageFilter: public ImageToImageFilter< TInputImage, TOutputImage>
//         {
//         public:
//         ITK_DISALLOW_COPY_AND_ASSIGN(AdaptiveThresholdImageFilter);

//         /** Number of dimensions. */
//         static constexpr unsigned int NDimensions = TInputImage::ImageDimension;

//         /** Standard class type aliases. */
//         using Self = AdaptiveThresholdImageFilter;

//         /** Standard class type aliases. */
//         using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
//         using Pointer = SmartPointer<Self>;
//         using ConstPointer = SmartPointer<const Self>;

//         /** Method for creation through the object factory. */
//         itkNewMacro(Self);

//         /** Run-time type information (and related methods). */
//         itkTypeMacro(AdaptiveThresholdImageFilter, ImageToImageFilter);

//         /** Image size type alias. */
//         using SizeType = Size<Self::NDimensions>;

//         /**Typedef for input image*/
//         using InputImageType = TInputImage;
//         using InputImageRegionType = typename TInputImage::RegionType;

//         /** Image pixel value type alias. */
//         using InputPixelType = typename TInputImage::PixelType;

//         /** Typedef for integral image. */
//         using InternalPixelType = typename NumericTraits<InputPixelType>::AccumulateType;
//         using InternalImageType = Image<InternalPixelType, NDimensions>;

//         /**Typedef for output image*/
//         using OutputImageType = TOutputImage;
//         using OutputImageRegionType = typename OutputImageType::RegionType;
//         using OutputPixelType = typename OutputImageType::PixelType;
//         using ImageIndexType = typename OutputImageType::IndexType;

//         /**Type of internal integral image filter*/
//         using IntegralImageFilterType =
//                 IntegralImageFilter<InputImageType, InternalImageType>;

//         /** Getters and setters for parameters*/
//         itkGetConstMacro(Sensitivity, double);
//         itkSetMacro(Sensitivity, double);

//         itkGetConstMacro(WindowHalfWidth, SizeValueType);
//         itkSetMacro(WindowHalfWidth, SizeValueType);

//         itkGetConstMacro(InsideValue, OutputPixelType);
//         itkSetMacro(InsideValue, OutputPixelType);

//         itkGetConstMacro(OutsideValue, OutputPixelType);
//         itkSetMacro(OutsideValue, OutputPixelType);
//       protected:
//         AdaptiveThresholdImageFilter(); // = default;
//         ~AdaptiveThresholdImageFilter() override = default;
//         void
//         PrintSelf(std::ostream & os, Indent indent) const override;
//         void
//         GenerateData() override;

//         InternalPixelType
//         integrate(ImageIndexType, ImageIndexType);
//         InternalPixelType
//             integrate_raw2(ImageIndexType, ImageIndexType);

//         private:
//         typename IntegralImageFilterType::Pointer m_integralImageFilter;
//         SizeValueType m_WindowHalfWidth;
//         double m_Sensitivity;
//         OutputPixelType m_InsideValue;
//         OutputPixelType m_OutsideValue;
//         typename InternalImageType::Pointer m_integralImage;
//         };
// } // end namespace itk

// #ifndef ITK_MANUAL_INSTANTIATION
// #include "itkAdaptiveThresholdImageFilter.hxx"
// #endif

// #endif
