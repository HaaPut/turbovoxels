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
// Created by tabish on 2021-04-12.
//

#ifndef itkGrayToRGBColormapFunction_h
#define itkGrayToRGBColormapFunction_h

#include "itkColormapFunction.h"
#include <algorithm>

namespace itk
{
    namespace Function
    {
    /**
     * \class GrayToRGBColormapFunction
     * \brief Function object which maps a scalar value into an RGB colormap value.
     *
     */
        template <typename TScalar, typename TRGBPixel>
            class ITK_TEMPLATE_EXPORT GrayToRGBColormapFunction : public ColormapFunction<TScalar, TRGBPixel>
        {
        public:


            using Self = GrayToRGBColormapFunction;
            using Superclass = ColormapFunction<TScalar, TRGBPixel>;
            using Pointer = SmartPointer<Self>;
            using ConstPointer = SmartPointer<const Self>;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            using RGBPixelType = typename Superclass::RGBPixelType;
            using ScalarType = typename Superclass::ScalarType;
            using RealType = typename Superclass::RealType;

            RGBPixelType operator()( const TScalar &v) const override{
                RGBPixelType pixel;
                NumericTraits<TRGBPixel>::SetLength(pixel, 3);
                // Map the input scalar between [0, 1].
                RealType value = std::min(this->RescaleInputValue(v),1.0);
                value = std::max(0.0,value);
                // Set the rgb components after rescaling the values.
                pixel[0] = this->RescaleRGBComponentValue(value);
                pixel[1] = pixel[0];
                pixel[2] = pixel[0];
                return pixel;
            }

        protected:
            GrayToRGBColormapFunction() = default;
            ~GrayToRGBColormapFunction() override = default;

        };
    } // end namespace Function
}

#endif
