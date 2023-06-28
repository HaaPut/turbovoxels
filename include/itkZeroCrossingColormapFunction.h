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
// Created by tabish on 2021-04-08.
//
#ifndef itkZeroCrossingColormapFunction_h
#define itkZeroCrossingColormapFunction_h

#include "itkColormapFunction.h"
#include <algorithm>
#include <cstdlib>

namespace itk
{
    namespace Function
    {
    /**
     * \class ZeroCrossingColormapFunction
     * \brief Function object which maps a scalar value into an RGB colormap value.
     *
     */
        template <typename TScalar, typename TRGBPixel>
            class ITK_TEMPLATE_EXPORT ZeroCrossingColormapFunction : public ColormapFunction<TScalar, TRGBPixel>
        {
        public:


            using Self = ZeroCrossingColormapFunction;
            using Superclass = ColormapFunction<TScalar, TRGBPixel>;
            using Pointer = SmartPointer<Self>;
            using ConstPointer = SmartPointer<const Self>;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            using RGBPixelType = typename Superclass::RGBPixelType;
            using ScalarType = typename Superclass::ScalarType;
            using RealType = typename Superclass::RealType;

            RGBPixelType operator()( const TScalar &v) const {
                RGBPixelType pixel;
                NumericTraits<TRGBPixel>::SetLength(pixel, 3);

                // Map the input scalar between [0, 1].
                if (std::abs(v) <= 0.5) {
                    pixel[0] = this->RescaleRGBComponentValue(1.0);
                    pixel[1] = 0;
                    pixel[2] = 0;
                }else{
                    RealType value = std::min(this->RescaleInputValue(v),1.0);
                    value = std::max(0.0,value);
                    // Set the rgb components after rescaling the values.
                    pixel[0] = this->RescaleRGBComponentValue(value);
                    pixel[1] = pixel[0];
                    pixel[2] = pixel[0];
                }
                return pixel;
            }

        protected:
            ZeroCrossingColormapFunction() = default;
            ~ZeroCrossingColormapFunction() override = default;
        };
    } // end namespace Function
} // end namespace itk

#endif
