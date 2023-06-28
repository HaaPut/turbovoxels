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
#ifndef itkBradleyAdaptiveThresholdImageCalculator_hxx
#define itkBradleyAdaptiveThresholdImageCalculator_hxx

#include "itkIntegralImageFilter.h"
#include "itkBradleyAdaptiveThresholdImageCalculator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <climits>
#include <itkImageRegionConstIterator.h>
#include <itkImageLinearIteratorWithIndex.h>
#include "itkNeighborhoodAlgorithm.h"
//#include <itkTotalProgressReporter.h>
#include <itkIntTypes.h>
#include <itkMacro.h>
#include <itkNumericTraits.h>
#include <itkProgressReporter.h>
#include <algorithm>
#include <bitset>
#include <itkMinimumMaximumImageCalculator.h>
#include <cmath>
#include <string>


namespace itk
{
    template <typename TInputImage>
    BradleyAdaptiveThresholdImageCalculator<TInputImage>::BradleyAdaptiveThresholdImageCalculator(){
        this->DynamicMultiThreadingOff();
        m_integralImageFilter = IntegralImageFilterType::New();
        m_WindowHalfWidth = 0;
        m_Sensitivity = 0.5;
    }

    template <typename TInputImage>
    typename BradleyAdaptiveThresholdImageCalculator<TInputImage>::InternalPixelType
    BradleyAdaptiveThresholdImageCalculator<TInputImage>::integrate(
        ImageIndexType start, ImageIndexType end) {
        InternalPixelType rectSum = 0;
        unsigned terms = (1 << NDimensions); //2^D permutation
        ImageIndexType idx;
        bool cornerOutsideImage;
        for(unsigned p = 0; p < terms; ++p){
            std::bitset<NDimensions> perm(p);
            int sign = perm.count() & 1 ? -1 : 1; //check if permutation is even or odd
            cornerOutsideImage = false;
            for(size_t d = 0 ; d < NDimensions; ++d){
                if(perm[d]){ // choose from start index
                    if (start[d] == 0){
                        cornerOutsideImage = true;
                        break;
                    }
                    idx[d] = start[d] - 1;
                }else{
                    idx[d] = end[d];
                }
            }
            if (!cornerOutsideImage)
                rectSum += sign * m_integralImage->GetPixel(idx);
        }
        return rectSum;
    }

    template <typename TInputImage>
    void
    BradleyAdaptiveThresholdImageCalculator<TInputImage>::GenerateData(){
        OutputImageType* output = this->GetOutput();
        const InputImageType* input  = this->GetInput();

        m_integralImageFilter->SetInput(input);
        m_integralImageFilter->Update();
        m_integralImage = m_integralImageFilter->GetOutput();

        using ImageCalculatorFilterType = MinimumMaximumImageCalculator<InputImageType>;

        typename ImageCalculatorFilterType::Pointer imageCalculatorFilter =
            ImageCalculatorFilterType::New();
        imageCalculatorFilter->SetImage(input);
        imageCalculatorFilter->Compute();

        double normalizer =
            static_cast<double>(imageCalculatorFilter->GetMaximum() -
                                imageCalculatorFilter->GetMinimum());
        itkDebugMacro("Set normalizer to : " + std::to_string(normalizer));
        output->SetLargestPossibleRegion(this->GetInput()->GetLargestPossibleRegion());
        output->SetBufferedRegion(this->GetInput()->GetLargestPossibleRegion());
        output->Allocate();

        InputImageRegionType  inputRegion = input->GetLargestPossibleRegion();
        OutputImageRegionType outputRegion = output->GetLargestPossibleRegion();
        SizeType outputSize = outputRegion.GetSize();
        if (m_WindowHalfWidth == 0) {
            //Default window width = max_dimension/8 + 1
            for (size_t d = 0; d < NDimensions; ++d) {
                m_WindowHalfWidth =
                    std::max(static_cast<SizeValueType>(std::floor(outputSize[d] / 16)),
                             m_WindowHalfWidth);
            }
        }

        using InputRegionIteratorType = ImageRegionConstIterator<InputImageType>;
        using OutputRegionIteratorType = ImageRegionIteratorWithIndex<OutputImageType>;
        InputRegionIteratorType ipIt(input, inputRegion);
        OutputRegionIteratorType outIt(output, outputRegion);
        ImageIndexType upperLeft, lowerRight, idx;
        SizeValueType rectArea;
        InternalPixelType rectSum;
        ipIt.GoToBegin();
        outIt.GoToBegin();
        //TotalProgressReporter progress(this, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels()*(NDimensions+1));
        ProgressReporter progress(this, 0, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels());
        while(!outIt.IsAtEnd()){
            idx = outIt.GetIndex();
            rectArea = 1;
            for(size_t d = 0; d < NDimensions; ++d){
                upperLeft[d] = std::max((int)idx[d] - (int)m_WindowHalfWidth, 0);
                lowerRight[d] = std::min(idx[d] +  m_WindowHalfWidth, outputSize[d]-1);
                rectArea *= (lowerRight[d] - upperLeft[d] + 1);
            }
            rectSum = integrate(upperLeft, lowerRight);
            outIt.Set((rectSum * std::exp(0.5 - m_Sensitivity)) /(rectArea * normalizer));
            ++outIt;
            ++ipIt;
            progress.CompletedPixel();
        }
    }

    template <typename TInputImage>
    void
    BradleyAdaptiveThresholdImageCalculator<TInputImage>::PrintSelf(std::ostream & os, Indent indent) const{
        Superclass::PrintSelf(os, indent);

        os << indent
           << "Bradley Adaptive Threshold Image Calculator: "
           << std::endl;
    }
} // namespace itk

#endif
