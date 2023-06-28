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
#ifndef itkIntegralImageFilter_hxx
#define itkIntegralImageFilter_hxx

#include "itkIntegralImageFilter.h"
#include "itkImageRegionIterator.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageLinearIteratorWithIndex.h>
#include "itkNeighborhoodAlgorithm.h"
//#include <itkTotalProgressReporter.h>
#include <itkProgressReporter.h>


namespace itk
{
    template <typename TInputImage, typename TOutputImage>
    IntegralImageFilter<TInputImage, TOutputImage>::IntegralImageFilter(){
        this->DynamicMultiThreadingOff();
    }
    template <typename TInputImage, typename TOutputImage>
    void
    IntegralImageFilter<TInputImage, TOutputImage>::GenerateData(){
        OutputImageType* output = this->GetOutput();
        const InputImageType* input  = this->GetInput();

        output->SetLargestPossibleRegion(this->GetInput()->GetLargestPossibleRegion());
        output->SetBufferedRegion(this->GetInput()->GetLargestPossibleRegion());
        output->Allocate();

        InputImageRegionType  inputRegion = input->GetLargestPossibleRegion();
        OutputImageRegionType outputRegion = output->GetLargestPossibleRegion();

        using InputRegionIteratorType = ImageRegionConstIterator<InputImageType>;
        using OutputRegionIteratorType = ImageRegionIterator<OutputImageType>;
        InputRegionIteratorType ipIt(input, inputRegion);
        OutputRegionIteratorType outCopier(output, outputRegion);
        ipIt.GoToBegin();
        outCopier.GoToBegin();
        //TotalProgressReporter progress(this, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels()*(NDimensions+1));
        ProgressReporter progress(this, 0, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels()*(NDimensions+1));
        while(!ipIt.IsAtEnd()){
            outCopier.Set(ipIt.Get());
            ++ipIt;
            ++outCopier;
            progress.CompletedPixel();
        }
        using OutputIteratorType = ImageLinearIteratorWithIndex<OutputImageType>;
        OutputIteratorType outIt(output, outputRegion);
        using SumType = double;

        for(size_t d = 0; d < NDimensions; ++d){
            outIt.SetDirection(d);
            for(outIt.GoToBegin(); !outIt.IsAtEnd(); outIt.NextLine()){
                SumType sum = NumericTraits<SumType>::ZeroValue();
                outIt.GoToBeginOfLine();
                while(!outIt.IsAtEndOfLine()){
                    sum += outIt.Get();
                    outIt.Set(sum);
                    ++outIt;
                    progress.CompletedPixel();
                }
            }
        }
    }

    template <typename TInputImage, typename TOutputImage>
    void
    IntegralImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const{
        Superclass::PrintSelf(os, indent);

        os << indent
           << "Integral Image Filter: "
           << std::endl;
    }
} // namespace itk

#endif
