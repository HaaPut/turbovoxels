//**********************************************************
//Copyright 2021 Tabish Syed
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//**********************************************************
//
// Created by tabish on 2021-06-01.
//


#ifndef TURBO_LEVELSETEVOLUTION_H
#define TURBO_LEVELSETEVOLUTION_H

#include "itkCommandLineArgumentParser.h"

#include <itkLogger.h>
#include <itkImage.h>

using FloatImageType3 = itk::Image<float,3>;
using FloatImageType2 = itk::Image<float,2>;

template <typename TRealImageType> //<typename InternalPixelType, unsigned int Dimension>
typename TRealImageType::Pointer
maskedSparseLevelSetEvolution( typename TRealImageType::Pointer input,
                               const itk::CommandLineArgumentParser::Pointer &parser,
                               const itk::Logger::Pointer &logger);

#include "levelSetEvolution.hxx"

#endif //TURBO_LEVELSETEVOLUTION_H
