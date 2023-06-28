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
// Created by tabish on 2021-06-12.
//

#ifndef TURBO_BOUNDARIES_H
#define TURBO_BOUNDARIES_H

#include "itkCommandLineArgumentParser.h"

#include <itkLogger.h>
#include <itkImage.h>

using FloatImageType3 = itk::Image<float,3>;
using FloatImageType2 = itk::Image<float,2>;

using LabelImageType3 = itk::Image<uint16_t ,3>;
using LabelImageType2 = itk::Image<uint16_t ,2>;

template<typename LevelSetImageType, typename LabelImageType>
typename LabelImageType::Pointer
levelSetImageToBoundary( typename LevelSetImageType::Pointer levelSetImage,
                         const itk::CommandLineArgumentParser::Pointer &parser,
                         const itk::Logger::Pointer &logger);

template<typename LevelSetImageType, typename LabelImageType>
typename LabelImageType::Pointer
levelSetImageToTurboVoxels( typename LevelSetImageType::Pointer levelSetImage,
                            const itk::CommandLineArgumentParser::Pointer &parser,
                            const itk::Logger::Pointer &logger);

#include "boundaries.hxx"


#endif //TURBO_BOUNDARIES
