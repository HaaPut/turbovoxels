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
// Created by tabish on 2021-06-14.
//

#include "boundaries.h"

template
typename LabelImageType2::Pointer
levelSetImageToBoundary<FloatImageType2,LabelImageType2>( typename FloatImageType2::Pointer levelSetImage,
                                                          const itk::CommandLineArgumentParser::Pointer &parser,
                                                          const itk::Logger::Pointer &logger);


template
typename LabelImageType3::Pointer
levelSetImageToBoundary<FloatImageType3,LabelImageType3>( typename FloatImageType3::Pointer levelSetImage,
                                                          const itk::CommandLineArgumentParser::Pointer &parser,
                                                          const itk::Logger::Pointer &logger);
