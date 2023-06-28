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
// Created by tabish on 2021-04-06.
//

#ifndef TURBO_SLIC_H
#define TURBO_SLIC_H

#include <itkLogger.h>
#include "itkCommandLineArgumentParser.h"

int slic(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger);


#endif //TURBO_SLIC_H
