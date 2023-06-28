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
// Created by tabish on 2021-04-07.
//

#ifndef TURBO_UTILS_H
#define TURBO_UTILS_H

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkLogger.h>
#include <string>
#include "itkCommandLineArgumentParser.h"
#include <itkWhitakerSparseLevelSetImage.h>
#include <itkLevelSetDenseImage.h>

using SparseLevelSetType2 = itk::WhitakerSparseLevelSetImage<float,2>;
using SparseLevelSetType3 = itk::WhitakerSparseLevelSetImage<float,3>;
using LevelSetType2 = itk::LevelSetDenseImage<itk::Image<float,2>>;
using LevelSetType3 = itk::LevelSetDenseImage<itk::Image<float,3>>;

using RGBPixelType = itk::RGBPixel<unsigned char>;
using Color = std::vector<double>;



template <typename OutputImageType>
int writeImage(typename OutputImageType::Pointer image,
               const std::string outputFileName,
               const itk::Logger::Pointer &logger,
               const std::string message = "Writing Image To");

template <typename InputImageType>
void generateSeeds(std::vector<typename InputImageType::IndexType> &seeds,
                   size_t d, typename InputImageType::IndexType seed,
                   typename InputImageType::SizeType cellSize,
                   typename InputImageType::Pointer &image);

template <typename LevelSetType, typename InputImageType>
typename LevelSetType::Pointer
initializeLevelSet(const itk::CommandLineArgumentParser::Pointer &parser,
                   const itk::Logger::Pointer &logger,
                   typename InputImageType::Pointer &input);

template <typename InternalImageType>
typename InternalImageType::Pointer
compute_edge_speed(const itk::CommandLineArgumentParser::Pointer &parser,
                   const itk::Logger::Pointer &logger,
                   typename InternalImageType::Pointer &input);

template<unsigned Dimension>
void
getImage(typename itk::WhitakerSparseLevelSetImage<float, Dimension>::Pointer& levelSet,
         typename itk::Image<float,Dimension>::Pointer& levelSetImage);

template<unsigned Dimension>
void
getImage(typename itk::LevelSetDenseImage<itk::Image<float, Dimension>>::Pointer& levelSet,
         typename itk::Image<float,Dimension>::Pointer& levelSetImage);

template<typename RGBImageType, typename MaskImageType>
void paintImage(typename RGBImageType::Pointer& input,
                typename MaskImageType::Pointer& mask,
                Color c);

template<typename InputImageType,typename RGBImageType>
typename RGBImageType::Pointer
paintCountours(typename InputImageType::Pointer& input,
               std::vector<typename InputImageType::PixelType>& values);

template<typename PixelType, unsigned Dimension>
void
set_boundary(typename itk::Image<PixelType,Dimension>::Pointer& input, PixelType value);

#include "utils.hxx"
#include <vector>

#endif //TURBO_UTILS_H
