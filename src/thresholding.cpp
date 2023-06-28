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
// Created by tabish on 2021-04-10.
//

#include "thresholding.h"
#include "utils.h"

#include <algorithm>
#include <cstdlib>
#include "filesystem.hpp"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkMinMaxCurvatureFlowImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include "itkBradleyAdaptiveThresholdImageCalculator.h"
#include "itkAdaptiveThresholdImageFilter.h"
#include <string>
#include <typeinfo>
#include <vector>
#include <iostream>
#include <cmath>

namespace fs = std::filesystem;

template <typename InputPixelType, unsigned int Dimension>
int bradley_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                 const itk::Logger::Pointer &logger) {
    using OutputPixelType = InputPixelType;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);
    logger->Info("Input image file " + reader->GetFileName() + "\n");

    using ThresholdCalculatorType =
        itk::BradleyAdaptiveThresholdImageCalculator<InputImageType>;
    typename ThresholdCalculatorType::Pointer adaptThresh =
    ThresholdCalculatorType::New();
    adaptThresh->SetInput(reader->GetOutput());
    using ThresholdImageType = typename ThresholdCalculatorType::OutputImageType;

    float foregroundSensitivity = 0.5;
    parser->GetCommandLineArgument("-sensitivity",foregroundSensitivity);
    if(foregroundSensitivity > 1 || foregroundSensitivity < 0){
        logger->Warning("Foreground sensitivity should be a number between 0 and 1\n");
    }
    foregroundSensitivity = std::max(std::min(foregroundSensitivity, 1.0f), 0.0f);
    adaptThresh->SetSensitivity(foregroundSensitivity);
    logger->Info("Set foreground Sensitivity to: " + std::to_string(foregroundSensitivity) + "\n");
    reader->Update();
    typename InputImageType::Pointer image = reader->GetOutput();
    typename InputImageType::RegionType region = image->GetLargestPossibleRegion();
    typename InputImageType::SizeType size = region.GetSize();

    unsigned windowHalfWidth = 0;
    for(size_t d = 0; d < Dimension; ++d){
        windowHalfWidth = std::max(static_cast<unsigned>(std::floor(size[d]/16)), windowHalfWidth);
    }
    parser->GetCommandLineArgument("-window", windowHalfWidth);
    adaptThresh->SetWindowHalfWidth(windowHalfWidth);
    logger->Info("Set Window width to 2 x " + std::to_string(windowHalfWidth) + "\n");

    using scaleAdaptiveThresholdFilterType = itk::MultiplyImageFilter<ThresholdImageType>;
    typename scaleAdaptiveThresholdFilterType::Pointer scaleFilter = scaleAdaptiveThresholdFilterType::New();
    scaleFilter->SetInput1(adaptThresh->GetOutput());
    InputPixelType scallingConstant = itk::NumericTraits<InputPixelType>::max();
    if (typeid(scallingConstant) == typeid(float) ||
        typeid(scallingConstant) == typeid(double)){
        std::cout << "Input is float \n";
        scallingConstant = itk::NumericTraits<InputPixelType>::OneValue();
    }
    scaleFilter->SetInput2(scallingConstant);
    logger->Info("scalling normalized threshold by " + std::to_string(scallingConstant) + "\n");

    using ThresholdFilterType =
        itk::AdaptiveThresholdImageFilter<InputImageType, ThresholdImageType,
                                          OutputImageType>;
    typename ThresholdFilterType::Pointer binarize = ThresholdFilterType::New();
    binarize->SetInput1(reader->GetOutput());
    binarize->SetInput2(scaleFilter->GetOutput());
    if(parser->ArgumentExists("-dark")){
        binarize->SetOutsideValue(static_cast<OutputPixelType>(scallingConstant));
        binarize->SetInsideValue(itk::NumericTraits<OutputPixelType>::ZeroValue());
    }else{
        binarize->SetOutsideValue(itk::NumericTraits<OutputPixelType>::ZeroValue());
        binarize->SetInsideValue(static_cast<OutputPixelType>(scallingConstant));
    }
    binarize->Update();

    std::string outputFileName;
    parser->GetCommandLineArgument("-output", outputFileName);
    logger->Info("output image file" + outputFileName + "\n");
    typename OutputImageType::Pointer result = binarize->GetOutput();
    writeImage<OutputImageType>(result, outputFileName, logger, "Wrote Adaptive Thresholded Image ");
    if(parser->ArgumentExists("-writeThreshold")){
        typename ThresholdImageType::Pointer threshold = scaleFilter->GetOutput();
        fs::path inputFilepath(inputFileName);
        fs::path threshFilepath = inputFilepath.parent_path() / inputFilepath.stem();
        threshFilepath += "_thresh.tif";
        writeImage<ThresholdImageType>(threshold, threshFilepath.string(), logger,
                                       "Wrote Adaptive Threshold to ");
    }
    return EXIT_SUCCESS;
}


template <typename InputPixelType, unsigned int Dimension>
int otsu_impl(const itk::CommandLineArgumentParser::Pointer &parser,
              const itk::Logger::Pointer &logger) {
    using OutputPixelType = InputPixelType;
    using InternalPixelType = float;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;
    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);

    using SmoothingFilterType =itk::MinMaxCurvatureFlowImageFilter<InputImageType, InternalImageType>;
    typename SmoothingFilterType::Pointer smoothingFilter =
        SmoothingFilterType::New();
    logger->Info("Using min-max Curvature smoothing for preprocessing\n");
    unsigned stencilRadius = 2;
    parser->GetCommandLineArgument("-radius", stencilRadius);
    logger->Info("Set radius to " + std::to_string(stencilRadius) +"\n");
    smoothingFilter->SetStencilRadius(stencilRadius);
    smoothingFilter->SetInput(reader->GetOutput());
    unsigned int numberOfIterations = 0;
    double timeStep = 0.125;
    parser->GetCommandLineArgument("-iterations", numberOfIterations);
    logger->Info("Number of Iteration :: " + std::to_string(numberOfIterations) + "\n");
    parser->GetCommandLineArgument("-dt", timeStep);
    logger->Info("Time Step :: " + std::to_string(timeStep) + "\n");

    smoothingFilter->SetNumberOfIterations(numberOfIterations);
    smoothingFilter->SetTimeStep(timeStep);

    using ThresholdFilterType =
        itk::OtsuThresholdImageFilter<InternalImageType, OutputImageType>;
    typename ThresholdFilterType::Pointer otsu = ThresholdFilterType::New();
    otsu->SetInput(smoothingFilter->GetOutput());

    OutputPixelType bright = itk::NumericTraits<OutputPixelType>::max();
    if (typeid(bright) == typeid(float) ||
        typeid(bright) == typeid(double)){
        bright = itk::NumericTraits<OutputPixelType>::OneValue();
    }
    OutputPixelType dark = itk::NumericTraits<OutputPixelType>::ZeroValue();

    if(parser->ArgumentExists("-dark")){
        std::swap(bright,dark);
    }
    otsu->SetOutsideValue(bright);
    otsu->SetInsideValue(dark);
    otsu->Update();

    std::string outputFileName;
    parser->GetCommandLineArgument("-output", outputFileName);

    typename OutputImageType::Pointer result = otsu->GetOutput();
    writeImage<OutputImageType>(result, outputFileName, logger, "Wrote Otsu Thresholded Image ");

    return EXIT_SUCCESS;
}

int thresholding(const itk::CommandLineArgumentParser::Pointer &parser,
                      const itk::Logger::Pointer &logger) {

    logger->Info("Starting Thresholding module\n");
    std::string inputFileName, outputFileName;
    fs::path inputFilePath, outputFilePath;
    bool hasInputFilePath =
        parser->GetCommandLineArgument("-input", inputFileName);
    bool hasOutputFilePath =
        parser->GetCommandLineArgument("-output", outputFileName);
    // Read the input image.
    if (!hasInputFilePath) {
        logger->Error("Input file not specified\n");
        return EXIT_FAILURE;
    }else{
        inputFilePath = inputFileName;
    }
    if (!hasOutputFilePath) {
        outputFilePath = inputFilePath.parent_path()/inputFilePath.stem();
        outputFilePath += "_seg";
        outputFilePath += inputFilePath.extension();
        std::vector<std::string> newargs;
        newargs.push_back("-output");
        newargs.push_back(outputFilePath.string());
        parser->AddCommandLineArguments(newargs);
    }else{
        outputFilePath = outputFileName;
    }
    logger->Info("Input Filename : " + inputFilePath.string() + "\n");
    logger->Info("Output Filename : "  + outputFilePath.string() + "\n");
    if(!fs::exists(inputFileName)){
        logger->Critical("File " + inputFileName + " not found\n");
        return EXIT_FAILURE;
    }
    itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(inputFileName.c_str(), itk::ImageIOFactory::FileModeType::ReadMode);
    imageIO->SetFileName(inputFileName);
    imageIO->ReadImageInformation();
    itk::ImageIOBase::IOPixelType pixelType = imageIO->GetPixelType();
    itk::ImageIOBase::IOComponentType componentType = imageIO->GetComponentType();
    unsigned int dimensions = imageIO->GetNumberOfDimensions();

    logger->Info("Component Type  : " +
                 imageIO->GetComponentTypeAsString(componentType) + "\n");
    logger->Info("Pixel Type      : " +
                 itk::ImageIOBase::GetPixelTypeAsString(pixelType) + "\n");
    logger->Info("Image Dimension : " + std::to_string(dimensions) + "\n");

    switch (dimensions){
    case 2:
        switch (componentType){
        case itk::ImageIOBase::UCHAR:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<unsigned char,2>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<unsigned short, 2>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<unsigned short,2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<float, 2>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<float,2>(parser, logger);
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          logger->Critical("Unknown ComponentType :: " +
                           std::string(typeid(componentType).name()) +
                           " in 2 dimensions\n");
        }
        break;
    case 3:
        switch (componentType) {
        case itk::ImageIOBase::UCHAR:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<unsigned char, 3>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<unsigned char,3>(parser, logger);
        case itk::ImageIOBase::USHORT:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<unsigned short, 3>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<unsigned short,3>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            if (parser->ArgumentExists("-otsu"))
                return otsu_impl<float, 3>(parser, logger);
            else if(parser->ArgumentExists("-bradley"))
                return bradley_impl<float,3>(parser, logger);
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          logger->Critical("Unknown ComponentType :: " +
                           std::string(typeid(componentType).name()) +
                           " in 3 dimensions\n");
        }
        break;
    default:
        logger->Critical("Only 2 and 3D images are supported\n");
    }
    return EXIT_FAILURE;
}
