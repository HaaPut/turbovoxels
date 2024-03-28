
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
// Created by tabish on 2021-03-26.
//

#include "smoothing.h"
#include "filesystem.hpp"
#include <cstdlib>
#include <memory>
#include <string>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCurvatureFlowImageFilter.h>
#include <itkMinMaxCurvatureFlowImageFilter.h>
#include <itkBinaryMinMaxCurvatureFlowImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkClampImageFilter.h>
#include <itkCastImageFilter.h>
#include <typeinfo>
#include <vector>

namespace fs = std::filesystem;

template <typename InputPixelType, unsigned int Dimension>
int curvature_smoothing_impl(
//# Bugs: doesn't work with rgb images. Converts rgb images to grayscale implicitly
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger) {
    using OutputPixelType = float;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);

    using SmoothingFilterType =
    itk::CurvatureFlowImageFilter<InputImageType, OutputImageType>;

    typename SmoothingFilterType::Pointer smoothingFilter =
        SmoothingFilterType::New();

    if (parser->ArgumentExists("-mmcs")){
        using MinMaxCurvatureFlowFilterType = itk::MinMaxCurvatureFlowImageFilter<InputImageType, OutputImageType>;
        typename MinMaxCurvatureFlowFilterType::Pointer mmSmoothingFilter = MinMaxCurvatureFlowFilterType::New();
        logger->Info("Using min-max Curvature smoothing\n");
        unsigned stencilRadius = 2;
        parser->GetCommandLineArgument("-radius", stencilRadius);
        logger->Info("Set radius to " + std::to_string(stencilRadius) +"\n");
        mmSmoothingFilter->SetStencilRadius(stencilRadius);
        smoothingFilter = mmSmoothingFilter;
    } else if (parser->ArgumentExists("-bmmcs")) {
      using BinaryMinMaxCurvatureFlowFilterType =
          itk::BinaryMinMaxCurvatureFlowImageFilter<InputImageType, OutputImageType>;
      typename BinaryMinMaxCurvatureFlowFilterType::Pointer bmmSmoothingFilter = BinaryMinMaxCurvatureFlowFilterType::New();
      logger->Info("Using binary min-max Curvature smoothing\n");
      unsigned stencilRadius = 2;
      parser->GetCommandLineArgument("-radius", stencilRadius);
      logger->Info("Set radius to " + std::to_string(stencilRadius) + "\n");
      bmmSmoothingFilter->SetStencilRadius(stencilRadius);
      double threshold = 0.5;
      parser->GetCommandLineArgument("-thresh", threshold);
      logger->Info("Set binary image threshold @ " + std::to_string(threshold)+"\n");
      bmmSmoothingFilter->SetThreshold(threshold);
      smoothingFilter = bmmSmoothingFilter;
    } else {
      logger->Info("Using mean Curvature smoothing\n");
    }
    smoothingFilter->SetInput(reader->GetOutput());
    unsigned int numberOfIterations = 1;
    double timeStep = 0.125;
    parser->GetCommandLineArgument("-iterations", numberOfIterations);
    logger->Info("Number of Iteration :: " + std::to_string(numberOfIterations) + "\n");

    parser->GetCommandLineArgument("-dt", timeStep);
    logger->Info("Time Step :: " + std::to_string(timeStep) + "\n");

    smoothingFilter->SetNumberOfIterations(numberOfIterations);
    smoothingFilter->SetTimeStep(timeStep);

    using ClampFilterType = itk::ClampImageFilter<OutputImageType, OutputImageType>;
    auto clampFilter = ClampFilterType::New();
    clampFilter->SetInput(smoothingFilter->GetOutput());

    if(parser->ArgumentExists("-clamp")){
	logger->Debug("Switched output to clamp to input range\n");
        using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<InputImageType>;
        auto imageCalculatorFilter = ImageCalculatorFilterType::New();
        imageCalculatorFilter->SetImage(reader->GetOutput());
	reader->Update();
        imageCalculatorFilter->Compute();
        OutputPixelType lowerBound = static_cast<OutputPixelType>(imageCalculatorFilter->GetMinimum());
	parser->GetCommandLineArgument("-lthresh",lowerBound);
	logger->Debug("Set lower clamp bound = " + std::to_string(lowerBound) + "\n");

	OutputPixelType upperBound = static_cast<OutputPixelType>(imageCalculatorFilter->GetMaximum());
	parser->GetCommandLineArgument("-uthresh",upperBound);
	logger->Debug("Set uppder clamp bound = " + std::to_string(upperBound) + "\n");
	clampFilter->SetBounds(lowerBound, upperBound);
    }
    
    bool rescale = parser->ArgumentExists("-rescale");
    std::string outputFileName;
    parser->GetCommandLineArgument("-output", outputFileName);

    if(rescale){
        using RescaleFilterType =
            itk::RescaleIntensityImageFilter<OutputImageType, OutputImageType>;
        typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
        rescaler->SetOutputMinimum(0);
        rescaler->SetOutputMaximum(1);
        rescaler->SetInput(clampFilter->GetOutput());

        using WriterImageType = OutputImageType;
        using WriterType = itk::ImageFileWriter<WriterImageType>;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(outputFileName);
        writer->SetInput(rescaler->GetOutput());
        logger->Info("Rescale ON\n");
        try {
          writer->Update();
          logger->Info("Wrote Smoothed File to " +
                       std::string(writer->GetFileName()) + "\n");
        }catch(const itk::ExceptionObject &noSuppEx){
            fs::path outpath(outputFileName);
            outpath = outpath.parent_path() / outpath.stem();
            outpath += ".tif";
            writer->SetFileName(outpath.string().c_str());
            try{
                writer->Update();
                logger->Warning("Failed Writing to " + outputFileName + "\n");
                logger->Info("Wrote to" + outpath.string() + " instead\n");
            }catch(const itk::ExceptionObject &nestedEx){
                logger->Warning("Failed Writing to " + outputFileName + "\n");
                logger->Warning("Failed Writing to" + outpath.string() + "\n");
                logger->Critical(nestedEx.what() + std::string("\n"));
                return EXIT_FAILURE;
            }
        }
    }else {
        using WriterImageType = InputImageType;
        using WriterType = itk::ImageFileWriter<WriterImageType>;
        using CastFilterType =
            itk::CastImageFilter<OutputImageType, WriterImageType>;
        typename CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(clampFilter->GetOutput());
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(outputFileName);
        writer->SetInput(castFilter->GetOutput());
        try {
            writer->Update();
            logger->Info("Wrote Smoothed File to " +
                         std::string(writer->GetFileName()) + "\n");
        } catch (const itk::ExceptionObject &excep) {
            logger->Critical("Exception caught !\n");
            logger->Critical(excep.what() + std::string("\n"));
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int curvature_smoothing(const itk::CommandLineArgumentParser::Pointer &parser,
                        const itk::Logger::Pointer &logger) {

    logger->Info("Starting mean curvature smoothing\n");
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
        outputFilePath += "_smooth";
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
            return curvature_smoothing_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return curvature_smoothing_impl<unsigned short, 2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return curvature_smoothing_impl<float, 2>(parser, logger);
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
            return curvature_smoothing_impl<unsigned char, 3>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return curvature_smoothing_impl<unsigned short, 3>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return curvature_smoothing_impl<float, 3>(parser, logger);
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
