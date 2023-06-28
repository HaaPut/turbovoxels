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

#include "slic.h"
#include "filesystem.hpp"
#include <cstdlib>
#include <memory>
#include <string>

#include "utils.h"
#include <algorithm>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkMinMaxCurvatureFlowImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSLICImageFilter.h>
#include <itkScalarToRGBColormapImageFilter.h>
#include <typeinfo>
#include <vector>

namespace fs = std::filesystem;


template <typename InputPixelType, unsigned int Dimension>
int slic_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                    const itk::Logger::Pointer &logger) {
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using RealImageType = itk::Image<float, Dimension>;
    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);

    using slicImageFilterType = itk::SLICImageFilter<InputImageType, InputImageType>;
    typename slicImageFilterType::Pointer slic = slicImageFilterType::New();
    slic->SetInput(reader->GetOutput());
    unsigned iters=10;
    parser->GetCommandLineArgument("-iterations",iters);
    logger->Info("Set Maximum iteration count = " + std::to_string(iters) + "\n");
    slic->SetMaximumNumberOfIterations(iters);
    double weight = 1;
    parser->GetCommandLineArgument("-weight",weight);
    logger->Info("Proximity weight = " + std::to_string(weight) + "\n");
    slic->SetSpatialProximityWeight(weight);


    unsigned seedcount = 100;
    parser->GetCommandLineArgument("-seedcount",seedcount);
    logger->Info("Set Number of regions = " + std::to_string(seedcount) + "\n");

    typename  InputImageType::Pointer input = reader->GetOutput();
    input->Update();
    typename InputImageType::SizeType regionSize = input->GetLargestPossibleRegion().GetSize();
    //typename InputImageType::SpacingType voxelSize = input->GetSpacing();

     typename slicImageFilterType::SuperGridSizeType superpixelSize;
     for (size_t d = 0; d < Dimension; ++d) {
         superpixelSize[d] =
                 std::llround(regionSize[d] / std::pow(seedcount, 1.0 / Dimension));
     }
     slic->SetSuperGridSize(superpixelSize);
     std::stringstream ss;
     ss << superpixelSize;
     logger->Info("Set superpixel grid size = " + ss.str()+ "\n");

    using GradientMagnitudeFilterType =
        itk::GradientMagnitudeImageFilter<InputImageType, RealImageType>;
    typename GradientMagnitudeFilterType::Pointer gradientMag =
        GradientMagnitudeFilterType::New();
    gradientMag->SetInput(slic->GetOutput());

    // Generate initial mask/ seeds for segmentation.
    using BoundaryFilterType =
        itk::BinaryThresholdImageFilter<RealImageType, InputImageType>;
    typename BoundaryFilterType::Pointer boundary =
        BoundaryFilterType::New();
    boundary->SetInput(gradientMag->GetOutput());
    InputPixelType outsideValue = itk::NumericTraits<InputPixelType>::One;
    InputPixelType insideValue = itk::NumericTraits<InputPixelType>::Zero;
    boundary->SetOutsideValue(outsideValue);
    boundary->SetInsideValue(insideValue);
    boundary->SetLowerThreshold(0.01);
    boundary->SetUpperThreshold(itk::NumericTraits<float>::max());

    using MultiplicationFilterType =
        itk::MultiplyImageFilter<InputImageType, InputImageType, InputImageType>;
    typename MultiplicationFilterType::Pointer multiplier = MultiplicationFilterType::New();
    multiplier->SetInput1(reader->GetOutput());
    multiplier->SetInput2(boundary->GetOutput());

    std::string outputFileName;
    parser->GetCommandLineArgument("-output", outputFileName);

    using WriterImageType = InputImageType;
    using WriterType = itk::ImageFileWriter<WriterImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFileName);
    writer->SetInput(multiplier->GetOutput());
    try {
        writer->Update();
        logger->Info("Wrote File to " +
                     std::string(writer->GetFileName()) + "\n");
    } catch (const itk::ExceptionObject &excep) {
        logger->Critical("Error Writing  File!\n");
        logger->Critical(excep.what() + std::string("\n"));
        return EXIT_FAILURE;
    }
    unsigned labelmap = 0;
    parser->GetCommandLineArgument("-writelabels",labelmap);
    if(labelmap > 0){
        logger->Info("Write superpixel label image option on\n");
        fs::path inputFilePath, outputFilePath;
        inputFilePath = inputFileName;
        outputFilePath = inputFilePath.parent_path() / inputFilePath.stem();
        outputFilePath += "_slic_super_pixel_labels";
        outputFilePath += inputFilePath.extension();

        if(labelmap == 1){
            slic->Update();
            writeImage<InputImageType>(slic->GetOutput(),outputFilePath.string(),logger,"Wrote SuperPixel Labels to ");
        }else if(labelmap > 1){
          using RGBPixelType = itk::RGBPixel<unsigned char>;
          using RGBImageType = itk::Image<RGBPixelType, 2>;
          using RGBFilterType =
              itk::ScalarToRGBColormapImageFilter<InputImageType,
                                                  RGBImageType>;
          typename RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
          colormapImageFilter->SetInput(slic->GetOutput());
          if(labelmap == 2)
              colormapImageFilter->SetColormap(RGBFilterType::ColormapEnumType::Hot);
          else if(labelmap == 3)
              colormapImageFilter->SetColormap(RGBFilterType::ColormapEnumType::Cool);
          else if(labelmap == 4)
              colormapImageFilter->SetColormap(RGBFilterType::ColormapEnumType::Autumn);
          else
              colormapImageFilter->SetColormap(RGBFilterType::ColormapEnumType::Jet);

          colormapImageFilter->Update();
          writeImage<RGBImageType>(colormapImageFilter->GetOutput(), outputFilePath.string(),
                                   logger, "Wrote colored Super pixel labels to ");
        }
    }
    return EXIT_SUCCESS;
}

int slic(const itk::CommandLineArgumentParser::Pointer &parser,
                 const itk::Logger::Pointer &logger) {

    logger->Info("Experiment\n");
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
        outputFilePath += "_slic";
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
    if (!fs::exists(inputFileName)) {
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
            return slic_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return slic_impl<unsigned short, 2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return slic_impl<float, 2>(parser, logger);
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
          return slic_impl<unsigned char, 3>(parser, logger);
      case itk::ImageIOBase::USHORT:
          return slic_impl<unsigned short, 3>(parser, logger);
      case itk::ImageIOBase::FLOAT:
          return slic_impl<float, 3>(parser, logger);
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
