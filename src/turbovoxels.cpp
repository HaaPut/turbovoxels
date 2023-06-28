//**********************************************************
//Copyright 2022 Tabish Syed
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
// Created by tabish on 07/10/22.
//

#include <itkLogger.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkScalarToRGBColormapImageFilter.h>
#include <itkRGBToLuminanceImageFilter.h>
#include <itkRGBPixel.h>

#include "turbovoxels.h"
#include "itkCommandLineArgumentParser.h"
#include "itkGrayToRGBColormapFunction.h"
#include "levelSetEvolution.h"
#include "boundaries.h"

template class turboVoxel<2>;
template class turboVoxel<3>;


template <typename InputPixelType, unsigned int Dimension>
int turbovoxels_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                    const itk::Logger::Pointer &logger) {
    logger->Info("Running Grayscale turbovoxels\n");
    using RealPixelType = float;
    using LabelPixelType = unsigned char;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using RealImageType = itk::Image<RealPixelType, Dimension>;
    using LabelImageType = itk::Image<LabelPixelType, Dimension>;

    using ReaderType = itk::ImageFileReader<InputImageType>;

    using LevelSetPixelType = float;
    using LevelSetImageType = itk::Image<LevelSetPixelType, Dimension>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);

    reader->SetFileName(inputFileName);
    logger->Info("Input image file " + reader->GetFileName() + "\n");
    using RescaleFilterType = itk::RescaleIntensityImageFilter<InputImageType, RealImageType>;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(reader->GetOutput());
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(1);
    rescaleFilter->Update();
    typename RealImageType::Pointer input = rescaleFilter->GetOutput();

    if (parser->ArgumentExists("-debug")) {
        fs::path inputSavepath(inputFileName);
        inputSavepath = inputSavepath.parent_path() / inputSavepath.stem();
        inputSavepath += "_turbo_input.tif";
        writeImage<LevelSetImageType>(input, inputSavepath.string(), logger,
                                      "Wrote Final Skeleton Image to " +
                                      inputSavepath.string() + "\n");
    }

    typename LevelSetImageType::Pointer levelSetImage;

    logger->Info("Started LevelSet Evolution\n");
    levelSetImage = maskedSparseLevelSetEvolution<RealImageType>(input, parser, logger);

    if (parser->ArgumentExists("-debug")) {
        fs::path finalLsSavepath(inputFileName);
        finalLsSavepath = finalLsSavepath.parent_path() / finalLsSavepath.stem();
        finalLsSavepath += "_turbo_ls.tif";
        writeImage<LevelSetImageType>(levelSetImage, finalLsSavepath.string(), logger,
                                      "Wrote post evolution Level set to " +
                                      finalLsSavepath.string() + "\n");
    }


    typename LabelImageType::Pointer turboVoxelBoundaries =
            levelSetImageToBoundary<LevelSetImageType, LabelImageType>(levelSetImage, parser, logger);

    Color red(3,0);
    red[0] = 1;
    using RGBImageType = itk::Image<RGBPixelType,Dimension>;
    using GrayToRGBMapType = itk::Function::GrayToRGBColormapFunction<InputPixelType, RGBPixelType>;
    typename GrayToRGBMapType::Pointer grayToRGBMap = GrayToRGBMapType::New();
    //grayToRGBMap->SetMinimumInputValue(0);
    //grayToRGBMap->SetMaximumInputValue(255);
    using ColormapFilterType =
            itk::ScalarToRGBColormapImageFilter<InputImageType, RGBImageType>;
    typename ColormapFilterType::Pointer grayToRGBFilter = ColormapFilterType::New();
    grayToRGBFilter->SetInput(reader->GetOutput());
    typename RGBImageType::Pointer coloredInput = grayToRGBFilter->GetOutput();
    coloredInput->Update();
    paintImage<RGBImageType, LabelImageType>(coloredInput, turboVoxelBoundaries, red);
    fs::path turboVoxelsSavepath(inputFileName);
    turboVoxelsSavepath = turboVoxelsSavepath.parent_path() / turboVoxelsSavepath.stem();
    turboVoxelsSavepath += "_turbo_voxels.tif";
    writeImage<RGBImageType>(coloredInput, turboVoxelsSavepath.string(), logger,
                                  "Wrote overlaid Turbo voxels to " +
                                  turboVoxelsSavepath.string() + "\n");
    return EXIT_SUCCESS;
}

template <typename InputPixelValueType, unsigned int Dimension>
int turbovoxels_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                     const itk::Logger::Pointer &logger, itk::RGBPixel<InputPixelValueType> p) {
    logger->Info("Running RGB turbovoxels\n");
    using InputPixelType = itk::RGBPixel<InputPixelValueType>;
    using RealPixelType = float;
    using LabelPixelType = unsigned char;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using RealImageType = itk::Image<RealPixelType, Dimension>;
    using LabelImageType = itk::Image<LabelPixelType, Dimension>;

    using ReaderType = itk::ImageFileReader<InputImageType>;

    using LevelSetPixelType = float;
    using LevelSetImageType = itk::Image<LevelSetPixelType, Dimension>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);

    reader->SetFileName(inputFileName);
    logger->Info("Input image file " + reader->GetFileName() + "\n");
    using RGB2GrayFilterType = itk::RGBToLuminanceImageFilter<InputImageType, RealImageType>;
    typename RGB2GrayFilterType::Pointer rgb2grayFilter = RGB2GrayFilterType::New();
    rgb2grayFilter->SetInput(reader->GetOutput());

    using RescaleFilterType = itk::RescaleIntensityImageFilter<RealImageType, RealImageType>;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(rgb2grayFilter->GetOutput());
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(1);
    rescaleFilter->Update();
    typename RealImageType::Pointer input = rescaleFilter->GetOutput();

    typename LevelSetImageType::Pointer levelSetImage;

    logger->Info("Started LevelSet Evolution\n");
    levelSetImage = maskedSparseLevelSetEvolution<RealImageType>(input, parser, logger);

    typename LabelImageType::Pointer turboVoxelBoundaries =
            levelSetImageToBoundary<LevelSetImageType, LabelImageType>(levelSetImage, parser, logger);

    Color red(3,0);
    red[0] = 1;
    typename InputImageType::Pointer coloredInput = reader->GetOutput();
    coloredInput->Update();
    paintImage<InputImageType, LabelImageType>(coloredInput, turboVoxelBoundaries, red);
    fs::path turboVoxelsSavepath(inputFileName);
    turboVoxelsSavepath = turboVoxelsSavepath.parent_path() / turboVoxelsSavepath.stem();
    turboVoxelsSavepath += "_turbo_voxels.tif";
    writeImage<InputImageType>(coloredInput, turboVoxelsSavepath.string(), logger,
                             "Wrote overlaid Turbo voxels to " +
                             turboVoxelsSavepath.string() + "\n");

    return EXIT_SUCCESS;
}

int turbovoxels(const itk::CommandLineArgumentParser::Pointer &parser,
               const itk::Logger::Pointer &logger) {

    logger->Info("Turbo Voxels\n");
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
        outputFilePath += "_turbovoxels";
        outputFilePath += inputFilePath.extension();
        std::vector<std::string> newargs;
        newargs.emplace_back("-output");
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
            if(pixelType == itk::ImageIOBase::RGB){
                switch (componentType) {
                    case itk::ImageIOBase::UCHAR:
                        return turbovoxels_impl<unsigned char, 2>(parser, logger, itk::RGBPixel<unsigned char>());
                    default:
                        logger->Critical("For RGB only unsigned char is supported " +
                                         std::string(typeid(componentType).name()) +
                                         " in 2 dimensions\n");
                }
            }else{
                switch (componentType) {
                    case itk::ImageIOBase::UCHAR:
                        return turbovoxels_impl<unsigned char, 2>(parser, logger);
                    case itk::ImageIOBase::USHORT:
                        return turbovoxels_impl<unsigned short, 2>(parser, logger);
                    case itk::ImageIOBase::FLOAT:
                        return turbovoxels_impl<float, 2>(parser, logger);
                    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                    default:
                        logger->Critical("Unknown ComponentType :: " +
                                         std::string(typeid(componentType).name()) +
                                         " in 2 dimensions\n");
                }
            }
            break;
        case 3:
            switch (componentType) {
                case itk::ImageIOBase::UCHAR:
                    return turbovoxels_impl<unsigned char, 3>(parser, logger);
                case itk::ImageIOBase::USHORT:
                    return turbovoxels_impl<unsigned short, 3>(parser, logger);
                case itk::ImageIOBase::FLOAT:
                    return turbovoxels_impl<float, 3>(parser, logger);
                case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                default:
                    logger->Critical("Unknown ComponentType :: " +
                                     std::string(typeid(componentType).name()) +
                                     " in 3 dimensions\n");
            }
            break;
        default:
            logger->Critical("File not supported..\n");
    }
    return EXIT_FAILURE;
}
