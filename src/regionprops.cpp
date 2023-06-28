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
// Created by tabish on 2021-04-27
//

#include "itkLabelGeometryImageFilter.h"
#include "utils.h"
#include "regionprops.h"
#include <climits>
#include <cstdlib>
#include "filesystem.hpp"
#include <itkShapeRelabelLabelMapFilter.h>
#include <itkCastImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkLabelShapeOpeningImageFilter.h>
#include <itkNumericTraits.h>
#include <itkShapeOpeningLabelMapFilter.h>
#include <string>
#include <algorithm>
#include <vector>
#include <typeinfo>

namespace fs = std::filesystem;

template <typename InputPixelType, unsigned int Dimension>
int regionprops_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                    const itk::Logger::Pointer &logger) {
    //using OutputPixelType = InputPixelType;
    using InternalPixelType = uint64_t;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    //using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;
    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);
    logger->Info("Input image file " + reader->GetFileName() + "\n");

    using CastFilterType = itk::CastImageFilter<InputImageType, InternalImageType>;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(reader->GetOutput());

    using ConnectedComponentImageFilterType =
        itk::ConnectedComponentImageFilter<InternalImageType, InternalImageType>;
    typename ConnectedComponentImageFilterType::Pointer connectedComponentImageFilter =
        ConnectedComponentImageFilterType::New();
    connectedComponentImageFilter->SetInput(castFilter->GetOutput());
    if (parser->ArgumentExists("-dark")){
      connectedComponentImageFilter->SetBackgroundValue(
          itk::NumericTraits<InternalPixelType>::OneValue());//change to max value in image...
    }else{
      connectedComponentImageFilter->SetBackgroundValue(
          itk::NumericTraits<InternalPixelType>::ZeroValue());
    }
    logger->Info(
        "Set background value to " +
        std::to_string(connectedComponentImageFilter->GetBackgroundValue()) +
        "\n");
    connectedComponentImageFilter->Update();

    using LabelImageToLabelMapFilterType =
        itk::LabelImageToShapeLabelMapFilter<InternalImageType>;
    typename LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter =
        LabelImageToLabelMapFilterType::New();
    labelImageToLabelMapFilter->SetInput(
        connectedComponentImageFilter->GetOutput());
    labelImageToLabelMapFilter->Update();

    using LabelShapeOpeningFilterType = itk::ShapeOpeningLabelMapFilter<
        typename LabelImageToLabelMapFilterType::OutputImageType>;
    typename LabelShapeOpeningFilterType::Pointer shapeFilter = LabelShapeOpeningFilterType::New();
    shapeFilter->SetInput(labelImageToLabelMapFilter->GetOutput());
    shapeFilter->SetAttribute("NumberOfPixels");
    int objectSize = 100;
    parser->GetCommandLineArgument("-size", objectSize);
    objectSize = std::max(objectSize,0);
    logger->Info("Set Object size Threshold = " + std::to_string(objectSize)+"\n");
    shapeFilter->SetLambda(objectSize);
    shapeFilter->Update();
    logger->Info(
        std::to_string(shapeFilter->GetOutput()->GetNumberOfLabelObjects()) + "/"+
        std::to_string(labelImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()) +
        " objects filtered\n");

    using RelabelFilterType = itk::ShapeRelabelLabelMapFilter<typename LabelShapeOpeningFilterType::OutputImageType>;
    typename RelabelFilterType::Pointer  labelFilter = RelabelFilterType::New();
    labelFilter->SetInput(shapeFilter->GetOutput());

    using LabelMapToLabelImageFilterType = itk::LabelMapToLabelImageFilter<
        typename RelabelFilterType::OutputImageType, InternalImageType>;
    typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter =
        LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput(labelFilter->GetOutput());
    labelMapToLabelImageFilter->Update();

    using LabelGeometryImageFilterType =
        itk::LabelGeometryImageFilter<typename LabelMapToLabelImageFilterType::OutputImageType>;
    typename LabelGeometryImageFilterType::Pointer labelGeometryImageFilter =
        LabelGeometryImageFilterType::New();
    labelGeometryImageFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
    labelGeometryImageFilter->SetCalculateOrientedLabelRegions(true);
    labelGeometryImageFilter->Update();
    typename LabelGeometryImageFilterType::LabelsType allLabels =
        labelGeometryImageFilter->GetLabels();
    typename LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
    logger->Info("Number of Objects Analyzed: " +
                 std::to_string(labelGeometryImageFilter->GetNumberOfLabels()) +
                 "\n");

    // Start Writing properties to json file....
    std::string dataFileName;
    parser->GetCommandLineArgument("-output", dataFileName);
    std::ofstream dataFile;
    dataFile.open(dataFileName);
    dataFile <<"{\n\"location\":[\n";
    for (allLabelsIt = allLabels.begin();
         allLabelsIt != allLabels.end();
         ++allLabelsIt) {
        typename LabelGeometryImageFilterType::LabelPixelType labelValue =
            *allLabelsIt;
        auto location = labelGeometryImageFilter->GetCentroid(labelValue);
        dataFile << "  [\n";
        dataFile << "    "<<location[0]<< ",\n";
        dataFile << "    "<<location[1]<< ",\n";
        if (Dimension == 2) dataFile << "    0\n  ]";
        else dataFile << "    " << location[2] << "\n  ]";
        if(allLabelsIt != allLabels.end()-1) dataFile << ",\n";
        else dataFile << "\n";
    }
    dataFile << "  ],\n"; //finish location

    dataFile << "\"length\":[\n";
    for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end();
         ++allLabelsIt) {
      typename LabelGeometryImageFilterType::LabelPixelType labelValue =
          *allLabelsIt;
      auto length = labelGeometryImageFilter->GetMajorAxisLength(labelValue);
      dataFile << "  [\n";
      dataFile << "    " << length << "\n  ]";
      if (allLabelsIt != allLabels.end() - 1)
        dataFile << ",\n";
      else
        dataFile << "\n";
    }
    dataFile << "  ],\n"; // finsh major axis length

    dataFile << "\"rotation\":[\n";
    for (allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end();
         ++allLabelsIt) {
      typename LabelGeometryImageFilterType::LabelPixelType labelValue =
          *allLabelsIt;
      auto rot = labelGeometryImageFilter->GetRotationMatrix(labelValue);
      dataFile << "  [\n"; //begin rot matrix
      if(Dimension == 3){
          for(size_t r = 0; r < 3; ++r){
              dataFile << "    [\n";
              for(size_t c = 0 ; c < 3; ++c){
                  dataFile << "      " << rot[r][c] << ",\n";
              }
              dataFile << "      " << "0\n";
              dataFile << "    ],\n";
          }
      }
      else{
          dataFile << "    [\n";
          dataFile << "      "<< rot[0][0]<<",\n";
          dataFile << "      "<<rot[0][1] <<",\n";
          dataFile << "      0,\n";
          dataFile << "      0\n";
          dataFile << "    ],\n";

          dataFile << "    [\n";
          dataFile << "      " << rot[1][0] << ",\n";
          dataFile << "      " << rot[1][1] << ",\n";
          dataFile << "      0,\n";
          dataFile << "      0\n";
          dataFile << "    ],\n";

          dataFile << "    [\n";
          dataFile << "      0,\n";
          dataFile << "      0,\n";
          dataFile << "      1,\n";
          dataFile << "      0\n";
          dataFile << "    ],\n";
      }
      dataFile << "    [\n";
      dataFile << "      0,\n";
      dataFile << "      0,\n";
      dataFile << "      0,\n";
      dataFile << "      1\n";
      dataFile << "    ]\n";
      dataFile << "  ]";
      if (allLabelsIt != allLabels.end() - 1)
        dataFile << ",\n";
      else
        dataFile << "\n";
      //end rot matrix
    }
    dataFile << "  ]\n"; //finish rotation

    dataFile << "}"; //finsh file

    dataFile.close();
    logger->Info("Wrote Object Properties File to " + dataFileName + "\n");

    return EXIT_SUCCESS;
}

int regionprops(const itk::CommandLineArgumentParser::Pointer &parser,
                const itk::Logger::Pointer &logger) {

    logger->Info("Starting Region Props\n");
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
        outputFilePath += "_props";
        outputFilePath += ".json";
        std::vector<std::string> newargs;
        newargs.push_back("-output");
        newargs.push_back(outputFilePath.string());
        parser->AddCommandLineArguments(newargs);
    }else{
        outputFilePath = outputFileName;
    }
    logger->Info("Input Filename : " + inputFilePath.string() + "\n");
    logger->Info("Output Filename : "  + outputFilePath.string() + "\n");

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
            return regionprops_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return regionprops_impl<unsigned short, 2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return regionprops_impl<float, 2>(parser, logger);
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
        return regionprops_impl<unsigned char, 3>(parser, logger);
      case itk::ImageIOBase::USHORT:
        return regionprops_impl<unsigned short, 3>(parser, logger);
      case itk::ImageIOBase::FLOAT:
        return regionprops_impl<float, 3>(parser, logger);
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
