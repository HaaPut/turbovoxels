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
// Created by tabish on 2021-03-27.
//

#include "chanvese.h"

#include <cstdlib>
#include "filesystem.hpp"
#include <itkNumericTraits.h>
#include <string>

#include <itkCastImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinMaxCurvatureFlowImageFilter.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLevelSetDenseImage.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryImageToLevelSetImageAdaptor.h>
#include <itkLevelSetEquationChanAndVeseInternalTerm.h>
#include <itkLevelSetEquationChanAndVeseExternalTerm.h>
#include <itkLevelSetEquationContainer.h>
#include <itkLevelSetEquationTermContainer.h>
#include <itkLevelSetEvolution.h>
#include <itkLevelSetContainer.h>
#include <itkLevelSetEvolutionNumberOfIterationsStoppingCriterion.h>
#include <itkSinRegularizedHeavisideStepFunction.h>
#include <typeinfo>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

template <typename InputPixelType, unsigned int Dimension>
int chanvese_segmentation_impl(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger) {
    using OutputPixelType = uint8_t;
    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using InternalImageType = itk::Image<float, Dimension>;
    using LevelSetPixelType = float;
    using LevelSetImageType = itk::Image<LevelSetPixelType, Dimension>;
    using LevelSetType = itk::LevelSetDenseImage<LevelSetImageType>;
    using LevelSetOutputType = typename LevelSetType::OutputType;
    using LevelSetRealType = typename LevelSetType::OutputRealType;

    using ReaderType = itk::ImageFileReader<InputImageType>;

    typename ReaderType::Pointer reader = ReaderType::New();
    std::string inputFileName;
    parser->GetCommandLineArgument("-input", inputFileName);
    reader->SetFileName(inputFileName);
    typename InputImageType::Pointer inputImage;

    using SmoothingFilterType =
    itk::MinMaxCurvatureFlowImageFilter<InputImageType, InternalImageType>;

    bool smooth = parser->ArgumentExists("-smooth");

    if(smooth){
        // Smoothing input image..
        typename SmoothingFilterType::Pointer smoothingFilter =
            SmoothingFilterType::New();
        logger->Info("Pre smoothing using min-max Curvature smoothing\n");
        unsigned stencilRadius = 2;
        parser->GetCommandLineArgument("-radius", stencilRadius);
        logger->Info("Set radius to " + std::to_string(stencilRadius) +"\n");
        smoothingFilter->SetStencilRadius(stencilRadius);
        smoothingFilter->SetInput(reader->GetOutput());
        unsigned int numberOfSmoothingIterations = 10;
        parser->GetCommandLineArgument("-smooth", numberOfSmoothingIterations);
        logger->Info("Number of Smoothing  Iteration :: " + std::to_string(numberOfSmoothingIterations) + "\n");

        double timeStep = 0.125;
        parser->GetCommandLineArgument("-dt", timeStep);
        logger->Info("Time Step :: " + std::to_string(timeStep) + "\n");

        smoothingFilter->SetNumberOfIterations(numberOfSmoothingIterations);
        smoothingFilter->SetTimeStep(timeStep);
        using CastFilterType =
            itk::CastImageFilter<InternalImageType, InputImageType>;
        typename CastFilterType::Pointer castFilter = CastFilterType::New();

        castFilter->SetInput(smoothingFilter->GetOutput());
        // Update to prevent down "PartitionDomainreturned more subdomains than were requested
        // declaring castfilter/smoothing filter outside if block doesn't help
        castFilter->Update();
        inputImage = castFilter->GetOutput();
    }else{
        inputImage = reader->GetOutput();
    }

    // Generate initial mask/ seeds for segmentation.
    using InitLevelSetImageType = itk::Image<LevelSetOutputType, Dimension>;
    using InitLevelSetFilterType =
        itk::BinaryThresholdImageFilter<InputImageType, InitLevelSetImageType>;
    typename InitLevelSetFilterType::Pointer initLevelSet = InitLevelSetFilterType::New();
    initLevelSet->SetInput(reader->GetOutput());
    LevelSetOutputType outsideValue =
        itk::NumericTraits<LevelSetOutputType>::One;
    LevelSetOutputType insideValue =
        itk::NumericTraits<LevelSetOutputType>::Zero;
    if(parser->ArgumentExists("-background")){
        logger->Info("Running in background segmentation mode\n");
        std::swap(insideValue, outsideValue);
    }
    initLevelSet->SetOutsideValue(outsideValue);
    initLevelSet->SetInsideValue(insideValue);

    InputPixelType lowerThreshold = itk::NumericTraits<InputPixelType>::NonpositiveMin();
    InputPixelType upperThreshold = itk::NumericTraits<InputPixelType>::max();
    parser->GetCommandLineArgument("-lthresh", lowerThreshold);
    parser->GetCommandLineArgument("-uthresh", upperThreshold);
    logger->Info("Lower init mask threshold = " + std::to_string(lowerThreshold) + "\n");
    logger->Info("Upper init mask threshold = " + std::to_string(upperThreshold) + "\n");

    lowerThreshold = static_cast<InputPixelType>(lowerThreshold);
    upperThreshold = static_cast<InputPixelType>(upperThreshold);
    initLevelSet->SetLowerThreshold(lowerThreshold);
    initLevelSet->SetUpperThreshold(upperThreshold);

    //Generate LevelSet from initLevelSet mask
    using InitToLevelSetType =
        itk::BinaryImageToLevelSetImageAdaptor<InitLevelSetImageType,
                                               LevelSetType>;
    typename InitToLevelSetType::Pointer adaptor = InitToLevelSetType::New();
    adaptor->SetInputImage(initLevelSet->GetOutput());
    adaptor->Initialize();
    typename LevelSetType::Pointer levelSet = adaptor->GetModifiableLevelSet();


    // The Heaviside function
    using HeavisideFunctionType = itk::SinRegularizedHeavisideStepFunction<LevelSetRealType, LevelSetRealType>;
    typename HeavisideFunctionType::Pointer heaviside = HeavisideFunctionType::New();
    heaviside->SetEpsilon(1.5);

    // Create the level set container
    using LevelSetContainerType = itk::LevelSetContainer<itk::IdentifierType, LevelSetType>;
    typename LevelSetContainerType::Pointer levelSetContainer = LevelSetContainerType::New();
    levelSetContainer->SetHeaviside(heaviside);
    typename LevelSetContainerType::LevelSetIdentifierType levelSetId = 0;
    levelSetContainer->AddLevelSet(levelSetId, levelSet);

    // Create the terms.
    //
    // // Chan and Vese internal term
    using ChanAndVeseInternalTermType =
        itk::LevelSetEquationChanAndVeseInternalTerm<InputImageType, LevelSetContainerType>;
    typename ChanAndVeseInternalTermType::Pointer cvInternalTerm =
        ChanAndVeseInternalTermType::New();
    cvInternalTerm->SetInput(inputImage);
    float lambda = 0.5;
    parser->GetCommandLineArgument("-lambda",lambda);
    logger->Info("Set Weight of internal term(lambda) = " + std::to_string(lambda) + "\n");
    cvInternalTerm->SetCoefficient(lambda);

    // // Chan and Vese external term
    using ChanAndVeseExternalTermType =
        itk::LevelSetEquationChanAndVeseExternalTerm<InputImageType, LevelSetContainerType>;
    typename ChanAndVeseExternalTermType::Pointer cvExternalTerm =
        ChanAndVeseExternalTermType::New();
    cvExternalTerm->SetInput(inputImage);

    // Create term container (equation rhs)
    using TermContainerType = itk::LevelSetEquationTermContainer<InputImageType, LevelSetContainerType>;
    typename TermContainerType::Pointer termContainer = TermContainerType::New();
    termContainer->SetLevelSetContainer(levelSetContainer);
    termContainer->SetInput(inputImage);
    termContainer->AddTerm(0, cvInternalTerm);
    termContainer->AddTerm(1, cvExternalTerm);

    // Create equation container
    using EquationContainerType = itk::LevelSetEquationContainer<TermContainerType>;
    typename EquationContainerType::Pointer equationContainer =
        EquationContainerType::New();
    equationContainer->SetLevelSetContainer(levelSetContainer);
    equationContainer->AddEquation(0, termContainer);

    // Create stopping criteria
    using StoppingCriterionType = itk::LevelSetEvolutionNumberOfIterationsStoppingCriterion<LevelSetContainerType>;
    typename StoppingCriterionType::Pointer criterion =
        StoppingCriterionType::New();
    unsigned numberOfIterations = 1;
    parser->GetCommandLineArgument("-iterations", numberOfIterations);
    logger->Info("Running For " + std::to_string(numberOfIterations)+
                 " iterations\n");
    criterion->SetNumberOfIterations(numberOfIterations);

    // Create evolution class
    using LevelSetEvolutionType = itk::LevelSetEvolution<EquationContainerType, LevelSetType>;
    typename LevelSetEvolutionType::Pointer evolution = LevelSetEvolutionType::New();
    evolution->SetEquationContainer(equationContainer);
    evolution->SetStoppingCriterion(criterion);
    evolution->SetLevelSetContainer(levelSetContainer);
    evolution->Update();
    typename LevelSetType::Pointer finalLevelSet = evolution->GetLevelSetContainer()->GetLevelSet(levelSetId);
    typename LevelSetImageType::Pointer levelSetImage = finalLevelSet->GetModifiableImage();

    using LevelSetThresholdFilter = itk::BinaryThresholdImageFilter<LevelSetImageType, OutputImageType>;
    typename LevelSetThresholdFilter::Pointer segmentationFilter = LevelSetThresholdFilter::New();
    segmentationFilter->SetInput(levelSetImage);
    segmentationFilter->SetInsideValue(itk::NumericTraits<OutputPixelType>::One);
    segmentationFilter->SetOutsideValue(itk::NumericTraits<OutputPixelType>::Zero);
    segmentationFilter->SetLowerThreshold(0);
    segmentationFilter->SetUpperThreshold(itk::NumericTraits<LevelSetOutputType>::max());
    //could be removed...
    using RescaleFilterType =
        itk::RescaleIntensityImageFilter<OutputImageType, OutputImageType>;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(itk::NumericTraits<OutputPixelType>::Zero);
    rescaler->SetOutputMaximum(itk::NumericTraits<OutputPixelType>::max());
    rescaler->SetInput(segmentationFilter->GetOutput());

    std::string outputFileName;
    parser->GetCommandLineArgument("-output", outputFileName);

    using WriterImageType = OutputImageType;
    using WriterType = itk::ImageFileWriter<WriterImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFileName);
    writer->SetInput(rescaler->GetOutput());
    try {
        writer->Update();
        logger->Info("Wrote Segmentation File to " +
                     std::string(writer->GetFileName()) + "\n");
    } catch (const itk::ExceptionObject &excep) {
        logger->Critical("Error Writing Segmentation File!\n");
        logger->Critical(excep.what() + std::string("\n"));
        return EXIT_FAILURE;
    }

    if (parser->ArgumentExists("-debug")){
        using WriterTypeInitSeg = itk::ImageFileWriter<InitLevelSetImageType>;
        typename WriterTypeInitSeg::Pointer writerInitSet =  WriterTypeInitSeg::New();
        fs::path initsavepath(inputFileName);
        initsavepath = initsavepath.parent_path()/initsavepath.stem();
        initsavepath += "_init.tif";
        writerInitSet->SetFileName(initsavepath.string());
        writerInitSet->SetInput(initLevelSet->GetOutput());
        try {
            writerInitSet->Update();
            logger->Info("Wrote Init levelset File to " +
                         std::string(writerInitSet->GetFileName()) + "\n");
        } catch (const itk::ExceptionObject &excep) {
            logger->Critical("Error Writing init File!\n");
            logger->Critical(excep.what() + std::string("\n"));
            return EXIT_FAILURE;
        }

        using WriterTypeLs = itk::ImageFileWriter<LevelSetImageType>;
        typename WriterTypeLs::Pointer writerFinalLS = WriterTypeLs::New();
        fs::path lssavepath(inputFileName);
        lssavepath = lssavepath.parent_path() / lssavepath.stem();
        lssavepath += "_ls.tif";
        writerFinalLS->SetFileName(lssavepath.string());
        writerFinalLS->SetInput(levelSetImage);
        try {
            writerFinalLS->Update();
            logger->Info("Wrote Final levelset File to " +
                         std::string(writerFinalLS->GetFileName()) + "\n");
        } catch (const itk::ExceptionObject &excep) {
            logger->Critical("Error Writing Final File!\n");
            logger->Critical(excep.what() + std::string("\n"));
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}


int chanvese_segmentation(const itk::CommandLineArgumentParser::Pointer &parser,
                           const itk::Logger::Pointer &logger) {

    logger->Info("Starting Chan-vese\n");
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
            return chanvese_segmentation_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return chanvese_segmentation_impl<unsigned short, 2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return chanvese_segmentation_impl<float, 2>(parser, logger);
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
            return chanvese_segmentation_impl<unsigned char, 3>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return chanvese_segmentation_impl<unsigned short, 3>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return chanvese_segmentation_impl<float, 3>(parser, logger);
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
