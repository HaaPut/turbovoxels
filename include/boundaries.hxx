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

#ifndef TURBO_BOUNDARIES_HXX
#define TURBO_BOUNDARIES_HXX

#include "boundaries.h"
#include "utils.h"

#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryPruningImageFilter.h>
#include <itkBinaryPruningImageFilterMy.h>
#include <itkZeroCrossingImageFilter.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkSkeletonImageFilter.h>

#include <algorithm>
#include <limits>
#include <string>
#include <itkOrImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include "filesystem.hpp"

namespace fs = std::filesystem;

template<typename LevelSetImageType, typename LabelImageType>
typename LabelImageType::Pointer
levelSetImageToBoundary( typename LevelSetImageType::Pointer levelSetImage,
                         const itk::CommandLineArgumentParser::Pointer &parser,
                         const itk::Logger::Pointer &logger){
    logger->Info("Starting LevelSet image to boundary skeleton computation\n");
    using LevelSetPixelType = typename LevelSetImageType::PixelType;
    using LabelPixelType = typename LabelImageType::PixelType;
    constexpr unsigned Dimension  = LevelSetImageType::ImageDimension;
    using FloatToLabelThresholdFilterType = itk::BinaryThresholdImageFilter<LevelSetImageType, LabelImageType>;
    typename FloatToLabelThresholdFilterType::Pointer thresholdFilter = FloatToLabelThresholdFilterType::New();
    thresholdFilter->SetInput(levelSetImage);
    thresholdFilter
        ->SetLowerThreshold(itk::NumericTraits<LevelSetPixelType>::Zero);
    thresholdFilter->SetUpperThreshold(
        itk::NumericTraits<LevelSetPixelType>::max());
    thresholdFilter->SetInsideValue(
        itk::NumericTraits<LabelPixelType>::One);
    thresholdFilter->SetOutsideValue(
        itk::NumericTraits<LabelPixelType>::Zero);
    // typename InternalImageType::Pointer backgroundb =
    // thresholdFilter->GetOutput();

    using ZeroCrossingFilterType =
        itk::ZeroCrossingImageFilter<LevelSetImageType,
                                     LabelImageType>;
    typename ZeroCrossingFilterType::Pointer contourFilter =
        ZeroCrossingFilterType::New();
    contourFilter->SetInput(levelSetImage);
    contourFilter->SetBackgroundValue(
        itk::NumericTraits<LabelPixelType>::Zero);
    contourFilter->SetForegroundValue(
        itk::NumericTraits<LabelPixelType>::One);
    typename LabelImageType::Pointer contour =
        contourFilter->GetOutput();

    //Performs bitwise or of two inputs
    using OrImageFilterType = itk::OrImageFilter<LabelImageType>;
    typename OrImageFilterType::Pointer backgroundPreFilter =
            OrImageFilterType::New();
    backgroundPreFilter->SetInput1(contourFilter->GetOutput());
    backgroundPreFilter->SetInput2(thresholdFilter->GetOutput());

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
    typename BinaryThresholdFilterType::Pointer backgroundFilter =
            BinaryThresholdFilterType::New();
    backgroundFilter->SetInput(backgroundPreFilter->GetOutput());
    backgroundFilter->SetLowerThreshold( itk::NumericTraits<LabelPixelType>::OneValue());
    backgroundFilter->SetUpperThreshold(itk::NumericTraits<LabelPixelType>::max());
    backgroundFilter->SetInsideValue(itk::NumericTraits<LabelPixelType>::max());//itk::NumericTraits<InternalPixelType>::OneValue());
    backgroundFilter->SetOutsideValue(itk::NumericTraits<LabelPixelType>::ZeroValue());

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path backgroundSavepath(inputFileName);
        backgroundSavepath =
                backgroundSavepath.parent_path() / backgroundSavepath.stem();
        backgroundSavepath += "_background.tif";
        typename LabelImageType::Pointer background = backgroundFilter->GetOutput();
        background->Update();
        writeImage<LabelImageType>(
                background, backgroundSavepath.string(), logger,
                "Wrote background File to " +
                std::string(backgroundSavepath.string()) + "\n");
    }
    // This may be unecessary..
    using HoleFillFilterType =
        itk::VotingBinaryIterativeHoleFillingImageFilter<LabelImageType>;
    typename HoleFillFilterType::InputSizeType radius;
    radius.Fill(1);
    typename HoleFillFilterType::Pointer fillInsideFilter =
        HoleFillFilterType::New();
    fillInsideFilter->SetInput(backgroundFilter->GetOutput());
    fillInsideFilter->SetRadius(radius);
    // filter->SetMajorityThreshold(majorityThreshold);
    fillInsideFilter->SetBackgroundValue(backgroundFilter->GetOutsideValue());
    fillInsideFilter->SetForegroundValue(backgroundFilter->GetInsideValue());
    fillInsideFilter->SetMaximumNumberOfIterations(1);

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);
        filledBackgroundSavepath =
                filledBackgroundSavepath.parent_path() / filledBackgroundSavepath.stem();
        filledBackgroundSavepath += "_back_fill_holes.tif";
        fillInsideFilter->Update();
        writeImage<LabelImageType>(
                fillInsideFilter->GetOutput(), filledBackgroundSavepath.string(), logger,
                "Wrote background hole filled File to " +
                std::string(filledBackgroundSavepath.string()) + "\n");
    }

    using CleanBinaryImageFilterType = itk::BinaryShapeOpeningImageFilter<LabelImageType>;
    typename CleanBinaryImageFilterType::Pointer clearBackgroundFilter = CleanBinaryImageFilterType::New();
    clearBackgroundFilter->SetAttribute("NumberOfPixels");
    clearBackgroundFilter->SetLambda(8);
    clearBackgroundFilter->FullyConnectedOn();
    clearBackgroundFilter->SetInput(fillInsideFilter->GetOutput());//backgroundFilter->GetOutput());
    clearBackgroundFilter->SetForegroundValue(fillInsideFilter->GetForegroundValue());
    clearBackgroundFilter->SetBackgroundValue(fillInsideFilter->GetBackgroundValue());

    auto shape = levelSetImage->GetLargestPossibleRegion().GetSize();
    float size = 1;
    float sizeThreshold = 2e-4;
    parser->GetCommandLineArgument("-sizeThreshold", sizeThreshold);

    logger->Info("Set super voxel fraction threshold = " + std::to_string(sizeThreshold) + "\n");
    if(sizeThreshold >= 1){
        logger->Warning("Size Threshold is the fraction of size of the image and meaningful values should be less than 1\n");
    }

    for(size_t d = 0; d < Dimension; ++d) size *= shape[d];
    size *= sizeThreshold;

    logger->Info("Set super voxel size threshold = " + std::to_string(size) + "\n");

    using RemoveSmallSuperPixelType = itk::BinaryShapeOpeningImageFilter<LabelImageType>;
    typename RemoveSmallSuperPixelType::Pointer badSuperPixelFilter = RemoveSmallSuperPixelType::New();
    badSuperPixelFilter->SetAttribute("NumberOfPixels");
    badSuperPixelFilter->SetLambda(size);
    badSuperPixelFilter->FullyConnectedOff();
    badSuperPixelFilter->SetInput(clearBackgroundFilter->GetOutput());//backgroundFilter->GetOutput());
    badSuperPixelFilter->SetForegroundValue(itk::NumericTraits<LabelPixelType>::ZeroValue());
    badSuperPixelFilter->SetBackgroundValue(itk::NumericTraits<LabelPixelType>::OneValue());

    auto background_clean = badSuperPixelFilter->GetOutput();

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);

        fs::path backgroundCleanSavepath(inputFileName);
        backgroundCleanSavepath = backgroundCleanSavepath.parent_path() / backgroundCleanSavepath.stem();
        backgroundCleanSavepath += "_background_clean.tif";
        background_clean->Update();
        writeImage<LabelImageType>(background_clean, backgroundCleanSavepath.string(), logger,
                                      "Wrote background_clean File to " +
                                      std::string(backgroundCleanSavepath.string()) +
                                      "\n");
    }

    using ThinningFilterType = itk::SkeletonImageFilter<LabelPixelType, Dimension>;//itk::BinaryThinningImageFilter<InternalImageType, InternalImageType>;
    typename ThinningFilterType::Pointer skeletonizeBackgroundFilter = ThinningFilterType::New();
    skeletonizeBackgroundFilter->SetInput(badSuperPixelFilter->GetOutput());
    skeletonizeBackgroundFilter->SetLowerThreshold(backgroundFilter->GetInsideValue());
    skeletonizeBackgroundFilter->SetMaxIterations(std::numeric_limits<unsigned>::max());
    auto skeleton_unclean = skeletonizeBackgroundFilter->GetOutput();

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);
        fs::path skeletonUncleanSavepath(inputFileName);
        skeletonUncleanSavepath =
                skeletonUncleanSavepath.parent_path() / skeletonUncleanSavepath.stem();
        skeletonUncleanSavepath += "_preskelton.tif";
        skeleton_unclean->Update();
        writeImage<typename ThinningFilterType::OutputImageType>(
                skeleton_unclean, skeletonUncleanSavepath.string(), logger,
                "Wrote unclean Skeleton File to " + std::string(skeletonUncleanSavepath.string()) +
                "\n");

    }

    using PrunningFilterType =
        itk::BinaryPruningImageFilterMy<LabelImageType, LabelImageType>;
    typename PrunningFilterType::Pointer pruneSkeletonFilter =
        PrunningFilterType::New();
    pruneSkeletonFilter->SetInput(skeletonizeBackgroundFilter->GetOutput());
    pruneSkeletonFilter->SetForegroundValue(skeletonizeBackgroundFilter->GetInsideValue());
    pruneSkeletonFilter->SetBackgroundValue(skeletonizeBackgroundFilter->GetOutsideValue());
    typename PrunningFilterType::OutputImageType::Pointer skeleton = pruneSkeletonFilter->GetOutput();
    skeleton->Update();
    skeleton->DisconnectPipeline();

    if(parser->ArgumentExists("-debug")){
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);
        fs::path skelSavepath(inputFileName);
        skelSavepath = skelSavepath.parent_path() / skelSavepath.stem();
        skelSavepath += "_pruned_skeleton.tif";
        writeImage<LabelImageType>(skeleton, skelSavepath.string(), logger,
                                      "Wrote Final Skeleton Image to " +
                                      skelSavepath.string() + "\n");

    }
    return skeleton;
}

template<typename LevelSetImageType, typename LabelImageType>
typename LabelImageType::Pointer
levelSetImageToTurboVoxels( typename LevelSetImageType::Pointer levelSetImage,
                         const itk::CommandLineArgumentParser::Pointer &parser,
                         const itk::Logger::Pointer &logger){
    logger->Info("Starting LevelSet image to boundary skeleton computation\n");
    using LevelSetPixelType = typename LevelSetImageType::PixelType;
    using LabelPixelType = typename LabelImageType::PixelType;
    constexpr unsigned Dimension  = LevelSetImageType::ImageDimension;
    using FloatToLabelThresholdFilterType = itk::BinaryThresholdImageFilter<LevelSetImageType, LabelImageType>;
    typename FloatToLabelThresholdFilterType::Pointer thresholdFilter = FloatToLabelThresholdFilterType::New();
    thresholdFilter->SetInput(levelSetImage);
    thresholdFilter
            ->SetLowerThreshold(itk::NumericTraits<LevelSetPixelType>::Zero);
    thresholdFilter->SetUpperThreshold(
            itk::NumericTraits<LevelSetPixelType>::max());
    thresholdFilter->SetInsideValue(
            itk::NumericTraits<LabelPixelType>::One);
    thresholdFilter->SetOutsideValue(
            itk::NumericTraits<LabelPixelType>::Zero);
    // typename InternalImageType::Pointer backgroundb =
    // thresholdFilter->GetOutput();

    using ZeroCrossingFilterType =
            itk::ZeroCrossingImageFilter<LevelSetImageType,
                    LabelImageType>;
    typename ZeroCrossingFilterType::Pointer contourFilter =
            ZeroCrossingFilterType::New();
    contourFilter->SetInput(levelSetImage);
    contourFilter->SetBackgroundValue(
            itk::NumericTraits<LabelPixelType>::Zero);
    contourFilter->SetForegroundValue(
            itk::NumericTraits<LabelPixelType>::One);
    typename LabelImageType::Pointer contour =
            contourFilter->GetOutput();

    //Performs bitwise or of two inputs
    using OrImageFilterType = itk::OrImageFilter<LabelImageType>;
    typename OrImageFilterType::Pointer backgroundPreFilter =
            OrImageFilterType::New();
    backgroundPreFilter->SetInput1(contourFilter->GetOutput());
    backgroundPreFilter->SetInput2(thresholdFilter->GetOutput());

    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
    typename BinaryThresholdFilterType::Pointer backgroundFilter =
            BinaryThresholdFilterType::New();
    backgroundFilter->SetInput(backgroundPreFilter->GetOutput());
    backgroundFilter->SetLowerThreshold( itk::NumericTraits<LabelPixelType>::OneValue());
    backgroundFilter->SetUpperThreshold(itk::NumericTraits<LabelPixelType>::max());
    backgroundFilter->SetInsideValue(itk::NumericTraits<LabelPixelType>::max());//itk::NumericTraits<InternalPixelType>::OneValue());
    backgroundFilter->SetOutsideValue(itk::NumericTraits<LabelPixelType>::ZeroValue());

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path backgroundSavepath(inputFileName);
        backgroundSavepath =
                backgroundSavepath.parent_path() / backgroundSavepath.stem();
        backgroundSavepath += "_background.tif";
        typename LabelImageType::Pointer background = backgroundFilter->GetOutput();
        background->Update();
        writeImage<LabelImageType>(
                background, backgroundSavepath.string(), logger,
                "Wrote background File to " +
                std::string(backgroundSavepath.string()) + "\n");
    }
    // This may be unecessary..
    using HoleFillFilterType =
            itk::VotingBinaryIterativeHoleFillingImageFilter<LabelImageType>;
    typename HoleFillFilterType::InputSizeType radius;
    radius.Fill(1);
    typename HoleFillFilterType::Pointer fillInsideFilter =
            HoleFillFilterType::New();
    fillInsideFilter->SetInput(backgroundFilter->GetOutput());
    fillInsideFilter->SetRadius(radius);
    // filter->SetMajorityThreshold(majorityThreshold);
    fillInsideFilter->SetBackgroundValue(backgroundFilter->GetOutsideValue());
    fillInsideFilter->SetForegroundValue(backgroundFilter->GetInsideValue());
    fillInsideFilter->SetMaximumNumberOfIterations(1);

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);
        filledBackgroundSavepath =
                filledBackgroundSavepath.parent_path() / filledBackgroundSavepath.stem();
        filledBackgroundSavepath += "_back_fill_holes.tif";
        fillInsideFilter->Update();
        writeImage<LabelImageType>(
                fillInsideFilter->GetOutput(), filledBackgroundSavepath.string(), logger,
                "Wrote background hole filled File to " +
                std::string(filledBackgroundSavepath.string()) + "\n");
    }

    using CleanBinaryImageFilterType = itk::BinaryShapeOpeningImageFilter<LabelImageType>;
    typename CleanBinaryImageFilterType::Pointer clearBackgroundFilter = CleanBinaryImageFilterType::New();
    clearBackgroundFilter->SetAttribute("NumberOfPixels");
    clearBackgroundFilter->SetLambda(8);
    clearBackgroundFilter->FullyConnectedOn();
    clearBackgroundFilter->SetInput(fillInsideFilter->GetOutput());//backgroundFilter->GetOutput());
    clearBackgroundFilter->SetForegroundValue(fillInsideFilter->GetForegroundValue());
    clearBackgroundFilter->SetBackgroundValue(fillInsideFilter->GetBackgroundValue());

    auto shape = levelSetImage->GetLargestPossibleRegion().GetSize();
    float size = 1;
    float sizeThreshold = 2e-4;
    parser->GetCommandLineArgument("-sizeThreshold", sizeThreshold);

    logger->Info("Set super voxel fraction threshold = " + std::to_string(sizeThreshold) + "\n");
    if(sizeThreshold >= 1){
        logger->Warning("Size Threshold is the fraction of size of the image and meaningful values should be less than 1\n");
    }

    for(size_t d = 0; d < Dimension; ++d) size *= shape[d];
    size *= sizeThreshold;

    logger->Info("Set super voxel size threshold = " + std::to_string(size) + "\n");

    using RemoveSmallSuperPixelType = itk::BinaryShapeOpeningImageFilter<LabelImageType>;
    typename RemoveSmallSuperPixelType::Pointer badSuperPixelFilter = RemoveSmallSuperPixelType::New();
    badSuperPixelFilter->SetAttribute("NumberOfPixels");
    badSuperPixelFilter->SetLambda(size);
    badSuperPixelFilter->FullyConnectedOff();
    badSuperPixelFilter->SetInput(clearBackgroundFilter->GetOutput());//backgroundFilter->GetOutput());
    badSuperPixelFilter->SetForegroundValue(itk::NumericTraits<LabelPixelType>::ZeroValue());
    badSuperPixelFilter->SetBackgroundValue(itk::NumericTraits<LabelPixelType>::OneValue());

    auto background_clean = badSuperPixelFilter->GetOutput();

    if(parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path filledBackgroundSavepath(inputFileName);

        fs::path backgroundCleanSavepath(inputFileName);
        backgroundCleanSavepath = backgroundCleanSavepath.parent_path() / backgroundCleanSavepath.stem();
        backgroundCleanSavepath += "_background_clean.tif";
        background_clean->Update();
        writeImage<LabelImageType>(background_clean, backgroundCleanSavepath.string(), logger,
                                   "Wrote background_clean File to " +
                                   std::string(backgroundCleanSavepath.string()) +
                                   "\n");
    }

    using ThinningFilterType = itk::SkeletonImageFilter<LabelPixelType, Dimension>;//itk::BinaryThinningImageFilter<InternalImageType, InternalImageType>;
    typename ThinningFilterType::Pointer skeletonizeBackgroundFilter = ThinningFilterType::New();
    skeletonizeBackgroundFilter->SetInput(badSuperPixelFilter->GetOutput());
    skeletonizeBackgroundFilter->SetLowerThreshold(backgroundFilter->GetInsideValue());
    skeletonizeBackgroundFilter->SetMaxIterations(std::numeric_limits<unsigned>::max());
    //auto skeleton_unclean = skeletonizeBackgroundFilter->GetOutput();


//    auto invertIntensityFilter = InvertIntensityImageFilterType::New();
//    invertIntensityFilter->SetInput(skeleton_unclean);
//    invertIntensityFilter->SetMaximum(skeletonizeBackgroundFilter->GetInsideValue());
//
//    using CCImageFilterType = itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>;
//    typename CCImageFilterType::Pointer ccImageFilter = CCImageFilterType::New();
//    ccImageFilter->SetInput(invertIntensityFilter->GetOutput());
//    ccImageFilter->FullyConnectedOff();
//    ccImageFilter->Update();

    //This won't work for 3D spur removal
    using PrunningFilterType =
            itk::BinaryPruningImageFilterMy<LabelImageType, LabelImageType>;
    typename PrunningFilterType::Pointer pruneSkeletonFilter =
            PrunningFilterType::New();
    pruneSkeletonFilter->SetInput(skeletonizeBackgroundFilter->GetOutput());
    pruneSkeletonFilter->SetForegroundValue(skeletonizeBackgroundFilter->GetInsideValue());
    pruneSkeletonFilter->SetBackgroundValue(skeletonizeBackgroundFilter->GetOutsideValue());
    //typename PrunningFilterType::OutputImageType::Pointer skeleton = pruneSkeletonFilter->GetOutput();
    //skeleton->Update();
    //skeleton->DisconnectPipeline();

    using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<LabelImageType>;
    auto invertIntensityFilter = InvertIntensityImageFilterType::New();
    invertIntensityFilter->SetInput(pruneSkeletonFilter->GetOutput());
    invertIntensityFilter->SetMaximum(pruneSkeletonFilter->GetForegroundValue());

    using CCImageFilterType = itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>;
    typename CCImageFilterType::Pointer ccImageFilter = CCImageFilterType::New();
    ccImageFilter->SetInput(invertIntensityFilter->GetOutput());
    ccImageFilter->FullyConnectedOff();
    ccImageFilter->Update();
    typename LabelImageType::Pointer turboVoxels = ccImageFilter->GetOutput();

    return turboVoxels;
}

#endif // TURBO_BOUNDARIES_HXX
