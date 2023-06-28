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
#ifndef TURBO_UTILS_HXX
#define TURBO_UTILS_HXX

#include "itkCommandLineArgumentParser.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "filesystem.hpp"
#include <itkAddImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryImageToLevelSetImageAdaptor.h>
#include <itkCurvatureFlowImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkExpImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkLogger.h>
#include <itkMultiplyImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <sstream>
#include <string>
#include <vector>
#include <itkLabelMapToLabelImageFilter.h>
#include <limits>

namespace fs = std::filesystem;

template <typename OutputImageType>
int writeImage(typename OutputImageType::Pointer image,
               const std::string outputFileName,
               const itk::Logger::Pointer &logger, const std::string message) {
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  writer->SetInput(image);
  writer->ReleaseDataFlagOff();
  //writer->PushBackInput(image);
  try {
    writer->Update();
    logger->Info(message + std::string(writer->GetFileName()) + "\n");
  } catch (const itk::ExceptionObject &excep) {
    logger->Critical("Error Writing  File!\n");
    logger->Critical(excep.what() + std::string("\n"));
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template<typename InputImageType>
void generateCriterionSeeds(std::vector<typename InputImageType::IndexType>& seeds,
                            size_t d, typename InputImageType::IndexType seed,
                            typename InputImageType::SizeType cellSize,
                            typename InputImageType::Pointer& image,
                            double criterionThreshold) {
    /*Generate seed which only in low average intensity regions*/
    typename InputImageType::RegionType imageRegion =
            image->GetLargestPossibleRegion();
    if (d == InputImageType::ImageDimension){
        typename InputImageType::PixelType currentValue, minValue;
        double meanValue = 0;
        typename InputImageType::RegionType seedRegion;
        using InputImageIterator =
                itk::ImageRegionConstIteratorWithIndex<InputImageType>;
        seedRegion.SetIndex(seed);
        typename InputImageType::SizeType searchSize;
        for(size_t dd = 0; dd < InputImageType::ImageDimension; ++dd)
            searchSize[dd] = cellSize[dd]/2;
        seedRegion.SetSize(searchSize);
        seedRegion.Crop(imageRegion);
        InputImageIterator it(image, seedRegion);
        it.GoToBegin();
        minValue = itk::NumericTraits<typename InputImageType::PixelType>::max();
        double alpha = 0;
        //alpha = n-1/n
        while (!it.IsAtEnd()) {
            currentValue = it.Get();
            // alpha = n-1/n
            //mean_n = alpha*mean_n + (1-alpha)*a_n
            meanValue = meanValue * alpha + (1 - alpha) * currentValue;
            //2 - alpha_n = (n+1)/n
            // 1/(2 - alpha_n) = n/n+1
            alpha = 1/(2-alpha);
            if (minValue > currentValue) {
                minValue = currentValue;
                seed = it.GetIndex();
            }
            ++it;
        }
        if(meanValue < criterionThreshold) {
            seeds.push_back(seed);
        }
        return;
    }
    typename InputImageType::SizeType imageSize = imageRegion.GetSize();
    for (size_t x = cellSize[d]/4; x < imageSize[d]-cellSize[d]/4; x += cellSize[d]) {
        seed[d] = x;
        generateCriterionSeeds<InputImageType>(seeds, d+1, seed, cellSize, image, criterionThreshold);
    }
}

template<typename InputImageType>
void generateSeeds(std::vector<typename InputImageType::IndexType>& seeds,
                   size_t d, typename InputImageType::IndexType seed,
                   typename InputImageType::SizeType cellSize,
                   typename InputImageType::Pointer& image) {
    typename InputImageType::RegionType imageRegion =
        image->GetLargestPossibleRegion();
    if (d == InputImageType::ImageDimension){
        typename InputImageType::PixelType currentValue, minValue;
        typename InputImageType::RegionType seedRegion;
        using InputImageIterator =
                itk::ImageRegionConstIteratorWithIndex<InputImageType>;
        seedRegion.SetIndex(seed);
        typename InputImageType::SizeType searchSize;
        for(size_t dd = 0; dd < InputImageType::ImageDimension; ++dd)
            searchSize[dd] = cellSize[dd]/2;
        seedRegion.SetSize(searchSize);
        seedRegion.Crop(imageRegion);
        InputImageIterator it(image, seedRegion);
        it.GoToBegin();
        minValue = itk::NumericTraits<typename InputImageType::PixelType>::max();
        while (!it.IsAtEnd()) {
            currentValue = it.Get();
            if (minValue > currentValue) {
                minValue = currentValue;
                seed = it.GetIndex();
            }
            ++it;
        }
        seeds.push_back(seed);
        return;
    }
    typename InputImageType::SizeType imageSize = imageRegion.GetSize();
    for (size_t x = cellSize[d]/4; x < imageSize[d]-cellSize[d]/4; x += cellSize[d]) {
        seed[d] = x;
        generateSeeds<InputImageType>(seeds, d+1, seed, cellSize, image);
    }
}


template <typename LevelSetType, typename InputImageType>
typename LevelSetType::Pointer
initializeLevelSet(const itk::CommandLineArgumentParser::Pointer &parser,
                  const itk::Logger::Pointer &logger,
                  typename InputImageType::Pointer &input) {
    using LevelSetOutputType = typename LevelSetType::OutputType;
    static constexpr unsigned Dimension = LevelSetType::Dimension;
    //using LevelSetRealType = typename LevelSetType::OutputRealType;
    // Generate a seedImage mask that will be used as initialization for the level
    // set evolution.
    std::stringstream ss;
    using SeedImageType = itk::Image<LevelSetOutputType, Dimension>;
    typename SeedImageType::Pointer seedImage = SeedImageType::New();
    seedImage->SetRegions(input->GetLargestPossibleRegion());
    seedImage->CopyInformation(input);
    seedImage->Allocate();
    seedImage->FillBuffer(itk::NumericTraits<LevelSetOutputType>::Zero);

    using IndexType = typename InputImageType::IndexType;
    int seedCount = 150;
    parser->GetCommandLineArgument("-seedcount", seedCount);
    seedCount = std::max(seedCount, 0);
    std::vector<std::string> args = {"-seedcount",std::to_string(seedCount)};
    parser->AddCommandLineArguments(args);
    logger->Info("Set initial Seed Count = " + std::to_string(seedCount) + "\n");
    using SizeType = typename InputImageType::SizeType;
    std::vector<IndexType> seedPositions;
    SizeType regionSize = input->GetLargestPossibleRegion().GetSize();

    SizeType gridSize;
    for (size_t d = 0; d < Dimension; ++d) {
//TODO:Fix this
// let l(b,h) = number of patches/segments in direction l(b,h)
// L(B,H) = Length of volume in l(b,h) direction
// distribute such that l is proportional  to L/L+B+H (longer size has larger segments)
// then l = a*L/(L+B+H)
// => l*b*h = c (seed count)
// a^3 * (LBH)/(L+B+H)^3 = c
// => a = cube_root(c* sum(lengths)^Dim/Prod(lengths))

    gridSize[d] =
        std::llround(regionSize[d] / std::pow(seedCount, 1.0 / Dimension));
    }
    ss << gridSize;
    logger->Debug("Set cell Size of seeding = " + ss.str() + "\n");
    ss.clear();

    IndexType seed;
    seed.Fill(0);
    if(parser->ArgumentExists("-meanseedcriterion")){
        double meanCriterion = 0.3;
        parser->GetCommandLineArgument("-meanseedcriterion",meanCriterion);
        logger->Info("Set mean intensity threshold for seed placing to " + std::to_string(meanCriterion) + "\n");
        generateCriterionSeeds<InputImageType>(seedPositions, 0, seed, gridSize, input, meanCriterion);
    }else {
        generateSeeds<InputImageType>(seedPositions, 0, seed, gridSize, input);
    }
    logger->Info("Placed " + std::to_string(seedPositions.size()) + "seeds\n");
    for (size_t id = 0; id < seedPositions.size(); ++id) {
        //ss << "Added fastmarching starting Node  @" << seedPositions[id]
        //   << std::endl;
        //logger->Debug(ss.str());
        //ss.clear();
        seedImage->SetPixel(seedPositions[id],
                            itk::NumericTraits<LevelSetOutputType>::One);
    }

    using StructuringElementType =
        itk::BinaryBallStructuringElement<LevelSetOutputType, Dimension>;
    StructuringElementType structuringElement;
    unsigned radiusValue = 3;
    parser->GetCommandLineArgument("-seedsize", radiusValue);
    logger->Info("Set seed size = " + std::to_string(radiusValue) + "\n");

    structuringElement.SetRadius(radiusValue);
    structuringElement.CreateStructuringElement();
    using BinaryDilateImageFilterType =
        itk::BinaryDilateImageFilter<SeedImageType, SeedImageType,
                                     StructuringElementType>;
    typename BinaryDilateImageFilterType::Pointer dilateFilter =
        BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(seedImage);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetForegroundValue(itk::NumericTraits<LevelSetOutputType>::One);
    dilateFilter->Update();

    // convert a seedImage mask to a level-set function
    using SeedImageToLevelSetType =
        itk::BinaryImageToLevelSetImageAdaptor<SeedImageType, LevelSetType>;
    typename SeedImageToLevelSetType::Pointer adaptor =
        SeedImageToLevelSetType::New();
    adaptor->SetInputImage(dilateFilter->GetOutput());
    adaptor->Initialize();
    // Returns signed distance function
    typename LevelSetType::Pointer levelSet = adaptor->GetModifiableLevelSet();
    if(parser->ArgumentExists("-debug")){
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path initsavepath(inputFileName);
        initsavepath = initsavepath.parent_path() / initsavepath.stem();
        initsavepath += "_seed.tif";
        writeImage<SeedImageType>(seedImage, initsavepath.string(), logger,
                                  "Wrote init seed image to " +
                                  initsavepath.string() + "\n");
        fs::path ilsSavepath(inputFileName);
        ilsSavepath = ilsSavepath.parent_path() / ilsSavepath.stem();
        ilsSavepath += "_ls_init.tif";

        typename InputImageType::Pointer ils;
        getImage<Dimension>(levelSet, ils);
        writeImage<InputImageType>(ils, ilsSavepath.string(), logger,
                                      "Wrote Initial levelset to " + ilsSavepath.string() + "\n");

    }
    return levelSet;
}


template<typename PixelType, unsigned Dimension>
void
set_boundary(typename itk::Image<PixelType,Dimension>::Pointer& input, PixelType value){
    /*Set boundary pixels of `input` image to value*/
    using ImageType = itk::Image<PixelType ,Dimension>;

    using FaceCalculatorType =
            itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<ImageType>;

    FaceCalculatorType  faceCalculator;
    typename FaceCalculatorType::FaceListType faceList;
    itk::Size<Dimension> radius;
    radius.Fill(1);
    faceList = faceCalculator(input, input->GetLargestPossibleRegion(), radius);

    typename FaceCalculatorType::FaceListType::iterator fit;
    fit = faceList.begin();
    ++fit;//skip the interior

    using ImageIteratorType = itk::ImageRegionIterator<ImageType>;
    ImageIteratorType it;

    for (; fit != faceList.end(); ++fit){
        it = ImageIteratorType(input, *fit);
        for (it.GoToBegin(); !it.IsAtEnd(); ++it){
            it.Set(value);
        }
    }
}


template <typename InternalImageType>
typename InternalImageType::Pointer
compute_edge_speed(const itk::CommandLineArgumentParser::Pointer &parser,
                   const itk::Logger::Pointer &logger,
                   typename InternalImageType::Pointer &input) {
    using InternalImagePointer = typename InternalImageType::Pointer;
    constexpr unsigned Dimension = InternalImageType::ImageDimension;

    using DiffusionFilterType =
        itk::CurvatureFlowImageFilter<InternalImageType, InternalImageType>;

    typename DiffusionFilterType::Pointer diffusionFilter =
        DiffusionFilterType::New();
    diffusionFilter->SetInput(input);
    int smoothingIterations = 10;
    parser->GetCommandLineArgument("-smoothing", smoothingIterations);
    logger->Info("Diffusion Smoothing iterations = " +
                 std::to_string(smoothingIterations) + "\n");
    diffusionFilter->SetNumberOfIterations(smoothingIterations);
    float timeStep = 0.1;
    parser->GetCommandLineArgument("-diffusionstep", timeStep);
    logger->Info("Diffusion Smoothing timeStep = " + std::to_string(timeStep) +
                 "\n");
    diffusionFilter->SetTimeStep(timeStep);
    diffusionFilter->UseImageSpacingOn();

    using MultFilterType = itk::MultiplyImageFilter<InternalImageType>;
    typename MultFilterType::Pointer multFilter = MultFilterType::New();
    multFilter->SetInput(diffusionFilter->GetOutput());
    multFilter->SetConstant(255);

    using GradientMagnitudeFilterType =
        itk::GradientMagnitudeImageFilter<InternalImageType, InternalImageType>;
    typename GradientMagnitudeFilterType::Pointer gradientMagnitudeFilter =
        GradientMagnitudeFilterType::New();
    gradientMagnitudeFilter->SetInput(multFilter->GetOutput());
    typename InternalImageType::Pointer mag =
        gradientMagnitudeFilter->GetOutput();
    if (parser->ArgumentExists("-debug")) {
        mag->Update();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path savepath(inputFileName);
        savepath = savepath.parent_path() / savepath.stem();
        savepath += "_mag.tif";
        writeImage<InternalImageType>(mag, savepath.string(), logger,
                                      "Wrote to " + savepath.string() +
                                      "\n");
    }

    using ResampleImageFilter =
        itk::ResampleImageFilter<InternalImageType, InternalImageType>;
    typename ResampleImageFilter::Pointer downSampler =
        ResampleImageFilter::New();
    using InterpolatorType =
        itk::LinearInterpolateImageFunction<InternalImageType, double>;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

    downSampler->SetInterpolator(interpolator);

    using TransformType = itk::IdentityTransform<double, Dimension>;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    downSampler->SetTransform(transform);

    using SpacingType = typename InternalImageType::SpacingType;
    SpacingType spacing;
    SpacingType inputSpacing = input->GetSpacing();
    for (size_t d = 0; d < Dimension; ++d) {
        spacing[d] = inputSpacing[d] * 2;
    }
    downSampler->SetOutputSpacing(spacing);

    input->Update();
    using SizeType = typename InternalImageType::SizeType;
    SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
    using SizeValueType = typename InternalImageType::SizeType::SizeValueType;
    SizeType size;
    double voxelCount = 1;
    for (size_t d = 0; d < Dimension; ++d) {
        size[d] = static_cast<SizeValueType>(inputSize[d] / 2);
        voxelCount *= inputSize[d];
    }
    downSampler->SetSize(size);

    downSampler->SetInput(mag);

    if (parser->ArgumentExists("-debug")) {
        InternalImagePointer img = downSampler->GetOutput();
        img->Update();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path savepath(inputFileName);
        savepath = savepath.parent_path() / savepath.stem();
        savepath += "_ss_mag.tif";
        writeImage<InternalImageType>(img, savepath.string(), logger,
                                      "Wrote to " + savepath.string() +
                                      "\n");
    }
    using SmoothingFilterType =
        itk::SmoothingRecursiveGaussianImageFilter<InternalImageType,
                                                   InternalImageType>;
    typename SmoothingFilterType::Pointer smoothFilter =
        SmoothingFilterType::New();
    unsigned seedCount = 10;
    parser->GetCommandLineArgument("-seedcount", seedCount);
    std::vector<std::string> args = {"-seedcount", std::to_string(seedCount)};
    parser->AddCommandLineArguments(args);
    logger->Info("Updated seedcount = " + std::to_string(seedCount) + "\n");
    //std::cout <<"NumofSuperpixels: " <<seedCount <<"\n";
    //std::cout <<"voxelcount: " <<voxelCount <<"\n";
    double expectedDistance = std::pow(voxelCount / seedCount, 1.0/Dimension);
    //std::cout <<"expected Distance: " <<expectedDistance <<"\n";
    double gradientSigma = std::floor(expectedDistance/1.0);
    gradientSigma = std::max(gradientSigma, std::numeric_limits<double>::epsilon());
    parser->GetCommandLineArgument("-gradientsigma", gradientSigma);
    logger->Info("Gradient Sigma = " + std::to_string(gradientSigma) + "\n");
    smoothFilter->SetSigma(gradientSigma);



    smoothFilter->SetInput(downSampler->GetOutput());

    if (parser->ArgumentExists("-debug")) {
        InternalImagePointer img = smoothFilter->GetOutput();
        img->Update();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path savepath(inputFileName);
        savepath = savepath.parent_path() / savepath.stem();
        savepath += "_smooth_ssmag.tif";
        writeImage<InternalImageType>(img, savepath.string(), logger,
                                      "Wrote to " + savepath.string() +
                                      "\n");
    }

    typename ResampleImageFilter::Pointer upSampler = ResampleImageFilter::New();
    upSampler->SetInterpolator(interpolator);
    upSampler->SetTransform(transform);
    upSampler->SetOutputSpacing(inputSpacing);
    upSampler->SetSize(inputSize);
    upSampler->SetInput(smoothFilter->GetOutput());

    ////////////////////////
    typename MultFilterType::Pointer normalizeFilter = MultFilterType::New();
    normalizeFilter->SetInput(upSampler->GetOutput());
    double peak = 1.0/std::pow(2*vnl_math::pi,1.0/Dimension)*gradientSigma;
    normalizeFilter->SetConstant(2*peak);

    ///////////////////////
    if (parser->ArgumentExists("-debug")) {
        InternalImagePointer img = normalizeFilter->GetOutput();
        img->Update();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path savepath(inputFileName);
        savepath = savepath.parent_path() / savepath.stem();
        savepath += "_smooth_mag.tif";
        writeImage<InternalImageType>(img, savepath.string(), logger,
                                      "Wrote to " + savepath.string() +
                                      "\n");
    }

    InternalImagePointer smoothMag =
        normalizeFilter->GetOutput(); // upSampler->GetOutput();

    float magHalfHeight = 10.0f; // 0.04;
    using AdditionFilterType =
        itk::AddImageFilter<InternalImageType, InternalImageType,
                            InternalImageType>;
    typename AdditionFilterType::Pointer additionFilter =
        AdditionFilterType::New();
    additionFilter->SetInput(smoothMag);
    additionFilter->SetConstant(magHalfHeight);
    InternalImagePointer denominator = additionFilter->GetOutput();

    using DivisionFilterType =
        itk::DivideImageFilter<InternalImageType, InternalImageType,
                               InternalImageType>;
    typename DivisionFilterType::Pointer divisionFilter =
        DivisionFilterType::New();
    divisionFilter->SetInput1(mag);
    divisionFilter->SetInput2(denominator);
    InternalImagePointer normGradMag = divisionFilter->GetOutput();
    if (parser->ArgumentExists("-debug")) {
        InternalImagePointer img = normGradMag;
        img->Update();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path savepath(inputFileName);
        savepath = savepath.parent_path() / savepath.stem();
        savepath += "_normGradMag.tif";
        writeImage<InternalImageType>(img, savepath.string(), logger,
                                      "Wrote to " + savepath.string() +
                                      "\n");
    }

    using ScaleFilterType = itk::MultiplyImageFilter<InternalImageType>;

    typename ScaleFilterType::Pointer scaleFilter = ScaleFilterType::New();
    scaleFilter->SetInput(normGradMag);
    double factor = 12.7; // 1/K value in exp(-grad(img)/K)
    parser->GetCommandLineArgument("-edgescale", factor);
    logger->Info("Exp edge scale = " + std::to_string(factor) + "\n");
    scaleFilter->SetConstant(-1 * factor);

    using ExponentialImageFilter =
        itk::ExpImageFilter<InternalImageType, InternalImageType>;
    typename ExponentialImageFilter::Pointer expFilter =
        ExponentialImageFilter::New();
    expFilter->SetInput(scaleFilter->GetOutput());

    typename InternalImageType::Pointer speed = expFilter->GetOutput();
    speed->Update();

    if (parser->ArgumentExists("-debug")) {
        std::string inputFileName;
        parser->GetCommandLineArgument("-input",inputFileName);
        fs::path speedSavepath(inputFileName);
        speedSavepath = speedSavepath.parent_path() / speedSavepath.stem();
        speedSavepath += "_speed.tif";
        writeImage<InternalImageType>(speed, speedSavepath.string(), logger,
                                      "Wrote speed to " + speedSavepath.string() +
                                      "\n");
    }
    return speed;
}

template<unsigned Dimension>
void
getImage(typename itk::WhitakerSparseLevelSetImage<float, Dimension>::Pointer& levelSet,
         typename itk::Image<float,Dimension>::Pointer& levelSetImage){
    using LevelSetType = itk::WhitakerSparseLevelSetImage<float, Dimension>;
    using LevelSetImageType = itk::Image<float,Dimension>;
    using LabelMapToMaskImageFilterType = itk::LabelMapToLabelImageFilter<typename LevelSetType::LabelMapType, LevelSetImageType>;
    typename LabelMapToMaskImageFilterType::Pointer labelToImageFilter;
    labelToImageFilter = LabelMapToMaskImageFilterType::New();
    labelToImageFilter->SetInput(levelSet->GetLabelMap());
    labelToImageFilter->Update();
    levelSetImage = labelToImageFilter->GetOutput();
}

template<unsigned Dimension>
void
getImage(typename itk::LevelSetDenseImage<itk::Image<float, Dimension>>::Pointer& levelSet,
         typename itk::Image<float,Dimension>::Pointer& levelSetImage){
    levelSetImage = levelSet->GetImage();
}

template<typename RGBImageType, typename MaskImageType>
void paintImage(typename RGBImageType::Pointer& input,
                typename MaskImageType::Pointer& mask,
                Color c){
    auto ipIt = itk::ImageRegionIterator<RGBImageType>(input, input->GetLargestPossibleRegion());
    auto maskIt = itk::ImageRegionConstIterator<MaskImageType>(mask, mask->GetLargestPossibleRegion());
    ipIt.GoToBegin();
    maskIt.GoToBegin();
    while(!ipIt.IsAtEnd()){
        if(maskIt.Get() > 0){
            RGBPixelType p;
            p[0] = c[0]*255; p[1] = c[1]*255; p[2] = c[2]*255;
            ipIt.Set(p);
        }
        ++ipIt;
        ++maskIt;
    }
}

template<typename InputImageType,typename RGBImageType>
typename RGBImageType::Pointer
paintCountours(typename InputImageType::Pointer& input,
               std::vector<typename InputImageType::PixelType>& values){
    std::vector<Color> cmap;
    double tol = 1;
    Color from(3,1);
    //from[2] = 1;
    Color to(3,0);
    to[1] = 0.7;
    for(size_t c = 0 ; c  < values.size(); ++c){
        Color current(3,0);
        double alpha = c/(values.size()-1.0);
        std::cout << "alpha: "<< alpha << std::endl;
        for(size_t i =0; i < 3; ++i) current[i] = (1-alpha)*from[i] + alpha*to[i];
        cmap.push_back(current);
    }
    std::cout << "Colormap:\n";
    for(auto& c: cmap){
        std::cout <<"("<< c[0] << "," << c[1] << ","<< c[2] << ")\n";
    }
    std::cout << "Input image size:\n" << input->GetLargestPossibleRegion() << std::endl;
    //using RGBPixelType = itk::RGBPixel<unsigned char>;
    //using RGBImageType = itk::Image<RGBPixelType, InputImageType::Dimension>;
    typename RGBImageType::Pointer coloredImage = RGBImageType::New();
    coloredImage->SetRegions(input->GetLargestPossibleRegion());
    coloredImage->Allocate(true);
    std::cout << "Allocated colored image:\n" << coloredImage->GetLargestPossibleRegion() << std::endl;

    auto ipIt = itk::ImageRegionIterator<InputImageType>(input, input->GetLargestPossibleRegion());
    auto colIt = itk::ImageRegionIterator<RGBImageType>(coloredImage, coloredImage->GetLargestPossibleRegion());
    ipIt.GoToBegin();
    colIt.GoToBegin();
    while(!ipIt.IsAtEnd()){
        auto ip = ipIt.Get();
        for(size_t c = 0 ; c  < values.size() ; ++c){
            if(std::abs(static_cast<double>(ip - values[c])) <= tol) {
                typename RGBImageType::PixelType p;
                p[0] = cmap[c][0] * 255;
                p[1] = cmap[c][1] * 255;
                p[2] = cmap[c][2] * 255;
                colIt.Set(p);
                break;
            }
        }
        ++ipIt;
        ++colIt;
    }
    return coloredImage;
}

#endif // TURBO_UTILS_HXX
