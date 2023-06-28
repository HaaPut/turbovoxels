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
// Created by tabish on 2021-06-01.
//


#ifndef TURBO_LEVELSETEVOLUTION_HXX
#define TURBO_LEVELSETEVOLUTION_HXX

#include "itkCommandLineArgumentParser.h"
#include "utils.h"
#include "monitor.h"

#include <itkLogger.h>
#include <itkImage.h>
#include <itkGradientImageFilter.h>
#include <itkSinRegularizedHeavisideStepFunction.h>
#include <itkLevelSetContainer.h>
#include <itkLevelSetEquationPropagationTerm.h>
#include <itkLevelSetEquationCurvatureTerm.h>
#include <itkLevelSetEquationAdvectionTerm.h>
#include <itkLevelSetSparseMaskedEquationTermContainer.h>
#include <itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion.h>
#include <itkLevelSetEvolution.h>

#include <cstdlib>
#include <itkLevelSetEquationContainer.h>

template <typename TRealImageType>
typename TRealImageType::Pointer
maskedSparseLevelSetEvolution( typename TRealImageType::Pointer input,
                               const itk::CommandLineArgumentParser::Pointer &parser,
                               const itk::Logger::Pointer &logger) {
    using InternalImageType = TRealImageType;
    using InternalPixelType = typename InternalImageType::PixelType;
    constexpr unsigned Dimension = InternalImageType::ImageDimension;
    static_assert(std::is_floating_point<InternalPixelType>::value,
                  "Internal pixel type for level sets should be floating point");

    using InternalImagePointer = typename InternalImageType::Pointer;
    using LevelSetType = itk::WhitakerSparseLevelSetImage<InternalPixelType, Dimension>;
    using LevelSetImageType = itk::Image<InternalPixelType, Dimension>;
    using LevelSetRealType = typename LevelSetType::OutputRealType;


    InternalImagePointer speed =
            compute_edge_speed<InternalImageType>(parser, logger, input);

    using GradientFilterType =
    itk::GradientImageFilter<InternalImageType, LevelSetRealType,
            LevelSetRealType>;
    using GradientImageType = typename GradientFilterType::OutputImageType;
    typename GradientFilterType::Pointer speedGradientFilter =
            GradientFilterType::New();
    speedGradientFilter->SetInput(speed);
    typename GradientImageType::Pointer speed_grad =
            speedGradientFilter->GetOutput();
    speed_grad->Update();

    typename LevelSetType::Pointer levelSet;
    if(parser->ArgumentExists("-gradientseeding")) {
        logger->Info("Using Gradient Magnitude image for seeding\n");
        using GradientMagnitudeFilterType =
                itk::GradientMagnitudeImageFilter<TRealImageType, InternalImageType>;
        typename GradientMagnitudeFilterType::Pointer gradMagnitudeFilter =
                GradientMagnitudeFilterType::New();
        gradMagnitudeFilter->SetInput(input);
        typename InternalImageType::Pointer gradMagnitude =
                gradMagnitudeFilter->GetOutput();
        gradMagnitude->Update();

        levelSet = initializeLevelSet<LevelSetType, InternalImageType>(parser, logger,
                                                                    gradMagnitude);
    }else{
        logger->Info("Using smoothed raw image for seeding\n");
        using SmoothingFilterType =
                itk::SmoothingRecursiveGaussianImageFilter<InternalImageType,
                        InternalImageType>;
        typename SmoothingFilterType::Pointer smoothFilter =
                SmoothingFilterType::New();
        double sigma = 2;
        parser->GetCommandLineArgument("-seedsigma", sigma);
        logger->Info("Set sigma for seeding image = " + std::to_string(sigma) + "\n");
        smoothFilter->SetSigma(sigma);
        smoothFilter->SetInput(input);
        typename InternalImageType::Pointer smoothInput = smoothFilter->GetOutput();
        smoothInput->Update();
        levelSet = initializeLevelSet<LevelSetType, InternalImageType>(parser, logger, smoothInput);
    }

    // The Heaviside function
    using HeavisideFunctionType =
    itk::SinRegularizedHeavisideStepFunction<LevelSetRealType , LevelSetRealType>;//InternalPixelType,InternalPixelType>;
    typename HeavisideFunctionType::Pointer heaviside =
            HeavisideFunctionType::New();
    heaviside->SetEpsilon(1.5);

    // Create the level set container
    using LevelSetContainerType =
    itk::LevelSetContainer<itk::IdentifierType, LevelSetType>;
    typename LevelSetContainerType::Pointer levelSetContainer =
            LevelSetContainerType::New();
    levelSetContainer->SetHeaviside(heaviside);
    typename LevelSetContainerType::LevelSetIdentifierType levelSetId = 0;
    levelSetContainer->AddLevelSet(levelSetId, levelSet);

    // Create the terms.

    using PropagationTermType = itk::LevelSetEquationPropagationTerm<InternalImageType, LevelSetContainerType, InternalImageType>;
    typename PropagationTermType::Pointer propagationTerm =
            PropagationTermType::New();
    double propagationScaling = 1.0;
    parser->GetCommandLineArgument("-propagation", propagationScaling);
    propagationTerm->SetPropagationImage(speed);
    logger->Info("Propagation scale: " + std::to_string(propagationScaling) +
                 "\n");
    propagationScaling *= -1; // weight is -ve for outward flow
    propagationTerm->SetCoefficient(propagationScaling);

    using CurvatureTermType =
    itk::LevelSetEquationCurvatureTerm<InternalImageType,
            LevelSetContainerType>;
    typename CurvatureTermType::Pointer curvatureTerm =
            CurvatureTermType::New();
    curvatureTerm->SetUseCurvatureImage(true);
    curvatureTerm->SetCurvatureImage(speed);
    double curvatureScaling = 0.3;
    parser->GetCommandLineArgument("-curvature", curvatureScaling);
    logger->Info("Curvature scale: " + std::to_string(curvatureScaling) +
                 "\n");
    curvatureScaling *= -1; //weight is -ve for outward flow
    curvatureTerm->SetCoefficient(curvatureScaling);

    using AdvectionTermType =
    itk::LevelSetEquationAdvectionTerm<InternalImageType,
            LevelSetContainerType>;
    typename AdvectionTermType::Pointer advectionTerm =
            AdvectionTermType::New();
    double advectionScaling = 0.5;
    parser->GetCommandLineArgument("-advection", advectionScaling);
    logger->Info("Advection scale: " + std::to_string(advectionScaling) +
                 "\n");
    advectionScaling *= -1; // weight is -ve for outward flow

    advectionTerm->SetCoefficient(advectionScaling);
    // advectionTerm = grad(phi).advectionImage
    // default advectionImage = gradient of input image
    advectionTerm->SetAdvectionImage(speed_grad);

    // Create term container (equation rhs)
    // using TermContainerType =
    // itk::LevelSetMaskedEquationTermContainer<InternalImageType>;
    using TermContainerType =
    itk::LevelSetSparseMaskedEquationTermContainer<InternalImageType,
            LevelSetType>;
    typename TermContainerType::Pointer termContainer =
            TermContainerType::New();
    termContainer->SetLevelSetContainer(levelSetContainer);
    termContainer->SetInput(input);
    termContainer->AddTerm(0, curvatureTerm);
    termContainer->AddTerm(1, propagationTerm);
    termContainer->AddTerm(2, advectionTerm);

    // Create equation container
    using EquationContainerType =
    itk::LevelSetEquationContainer<TermContainerType>;
    typename EquationContainerType::Pointer equationContainer =
            EquationContainerType::New();
    equationContainer->SetLevelSetContainer(levelSetContainer);
    equationContainer->AddEquation(0, termContainer);

    // Create stopping criteria
    using StoppingCriterionType =
    itk::LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion<
            LevelSetContainerType>;
    typename StoppingCriterionType::Pointer criterion =
            StoppingCriterionType::New();

    int numberOfIterations = 500;
    parser->GetCommandLineArgument("-iterations", numberOfIterations);
    logger->Info("Set maximum number of iterations = " +
                 std::to_string(numberOfIterations) + "\n");
    criterion->SetNumberOfIterations(numberOfIterations);

    double rmsThreshold = 0.1;
    parser->GetCommandLineArgument("-rmsthreshold", rmsThreshold);
    criterion->SetRMSThreshold(rmsThreshold);
    criterion->SetRMSChangeAccumulator(criterion->GetRMSThreshold() * 2 + 1);//std::numeric_limits<double>::max());
    logger->Info("Set RMSE threshold = " +
                 std::to_string(criterion->GetRMSThreshold()) + "\n");

    // Create evolution class
    using LevelSetEvolutionType =
    itk::LevelSetEvolution<EquationContainerType, LevelSetType>;
    typename LevelSetEvolutionType::Pointer evolution =
            LevelSetEvolutionType::New();
    evolution->SetEquationContainer(equationContainer);
    evolution->SetStoppingCriterion(criterion);
    evolution->SetLevelSetContainer(levelSetContainer);
    double timestep;
    if (parser->GetCommandLineArgument("-dt", timestep)){
        evolution->SetTimeStep(timestep);
        logger->Info("Set timestep manually = " + std::to_string(timestep) +  "\n");
        logger->Warning("Make sure the time step value is in accordance with the CFL conditions\n");
    }

    if(parser->ArgumentExists("-observe")) {
        typename SegmentationIterationUpdate<LevelSetEvolutionType>::Pointer
                observer = SegmentationIterationUpdate<LevelSetEvolutionType>::New();
        std::string inputFileName;
        parser->GetCommandLineArgument("-input", inputFileName);
        fs::path inputfilepath(inputFileName);
        observer->SetOutputPath(inputfilepath);
        bool rgb = parser->ArgumentExists("-rgb");
        observer->SetRGB(rgb);
        evolution->AddObserver(itk::IterationEvent(), observer);
    }else if(parser->ArgumentExists("-slient")){
        // no observer...
    }else{
        typename IterationLog<LevelSetEvolutionType>::Pointer
                observer = IterationLog<LevelSetEvolutionType>::New();
        evolution->AddObserver(itk::IterationEvent(), observer);
    }
    evolution->Update();
    logger->Debug("Evolution Finshed:\n");
    logger->Debug("Completed : " +
                  std::to_string(evolution->GetNumberOfIterations()) + " Iterations\n");

    typename LevelSetType::Pointer finalLevelSet =
            evolution->GetLevelSetContainer()->GetLevelSet(levelSetId);
    typename LevelSetImageType::Pointer levelSetImage;
    getImage<Dimension>(finalLevelSet, levelSetImage);
    return levelSetImage;
}



#endif //TURBO_LEVELSETEVOLUTION_HXX
