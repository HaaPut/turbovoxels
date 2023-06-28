
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
// Created by tabish on 2021-05-14.
//

#ifndef TURBO_MONITOR_H
#define TURBO_MONITOR_H

#include "itkGrayToRGBColormapFunction.h"
#include "itkLevelSetSparseMaskedEquationTermContainer.h"
#include "itkZeroCrossingColormapFunction.h"
#include "utils.h"
#include <cstdio>
#include <ctime>
#include "filesystem.hpp"
#include <iomanip>
#include <iostream>
#include <itkCommand.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLevelSetContainer.h>
#include <itkLevelSetDenseImage.h>
#include <itkLevelSetEquationContainer.h>
#include <itkLevelSetEvolution.h>
#include <itkLevelSetMaskedEquationTermContainer.h>
#include <itkRGBPixel.h>
#include <itkScalarToRGBColormapImageFilter.h>
#include <itkZeroCrossingImageFilter.h>
#include <sstream>
#include <string>

namespace fs = std::filesystem;


using LevelSetImageType2 = itk::Image<float, 2>;
using LevelSetType2 = itk::LevelSetDenseImage<LevelSetImageType2>;
using LevelSetContainerType2 = itk::LevelSetContainer<itk::IdentifierType, LevelSetType2>;
using TermContainerType2 = itk::LevelSetMaskedEquationTermContainer<itk::Image<float,2>>;
using EquationContainerType2 = itk::LevelSetEquationContainer<TermContainerType2>;
using LevelSetEvolutionType2 = itk::LevelSetEvolution<EquationContainerType2, LevelSetType2>;

using LevelSetImageType3 = itk::Image<float, 3>;
using LevelSetType3 = itk::LevelSetDenseImage<LevelSetImageType3>;
using LevelSetContainerType3 = itk::LevelSetContainer<itk::IdentifierType, LevelSetType3>;
using TermContainerType3 = itk::LevelSetMaskedEquationTermContainer<itk::Image<float,3>>;
using EquationContainerType3 = itk::LevelSetEquationContainer<TermContainerType3>;
using LevelSetEvolutionType3 = itk::LevelSetEvolution<EquationContainerType3, LevelSetType3>;

using SparseLevelSetType2 = itk::WhitakerSparseLevelSetImage<float, 2>;
using SparseLevelSetContainerType2 = itk::LevelSetContainer<itk::IdentifierType, SparseLevelSetType2>;
using SparseTermContainerType2 = itk::LevelSetSparseMaskedEquationTermContainer<itk::Image<float,2>,SparseLevelSetType2>;
using SparseEquationContainerType2 = itk::LevelSetEquationContainer<SparseTermContainerType2>;
using SparseLevelSetEvolutionType2 = itk::LevelSetEvolution<SparseEquationContainerType2, SparseLevelSetType2>;

using SparseLevelSetType3 = itk::WhitakerSparseLevelSetImage<float, 3>;
using SparseLevelSetContainerType3 = itk::LevelSetContainer<itk::IdentifierType, SparseLevelSetType3>;
using SparseTermContainerType3 = itk::LevelSetSparseMaskedEquationTermContainer<itk::Image<float,3>,SparseLevelSetType3>;
using SparseEquationContainerType3 = itk::LevelSetEquationContainer<SparseTermContainerType3>;
using SparseLevelSetEvolutionType3 = itk::LevelSetEvolution<SparseEquationContainerType3, SparseLevelSetType3>;

template <typename LevelSetEvolutionType>
class SegmentationIterationUpdate : public itk::Command
{
public:
    using Self = SegmentationIterationUpdate;
    using Superclass = itk::Command;
    using Pointer = itk::SmartPointer<Self>;
    itkNewMacro(Self);

    using LevelSetType = typename LevelSetEvolutionType::LevelSetType;
    using TermContainerType = typename LevelSetEvolutionType::TermContainerType;
    using LevelSetImageType = typename TermContainerType::MaskImageType;
    using SkeletonImageType = typename TermContainerType::MaskImageType;

    using PixelType = typename LevelSetImageType::PixelType;
    //using RGBPixelType = itk::RGBPixel<unsigned char>;
    static constexpr unsigned Dimension = LevelSetImageType::ImageDimension;
    using RGBImageType = itk::Image<RGBPixelType, Dimension>;

    using ZeroCrossingColormapType =
        itk::Function::ZeroCrossingColormapFunction<PixelType, RGBPixelType>;

    using GrayToRGBColormapType = itk::Function::GrayToRGBColormapFunction<PixelType, RGBPixelType>;

    using ColormapFilterType =
        itk::ScalarToRGBColormapImageFilter<LevelSetImageType, RGBImageType>;
    //using Color = std::vector<typename ZeroCrossingColormapType::RealType>;

    using ZeroCrossingFilterType =
        itk::ZeroCrossingImageFilter<LevelSetImageType, LevelSetImageType>;

    using RGBLSWriterType = itk::ImageFileWriter<RGBImageType>;
    using LSWriterType = itk::ImageFileWriter<LevelSetImageType>;

  protected:
    SegmentationIterationUpdate() {
        //colormap = ZeroCrossingColormapType::New();
        colormap = GrayToRGBColormapType::New();
        colormap->SetMinimumInputValue(m_Minimum);
        colormap->SetMaximumInputValue(m_Maximum);
        contourFilter = ZeroCrossingFilterType::New();

        colormapFilter = ColormapFilterType::New();
        colormapFilter->SetColormap(colormap);
        m_lswriter = LSWriterType::New();
        m_rgblswriter = RGBLSWriterType::New();
        tStart = clock();

        //labelToImageFilter = LabelMapToMaskImageFilterType::New();
        m_everyNth = 10;
    }

public:
    void
    Execute(itk::Object *caller, const itk::EventObject &event) override {
      Execute((const itk::Object *)caller, event);
    }

    void
    Execute(const itk::Object * object, const itk::EventObject & event) override {
        using ConstEvolutionPointer = const LevelSetEvolutionType *;
        using StoppingCriterionType = typename LevelSetEvolutionType::StoppingCriterionType;
        auto evolution = static_cast<ConstEvolutionPointer>(object);
        if (!(itk::IterationEvent().CheckEvent(&event))){
            return;
        }

        StoppingCriterionType *stoppingCriterion =const_cast<StoppingCriterionType*>(evolution->GetStoppingCriterion());
        unsigned iteration = stoppingCriterion->GetCurrentIteration();
        double rmse = stoppingCriterion->GetRMSChangeAccumulator();
        std::cout << "Iteration:: " << iteration
                  << " RMSE :: " << rmse
                  << std::endl;
        if(iteration % m_everyNth == 0){
            typename LevelSetType::Pointer levelSet =
                evolution->GetLevelSetContainer()
                ->GetLevelSet(m_LevelSetId);
                typename LevelSetImageType::Pointer levelSetImage;
                getImage<Dimension>(levelSet, levelSetImage);

                typename LevelSetImageType::Pointer input = evolution->GetEquationContainer()->GetEquation(0)->GetInput();
                //colormapFilter->SetInput(levelSetImage);
                colormapFilter->SetInput(input);
                typename RGBImageType::Pointer coloredImage =
                    colormapFilter->GetOutput();
                coloredImage->Update();
                coloredImage->DisconnectPipeline();
                if(m_colorSkel) {
                    auto termContainer =
                            evolution->GetEquationContainer()->GetEquation(0);
                    typename SkeletonImageType::Pointer skeleton = termContainer->GetMask();
                    Color green(3, 0);
                    green[1] = 1;

                    paintImage<RGBImageType, SkeletonImageType>(coloredImage,
                                                                skeleton, green);
                }
                contourFilter->SetInput(levelSetImage);
                contourFilter->SetBackgroundValue(
                    itk::NumericTraits<float>::Zero);
                contourFilter->SetForegroundValue(
                    itk::NumericTraits<float>::One);
                typename LevelSetImageType::Pointer contour =
                    contourFilter->GetOutput();
                contour->Update();
                contour->DisconnectPipeline();
                Color red(3,0);
                red[0] = 1;
                paintImage<RGBImageType, LevelSetImageType>(coloredImage,
                                                            contour, red);
                fs::path lssavepath = m_OutputPath / m_FileNameStem;
                std::stringstream suffix;
                if (m_rgb) suffix << "-lsrgb";
                else suffix << "-ls";
                suffix << std::setfill('0')<< std::setw(5) << evolution->GetNumberOfIterations() << ".tif";
                lssavepath += suffix.str();
                m_lswriter->SetFileName(lssavepath.string());
                m_rgblswriter->SetFileName(lssavepath.string());
                m_rgblswriter->SetInput(coloredImage);
                m_lswriter->SetInput(levelSetImage);
                try {
                    if(m_rgb) m_rgblswriter->Update();
                    else m_lswriter->Update();
                    std::cout << "Wrote levelset File to " + std::string(m_lswriter->GetFileName()) + "\n";
                } catch (const itk::ExceptionObject &excep) {
                    std::cerr << "Error Writing Final File!\n";
                    std::cerr << excep.what() << "\n";
                    return;
                }
                printf("Time taken: %.2fs\n",
                       (double)(clock() - tStart) / (CLOCKS_PER_SEC*m_everyNth));
                tStart = clock();
        }
    }
    void
    SetOutputPath(fs::path& path){
        m_OutputPath = path;
        if(m_OutputPath.has_filename()){
            m_FileNameStem = m_OutputPath.stem();
        }else{
            m_FileNameStem = "frame";
        }
        m_OutputPath.remove_filename();
        m_OutputPath = m_OutputPath.parent_path() / "frames";
        //does nothing if directory exists
        fs::create_directories(m_OutputPath);
    }
    void SetLevelSetId(unsigned id){
        m_LevelSetId = id;
    }
    void SetMinumum(float min){
        m_Minimum = min;
    }
    void SetMaximum(float max){
        m_Maximum = max;
    }
    void SetRGB(bool state){
      m_rgb = state;
    }

    void ColorSkeleton(bool state){
        m_colorSkel = state;
    }

    void SetFrequency(int every){
        m_everyNth = every;
    }

    bool m_rgb = true;
    bool m_colorSkel = true;
    typename LSWriterType::Pointer m_lswriter;
    typename RGBLSWriterType::Pointer m_rgblswriter;

    //typename ZeroCrossingColormapType::Pointer colormap;
    typename GrayToRGBColormapType::Pointer colormap;

    typename ColormapFilterType::Pointer colormapFilter;
    //typename LabelMapToMaskImageFilterType::Pointer labelToImageFilter;
    typename ZeroCrossingFilterType::Pointer contourFilter;
    
    clock_t tStart;
    float m_Minimum = -2;
    float m_Maximum = 2;
    fs::path m_FileNameStem;
    fs::path m_OutputPath;
    unsigned m_LevelSetId = 0;
    unsigned m_everyNth;
};

template <typename LevelSetEvolutionType>
class IterationLog : public itk::Command
{
public:
    using Self = IterationLog;
    using Superclass = itk::Command;
    using Pointer = itk::SmartPointer<Self>;
    itkNewMacro(Self);

//    using LevelSetType = typename LevelSetEvolutionType::LevelSetType;
//    using TermContainerType = typename LevelSetEvolutionType::TermContainerType;
//    using LevelSetImageType = typename TermContainerType::MaskImageType;
//    using SkeletonImageType = typename TermContainerType::MaskImageType;
//
//    using PixelType = typename LevelSetImageType::PixelType;
//    //using RGBPixelType = itk::RGBPixel<unsigned char>;
//    static constexpr unsigned Dimension = LevelSetImageType::ImageDimension;
//    using RGBImageType = itk::Image<RGBPixelType, Dimension>;
//
//    using ZeroCrossingColormapType =
//            itk::Function::ZeroCrossingColormapFunction<PixelType, RGBPixelType>;
//
//    using GrayToRGBColormapType = itk::Function::GrayToRGBColormapFunction<PixelType, RGBPixelType>;
//
//    using ColormapFilterType =
//            itk::ScalarToRGBColormapImageFilter<LevelSetImageType, RGBImageType>;
//    //using Color = std::vector<typename ZeroCrossingColormapType::RealType>;
//
//    using ZeroCrossingFilterType =
//            itk::ZeroCrossingImageFilter<LevelSetImageType, LevelSetImageType>;
//
//    using RGBLSWriterType = itk::ImageFileWriter<RGBImageType>;
//    using LSWriterType = itk::ImageFileWriter<LevelSetImageType>;

protected:
    IterationLog() { }

public:
    void
    Execute(itk::Object *caller, const itk::EventObject &event) override {
        Execute((const itk::Object *)caller, event);
    }

    void
    Execute(const itk::Object * object, const itk::EventObject & event) override {
        using ConstEvolutionPointer = const LevelSetEvolutionType *;
        using StoppingCriterionType = typename LevelSetEvolutionType::StoppingCriterionType;
        auto evolution = static_cast<ConstEvolutionPointer>(object);
        if (!(itk::IterationEvent().CheckEvent(&event))){
            return;
        }

        auto *stoppingCriterion =const_cast<StoppingCriterionType*>(evolution->GetStoppingCriterion());
        unsigned iteration = stoppingCriterion->GetCurrentIteration();
        double rmse = stoppingCriterion->GetRMSChangeAccumulator();
        std::cout << "Iteration:: " << iteration
                  << " RMSE :: " << rmse
                  << std::endl;
    }
//    void
//    SetOutputPath(fs::path& path){
//        m_OutputPath = path;
//        if(m_OutputPath.has_filename()){
//            m_FileNameStem = m_OutputPath.stem();
//        }else{
//            m_FileNameStem = "frame";
//        }
//        m_OutputPath.remove_filename();
//        m_OutputPath = m_OutputPath.parent_path() / "frames";
//        //does nothing if directory exists
//        fs::create_directories(m_OutputPath);
//    }
//    void SetLevelSetId(unsigned id){
//        m_LevelSetId = id;
//    }
//    void SetMinumum(float min){
//        m_Minimum = min;
//    }
//    void SetMaximum(float max){
//        m_Maximum = max;
//    }
//    void SetRGB(bool state){
//        m_rgb = state;
//    }
//
//    void ColorSkeleton(bool state){
//        m_colorSkel = state;
//    }
//

//    bool m_rgb = true;
//    bool m_colorSkel = true;
//    typename LSWriterType::Pointer m_lswriter;
//    typename RGBLSWriterType::Pointer m_rgblswriter;
//
//    //typename ZeroCrossingColormapType::Pointer colormap;
//    typename GrayToRGBColormapType::Pointer colormap;
//
//    typename ColormapFilterType::Pointer colormapFilter;
//    //typename LabelMapToMaskImageFilterType::Pointer labelToImageFilter;
//    typename ZeroCrossingFilterType::Pointer contourFilter;
//
//    clock_t tStart;
//    float m_Minimum = -2;
//    float m_Maximum = 2;
//    fs::path m_FileNameStem;
//    fs::path m_OutputPath;
//    unsigned m_LevelSetId = 0;
//    unsigned m_everyNth;
};
#endif
