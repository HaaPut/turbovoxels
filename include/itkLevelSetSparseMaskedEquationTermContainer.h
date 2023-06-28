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

#ifndef itkLevelSetSparseMaskedEquationTermContainer_h
#define itkLevelSetSparseMaskedEquationTermContainer_h


#include <itkLevelSetEquationTermBase.h>
#include <itkObject.h>
#include <itkMacro.h>
#include "itkSkeletonImageFilter.h"
#include <itkLevelSetContainer.h>
#include <itkLabelMapToLabelImageFilter.h>
#include "itkLevelSetMaskedEquationTermContainerBase.h"

#include <unordered_map>


namespace itk
{
/**
 *  \class LevelSetSparseMaskedEquationTermContainer
 *  \brief Class for container holding the terms of a given level set update equation
 *
 *  \tparam TInputImage Input image or speed image or feature image for segmentation
 *  \tparam TLevelSetContainer Container holding the all the level set functions
 *
 *  \ingroup ITKLevelSetsv4
 */

    template <typename TInputImage, typename TLevelSetType>
    class ITK_TEMPLATE_EXPORT LevelSetSparseMaskedEquationTermContainer : public LevelSetMaskedEquationTermContainerBase<TInputImage, LevelSetContainer<IdentifierType,TLevelSetType>>
    {
    public:
        //ITK_DISALLOW_COPY_AND_MOVE(LevelSetSparseMaskedEquationTermContainer);
        using LevelSetContainerType = LevelSetContainer<IdentifierType,TLevelSetType>;

        using Self = LevelSetSparseMaskedEquationTermContainer;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;
        using Superclass = LevelSetMaskedEquationTermContainerBase< TInputImage, LevelSetContainerType>;

        /** Method for creation through object factory */
        itkNewMacro(Self);

        /** Run-time type information */
        itkTypeMacro(LevelSetSparseMaskedEquationTermContainer, LevelSetMaskedEquationTerContainerBase);

        using TermIdType = unsigned int;

        using InputImageType = TInputImage;
        using InputImagePointer = typename InputImageType::Pointer;
        using InputPixelType = typename InputImageType::PixelType;

        using LevelSetContainerPointer = typename Superclass::LevelSetContainerPointer;

        using LevelSetType = typename Superclass::LevelSetType;
        using LevelSetPointer = typename Superclass::LevelSetPointer;

        using LevelSetIdentifierType =
            typename Superclass::LevelSetIdentifierType;
        using LevelSetOutputPixelType =
            typename Superclass::LevelSetOutputPixelType;
        using LevelSetOutputRealType =
            typename Superclass::LevelSetOutputRealType;
        using LevelSetDataType = typename Superclass::LevelSetDataType;
        using LevelSetInputIndexType =
            typename Superclass::LevelSetInputIndexType;
        using LevelSetGradientType = typename Superclass::LevelSetGradientType;
        using LevelSetHessianType = typename Superclass::LevelSetHessianType;
        using LevelSetLabelMapType = typename LevelSetType::LabelMapType;
        using MaskImageType = typename Superclass::MaskImageType;
        using LabelMapToMaskImageFilterType = LabelMapToLabelImageFilter<LevelSetLabelMapType,MaskImageType>;
        using MaskImagePointer = typename Superclass::MaskImagePointer;
        using MaskPixelType = typename Superclass::MaskPixelType;

        using SkeletonFilterType =
            SkeletonImageFilter<MaskPixelType, InputImageType::ImageDimension>;
        using SkeletonOutputType = MaskPixelType;

        using TermType =
            LevelSetEquationTermBase<InputImageType, LevelSetContainerType>;
        using TermPointer = typename Superclass::TermPointer;

        /** Update background skeleton */
        void UpdateMask() override {
            auto tIt = this->m_Container.begin();
            LevelSetPointer levelSet =(tIt->second)->GetModifiableCurrentLevelSetPointer();
            this->m_LabelToImageFilter->SetInput(levelSet->GetModifiableLabelMap());
            this->m_SkeletonFilter->SetInput(m_LabelToImageFilter->GetOutput());
            this->m_SkeletonFilter->Update();
            this->m_Mask = this->m_SkeletonFilter->GetOutput();
            this->m_Mask->DisconnectPipeline();
        }


    protected:

        LevelSetSparseMaskedEquationTermContainer(){
            this->m_CurrentLevelSetId = LevelSetIdentifierType();
            this->m_Mask = nullptr;
            this->m_SkeletonFilter = SkeletonFilterType::New();
            this->m_SkeletonFilter->ReleaseDataFlagOn();
            this->m_LabelToImageFilter = LabelMapToMaskImageFilterType::New();
            this->m_LabelToImageFilter->ReleaseDataFlagOn();// this should perhaps go
        }
        typename LabelMapToMaskImageFilterType::Pointer m_LabelToImageFilter;
        unsigned m_idx = 0;
        ~LevelSetSparseMaskedEquationTermContainer() override = default;
    };

} // namespace itk

#endif // itkLevelSetSparseMaskedEquationTermContainer_h
