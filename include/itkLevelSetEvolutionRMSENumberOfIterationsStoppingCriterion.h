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
// Created by tabish on 2021-05-17.
//

#ifndef itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion_h
#define itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion_h

#include <itkObjectFactory.h>
#include <itkLevelSetEvolutionStoppingCriterion.h>

namespace itk
{
/**
 *\class LevelSetEvolutionStoppingCriterion
 * stop if rmse change is below a threshold upperbounded by Number
 * of iterations
 \ingroup ITKLevelSetsv4
*/
    template <typename TLevelSetContainer>
    class ITK_TEMPLATE_EXPORT LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion
        : public LevelSetEvolutionStoppingCriterion<TLevelSetContainer>
    {
    public:
        //ITK_DISALLOW_COPY_AND_MOVE(LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion);

        using Self = LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion;
        using Superclass = LevelSetEvolutionStoppingCriterion<TLevelSetContainer>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        /** Method for creation through object factory */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion, LevelSetEvolutionStoppingCriterion);

        using LevelSetContainerType = TLevelSetContainer;
        using LevelSetContainerPointer = typename LevelSetContainerType::Pointer;

        using LevelSetIdentifierType = typename LevelSetContainerType::LevelSetIdentifierType;
        using LevelSetType = typename LevelSetContainerType::LevelSetType;
        using LevelSetPointer = typename LevelSetContainerType::LevelSetPointer;

        using InputIndexType = typename LevelSetContainerType::InputIndexType;
        using OutputType = typename LevelSetContainerType::OutputType;
        using OutputRealType = typename LevelSetContainerType::OutputRealType;
        using GradientType = typename LevelSetContainerType::GradientType;
        using HessianType = typename LevelSetContainerType::HessianType;

        using HeavisideType = typename LevelSetContainerType::HeavisideType;
        using HeavisidePointer = typename LevelSetContainerType::HeavisideType;

        bool
            IsSatisfied() const override;

        std::string
            GetDescription() const override;

        itkSetMacro(RMSThreshold, OutputRealType);
        itkGetConstMacro(RMSThreshold, OutputRealType);

    protected:
        /** Constructor */
        LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion();

        /** Destructor */
        ~LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion() override = default;

        OutputRealType m_RMSThreshold;

    };
} // namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion.hxx"
#include <string>

#endif
#endif
