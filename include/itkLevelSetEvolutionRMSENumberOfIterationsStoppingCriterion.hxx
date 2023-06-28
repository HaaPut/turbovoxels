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

#ifndef itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion_hxx
#define itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion_hxx

#include "itkLevelSetEvolutionRMSENumberOfIterationsStoppingCriterion.h"
#include <string>
#include <limits>

namespace itk
{
    template <typename TLevelSetContainer>
    LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion<TLevelSetContainer>::LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion(){
        this->m_RMSThreshold = std::numeric_limits<OutputRealType>::epsilon();
    }

    template <typename TLevelSetContainer>
    bool
    LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion<TLevelSetContainer>::IsSatisfied() const{
        return (this->m_CurrentIteration >= this->m_NumberOfIterations || this->m_RMSChangeAccumulator < m_RMSThreshold);
    }

    template <typename TLevelSetContainer>
    std::string
    LevelSetEvolutionRMSENumberOfIterationsStoppingCriterion<TLevelSetContainer>::GetDescription() const{
        return "Current Iteration Number >= Number Of Iterations";
    }

} // namespace itk
#endif
