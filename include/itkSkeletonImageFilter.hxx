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
// Created by tabish on 10-05-2021.
//

#ifndef itkSkeletonImageFilter_hxx
#define itkSkeletonImageFilter_hxx

#include "itkSkeletonImageFilter.h"
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <itkImageFileWriter.h>
#include <itkApproximateSignedDistanceMapImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <string>
#include <functional>

namespace itk{

    template< class PixelType>
    void
    SkeletonImageFilter<PixelType , 3>::GenerateInputRequestedRegion()
    {
        // call the superclass' implementation of this method
        Superclass::GenerateInputRequestedRegion();

        // get pointers to the input and output
        InputImagePointer  inputPtr = const_cast<InputImageType *>(this->GetInput());
        OutputImagePointer outputPtr = this->GetOutput();

        if (!inputPtr || !outputPtr)
        {
            return;
        }

        // get a copy of the input requested region (should equal the output
        // requested region)
        typename InputImageType::RegionType inputRequestedRegion = inputPtr->GetRequestedRegion();

        // pad the input requested region by one, the size of neighbourhood for labelling
        inputRequestedRegion.PadByRadius(1);

        // crop the input requested region at the input's largest possible region
        if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
        {
            inputPtr->SetRequestedRegion(inputRequestedRegion);
            return;
        }
        else
        {
            // Couldn't crop the region (requested region is outside the largest
            // possible region).  Throw an exception.

            // store what we tried to request (prior to trying to crop)
            inputPtr->SetRequestedRegion(inputRequestedRegion);

            // build an exception
            InvalidRequestedRegionError e(__FILE__, __LINE__);
            e.SetLocation(ITK_LOCATION);
            e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
            e.SetDataObject(inputPtr);
            throw e;
        }
    }

    template< class PixelType>
    unsigned
    SkeletonImageFilter<PixelType , 3>::BackgroundLabelling(IndexType index){
        //OutputImageType * output = this->GetOutput();

        unsigned regions = 0;
        std::queue<size_t> Q;
        std::vector<bool> visited(18,false);
        for(size_t i = 0; i < m_Neighbors18.size(); ++i){ // starting point
            if( (m_n6[i]) && (! visited[i]) &&                       // not already visited
                (m_Accessor.GetPixel(index + m_Neighbors18[i], m_Output) <= NumericTraits<PixelType>::ZeroValue())){     // is outside object
                ++regions; // new component
                Q.push(i);
                visited[i] = true;
                while(Q.size() > 0){
                    size_t current = Q.front();
                    Q.pop();
                    visited[current] = true;
                    for(size_t neighbor : m_Graph18[current]){
                        --neighbor;
                        if((!visited[neighbor]) &&
                           (m_Accessor.GetPixel(index + m_Neighbors18[neighbor], m_Output) <= NumericTraits<PixelType>::ZeroValue())){
                            visited[neighbor] = true;
                            Q.push(neighbor);
                        }
                    }
                }

            }
        }
        return regions;
    }

    template< class PixelType>
    unsigned
    SkeletonImageFilter<PixelType , 3>::ForegroundLabelling(IndexType index){
        //OutputImageType * output = this->GetOutput();

        unsigned regions = 0;
        std::queue<size_t> Q;
        std::vector<bool> visited(26,false);
        for(size_t i = 0; i < m_Neighbors26.size(); ++i){ // starting point
            if( (! visited[i]) &&                       // not already visited
                (m_Accessor.GetPixel(index + m_Neighbors26[i], m_Output) > NumericTraits<PixelType>::ZeroValue())){     // is in object
                ++regions; // new component
                Q.push(i);
                visited[i] = true;
                while(Q.size() > 0){
                    size_t current = Q.front();
                    Q.pop();
                    visited[current] = true;
                    for(size_t neighbor : m_Graph26[current]){
                        --neighbor;
                        if((!visited[neighbor]) &&
                           (m_Accessor.GetPixel(index + m_Neighbors26[neighbor], m_Output) > NumericTraits<PixelType>::ZeroValue())){
                            visited[neighbor] = true;
                            Q.push(neighbor);
                        }
                    }
                }

            }
        }
        return regions;
    }

    template< class PixelType>
    unsigned
    SkeletonImageFilter<PixelType , 3>::TopologicalLabel(IndexType index){
        unsigned label;
        unsigned Cstar = this->ForegroundLabelling(index);
        unsigned Cbar = this->BackgroundLabelling(index);

        if (Cbar==0)
            label =   2 ; //interior point
        else if (Cstar==0)
            label =   3; //isolated point
        else if ((Cbar==1)&&(Cstar==1))
            label =   4; //simple point
        else if ((Cbar==1)&&(Cstar==2))
            label =   5; //candidate curve point
        else if ((Cbar==1)&&(Cstar>2))
            label =   6; //junction of curves
        else if ((Cbar==2)&&(Cstar==1))
            label =   7; //candidate surface point
        else if ((Cbar==2)&&(Cstar>=2))
            label =   8;  //junction between curve(s) and surface
        else if ((Cbar>2)&&(Cstar==1))
            label =   9; //junction of surfaces
        else if ((Cbar>2)&&(Cstar>=2))
            label =  10; //junction between surface(s) and curve
        else
            label = 0;

        return label;
    }

    template< class PixelType>
    void
    SkeletonImageFilter<PixelType , 3>::GenerateData() {
        this->AllocateOutputs();
        InputImagePointer input = const_cast<InputImageType *>(this->GetInput(0));
        this->m_Output = this->GetOutput(0);
        this->m_RemoveCount = 0;
        this->m_Count = 0;
        const InternalPixelType maximumDistance = (static_cast<InternalPixelType>(m_MaxIterations))*m_MinSpacing;

        input->SetRequestedRegionToLargestPossibleRegion();

        m_Output->ReleaseDataFlagOff();

        InputConstIteratorType skit = InputConstIteratorType(input, input->GetRequestedRegion());
        OutputIteratorType outIt = OutputIteratorType(m_Output, m_Output->GetRequestedRegion());

        PixelType value;
        skit.GoToBegin();
        outIt.GoToBegin();
        while(!skit.IsAtEnd()){
            value = skit.Get();
            if(value >= this->m_LowerThreshold)
                outIt.Set(m_InsideValue);
                //BUG: potential overflow for not float types
                //outIt.Set(maximumDistance+10*m_MinSpacing);
            else
                outIt.Set(m_OutsideValue);
            ++skit;
            ++outIt;
        }

        m_ThresholdFilter->SetInput(input);
        m_ThresholdFilter->SetLowerThreshold( this->m_LowerThreshold);
        m_ThresholdFilter->SetUpperThreshold(NumericTraits<PixelType>::max());
        m_ThresholdFilter->SetInsideValue(maximumDistance+10*m_MinSpacing);
        m_ThresholdFilter->SetOutsideValue(NumericTraits<InternalPixelType>::ZeroValue());

        m_ChamferFilter->SetInput(m_ThresholdFilter->GetOutput());
        m_ChamferFilter->SetMaximumDistance(maximumDistance+m_MinSpacing);
        m_ChamferFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        m_ChamferFilter->Update();
        auto distanceMap = m_ChamferFilter->GetOutput();

        using DistIteratorType = ImageRegionIteratorWithIndex< InternalImageType >;
        DistIteratorType dIt = DistIteratorType(distanceMap, distanceMap->GetRequestedRegion());

        using NodeType = std::pair<itk::Index<3>, float>;
        const auto cmp = [](NodeType const& left, NodeType const& right) { return left.second > right.second; };
        std::priority_queue<NodeType, std::vector<NodeType>, decltype(cmp)> q(cmp);

        dIt.GoToBegin();
        while(!dIt.IsAtEnd()) {
            if(dIt.Get() > 0){
                q.push({dIt.GetIndex(), std::abs(dIt.Get())});
            }
            ++dIt;
        }

        float current_distance = 0;
        IndexType current_index;
        unsigned label;
        while(current_distance < maximumDistance && !q.empty()){
            NodeType current_node = q.top();
            current_distance = current_node.second;
            current_index = current_node.first;
            q.pop();
            label = this->TopologicalLabel(current_index);
            if(label == 4){
                m_Output->SetPixel(current_index, m_OutsideValue);
                ++this->m_RemoveCount;
            }
            ++this->m_Count;
        }
    }

    template< class PixelType>
    bool
    SkeletonImageFilter<PixelType , 2>::isSimple2(IndexType index){
        std::array<int,8> nbrs = {{0,0,0,0,0,0,0,0}};
        int numNeighbors = 0;
        int numEdges = 0;
        for(size_t i = 0 ;i < m_Neighbors8.size(); ++i){
            nbrs[i] = m_Accessor.GetPixel(index + m_Neighbors8[i], m_Output) > 0;
            size_t j = (i+1)%8;
            nbrs[j] = m_Accessor.GetPixel(index + m_Neighbors8[j], m_Output) > 0;
            if(nbrs[i] == 1 && nbrs[j] == 1){
                numNeighbors+=2;
                ++numEdges;
            }else if(nbrs[i] == 1 || nbrs[j] == 1){
                ++numNeighbors;
            }
        }
        //remove double counted
        numNeighbors /=2;
        //add corner diagonals if corner is 0
        for (size_t i = 0; i < 8; i+=2){
            numEdges += (nbrs[(i+7)%8] == 1 && nbrs[i] == 0 && nbrs[(i+1)%8] == 1)? 1: 0;
        }
        return (numNeighbors - numEdges == 1);
    }

    template< class PixelType>
    void
    SkeletonImageFilter<PixelType , 2>::GenerateData() {
        this->AllocateOutputs();
        auto input = const_cast<InputImageType *>(this->GetInput(0));
        this->m_Output = this->GetOutput(0);
        this->m_RemoveCount = 0;
        this->m_Count = 0;
        const InternalPixelType maximumDistance = (static_cast<InternalPixelType>(m_MaxIterations))*m_MinSpacing;
        input->SetRequestedRegionToLargestPossibleRegion();

        //m_Output->ReleaseDataFlagOff();
        InputConstIteratorType skit = InputConstIteratorType(input, input->GetRequestedRegion());
        OutputIteratorType outIt = OutputIteratorType(m_Output, m_Output->GetRequestedRegion());

        PixelType value;
        skit.GoToBegin();
        outIt.GoToBegin();
        while(!skit.IsAtEnd()){
            value = skit.Get();
            if(value >= this->m_LowerThreshold)
                outIt.Set(m_InsideValue);
            else
                outIt.Set(m_OutsideValue);
            ++skit;
            ++outIt;
        }

        m_ThresholdFilter->SetInput(input);
        m_ThresholdFilter->SetLowerThreshold( this->m_LowerThreshold);
        m_ThresholdFilter->SetUpperThreshold(NumericTraits<PixelType>::max());
        m_ThresholdFilter->SetInsideValue(maximumDistance+10*m_MinSpacing);
        m_ThresholdFilter->SetOutsideValue(NumericTraits<InternalPixelType>::ZeroValue());

        m_ChamferFilter->SetInput(m_ThresholdFilter->GetOutput());
        m_ChamferFilter->SetMaximumDistance(maximumDistance+1);
        m_ChamferFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        // Create the distance map
        m_ChamferFilter->Update();
        auto distanceMap = m_ChamferFilter->GetOutput();

        using DistIteratorType = ImageRegionIteratorWithIndex< InternalImageType >;
        DistIteratorType dIt = DistIteratorType(distanceMap, distanceMap->GetRequestedRegion());

        using NodeType = std::pair<itk::Index<2>, float>;
        const auto cmp = [](NodeType left, NodeType right) { return left.second > right.second; };
        std::priority_queue<NodeType, std::vector<NodeType>, decltype(cmp)> q(cmp);

        dIt.GoToBegin();
        while(!dIt.IsAtEnd()) {
            if(dIt.Get() > 0){
                q.push({dIt.GetIndex(), dIt.Get()});
            }
            ++dIt;
        }

        float current_distance = 0;
        IndexType current_index;
        while(current_distance < maximumDistance && !q.empty()){
            NodeType current_node = q.top();
            current_distance = current_node.second;
            current_index = current_node.first;
            q.pop();
            if(isSimple2(current_index)){
                m_Output->SetPixel(current_index,m_OutsideValue);
                ++this->m_RemoveCount;
            }
            ++this->m_Count;
        }
    }
/**
 *  Print Self
 */
    template< class PixelType , unsigned Dimension>
    void
    SkeletonImageFilter<PixelType , Dimension>::PrintSelf(std::ostream& os, Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
        os << indent << "Skeleton Filter in dimension" << Dimension << std::endl;
    }
} // namespace itk

#endif //itkSkeletonImageFilter_hxx