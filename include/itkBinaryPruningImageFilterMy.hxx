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
// Created by tabish on 2021-05-24.
//
#ifndef itkBinaryPruningImageFilterMy_hxx
#define itkBinaryPruningImageFilterMy_hxx

#include "itkBinaryPruningImageFilterMy.h"
#include "itkImageRegionIterator.h"
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>

namespace itk
{
/**
 *    Constructor
 */
template <typename TInputImage, typename TOutputImage>
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::BinaryPruningImageFilterMy()
{
  this->SetNumberOfRequiredOutputs(1);

  OutputImagePointer pruneImage = OutputImageType::New();
  this->SetNthOutput(0, pruneImage.GetPointer());

  m_Iteration = NumericTraits<unsigned>::max();
  m_ForegroundValue = NumericTraits<PixelType>::One;
  m_BackgroundValue = NumericTraits<PixelType>::Zero;
}

/**
 *  Return the pruning Image pointer
 */
template <typename TInputImage, typename TOutputImage>
typename BinaryPruningImageFilterMy<TInputImage, TOutputImage>::OutputImageType *
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::GetPruning()
{
  return dynamic_cast<OutputImageType *>(this->ProcessObject::GetOutput(0));
}

/**
 *  Prepare data for computation
 */
template <typename TInputImage, typename TOutputImage>
void
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::PrepareData()
{
  itkDebugMacro(<< "PrepareData Start");
  OutputImagePointer pruneImage = GetPruning();

  InputImagePointer inputImage = dynamic_cast<const TInputImage *>(ProcessObject::GetInput(0));

  pruneImage->SetBufferedRegion(pruneImage->GetRequestedRegion());
  pruneImage->Allocate();

  typename OutputImageType::RegionType region = pruneImage->GetRequestedRegion();

  ImageRegionConstIterator<TInputImage> it(inputImage, region);
  ImageRegionIterator<TOutputImage>     ot(pruneImage, region);

  it.GoToBegin();
  ot.GoToBegin();

  itkDebugMacro(<< "PrepareData: Copy input to output");

  while (!ot.IsAtEnd())
  {
    ot.Set(static_cast<typename OutputImageType::PixelType>(it.Get()));
    ++it;
    ++ot;
  }


  m_NeighborhoodImage= TOutputImage::New();
  m_NeighborhoodImage->SetSpacing( pruneImage->GetSpacing() );
  RegionType nbdRegion;
  SizeType nbdSize;
  nbdSize.Fill(3);
  IndexType start;
  start.Fill(-1);
  nbdRegion.SetIndex(start);
  nbdRegion.SetSize(nbdSize);
  m_NeighborhoodImage->SetRegions(nbdRegion);
  m_NeighborhoodImage->Allocate();

  itkDebugMacro(<< "PrepareData End");
}
    template <typename TInputImage, typename TOutputImage>
    unsigned
    BinaryPruningImageFilterMy<TInputImage, TOutputImage>::find(unsigned idx, std::vector<unsigned>& parent_vec){
        if(parent_vec[idx] == idx)
            return idx;
        return parent_vec[idx] = find(parent_vec[idx],parent_vec);
    }


    template <typename TInputImage, typename TOutputImage>
    bool
    BinaryPruningImageFilterMy<TInputImage, TOutputImage>::IsRemovable(NeighborhoodIteratorType& ot)
    {
        typename NeighborhoodIteratorType::NeighborhoodType pixelNeighborhood = ot.GetNeighborhood();
        std::vector<unsigned int> parent_vec(pixelNeighborhood.Size());
        RegionType validRegion = m_NeighborhoodImage->GetLargestPossibleRegion();
        OffsetType  zeroOffset,oneOffset;
        zeroOffset.Fill(0);
        oneOffset.Fill(1);
        IndexType currentIndex, zeroIndex, oneIndex;
        zeroIndex.Fill(0);
        oneIndex.Fill(1);
        typename OutputImageType::PixelType value;
        int background_parent = -1;
        //Connected Components With center removed
        for(size_t i = 0; i < pixelNeighborhood.Size(); ++i)
        {
            currentIndex.Fill(0);
            currentIndex += pixelNeighborhood.GetOffset(i);
            if(currentIndex != zeroIndex){
                value = pixelNeighborhood[i];
            }else{
                value = m_BackgroundValue;
            }

            if(value == m_BackgroundValue){
                if(background_parent == -1){
                    background_parent = i;
                }
                parent_vec[i] = background_parent;
            }else{
                parent_vec[i] = i;
            }
            m_NeighborhoodImage->SetPixel(currentIndex, value);
        }
         std::vector<OffsetType> partial_nbd;
         for(size_t d = 0; d < InputImageType::ImageDimension; ++d){
             for(int direction = -1; direction <=1 ; direction += 2){
                 OffsetType offset;
                 offset.Fill(0);
                 offset[d] = direction;
                 partial_nbd.push_back(offset);
             }
         }
        for(size_t i = 0; i < pixelNeighborhood.Size(); ++i)
        {
            currentIndex.Fill(0);
            currentIndex += pixelNeighborhood.GetOffset(i);
            if(m_NeighborhoodImage->GetPixel(currentIndex) == m_BackgroundValue) continue;
            OffsetType nbrOffset;
            IndexType nbrIndex;
            for (unsigned n = 0; n < partial_nbd.size(); ++n)
            {
                nbrIndex.Fill(0);
                nbrIndex += pixelNeighborhood.GetOffset(i); // current Index
                nbrIndex += partial_nbd[n];// nbr of current Index.
                if(nbrIndex == currentIndex || !validRegion.IsInside(nbrIndex))
                    continue;
                nbrOffset = nbrIndex - zeroIndex;
                if(m_NeighborhoodImage->GetPixel(nbrIndex) == m_ForegroundValue)
                {
                    unsigned nbrId = pixelNeighborhood.GetNeighborhoodIndex(nbrOffset);
                    //std::cout << nbrIndex <<"(" << nbrId << ") is connected to " << currentIndex<<"(" << i <<")\n";
                    unsigned currentParent = find(i,parent_vec);
                    unsigned nbrParent = find(nbrId,parent_vec);
                    if(currentParent != nbrParent)
                        parent_vec[nbrParent] = currentParent;
                }
            }
        }
        std::unordered_set<unsigned> without_components;
        for(size_t i = 0 ;i < parent_vec.size(); ++i)
        {
            without_components.insert(find(i,parent_vec));
        }
        // Connected Components with center pixel intact
        background_parent = -1;
        for(size_t i = 0; i < pixelNeighborhood.Size(); ++i)
        {
            currentIndex.Fill(0);
            currentIndex += pixelNeighborhood.GetOffset(i);
            if(pixelNeighborhood[i] == m_BackgroundValue)
            {
                if(background_parent == -1)
                    background_parent = i;
                parent_vec[i] = background_parent;
            }else
            {
                parent_vec[i] = i;
            }
            m_NeighborhoodImage->SetPixel(currentIndex, pixelNeighborhood[i]);
        }
         for(size_t i = 0; i < pixelNeighborhood.Size(); ++i)
        {
            currentIndex.Fill(0);
            currentIndex += pixelNeighborhood.GetOffset(i);
            if(m_NeighborhoodImage->GetPixel(currentIndex) == m_BackgroundValue) continue;
            OffsetType nbrOffset;
            IndexType nbrIndex;
            for (unsigned n = 0; n < pixelNeighborhood.Size(); ++n)
            {
                nbrIndex.Fill(0);
                nbrIndex += pixelNeighborhood.GetOffset(i); // current Index
                nbrIndex += pixelNeighborhood.GetOffset(n);//partial_nbd[n];
                if(nbrIndex == currentIndex || !validRegion.IsInside(nbrIndex)){
                    continue;
                }
                nbrOffset = nbrIndex - zeroIndex;
                if(m_NeighborhoodImage->GetPixel(nbrIndex) == m_ForegroundValue)
                {
                    unsigned nbrId = pixelNeighborhood.GetNeighborhoodIndex(nbrOffset);
                    unsigned currentParent = find(i,parent_vec);
                    unsigned nbrParent = find(nbrId,parent_vec);
                    if(currentParent != nbrParent)
                        parent_vec[nbrParent] = currentParent;

                }
             }
        }
        std::unordered_set<unsigned> with_components;
        for(size_t i = 0 ;i < parent_vec.size(); ++i)
        {
            with_components.insert(find(i,parent_vec));
        }
        return without_components.size() == with_components.size();
    }
/**
 *  Post processing for computing thinning
 */
template <typename TInputImage, typename TOutputImage>
void
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::ComputePruneImage()
{
  itkDebugMacro(<< "ComputeThinImage Start");
  OutputImagePointer pruneImage = GetPruning();

  typename OutputImageType::RegionType region = pruneImage->GetRequestedRegion();

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  NeighborhoodIteratorType ot(radius, pruneImage, region);

  unsigned int count = 0;
  bool changed = true;
  SizeType sz;
  sz.Fill(3);
  region.SetSize(sz);
  std::vector<IndexType> pruneList;
  if(m_Iteration < NumericTraits<unsigned>::max()){
  while(count < m_Iteration && changed){
      changed = false;
      ot.GoToBegin();
      while(!ot.IsAtEnd()){
          if(ot.GetCenterPixel() == m_ForegroundValue){
              if(IsRemovable(ot)){
                  pruneList.push_back(ot.GetIndex());
                  changed = true;
              }
          }
          ++ot;
      }
      for(auto& idx : pruneList){
          pruneImage->SetPixel(idx, static_cast<typename OutputImageType::PixelType>(m_BackgroundValue));
      }
      pruneList.clear();
      ++count;
  }
  }else{
    while (changed) {
      changed = false;
      ot.GoToBegin();
      while (!ot.IsAtEnd()) {
        if (ot.GetCenterPixel() == m_ForegroundValue) {
          if (IsRemovable(ot)) {
              pruneImage->SetPixel(ot.GetIndex(),static_cast<typename OutputImageType::PixelType>(m_BackgroundValue));
              changed = true;
          }
        }
        ++ot;
      }
      while (!ot.IsAtBegin()) {
        if (ot.GetCenterPixel() == m_ForegroundValue) {
          if (IsRemovable(ot)) {
            pruneImage->SetPixel(
                ot.GetIndex(),
                static_cast<typename OutputImageType::PixelType>(m_BackgroundValue));
            changed = true;
          }
        }
        --ot;
      }
      ++count;
    }
  }
  itkDebugMacro(<< "ComputeThinImage End after " + std::to_string(count) +" iterations");
}

/**
 *  Generate PruneImage
 */
template <typename TInputImage, typename TOutputImage>
void
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::GenerateData()
{
  this->PrepareData();

  itkDebugMacro(<< "GenerateData: Computing Thinning Image");
  this->ComputePruneImage();
} // end GenerateData()

/**
 *  Print Self
 */
template <typename TInputImage, typename TOutputImage>
void
BinaryPruningImageFilterMy<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Pruning image: " << std::endl;
  os << indent << "Iteration: " << m_Iteration << std::endl;
}
} // end namespace itk

#endif
