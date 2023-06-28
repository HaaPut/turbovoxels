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

#ifndef itkSkeletonImageFilter_h
#define itkSkeletonImageFilter_h

#include "itkBinaryThinningIterationImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkConstantBoundaryCondition.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageToImageFilter.h>
#include <itkProgressAccumulator.h>
#include <itkFastChamferDistanceImageFilter.h>
#include <vector>
#include <map>
#include <itkMacro.h>

/** TODO: m_MaxIterations souldn't be repeated in specialization as it is now!*/

namespace itk
{

    template< typename PixelType, unsigned Dimension>
    class ITK_TEMPLATE_EXPORT  SkeletonImageFilter:
        public ImageToImageFilter<Image<PixelType, Dimension>, Image<PixelType, Dimension> >
    {
    public:
        /** Standard class typedefs. */
        using InputImageType = Image<PixelType, Dimension>;
        using Self = SkeletonImageFilter;
        using Superclass = ImageToImageFilter<InputImageType, InputImageType>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        /** Method for creation through the object factory */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro( SkeletonImageFilter, ImageToImageFilter );

        using InputPointerType = typename InputImageType::ConstPointer;

        using OutputImageType = InputImageType;
        using OutputPointerType = typename InputImageType::Pointer;

        using InputConstIteratorType = ImageRegionConstIterator< InputImageType >;
        using OutputIteratorType = ImageRegionIteratorWithIndex< InputImageType >;
        using IndexType = typename InputImageType::IndexType;


    protected:
        SkeletonImageFilter() = default;
        ~SkeletonImageFilter() = default;


        void GenerateData() override;

        void PrintSelf(std::ostream& os, Indent indent) const;
    };

    //---------------------------------------------------------------------
    template< typename PixelType>
    class ITK_TEMPLATE_EXPORT  SkeletonImageFilter<PixelType, 3>:
        public ImageToImageFilter<Image<PixelType, 3>, Image<PixelType, 3> >
    {
    public:
        /** Standard class typedefs. */
        using InputImageType = Image<PixelType, 3>;
        using Self = SkeletonImageFilter;
        using Superclass = ImageToImageFilter<InputImageType, InputImageType>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        /** Method for creation through the object factory */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro( SkeletonImageFilter, ImageToImageFilter );

        using InputImagePointer = typename Superclass::InputImagePointer;
        using OutputImageType = InputImageType;
        using OutputImagePointer = typename InputImageType::Pointer;

        using InputConstIteratorType = ImageRegionConstIterator< InputImageType >;
        using OutputIteratorType = ImageRegionIteratorWithIndex< OutputImageType >;
        using IndexType = typename InputImageType::IndexType;
        using BoundaryConditionType = ConstantBoundaryCondition<InputImageType>;

        using InternalPixelType = float;
        using InternalImageType = Image<InternalPixelType,3>;
        using ChamferThresholdFilterType = itk::BinaryThresholdImageFilter<InputImageType, InternalImageType>;
        using ChamferType = FastChamferDistanceImageFilter<InternalImageType, InternalImageType>;

        void
        GenerateInputRequestedRegion() override;

        itkSetMacro(MaxIterations, unsigned);
        itkGetConstMacro(MaxIterations, unsigned);

        itkSetMacro(LowerThreshold, PixelType);
        itkGetConstMacro(LowerThreshold, PixelType);

        itkGetConstMacro(RemoveCount, double);

        itkGetConstMacro(InsideValue, PixelType);

        itkGetConstMacro(OutsideValue, PixelType);
      protected:
        SkeletonImageFilter(){
            m_Accessor.SetConstant(NumericTraits<PixelType>::OneValue());
            //this->DebugOn();
            m_RemoveCount = 0;
            m_InsideValue = NumericTraits<PixelType>::OneValue();
            m_OutsideValue = NumericTraits<PixelType>::ZeroValue();
            m_ChamferFilter = ChamferType::New();
            m_ThresholdFilter = ChamferThresholdFilterType::New();
        }
        ~SkeletonImageFilter() = default;


        void GenerateData() override;

    private:
        unsigned m_MaxIterations = 4;
        PixelType m_LowerThreshold = NumericTraits<PixelType>::ZeroValue();
        float m_MinSpacing = 1;
        //bool m_SaveProgressImage = true;
        typename ChamferType::Pointer m_ChamferFilter;
        typename ChamferThresholdFilterType::Pointer m_ThresholdFilter;

        OutputImagePointer m_Output;
        PixelType m_InsideValue;
        PixelType m_OutsideValue;
        double m_RemoveCount;
        double m_Count;
        unsigned TopologicalLabel(IndexType index);
        unsigned ForegroundLabelling(IndexType index);
        unsigned BackgroundLabelling(IndexType index);
        const std::vector<std::string> named_labels = {"unknown", "unassigned value", "interior point", "isolated point",
                                                       "simple point","candidate curve point","junction of curves",
                                                       "candidate surface point", " junction between curve(s] and surface",
                                                       " junction of surfaces","junction between surface(s) and curve" };

        const std::vector<bool> m_n6 = {false, false, true,  false, false, false,
                                        true,  false, true,  true,  false, true,
                                        false, false, false, true,  false, false};
        const std::vector<itk::Offset<3>> m_Neighbors18 = {
                {{-1, -1, 0}}, {{-1, 0, -1}}, {{-1, 0, 0}}, {{-1, 0, 1}},
                {{-1, 1, 0}},  {{0, -1, -1}}, {{0, -1, 0}}, {{0, -1, 1}},
                {{0, 0, -1}},  {{0, 0, 1}},   {{0, 1, -1}}, {{0, 1, 0}},
                {{0, 1, 1}},   {{1, -1, 0}},  {{1, 0, -1}}, {{1, 0, 0}},
                {{1, 0, 1}},   {{1, 1, 0}}};

        const std::vector<itk::Offset<3>> m_Neighbors26 = {
                {{-1, -1, -1}}, {{-1, -1, 0}}, {{-1, -1, 1}}, {{-1, 0, -1}},
                {{-1, 0, 0}},   {{-1, 0, 1}},  {{-1, 1, -1}}, {{-1, 1, 0}},
                {{-1, 1, 1}},   {{0, -1, -1}}, {{0, -1, 0}},  {{0, -1, 1}},
                {{0, 0, -1}},   {{0, 0, 1}},   {{0, 1, -1}},  {{0, 1, 0}},
                {{0, 1, 1}},    {{1, -1, -1}}, {{1, -1, 0}},  {{1, -1, 1}},
                {{1, 0, -1}},   {{1, 0, 0}},   {{1, 0, 1}},   {{1, 1, -1}},
                {{1, 1, 0}},    {{1, 1, 1}}};
        const std::vector<std::vector<size_t>> m_Graph26 = {
                {1, 2, 4, 5, 10, 11, 13},
                {1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14},
                {2, 3, 5, 6, 11, 12, 14},
                {1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 16},
                {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
                {2, 3, 5, 6, 8, 9, 11, 12, 14, 16, 17},
                {4, 5, 7, 8, 13, 15, 16},
                {4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17},
                {5, 6, 8, 9, 14, 16, 17},
                {1, 2, 4, 5, 10, 11, 13, 18, 19, 21, 22},
                {1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23},
                {2, 3, 5, 6, 11, 12, 14, 19, 20, 22, 23},
                {1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 16, 18, 19, 21, 22, 24, 25},
                {2, 3, 5, 6, 8, 9, 11, 12, 14, 16, 17, 19, 20, 22, 23, 25, 26},
                {4, 5, 7, 8, 13, 15, 16, 21, 22, 24, 25},
                {4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 21, 22, 23, 24, 25, 26},
                {5, 6, 8, 9, 14, 16, 17, 22, 23, 25, 26},
                {10, 11, 13, 18, 19, 21, 22},
                {10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23},
                {11, 12, 14, 19, 20, 22, 23},
                {10, 11, 13, 15, 16, 18, 19, 21, 22, 24, 25},
                {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
                {11, 12, 14, 16, 17, 19, 20, 22, 23, 25, 26},
                {13, 15, 16, 21, 22, 24, 25},
                {13, 14, 15, 16, 17, 21, 22, 23, 24, 25, 26},
                {14, 16, 17, 22, 23, 25, 26}};
        const std::vector<std::vector<size_t>> m_Graph18 = {
                {1, 3, 7},         {2, 3, 9},
                {1, 2, 3, 4, 5},   {3, 4, 10},
                {3, 5, 12},        {6, 7, 9},
                {1, 6, 7, 8, 14},  {7, 8, 10},
                {2, 6, 9, 11, 15}, {4, 8, 10, 13, 17},
                {9, 11, 12},       {5, 11, 12, 13, 18},
                {10, 12, 13},      {7, 14, 16},
                {9, 15, 16},       {14, 15, 16, 17, 18},
                {10, 16, 17},      {12, 16, 18}};
        BoundaryConditionType  m_Accessor;
    };

//--------------------------------------------------------------------------------------
    template <typename PixelType>
    class ITK_TEMPLATE_EXPORT SkeletonImageFilter<PixelType, 2>
            : public ImageToImageFilter<Image<PixelType, 2>, Image<PixelType, 2>> {
    public:
        /** Standard class typedefs. */
        using InputImageType = Image<PixelType, 2>;
        using Self = SkeletonImageFilter;
        using Superclass = ImageToImageFilter<InputImageType, InputImageType>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        /** Method for creation through the object factory */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro( SkeletonImageFilter, ImageToImageFilter );

        //using InputPointerType = typename InputImageType::ConstPointer;

        using OutputImageType = InputImageType;
        using OutputImagePointer = typename InputImageType::Pointer;

        using InputConstIteratorType = ImageRegionConstIterator< InputImageType >;
        using OutputIteratorType = ImageRegionIteratorWithIndex< InputImageType >;
        using IndexType = typename InputImageType::IndexType;
        using BoundaryConditionType = ConstantBoundaryCondition<InputImageType>;

        using InternalPixelType = float;
        using InternalImageType = Image<InternalPixelType,2>;
        using ChamferThresholdFilterType = itk::BinaryThresholdImageFilter<InputImageType, InternalImageType>;
        using ChamferType = FastChamferDistanceImageFilter<InternalImageType, InternalImageType>;

        itkSetMacro(MaxIterations, unsigned);
        itkGetConstMacro(MaxIterations, unsigned);

        itkSetMacro(LowerThreshold, PixelType);
        itkGetConstMacro(LowerThreshold, PixelType);

        itkGetConstMacro(RemoveCount, double);
        itkGetConstMacro(Count, double);

        itkGetConstMacro(InsideValue, PixelType);

        itkGetConstMacro(OutsideValue, PixelType);

    protected:
        SkeletonImageFilter() {
            m_Accessor.SetConstant(NumericTraits<PixelType>::OneValue());

            m_InsideValue = NumericTraits<PixelType>::OneValue();
            m_OutsideValue = NumericTraits<PixelType>::ZeroValue();
            m_ChamferFilter = ChamferType::New();
            m_ThresholdFilter = ChamferThresholdFilterType::New();
        }

        ~SkeletonImageFilter() = default;

        void GenerateData() override;
        bool isSimple2(IndexType index);

    private:
        //bool m_SaveProgressImage = true;
        typename ChamferType::Pointer m_ChamferFilter;
        typename ChamferThresholdFilterType::Pointer m_ThresholdFilter;
        PixelType m_LowerThreshold =
                NumericTraits<PixelType>::ZeroValue();
        unsigned m_MaxIterations = 4;
        float m_MinSpacing = 1;
        double m_RemoveCount = 0;
        double m_Count = 0;
        PixelType m_InsideValue;
        PixelType m_OutsideValue;

        OutputImagePointer m_Output;

        //order is important here
        //0 1 2
        //7 x 3
        //6 5 4
        const std::vector<itk::Offset<2>> m_Neighbors8 = {
                {{-1, -1}},
                {{-1, 0}},
                {{-1, 1}},
                {{0,  1}},
                //{{0,  0}},
                {{1,  1}},
                {{1,  0}},
                {{1,  -1}},
                {{0,  -1}}
        };
        BoundaryConditionType m_Accessor;
    };
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSkeletonImageFilter.hxx"
#endif

#endif //itkSkeletonImageFilter_h