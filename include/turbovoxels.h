//**********************************************************
//Copyright 2022 Tabish Syed
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
// Created by tabish on 07/10/22.
//

#ifndef TURBO_TURBOVOXELS_H
#define TURBO_TURBOVOXELS_H

#include <set>

#include <itkLogger.h>
#include <itkIndex.h>

#include "filesystem.hpp"
#include "itkCommandLineArgumentParser.h"

namespace fs = std::filesystem;

int turbovoxels(const itk::CommandLineArgumentParser::Pointer &parser,
                const itk::Logger::Pointer &logger);


template<unsigned Dimension>
class turboVoxel{
    struct indexLessCmp
    {
        bool operator() (const itk::Index<Dimension>& i1, const itk::Index<Dimension>& i2) const{
            for(size_t d = 0 ; d < Dimension; ++d){
                if(i1[d] != i2[d])
                    return i1[d] < i2[d];
            }
            return false;
        }
    };
public:
    turboVoxel(){
        id = 0;
    }
    turboVoxel(size_t _id): id{_id}{ };
    size_t id;
    std::set<itk::Index<Dimension>, indexLessCmp> pixels;
    std::set<size_t> neighbors;
    void merge(turboVoxel<Dimension>& other){
        std::merge(pixels.begin(),pixels.end(),other.pixels.begin(),other.pixels.end(),
                   std::inserter(pixels,pixels.begin()));
        std::merge(neighbors.begin(),neighbors.end(),other.neighbors.begin(),other.neighbors.end(),
                   std::inserter(neighbors,neighbors.begin()));
        neighbors.erase(other.id);
    }
    /*
     std::set<Index,indexLessCmp> boundary
     merge -> union(b1,b2) - intersection(b1,b2);
     */
    size_t size(){return pixels.size();}
};


#endif //TURBO_TURBOVOXELS_H
