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
// Created by tabish on 2021-03-26.
//

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <itkLogger.h>
#include <itkStdStreamLogOutput.h>

#include "itkCommandLineArgumentParser.h"
#include "smoothing.h"
#include "slic.h"
#include "chanvese.h"
#include "experiment.h"
#include "thresholding.h"
#include "regionprops.h"
#include "turbovoxels.h"

#define DEBUG_MODE

std::string helpstring() {
    std::ostringstream ss;
    ss << "Turbo Voxels\n";
    ss << "$ turbo -module -input <input_file_path> -output <output_file_name>\n\n";
    ss << "Options:: \n";
    ss << "===========================================\n";

    ss << "\t -input\n";
    ss << "\t\t path to input file\n";

    ss << "\t -output\n";
    ss << "\t\t path to output file (default <inputFileName>_smooth.ext)\n";

    ss << "\t -h, --help\n";
    ss << "\t\t display this help\n";

    ss << "\t -mcs       ::(mean Curvature Smoothing)\n";
    ss << "\t -mmcs      ::(min-max Curvature Smoothing)\n";
    ss << "\t -bmmcs     ::(binary min-max Curvature Smoothing)\n";
    ss << "\t -slic      ::(SLIC superpixels)\n";
    ss << "\t -otsu      ::(Otsu Thresholding)\n";
    ss << "\t -bradley   ::(Bradley Adaptive Thresholding)\n";
    ss << "\t -regionprops::(Region props of objects in image)\n";
    ss << "\t -chanvese  ::(chan-Vese segmentation)\n";
    ss << "\t -turbovoxels ::(Turbovoxel superpixels)\n";

    ss << "Mean/min-max Curvature Smoothing:\n";
    ss << "--------------------------------\n";
    ss << "\t\t -iterations N : (default = 1) number of iterations to run\n";
    ss << "\t\t -dt t         : (default = 0.125) time step for smoothing\n";
    ss << "\t\t -rescale      : (optional,false) weather to rescale output to [0,1]\n";
    ss << "\t\t                 rescale changes extension of output to tif\n";
    ss << "\t\t                 if float is not supported by output file\n";
    ss << "\t\t -radius R     : (unsigned) stencil radius for min-max flow (~noise size)\n";
    ss << "\t\t -thresh th    : (double) threshold value of binary image\n";
    // -------------------------------------------------------------------
    ss << "\n\n";
    ss << "Chan-Vese Segmentation:\n";
    ss << "----------------------\n";
    ss << "\t\t -iterations N : (default = 1) number of iterations to run\n";
    ss << "\t\t -smooth N     : (default off) smooth input using min-max\n";
    ss << "\t\t                 curvature smoothing, default N = 10\n";
    ss << "\t\t -dt t         : (default = 0.125) time step for smoothing\n";
    ss << "\t\t -radius R     : (default = 2) stencil radius for min-max flow ";
    ss << "\t\t                 choose proportional to noise size\n";
    ss << "\t\t -lthresh lt   : (default = non +ve min) lower thresh for init mask\n";
    ss << "\t\t -uthresh ut   : (default = max) upper thresh for init mask\n";
    ss << "\t\t -background   : flip initial mask\n";
    ss << "\t\t -lambda l     : (default = 0.5) weight of internal term.\n";
    ss << "\t\t                  Weight of external term is fixed at 1\n";
    ss << "\t\t -debug        : saves initial mask and final levelset image\n";
    ss << "\n\n";
    // -------------------------------------------------------------------
    ss << "\n\n";
    ss << "Slic SuperPixel:\n";
    ss << "----------------------\n";
    ss << "\t\t -iterations N : (default = 10) number of iterations to run\n";
    ss << "\t\t -weight w     : (default = 1) Proximity weight for superpixels\n";
    ss << "\t\t -writelabels N: (optional, 0)Write output image labelled by "
        "superpixel id\n";
    ss << "\t\t                  0(off), 1(labelId), 2(Hot), 3(Cool), 4(Autumn), "
        "(>5)Jet\n";
    ss << "\t\t -seedcount N  : (optional, 100) approx. number of superpixels\n";

    // -------------------------------------------------------------------
    ss << "\n\n";
    ss << "Otsu Thresholding:\n";
    ss << "----------------------\n";
    ss << "\t\t -iterations N : (optional, 0) number of iterations of pre-process smoothing\n";
    ss << "\t\t -radius w     : (optional, 2) stencil radius for smoothing\n";
    ss << "\t\t -dt t         : (optional, 0.125) time step for pre-smoothing\n";
    ss << "\t\t -dark         : (optional,false) flip polarity of output\n";
    // -------------------------------------------------------------------
    ss << "\n\n";
    ss << "Bradley Adaptive Thresholding:\n";
    ss << "----------------------\n";
    ss << "\t\t -sensitivity s : (optional, 0.5) foreground sensitivity\n";
    ss << "\t\t -window w      : (optional, max_image_dim/8) size of window for "
        "estimating local threshold\n";
    ss << "\t\t -dark         : (optional,false) flip polarity of output\n";
    ss << "\t\t -writeThreshold: (optional, false) write the threshold image used for thresholding\n";
    // -------------------------------------------------------------------
    ss << "\n\n";
    ss << "Region props:\n";
    ss << "----------------------\n";
    ss << "\t\t -size s : (optional, 100) size(pixels) of smallest object to consider\n";
    ss << "\t\t -dark   : (optional, false) flip image intensity\n(untested)\n";
    //------------------------------------------------------------------------
    ss << "\n\n";
    ss << "Turbo Voxels:\n";
    ss << "----------------------\n";
    ss << "\t\t -seedcount      : (optional, 150) Number of turbovoxels to seed."
          "Can be overriden if -size is specified\n";
    ss << "\t\t -smoothing      : (optional, 10) smoothing iteration for speed computation\n";
    ss << "\t\t -propagation    : (optional, 1) weight of the constant prop. term in level set evoluiton\n";
    ss << "\t\t -curvature      : (optional, 0.3) weight of the curvature term in level set evoluiton\n";
    ss << "\t\t -advection      : (optional, 5) weight of the doublet (advection) term in level set evoluiton\n";
    ss << "\t\t -iterations    : (optional, 1000) Maximum number of iteration to run for."
          " Will stop early if rms difference between updates is below -rmsthreshold\n";
    ss << "\t\t -rmsthreshold  : (optional, 1) Stopping threshold for RMS difference between updates. \n";
    ss << "\t\t -dt             : (optional, _) Value of timestep for evolution."
          " Should not be set manually unless you know what you are doing\n";
    ss << "\t\t -observe        : (optional[debug], false) Runs debug Evolution observer writing intermediate images and log\n";
    ss << "\t\t -rgb            : (optional[debug], false) Saves pseudo-color rgb images inside evolution observer \n";
    //------------------------------------------------------------------------

    ss << "\n\n";
    ss << "Examples:\n";
    ss << "=========\n";
    ss << "SLIC superpixel example:\n";
    ss << "$./turbo -slic -input lizard.jpg -writelabels 5 -seedcount 250\n";

    ss << "Chan-Vese segmentation example:\n";
    ss << "$ ./turbo -chanvese -lthresh 1 -uthresh 1000 -smooth -input ../data/dapi3d.tif\n";

    ss << "Perform mean curvature smoothing:\n";
    ss << "$./turbo -mcs -dt 0.05 -iterations 200  -input ../data/EM2D.tif  \n";

    ss << "Turbovoxel example:\n";
    ss << "$./turbo -turbovoxels -input ../data/wga3d.tif -seedcount 500 -iterations 500 -rmsthreshold 0.1 -edgescale 25 -curvature 0.3 -advection 0.2\n";
    ss << "$./turbo -turbovoxels -input ../data/lizard.jpg -seedcount 500 -gradientseeding -iterations 500 -edgescale 12.5 -curvature 0.3 -advection 0.2\n" ;
    ss << "\n------------------------------------------";
    ss << "\n\n";

    std::string helpString(ss.str());
    return helpString;
}

int main(int argc, char* argv[]){
    itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();
    itk::Logger::Pointer logger = itk::Logger::New();
    itk::StdStreamLogOutput::Pointer itkcout = itk::StdStreamLogOutput::New();
    itk::StdStreamLogOutput::Pointer itkfout = itk::StdStreamLogOutput::New();
    itkcout->SetStream(std::cout);
    logger->SetLevelForFlushing(itk::LoggerBaseEnums::PriorityLevel::DEBUG);

#ifdef DEBUG_MODE
    logger->SetPriorityLevel(itk::LoggerBaseEnums::PriorityLevel::DEBUG);
#endif

#ifdef RELEASE_MODE
    logger->SetPriorityLevel(itk::LoggerBaseEnums::PriorityLevel::CRITICAL);
#endif
    logger->AddLogOutput(itkcout);
    std::string humanReadableFormat = "[%b-%d-%Y, %H:%M:%S]";
    logger->SetHumanReadableFormat(humanReadableFormat);
    logger->SetTimeStampFormat(itk::LoggerBaseEnums::TimeStampFormat::HUMANREADABLE);

    // map to list module code, module name.
    std::map<int, std::string> modesMap;
    modesMap[0] = "Mean Curvature Smoothing\n";
    modesMap[1] = "Min-Max Curvature Smoothing\n";
    modesMap[2] = "Binary Min-Max Curvature Smoothing\n";
    modesMap[3] = "Chan Vese segmentation\n";
    modesMap[4] = "Slic SuperPixel\n";
    modesMap[5] = "Otsu Thresholding\n";
    modesMap[6] = "Bradley Adaptive Thresholding\n";
    modesMap[7] = "Region Props\n";
    modesMap[8] = "Turbo Voxels\n";
    modesMap[93] = "Undocumented/Experimental\n";

    std::vector<std::string> modules;
    modules.emplace_back("-mcs");
    modules.emplace_back("-mmcs");
    modules.emplace_back("-bmmcs");
    modules.emplace_back("-chanvese");
    modules.emplace_back("-slic");
    modules.emplace_back("-otsu");
    modules.emplace_back("-bradley");
    modules.emplace_back("-regionprops");
    modules.emplace_back("-turbovoxels");
    modules.emplace_back("-experiment");

    parser->SetCommandLineArguments(argc, argv);
    parser->SetProgramHelpText(helpstring());

    if(parser->CheckForRequiredArguments() == itk::CommandLineArgumentParser::HELPREQUESTED){
        return EXIT_SUCCESS;
    }
    if(parser->CheckForRequiredArguments() == itk::CommandLineArgumentParser::FAILED ||
        !parser->ExactlyOneExists(modules)){
        std::cout << parser->GetProgramHelpText() << std::endl;
        return EXIT_FAILURE;
    }
    std::string logFileName="/tmp/turbovoxels.log";

    if(parser->GetCommandLineArgument("-logfile", logFileName)){
        std::ofstream fileStream(logFileName);
        itkfout->SetStream(fileStream);
        logger->AddLogOutput(itkfout);
    }
    logger->Info( "Starting Turbo Voxels\n");

    int module = 0;
    if(parser->ArgumentExists("-mcs")){
        module = 0;
        logger->Info("Selected : " + modesMap[module] + "\n");
    }else if(parser->ArgumentExists("-mmcs")){
        module = 1;
        logger->Info("Selected : " + modesMap[module] + "\n");
    }else if (parser->ArgumentExists("-bmmcs")) {
      module = 2;
      logger->Info("Selected : " + modesMap[module] + "\n");
    } else if (parser->ArgumentExists("-chanvese")) {
      module = 3;
      logger->Info("Selected : " + modesMap[module] + "\n");
    }else if (parser->ArgumentExists("-slic")) {
        module = 4;
        logger->Info("Selected : " + modesMap[module] + "\n");
    } else if (parser->ArgumentExists("-otsu")) {
      module = 5;
      logger->Info("Selected : " + modesMap[module] + "\n");
    } else if (parser->ArgumentExists("-bradley")) {
      module = 6;
      logger->Info("Selected : " + modesMap[module] + "\n");
    } else if (parser->ArgumentExists("-regionprops")) {
      module = 7;
      logger->Info("Selected : " + modesMap[module] + "\n");
    } else if (parser->ArgumentExists("-turbovoxels")) {
        module = 8;
        logger->Info("Selected : " + modesMap[module] + "\n");
    }else if (parser->ArgumentExists("-experiment")) {
      module = 93;
      logger->Info("Selected : " + modesMap[module] + "\n");
    } else {
      logger->Info("Using default : " + modesMap[module] + "\n");
    }
    switch(module){
        case 0: // mean curvature smoothing
        case 1: // min-max curvature smoothing
        case 2: // binary min-max curvature smoothing
            curvature_smoothing(parser, logger);
          break;
        case 3:
            chanvese_segmentation(parser, logger);
            break;
        case 4:
            slic(parser, logger);
            break;
        case 5: // otsu thresholding
        case 6: // bradley adaptive thresholding
            thresholding(parser, logger);
            break;
        case 7:
            regionprops(parser, logger);
            break;
        case 8:
            turbovoxels(parser,logger);
            break;
        case 93:
            experiment(parser, logger);
            break;
        default:
            std::cout<<"Unknown Module"<<std::endl;
    }
    return EXIT_SUCCESS;
}
