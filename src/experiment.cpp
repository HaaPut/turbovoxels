#include "experiment.h"
#include "filesystem.hpp"
#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <string>
#include <vector>
#include <itkTimeProbe.h>


namespace fs = std::filesystem;

itk::TimeProbe itkclock;


template <typename InputPixelType, unsigned int Dimension>
int experiment_impl(const itk::CommandLineArgumentParser::Pointer &parser,
                    const itk::Logger::Pointer &logger) {

    return EXIT_SUCCESS;
}

int experiment(const itk::CommandLineArgumentParser::Pointer &parser,
                 const itk::Logger::Pointer &logger) {

    logger->Info("Experiment\n");
    std::string inputFileName, outputFileName;
    fs::path inputFilePath, outputFilePath;
    bool hasInputFilePath =
        parser->GetCommandLineArgument("-input", inputFileName);
    bool hasOutputFilePath =
        parser->GetCommandLineArgument("-output", outputFileName);
    // Read the input image.
    if (!hasInputFilePath) {
        logger->Error("Input file not specified\n");
        return EXIT_FAILURE;
    }else{
        inputFilePath = inputFileName;
    }
    if (!hasOutputFilePath) {
        outputFilePath = inputFilePath.parent_path()/inputFilePath.stem();
        outputFilePath += "_expt";
        outputFilePath += inputFilePath.extension();
        std::vector<std::string> newargs;
        newargs.push_back("-output");
        newargs.push_back(outputFilePath.string());
        parser->AddCommandLineArguments(newargs);
    }else{
        outputFilePath = outputFileName;
    }
    logger->Info("Input Filename : " + inputFilePath.string() + "\n");
    logger->Info("Output Filename : "  + outputFilePath.string() + "\n");
    if (!fs::exists(inputFileName)) {
      logger->Critical("File " + inputFileName + " not found\n");
      return EXIT_FAILURE;
    }

    itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(inputFileName.c_str(), itk::ImageIOFactory::FileModeType::ReadMode);
    imageIO->SetFileName(inputFileName);
    imageIO->ReadImageInformation();
    itk::ImageIOBase::IOPixelType pixelType = imageIO->GetPixelType();
    itk::ImageIOBase::IOComponentType componentType = imageIO->GetComponentType();
    unsigned int dimensions = imageIO->GetNumberOfDimensions();

    logger->Info("Component Type  : " +
                 imageIO->GetComponentTypeAsString(componentType) + "\n");
    logger->Info("Pixel Type      : " +
                 itk::ImageIOBase::GetPixelTypeAsString(pixelType) + "\n");
    logger->Info("Image Dimension : " + std::to_string(dimensions) + "\n");

    switch (dimensions){
    case 2:
        switch (componentType){
        case itk::ImageIOBase::UCHAR:
            return experiment_impl<unsigned char,2>(parser, logger);
        case itk::ImageIOBase::USHORT:
            return experiment_impl<unsigned short, 2>(parser, logger);
        case itk::ImageIOBase::FLOAT:
            return experiment_impl<float, 2>(parser, logger);
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          logger->Critical("Unknown ComponentType :: " +
                           std::string(typeid(componentType).name()) +
                           " in 2 dimensions\n");
        }
        break;
    case 3:
       switch (componentType) {
       case itk::ImageIOBase::UCHAR:
          return experiment_impl<unsigned char, 3>(parser, logger);
       case itk::ImageIOBase::USHORT:
         return experiment_impl<unsigned short, 3>(parser, logger);
       case itk::ImageIOBase::FLOAT:
          return experiment_impl<float, 3>(parser, logger);
       case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
       default:
         logger->Critical("Unknown ComponentType :: " +
                          std::string(typeid(componentType).name()) +
                          " in 3 dimensions\n");
       }
       break;
    default:
        logger->Critical("File not supported..\n");
    }
    return EXIT_FAILURE;
}
