#ifndef TURBO_EXPERIMENT_H
#define TURBO_EXPERIMENT_H


#include "itkCommandLineArgumentParser.h"
#include "itkLogger.h"

int experiment(const itk::CommandLineArgumentParser::Pointer &parser,
               const itk::Logger::Pointer &logger);

#endif
