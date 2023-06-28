#include "utils.h"

template itk::Image<float, 2>::Pointer compute_edge_speed<itk::Image<float, 2>>(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger, itk::Image<float, 2>::Pointer &input);


template itk::Image<float, 3>::Pointer compute_edge_speed<itk::Image<float, 3>>(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger, itk::Image<float, 3>::Pointer &input);

template SparseLevelSetType2::Pointer
initializeLevelSet<SparseLevelSetType2, itk::Image<float,2>>(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger, itk::Image<float, 2>::Pointer &input);

template SparseLevelSetType3::Pointer
initializeLevelSet<SparseLevelSetType3,
                   itk::Image<float, 3>>(
    const itk::CommandLineArgumentParser::Pointer &parser,
    const itk::Logger::Pointer &logger, itk::Image<float, 3>::Pointer &input);


//void drawContours(const itk::CommandLineArgumentParser::Pointer &parser,
//                  const itk::Logger::Pointer &logger){
//    constexpr unsigned Dimension = 3;
//    logger->Info("Running RGB turbovoxels\n");
//    using RealPixelType = float;
//    using RealImageType = itk::Image<RealPixelType, Dimension>;
//
//    using ReaderType = itk::ImageFileReader<RealImageType>;
//
//    typename ReaderType::Pointer reader = ReaderType::New();
//    std::string inputFileName = "../viz/sum.tif";
//    parser->GetCommandLineArgument("-input", inputFileName);
//    reader->SetFileName(inputFileName);
//    logger->Info("Input image file " + reader->GetFileName() + "\n");
//    reader->Update();
//    typename RealImageType::Pointer finalDenseLevelSet = reader->GetOutput();
//    using RGBPixelType = itk::RGBPixel<unsigned char>;
//    using RGBImageType = itk::Image<RGBPixelType, Dimension>;
//
//    std::vector<RealPixelType> values = {2,9,15};
//
//    std::stringstream ss;
//    ss << "Set countour values = {";
//    for (auto value: values) ss << std::to_string(value) << ", ";
//    ss.seekp(-2, ss.cur);
//    ss << "}\n";
//    logger->Info(ss.str());
//    ss.clear();
//
//    typename RGBImageType::Pointer colored = paintCountours<RealImageType, RGBImageType>(finalDenseLevelSet,values);
//
//    fs::path outputpath = fs::path(inputFileName).parent_path()/"coloredContours.tif";
//    writeImage<RGBImageType>(colored, outputpath.string(), logger,
//                             "Wrote overlaid Turbo voxels to " +
//                             outputpath.string() + "\n");
//}