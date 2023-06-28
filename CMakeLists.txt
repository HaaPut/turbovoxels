cmake_minimum_required(VERSION 3.13)
project(turbovoxels VERSION 1.0
        DESCRIPTION "Turbovoxels Segmentation"
        LANGUAGES CXX)
message(STATUS "Found CMAKE-${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
message(STATUS "----------------Running CMake---------------------------------")
cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)
message(STATUS "Invoking CMAKE on : ${HOSTNAME} running: ${CMAKE_HOST_SYSTEM}")

set(ITK_DIR $ENV{HOME}/opt/itk/lib/cmake/ITK-5.3)
message(STATUS "ITK Directory ${ITK_DIR}")
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -O3  -Wall -march=native")
set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(include)

# Tests
#enable_testing()
#add_subdirectory(tests)

# Source

set(SOURCES
  main.cpp
  src/itkCommandLineArgumentParser.cxx
  src/smoothing.cpp
  src/chanvese.cpp
  src/slic.cpp
  src/thresholding.cpp
  src/regionprops.cpp
  src/experiment.cpp
  src/utils.cpp
  src/monitor.cpp
  src/levelSetEvolution.cpp
  src/boundaries.cpp
  src/turbovoxels.cpp)

add_executable(turbo ${SOURCES})


if((${CMAKE_COMPILER_IS_GNUCXX}) AND (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS_EQUAL 8))
        message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}: adding stdc++fs library to link")
        target_link_libraries(turbo ${ITK_LIBRARIES} stdc++fs)
else()
        target_link_libraries(turbo ${ITK_LIBRARIES})
endif()
