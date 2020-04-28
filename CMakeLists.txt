cmake_minimum_required(VERSION 3.5)

project(GNI_predictors LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O9 -g -finline-functions -m64)
# Enable sanitizers - To take account run this: "cmake -DCMAKE_ASAN=ON .."
if(CMAKE_ASAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
endif()
include(ExternalProject)
ExternalProject_Add(
  NetworKit
  GIT_REPOSITORY      https://github.com/networkit/networkit.git
  GIT_TAG             origin/release-6.1
  UPDATE_DISCONNECTED ON
  INSTALL_COMMAND     ""
  TMP_DIR             /tmp
  DOWNLOAD_DIR        lib/src
  SOURCE_DIR          lib/src/NetworKit
  BINARY_DIR          lib/src/build
)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(GNI_predictors src/gnip.cpp
                              src/matrixPreparation.cpp
                              src/randomWalk.cpp)
target_include_directories(GNI_predictors PRIVATE src)
target_link_libraries(GNI_predictors PUBLIC "libnetworkit.so")
target_link_libraries(GNI_predictors PUBLIC ${ARMADILLO_LIBRARIES})
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()

target_link_libraries(GNI_predictors PUBLIC OpenMP::OpenMP_CXX)
set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.69.0 COMPONENTS date_time filesystem system program_options)

if(${Boost_FOUND})
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(GNI_predictors PUBLIC ${Boost_LIBRARIES})
endif()
