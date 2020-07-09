cmake_minimum_required(VERSION 3.9)

project(GNI_predictors LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O0")
if(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES
                                         "RelWithDebInfo")
  set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O3")
endif()

if(CMAKE_OPTIMIZATION_LEVEL)
  set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O${CMAKE_OPTIMIZATION_LEVEL}")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE}")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g")

set(CMAKE_CXX_FLAGS_DEBUG
    "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g -Wall -Wextra -Wpedantic")

# Enable sanitizers - To take account run this: "cmake -DCMAKE_ASAN=ON .."
if(CMAKE_ASAN)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined"
  )
endif()

include(ExternalProject)

# ##############################################################################
# External Libs
# ##############################################################################

# Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Boost
set(Boost_USE_STATIC_LIBS OFF) # only find static libs
set(Boost_USE_RELEASE_LIBS ON) # only find release libs
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS date_time filesystem system
                                       program_options)

# NetworKit we build ourselves
set(NETWORKIT_SRC_DIR "${CMAKE_BINARY_DIR}/lib/src/networkit")
set(NETWORKIT_BIN_DIR "${CMAKE_BINARY_DIR}/lib/src/networkit-build")
set(NETWORKIT_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/lib/src/networkit/include")
set(NETWORKIT_TLX_EXT "${CMAKE_BINARY_DIR}/lib/src/networkit/extlibs/tlx")
set(NETWORKIT_TTMATH_EXT "${CMAKE_BINARY_DIR}/lib/src/networkit/extlibs/ttmath")

ExternalProject_Add(
  NetworKit
  CMAKE_ARGS "-DNETWORKIT_STATIC=ON"
  GIT_REPOSITORY https://github.com/networkit/networkit.git
  GIT_TAG origin/release-6.1
  UPDATE_DISCONNECTED ON
  INSTALL_COMMAND ""
  TMP_DIR /tmp
  DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/lib/src
  SOURCE_DIR ${NETWORKIT_SRC_DIR}
  BINARY_DIR ${NETWORKIT_BIN_DIR})

# OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)
if(OPENMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Setup target
add_executable(GNI_predictors src/gnip.cpp src/matrixPreparation.cpp
                              src/randomWalk.cpp)

# Need to manually register NetworKit as a dep
add_dependencies(GNI_predictors NetworKit)


# Setup some more dirs the compiler needs to know about
target_include_directories(GNI_predictors PRIVATE src)
target_link_directories(GNI_predictors PUBLIC ${NETWORKIT_BIN_DIR})
include_directories(${NETWORKIT_INCLUDE_DIRS})
include_directories(${NETWORKIT_TLX_EXT})
include_directories(${NETWORKIT_TTMATH_EXT})

# Link stuff
target_link_libraries(GNI_predictors PUBLIC networkit)
target_link_libraries(GNI_predictors PUBLIC ${ARMADILLO_LIBRARIES})
target_link_libraries(GNI_predictors PUBLIC ${Boost_LIBRARIES})

# Install config
install(TARGETS GNI_predictors DESTINATION "bin")
