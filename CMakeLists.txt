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

find_package(Armadillo REQUIRED)
find_package(Boost COMPONENTS date_time filesystem system program_options)
if(UNIX)
    set(LOCAL_BOOST_DIR "${CMAKE_BINARY_DIR}/lib/src/boost")
    set(BOOST_SHA256_HASH 9995e192e68528793755692917f9eb6422f3052a53c5e13ba278a228af6c7acf)
    set(BOOST_DOWNLOAD_LINK "https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.gz")
else()
    set(LOCAL_BOOST_DIR "${CMAKE_BINARY_DIR}/lib/src/boost")
    set(BOOST_SHA256_HASH 0909a79524f857ef54570ceef8f397cc0629202532cc997785479c7c08bbc2a4)
    set(BOOST_DOWNLOAD_LINK "https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.zip")
endif()

include_directories(${ARMADILLO_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)


if(${Boost_FOUND})
    include_directories(${Boost_INCLUDE_DIRS})

else()
   # set(Boost_LIBRARIES "${CMAKE_BINARY_DIR}/lib/src/networkit-build/lib/libboost_date_time.so;${CMAKE_BINARY_DIR}/lib/src/networkit-build/lib/libboost_filesystem.so;${CMAKE_BINARY_DIR}/lib/src/networkit-build/lib/libboost_system.so;${CMAKE_BINARY_DIR}/lib/src/networkit-build/lib/program_options.so")
    set(Boost_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/lib/src/networkit/include")
    ExternalProject_Add(
      Boost
      URL         ${LOCAL_DOWNLOAD_LINK}
      URL_HASH    SHA256=${BOOST_SHA256_HASH}
      UPDATE_DISCONNECTED ON
      INSTALL_COMMAND     ""
      TMP_DIR             /tmp
      DOWNLOAD_DIR        ${CMAKE_BINARY_DIR}/lib/src
      SOURCE_DIR          ${CMAKE_BINARY_DIR}/lib/src/boost
      BINARY_DIR          ${CMAKE_BINARY_DIR}/lib/src/boost-build
    )
    find_package(Boost COMPONENTS date_time filesystem system program_options)
    message(STATUS "${Boost_LIBRARIES}")
endif()

set(NETWORKIT_SRC_DIR "${CMAKE_BINARY_DIR}/lib/src/networkit")
set(NETWORKIT_BIN_DIR "${CMAKE_BINARY_DIR}/lib/src/networkit-build")
set(NETWORKIT_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/lib/src/networkit/include")
set(NETWORKIT_TLX_EXT "${CMAKE_BINARY_DIR}/lib/src/networkit/extlibs/tlx")
set(NETWORKIT_TTMATH_EXT "${CMAKE_BINARY_DIR}/lib/src/networkit/extlibs/ttmath")


ExternalProject_Add(
  NetworKit
  GIT_REPOSITORY      https://github.com/networkit/networkit.git
  GIT_TAG             origin/release-6.1
  UPDATE_DISCONNECTED ON
  INSTALL_COMMAND     ""
  TMP_DIR             /tmp
  DOWNLOAD_DIR        ${CMAKE_BINARY_DIR}/lib/src
  SOURCE_DIR          ${NETWORKIT_SRC_DIR}
  BINARY_DIR          ${NETWORKIT_BIN_DIR}
)

include_directories(${NETWORKIT_INCLUDE_DIRS})
include_directories(${NETWORKIT_TLX_EXT})
include_directories(${NETWORKIT_TTMATH_EXT})

add_executable(GNI_predictors src/gnip.cpp
                              src/matrixPreparation.cpp
                              src/randomWalk.cpp)
target_include_directories(GNI_predictors PRIVATE src)

target_link_libraries(GNI_predictors PUBLIC "${NETWORKIT_BIN_DIR}/libnetworkit.so")
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

if(${Boost_FOUND})
    target_link_libraries(GNI_predictors PUBLIC ${Boost_LIBRARIES})

endif()

