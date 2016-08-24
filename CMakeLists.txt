cmake_minimum_required(VERSION 3.0)
project(fast-ibem)

# Add FindCMake files to module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# C++11 required to build this library
set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

# Add the boost library
find_package(Boost REQUIRED COMPONENTS filesystem iostreams system program_options)
if( Boost_FOUND )
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  add_definitions( "-DHAS_BOOST" )
  set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${Boost_LIBRARIES})
endif()

# Add the Hlib library required for fast dense matrix computations
set(hlibprodir "NOT-DEFINED" CACHE PATH "/opt")
if(hlibprodir)
    include_directories(${hlibprodir}/include)
    link_directories(${hlibprodir}/lib)
else()
    message("HlibPro directory not defined")
    return()
endif()

find_package(blas)
if(BLAS_FOUND)
    set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${BLAS_LIBRARIES})
endif()

# add Eigen library include fileds
include_directories(SYSTEM ${fast-ibem-SOURCE_DIR}/libs})

# Add conforming nurbs submodule
ADD_SUBDIRECTORY(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include/quadrature)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include/2D)
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} confnurbs)

# Set the include directories of this library
set(LIB_INCLUDE_DIR ${fast-ibem_SOURCE_DIR}/include)
include_directories(${LIB_INCLUDE_DIR})
set(LIBHEADERS  ${LIB_INCLUDE_DIR}/CollocationAssembly.h
                ${LIB_INCLUDE_DIR}/Functor.h
                ${LIB_INCLUDE_DIR}/GalerkinAssembly.h
                ${LIB_INCLUDE_DIR}/HAssembly.h
                ${LIB_INCLUDE_DIR}/Kernel.h
                ${LIB_INCLUDE_DIR}/MaterialParam.h)


# Add the .cpp source files to the library
set(LIB_SOURCE_DIR ${fast-ibem_SOURCE_DIR}/src)
set(LIBSOURCES 	${LIB_SOURCE_DIR}/HAssembly.cpp
                ${LIBHEADERS})
add_library(fastibem ${LIBSOURCES})
target_link_libraries(fastibem ${THIS_LINK_LIBS})

INCLUDE(CTest)
ADD_SUBDIRECTORY(tests)
ENABLE_TESTING()
