cmake_minimum_required(VERSION 3.0)
project(fast-ibem)

# C++11 required to build this library
set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COM PILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

find_package(Boost REQUIRED)
if( Boost_FOUND )
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions( "-DHAS_BOOST" )
endif()

set(hlibprodir "NOT-DEFINED" CACHE PATH "/opt")
include_directories(${hlibprodir}/include)
link_directories(${hlibprodir}/lib)

# add Eigen library include fileds
include_directories(${fast-ibem-SOURCE_DIR}/libs/Eigen})

# Add conforming nurbs submodule
ADD_SUBDIRECTORY(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include/quadrature)
INCLUDE_DIRECTORIES(${fast-ibem_SOURCE_DIR}/libs/conformingNurbs/include/2D)
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} confnurbs)

set(LIB_INCLUDE_DIR ${fast-ibem_SOURCE_DIR}/include)
include_directories(${LIB_INCLUDE_DIR})
set(LIBHEADERS  ${LIB_INCLUDE_DIR}/CollocationAssembly.h
                ${LIB_INCLUDE_DIR}/Kernel.h
                ${LIB_INCLUDE_DIR}/KernelIntegrator.h)

set(LIB_SOURCE_DIR ${fast-ibem_SOURCE_DIR}/src)
set(LIBSOURCES 	${LIB_SOURCE_DIR}/CollocationAssembly.cpp
                ${LIB_SOURCE_DIR}/KernelIntegrator.cpp
		${LIBHEADERS})
add_library(fastibem ${LIBSOURCES})
target_link_libraries(fastibem ${THIS_LINK_LIBS})
	
add_executable(driver main.cpp)
target_link_libraries(driver hpro)

INCLUDE(CTest)
ADD_SUBDIRECTORY(tests)
ENABLE_TESTING()
