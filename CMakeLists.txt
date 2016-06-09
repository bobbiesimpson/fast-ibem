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

add_executable(driver main.cpp)
target_link_libraries(driver hpro)
