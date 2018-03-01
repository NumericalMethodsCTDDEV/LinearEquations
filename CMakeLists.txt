cmake_minimum_required(VERSION 2.8)

project(LinearEquations)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-pthread")

add_library(googleTests STATIC
    gtest/gtest-all.cc
    gtest/gtest.h
    gtest/gtest_main.cc)

add_library(gaussMethod STATIC
    gauss.cpp)

#add_library(Gauss Method)

enable_testing()

add_executable(gauss_testing testing.cpp)
add_executable(gauss gauss_main.cpp)

target_link_libraries(gauss gaussMethod)
target_link_libraries(gauss_testing googleTests gaussMethod)


