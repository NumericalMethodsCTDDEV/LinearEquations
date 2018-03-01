cmake_minimum_required(VERSION 2.8)

project(LinearEquations)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -Wall -Wextra -Wshadow -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 -O2 -Wall -pthread")

add_library(googleTests STATIC
    gtest/gtest-all.cc
    gtest/gtest.h
    gtest/gtest_main.cc)

add_library(gaussMethod STATIC
    gauss.cpp)

enable_testing()

add_executable(gauss_testing testing.cpp)
add_executable(gauss gauss_main.cpp)

target_link_libraries(gauss gaussMethod)
target_link_libraries(gauss_testing googleTests gaussMethod)


