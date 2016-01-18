cmake_minimum_required(VERSION 3.3)
project(swalign)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/swalign.c
    src/swalign.h
    Makefile)

add_executable(swalign ${SOURCE_FILES})