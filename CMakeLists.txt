cmake_minimum_required(VERSION 3.3)
project(swalign)

#include_directories("/usr/include/python3.5m/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g -fPIC")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/home/andreas/workspace/swalign/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "/home/andreas/workspace/swalign/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "/home/andreas/workspace/swalign/lib")

set(SOURCE_FILES src/swalign.c src/swalign.h)

add_library(sharedswalign SHARED ${SOURCE_FILES})
add_library(staticswalign STATIC ${SOURCE_FILES})
add_executable(swalign ${SOURCE_FILES})

