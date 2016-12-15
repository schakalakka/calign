cmake_minimum_required(VERSION 3.3)
project(calign)

#include_directories("/usr/include/python3.5m/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g -fPIC")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/home/andreas/GDrive/workspace/calign/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "/home/andreas/GDrive/workspace/calign/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "/home/andreas/GDrive/workspace/calign/lib")

set(SOURCE_FILES src/calign.c src/calign.h)

add_library(calignshared SHARED ${SOURCE_FILES})
#add_library(staticcalign STATIC ${SOURCE_FILES})
add_executable(calign ${SOURCE_FILES})

