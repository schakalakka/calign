cmake_minimum_required(VERSION 3.3)
project(swalign)

include_directories("/usr/include/python3.5m/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g -shared")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/home/andreas/workspace/swalign/build")


set(SOURCE_FILES
    src/swalign.c
        src/swalign.h src/swalignmodule.c)

add_executable(swalign ${SOURCE_FILES})