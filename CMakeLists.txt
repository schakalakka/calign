cmake_minimum_required(VERSION 3.3)
project(swalign)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g")

set(SOURCE_FILES
    src/swalign.c
        src/swalign.h)

add_executable(swalign ${SOURCE_FILES})