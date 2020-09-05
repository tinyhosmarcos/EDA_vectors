project("course_vec")

cmake_minimum_required(VERSION 3.1.5)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSION OFF)

include_directories($(CMAKE_CURRENT_SOURCE_DIR)/source)

add_executable(course_vec.exe source/vectores.cpp)


