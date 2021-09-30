cmake_minimum_required(VERSION 2.9)

project(MC)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(include/fmt EXCLUDE_FROM_ALL)

add_executable(MC MC.cpp)

target_link_libraries(MC fmt)