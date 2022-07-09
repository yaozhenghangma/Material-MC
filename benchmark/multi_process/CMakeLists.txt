# project
cmake_minimum_required(VERSION 3.13.0)
project(2D_Heisenberg)
set(CMAKE_CXX_STANDARD 17)

find_package(pybind11 CONFIG)

# target
pybind11_add_module(heisenberg heisenberg.cpp)