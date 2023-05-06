# project
cmake_minimum_required(VERSION 3.13.0)
project(MMC VERSION "0.0.1")
set(CMAKE_CXX_STANDARD 17)

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific directory
  # here.
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

find_package(pybind11 CONFIG REQUIRED)

# target
pybind11_add_module(c_cell src/include/base/cell.cpp)

install(TARGETS c_cell DESTINATION .)
