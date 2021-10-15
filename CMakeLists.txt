cmake_minimum_required(VERSION 3.9)

project(MC CXX)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
add_subdirectory(include/scnlib EXCLUDE_FROM_ALL)
add_subdirectory(include/ctre EXCLUDE_FROM_ALL)

link_directories("/usr/local/Cellar/open-mpi/4.1.1_2/lib")
link_directories("/usr/local/Cellar/boost/1.76.0/lib")
link_directories("/usr/local/Cellar/boost-mpi/1.76.0/lib")

include_directories("/usr/local/Cellar/open-mpi/4.1.1_2/include")
include_directories("/usr/local/Cellar/boost/1.76.0/include")

add_executable(MC MC.cpp)

target_link_libraries(MC fmt)
target_link_libraries(MC scn::scn)
target_link_libraries(MC ctre)
find_package(MPI)
target_link_libraries(MC MPI::MPI_CXX)
target_link_libraries(MC boost_mpi)
target_link_libraries(MC boost_serialization)