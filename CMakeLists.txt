cmake_minimum_required(VERSION 3.9)

project(MC CXX)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
add_subdirectory(include/scnlib EXCLUDE_FROM_ALL)
add_subdirectory(include/ctre EXCLUDE_FROM_ALL)
include_directories("include/cereal/include/cereal")

add_executable(MC.out MC.cpp)

target_link_libraries(MC.out fmt)
target_link_libraries(MC.out scn::scn)
target_link_libraries(MC.out ctre)
find_package(MPI)
target_link_libraries(MC.out MPI::MPI_CXX)