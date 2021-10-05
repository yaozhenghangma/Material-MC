cmake_minimum_required(VERSION 2.9)

project(MC CXX)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
add_subdirectory(include/scnlib EXCLUDE_FROM_ALL)
add_subdirectory(include/ctre EXCLUDE_FROM_ALL)

add_executable(MC MC.cpp)

target_link_libraries(MC fmt)
target_link_libraries(MC scn::scn)
target_link_libraries(MC ctre)