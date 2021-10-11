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
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()
target_link_libraries(MC OpenMP::OpenMP_CXX)