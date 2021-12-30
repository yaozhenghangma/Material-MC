# project
cmake_minimum_required(VERSION 3.13.0)
project(Material_MC VERSION 0.3.4)
set(CMAKE_CXX_STANDARD 17)

# target
add_executable(MMC src/MMC.cpp)
set_target_properties(MMC PROPERTIES OUTPUT_NAME "MMC")
set_target_properties(MMC PROPERTIES RUNTIME_OUTPUT_DIRECTORY "build/")

# Boost
find_package(Boost COMPONENTS mpi serialization)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    
    MESSAGE( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    MESSAGE( STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    MESSAGE( STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")

    target_link_libraries (MMC PRIVATE ${Boost_LIBRARIES})
endif()

# MPI
find_package(MPI)
if (MPI_FOUND)
    include_directories(MMC ${MPI_CXX_INCLUDE_PATH})

    MESSAGE( STATUS "MPI_INCLUDE_DIRS = ${MPI_CXX_INCLUDE_DIRS}.")
    MESSAGE( STATUS "MPI_CXX_VERSION = ${MPI_CXX_VERSION}")

    target_link_libraries(MMC PRIVATE MPI::MPI_CXX)
endif ()

# scnlib
add_subdirectory(submodules/scnlib)
target_link_libraries(MMC PRIVATE scn::scn)

# fmt
add_subdirectory(submodules/fmt)
target_link_libraries(MMC PRIVATE fmt::fmt)

# ctre
add_subdirectory(submodules/ctre)
target_link_libraries(MMC PRIVATE ctre)

# toml++
add_subdirectory(submodules/toml++)
target_include_directories(MMC PRIVATE submodules/toml++/include)

target_sources(MMC PRIVATE
    src/log.cpp
    src/structure_in.cpp
    src/configure_in.cpp
    src/result_out.cpp
    src/spin_out.cpp
    src/MC_structure.cpp
    src/Hamiltonian.cpp
    src/rotation.cpp
    src/initialization.cpp
    src/random_function.cpp
    src/methods/classical.cpp
    src/methods/parallel_tempering.cpp
    src/methods/local_update.cpp
    custom/custom.cpp
)

