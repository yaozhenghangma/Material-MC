# project
cmake_minimum_required(VERSION 3.13.0)
project(Material_MC VERSION 0.3.4 LANGUAGES C CXX ASM)

# target
add_executable(MMC src/MMC.cpp)
set_target_properties(MMC PROPERTIES OUTPUT_NAME "MMC")
set_target_properties(MMC PROPERTIES RUNTIME_OUTPUT_DIRECTORY "build/")
target_include_directories(MMC PRIVATE
    /usr/local/Cellar/open-mpi/4.1.1_2/include
    /usr/local/Cellar/boost/1.76.0/include
    /usr/local/Cellar/boost-mpi/1.76.0/include
)

# scnlib
add_subdirectory(submodules/scnlib)
target_link_libraries(MMC PRIVATE scn::scn)

# fmt
add_subdirectory(submodules/fmt)
target_link_libraries(MMC PRIVATE fmt::fmt)

# ctre
add_subdirectory(submodules/ctre)
target_link_libraries(MMC PRIVATE ctre)

# spdlog
add_subdirectory(submodules/spdlog)
target_link_libraries(MMC PRIVATE spdlog::spdlog)

# toml++
add_subdirectory(submodules/toml++)
target_include_directories(MMC PRIVATE submodules/toml++/include)

target_link_libraries(MMC PRIVATE
    mpi
    boost_serialization
    boost_mpi
    z
)
target_link_directories(MMC PRIVATE
    /usr/local/Cellar/open-mpi/4.1.1_2/lib
    /usr/local/Cellar/boost/1.76.0/lib
    /usr/local/Cellar/boost-mpi/1.76.0/lib
)
target_link_options(MMC PRIVATE
    -stdlib=libc++
)
target_sources(MMC PRIVATE
    src/log.cpp
    src/structure_in.cpp
    src/configure_in.cpp
    src/result_out.cpp
    src/spin_out.cpp
    src/MC_structure.cpp
    src/Hamiltonian.cpp
    src/rotation.cpp
    custom/custom.cpp
)

