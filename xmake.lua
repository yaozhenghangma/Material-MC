set_project("MC")
set_version("0.0.1")

add_rules("mode.release", "mode.debug")

set_warnings("all", "error")
set_languages("c++17")

target("MC")
    set_kind("binary")
    add_files("MC.cpp")
    add_includedirs("/usr/local/Cellar/open-mpi/4.1.1_2/include")
    add_linkdirs("/usr/local/Cellar/open-mpi/4.1.1_2/lib")
    add_links("mpi")
    add_includedirs("/usr/local/Cellar/boost/1.76.0/include")
    add_linkdirs("/usr/local/Cellar/boost/1.76.0/lib")
    add_linkdirs("/usr/local/Cellar/boost-mpi/1.76.0/lib")
    add_links("boost_mpi", "boost_serialization")