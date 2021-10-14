set_project("MC")
set_version("0.0.1")

add_rules("mode.release", "mode.debug")

set_warnings("all", "error")
set_languages("c++17")

target("MC")
    set_kind("binary")
    add_files("MC.cpp")