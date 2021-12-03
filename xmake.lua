set_project("Material_MC")
set_version("0.1.0")

add_rules("mode.release", "mode.debug")

--set_warnings("all", "error")
set_languages("c++20")

package("mpi")
    on_fetch(function (package, opt)
        return package:find_package("cmake::MPI", opt)
    end)
package_end()

package("boost")
    on_fetch(function (package, opt)
        opt.presets = {Boost_USE_STATIC_LIB = true}
        opt.components = {}
        table.insert(opt.components, "serialization")
        table.insert(opt.components, "mpi")

        return package:find_package("cmake::Boost", opt)
    end)
package_end()

add_requires("cmake")
add_requires("scnlib")
add_requires("fmt")
add_requires("ctre")
add_requires("mpi")
add_requires("spdlog")
add_requires("toml++")

if is_plat("macosx") then 
    add_requires("brew::boost", {alias = "boost"})
    add_requires("brew::boost-mpi", {alias = "boost-mpi"})
else
    add_requires("boost")
end 

option("boost_mpi")
    set_showmenu(true)
    set_values("boost_mpi", "boost_mpi_mt")
    set_default("boost_mpi")
option_end()

target("MMC")
    set_kind("binary")
    add_files("src/*.cpp", "custom/custom.cpp")
    -- Boost
    if is_plat("macosx") then
        add_packages("mpi")
        add_packages("boost")
        add_packages("boost-mpi")
    else 
        add_packages("mpi")
        add_packages("boost")
        if is_config("boost_mpi", "boost_mpi") then
            add_links("boost_mpi")
        elseif is_config("boost_mpi", "boost_mpi_mt") then
            add_linkdirs("/home/linuxbrew/.linuxbrew/Cellar/boost-mpi/1.76.0/lib")
            add_linkdirs("/home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/lib")
            add_links("boost_mpi-mt")
        end
        add_links("boost_serialization")
    end

    -- Third part
    add_packages("scnlib")
    add_packages("fmt")
    add_packages("ctre")
    add_packages("spdlog")
    add_packages("toml++")
