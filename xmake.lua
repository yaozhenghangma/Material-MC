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

if is_plat("macosx") then 
    add_requires("brew::boost", {alias = "boost"})
    add_requires("brew::boost-mpi", {alias = "boost-mpi"})
else
    add_requires("boost")
    option("boost_mpi")
        set_default(true)
        add_links("boost_mpi")
    option("boost_mpi-mt")
        set_default(true)
        add_links("boost_mpi-mt")
end 

target("MMC")
    set_kind("binary")
    add_files("src/MMC.cpp")
    -- Boost
    if is_plat("macosx") then
        add_packages("mpi")
        add_packages("boost")
        add_packages("boost-mpi")
    else 
        add_packages("mpi")
        add_packages("boost")
        add_links("boost_serialization")
        add_options("boost_mpi")
        add_options("boost_mpi-mt")
    end

    -- Third part
    add_packages("scnlib")
    add_packages("fmt")
    add_packages("ctre")

    after_build(function (target)
        os.cp("$(buildir)/$(plat)/$(arch)/$(mode)/MMC", "$(curdir)/MMC")
    end)
