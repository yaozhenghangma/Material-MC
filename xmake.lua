set_project("Material_MC")
set_version("0.1.0")

add_rules("mode.release", "mode.debug")

--set_warnings("all", "error")
set_languages("c++17")

package("scn_local")
    add_deps("cmake")
    set_sourcedir(path.join(os.scriptdir(), "include/scnlib"))
    on_install(function (package)
        local configs = {"-DSCN_TESTS=OFF", "-DSCN_DOCS=OFF", "-DSCN_EXAMPLES=OFF", "-DSCN_BENCHMARKS=OFF", "-DSCN_PENDANTIC=OFF", "-DSCN_BUILD_FUZZING=OFF"}
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:debug() and "Debug" or "Release"))
        table.insert(configs, "-DBUILD_SHARED_LIBS=" .. (package:config("shared") and "ON" or "OFF"))
        import("package.tools.cmake").install(package, configs)
    end)

    on_test(function (package)
        assert(package:check_cxxsnippets({test = [[
            #include <scn/scn.h>
            #include <cstdio>

            static void test() {
                int i;
                scn::prompt("What's your favorite number? ", "{}", i);
                printf("Oh, cool, %d!", i);
            }
        ]]}, {configs = {languages = "c++17"}, includes = "scn/scn.h"}))
    end)
package_end()

package("fmt_local")
    set_sourcedir(path.join(os.scriptdir(), "include/fmt"))
    on_load(function (package)
        package:add("deps", "cmake")
        if package:config("shared") then
            package:add("defines", "FMT_EXPORT")
        end
    end)

    on_install(function (package)
        local configs = {"-DFMT_TEST=OFF", "-DFMT_DOC=OFF", "-DFMT_FUZZ=OFF"}
        table.insert(configs, "-DBUILD_SHARED_LIBS=" .. (package:config("shared") and "ON" or "OFF"))
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:debug() and "Debug" or "Release"))
        if package:config("pic") ~= false then
            table.insert(configs, "-DCMAKE_POSITION_INDEPENDENT_CODE=ON")
        end
        import("package.tools.cmake").install(package, configs)
    end)

    on_test(function (package)
        assert(package:check_cxxsnippets({test = [[
            #include <fmt/format.h>
            #include <string>
            #include <assert.h>
            static void test() {
                std::string s = fmt::format("{}", "hello");
                assert(s == "hello");
            }
        ]]}, {configs = {languages = "c++11"}, includes = "fmt/format.h"}))
    end)

package_end()

package("ctre_local")
    add_deps("cmake")
    set_sourcedir(path.join(os.scriptdir(), "include/ctre"))
    on_install(function (package)
        local configs = {"-DCTRE_BUILD_TESTS=OFF"}
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:debug() and "Debug" or "Release"))
        table.insert(configs, "-DBUILD_SHARED_LIBS=" .. (package:config("shared") and "ON" or "OFF"))
        import("package.tools.cmake").install(package, configs)
    end)
package_end()

add_requires("ctre_local")
add_requires("fmt_local")
add_requires("scn_local")
add_requires("brew::open-mpi/ompi-cxx", {alias = "mpi"})
--add_requires("brew::boost", {alias = "boost"})
--add_requires("brew::boost-mpi", {alias = "boost-mpi"})

target("MMC")
    set_kind("binary")
    add_files("src/MMC.cpp")
    -- MPI
    add_packages("mpi")
    -- Boost
    if is_plat("macos") then
        add_linkdirs("/usr/local/Cellar/boost/1.76.0/lib")
        add_includedirs("/usr/local/Cellar/boost/1.76.0/include")
        add_linkdirs("/usr/local/Cellar/boost-mpi/1.76.0/lib")
    elseif is_plat("linux") then
        add_linkdirs("/home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/lib")
        add_includedirs("/home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/include")
        add_linkdirs("/home/linuxbrew/.linuxbrew/Cellar/boost-mpi/1.76.0/lib")
    end
    add_links("boost_mpi-mt", "boost_serialization")
    --add_packages("boost")
    --add_packages("boost-mpi")
    -- Third part
    add_packages("scn_local")
    add_packages("fmt_local")
    add_packages("ctre_local")

    after_build(function (target)
        os.mv("$(buildir)/$(plat)/$(arch)/$(mode)/MMC", "$(curdir)/MMC")
    end)
