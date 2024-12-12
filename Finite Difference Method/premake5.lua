workspace "ASBsolver"
    configurations { "Debug", "Release" }
    platforms {"x64"}

    location "build"

project "FDM example"
    kind "ConsoleApp"
    language "C++"
    cppdialect "C++20"
    files {"source/FDM.cpp","source/example.cpp"}

    filter "configurations:Debug"
        symbols "On"

    filter "configurations:Release"
        optimize "On"
