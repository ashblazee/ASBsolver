@echo off
set /p userInput="1)GNUMake 2)VS2022"

if "%userInput%"=="1" (
    vendor\Windows\premake5.exe --file=premake5.lua gmake
) else if "%userInput%"=="2" (
    vendor\Windows\premake5.exe --file=premake5.lua vs2022
) else (
    echo Invalid input. Please enter either 1 or 2.
)