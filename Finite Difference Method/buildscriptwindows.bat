@echo off
set /p userInput="1)GNUMake 2)VS2022"

if "%userInput%"=="1" (
    premake5 gmake
) else if "%userInput%"=="2" (
    premake5 vs2022
) else (
    echo Invalid input. Please enter either 1 or 2.
)