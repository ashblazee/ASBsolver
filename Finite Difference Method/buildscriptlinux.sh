#!/bin/bash

echo -n "1)GNUMake 2)VS2022"
read userInput

if [ "$userInput" == "1" ]; then
    vendor/Linux/premake5 --file=premake5.lua gmake
elif [ "$userInput" == "2" ]; then
    vendor/Linux/premake5 --file=premake5.lua vs2022
else
    echo "Invalid input. Please enter either 1 or 2."
fi