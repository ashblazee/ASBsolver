#!/bin/bash

echo -n "1)GNUMake 2)VS2022"
read userInput

if [ "$userInput" == "1" ]; then
    premake5 gmake
elif [ "$userInput" == "2" ]; then
    premake5 vs2022
else
    echo "Invalid input. Please enter either 1 or 2."
fi