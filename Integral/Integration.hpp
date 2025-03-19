#pragma once
#include <functional>

namespace ASBsolver{

float integral(std::function<float(float)> function, std::array<float,2>& ranges, float& stepsize);

float integral(std::function<float(float,float)> function, std::array<std::array<float,2>,2>& ranges, std::array<float,2>& stepsize);

float integral(std::function<float(float,float,float)> function, std::array<std::array<float,2>,3>& ranges, std::array<float,3>& stepsize);

}