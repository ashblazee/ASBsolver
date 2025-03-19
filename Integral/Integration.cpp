#include "Integration.hpp"
#include "progressbar.hpp"
#include <array>
#include <cmath>
#include <omp.h>

#define allowprogressbar 1

namespace ASBsolver{
float integral(std::function<float(float)> function, std::array<float,2>& ranges, float& stepsize){
    float result = 0;
    long long totalsteps = round((ranges[1] - ranges[0]) / stepsize);
    #ifndef allowprogressbar
    #pragma omp parallel for
    #endif // don't allowprogressbar
    for(long long i = 0; i <= totalsteps; i += 1){
        #ifdef allowprogressbar
        showProgressBar(i,totalsteps,ranges[0] + (i * stepsize));
        #endif
        result += function(ranges[0] + (i * stepsize)) * stepsize;
    }
    return result;
}

float integral(std::function<float(float,float)> function, std::array<std::array<float,2>,2>& ranges, std::array<float,2>& stepsize){
    float result = 0;
    long long totalstepsx = round((ranges[0][1] - ranges[0][0]) / stepsize[0]);
    long long totalstepsy = round((ranges[1][1] - ranges[1][0]) / stepsize[1]);
    #ifndef allowprogressbar
    #pragma omp parallel for collapse(2)
    #endif // don't allowprogressbar
    for(long long i = 0; i < totalstepsx; i += 1){
        #ifdef allowprogressbar
        showProgressBar(i,totalstepsx,ranges[0][0] + (i * stepsize[0]));
        #pragma omp parallel for
        #endif
        for(long long j = 0; j < totalstepsy; j+= 1){
            result += function(ranges[0][0] + (i * stepsize[0]),ranges[1][0] + (j * stepsize[1])) * stepsize[0] * stepsize[1];
        }
    }
    return result;
}

float integral(std::function<float(float,float,float)> function, std::array<std::array<float,2>,3>& ranges, std::array<float,3>& stepsize){
    float result = 0;
    long long totalstepsx = round((ranges[0][1] - ranges[0][0]) / stepsize[0]);
    long long totalstepsy = round((ranges[1][1] - ranges[1][0]) / stepsize[1]);
    long long totalstepsz = round((ranges[2][1] - ranges[2][0]) / stepsize[2]);
    #ifndef allowprogressbar
    #pragma omp parallel for collapse(3)
    #endif // don't allowprogressbar
    for(long long i = 0; i < totalstepsx; i += 1){
        #ifdef allowprogressbar
        showProgressBar(i,totalstepsx,ranges[0][0] + (i * stepsize[0]));
        #pragma omp parallel for collapse(2)
        #endif
        for(long long j = 0; j < totalstepsy; j+= 1){
            for(long long k = 0; k < totalstepsz; k += 1){
                result += function(ranges[0][0] + (i * stepsize[0]),ranges[1][0] + (j * stepsize[1]),ranges[2][0] + (k * stepsize[2])) * stepsize[0] * stepsize[1] * stepsize[2];
            }
        }
    }
    return result;
}
}