#pragma once
#define LAPLCIAN1(x,i,h) (((x[i-1] + x[i+1] - (2 * x[i])) / h) / h)
#define LAPLCIAN2(x,i,j,h) (((x[i][j+1] + x[i][j-1] + x[i-1][j] + x[i-1][j-1] + x[i-1][j]\
        + x[i+1][j] + x[i+1][j+1] + x[i+1][j-1] - (8 * x[i][j])) / h) / h)
#include <functional>
#include <list>
#include <vector>
#include <array>
#define _IGNORE_WARNINGS (warnings == 0)
#define _WRITETOFILE (filewrite == 1)

extern bool warnings;
extern bool filewrite;
extern const char* filename;

namespace ASBsolver{
template <typename T>
class FiniteDifferenceMethod {
    public:
        int skips = 1;
        std::vector<float> stepsize;
        std::vector<std::array<float,2>> Rangelist;
        
        std::list<T> initialconditions;

        std::list<std::function<void(T&)>> boundaries;

    void InitialCondition(T initialcondition);

    void BoundaryCondition(std::function<void(T&)> boundarycondition);

    void AddRange(float start, float stop);

    void AddStepsize(float step);

    void AddSkips(int skip); // record the data every n times to save storage and memory

    std::vector<float> ODE1Solve(std::function<float(float,float)> ode); // takes in current value (ODE1)

    std::vector<float> ODE2Solve(std::function<float(float,float,float)> ode); // takes in past and current value (ODE2)

    std::vector<float> ODE3Solve(std::function<float(float,float,float,float)> ode);

    std::vector<float> ODE4Solve(std::function<float(float,float,float,float,float)> ode);

    std::vector<float> ODE5Solve(std::function<float(float,float,float,float,float,float)> ode);
    
                                                                            //tstep xstep currentlocation grid
    std::vector<std::vector<float>> PDE1Solve(std::function<void(float,float,long long,std::vector<std::vector<float>>&)> pde, int dim1, int dim2); 
    
    //* dim1 and dim2 refer to the order of the variables in the equation
    //*PDE1Solve the 1 refers to the fact that besides time it only has 1 other dimension

    std::vector<std::vector<std::vector<float>>> PDE2Solve(std::function<void(float,float,float,long long, long long,std::vector<std::vector<std::vector<float>>>&)> pde, int dim1, int dim2, int dim3);

    std::vector<std::vector<std::vector<std::vector<float>>>> PDE3Solve(std::function<void(float,float,float,float,long long, long long,long long,std::vector<std::vector<std::vector<std::vector<float>>>>&)> pde,
                                                                    int dim1, int dim2, int dim3,int dim4);
};
}