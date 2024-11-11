#include "FDM.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#define c 0.01

using namespace std;
using namespace ASBsolver;

void heatequation(float timestep,float xstep,float ystep, long long jpos,long long kpos, vector<vector<vector<float>>>& gridlist){
     // since xstep and ystep are the same for square grids
    gridlist[gridlist.size()-1][jpos][kpos] += LAPLCIAN2(gridlist[gridlist.size()-2],jpos,kpos,ystep) * timestep * c; // c is coefficient describing rate of heat transfer
}

void boundarycondi1(vector<vector<float>>& result){
    result[0] = result[1];
    result[result.size()-1] = result[result.size()-2];
}

void boundarycondi2(vector<vector<float>>& result){
    for(int i = 0; i < result.size()-1; i++){
        result[i][0] = result[i][1];
        result[i][result.size()-1] = result[i][result.size()-1];
    }
}

int main(){
    FiniteDifferenceMethod<vector<vector<float>>> heatsolver;
    heatsolver.AddRange(0,5); // time range
    heatsolver.AddRange(0,1); // gridx range
    heatsolver.AddRange(0,1); // gridy range
    
    function<void(vector<vector<float>>&)> boundarycond1 = boundarycondi1;
    function<void(vector<vector<float>>&)> boundarycond2 = boundarycondi2;

    heatsolver.BoundaryCondition(boundarycond1);
    heatsolver.BoundaryCondition(boundarycond2);
    
    heatsolver.AddSkips(500);

    heatsolver.AddStepsize(0.00001);
    heatsolver.AddStepsize(0.01);
    heatsolver.AddStepsize(0.01);
    std::cout << "Added step sizes\n";

    vector<vector<float>> initial;
    
    for(int i = 0; i < 100; i++){
        vector<float> stuff;
        initial.emplace_back(stuff);
        for(int j = 0; j < 100; j++){
            if(pow(i-50,2) + pow(j-50,2) <= pow(20,2)){
                initial[i].push_back(1);
            }
            else{
                initial[i].push_back(0);
            }
        }
    }

    std::cout << "Created initial condition\n";


    heatsolver.InitialCondition(initial);

    vector<vector<vector<float>>> sol;

    std::cout << "Added initial condition\n";

    sol = heatsolver.PDE2Solve(heatequation,1,2,2);
}