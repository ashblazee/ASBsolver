#include "FDM.hpp"
#include "progressbar.hpp"
#include <functional>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <omp.h>

bool warnings = 1;
bool filewrite = 1;
const char* filename = "result.csv";

using namespace ASBsolver;


void instabilitycheck(float initial, float changed){
    if (_IGNORE_WARNINGS){
        if(changed >= initial*2 || changed <= initial / 2){
        std::cout << "WARNING: Stiff system detected. There is a possibility of a singularity forming.";
        }
    } // if warnings are not ignored
}

template<typename T>
void FiniteDifferenceMethod<T>::AddRange(float start, float stop){
    std::array<float,2> range = {start, stop}; 
    this->Rangelist.emplace_back(range);
}

template<typename T>
void FiniteDifferenceMethod<T>::InitialCondition(T initialcondition){
    this->initialconditions.emplace_back(initialcondition);
}

template<typename T>
void FiniteDifferenceMethod<T>::AddStepsize(float step){
    this->stepsize.push_back(step); //* follows t, x, y, z... format
}

template<typename T>
void FiniteDifferenceMethod<T>::AddSkips(int skip){
    this->skips = skip;
}

//ODE1 function

template<>
std::vector<float> FiniteDifferenceMethod<float>::ODE1Solve(std::function<float(float,float)> ode){
    try{
        if (this->initialconditions.empty()){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    float current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> range = this->Rangelist[0];
    float step = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((range[1] - range[0]) / step);

    std::vector<float> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        float next = 0;
        next = ode(current,step);
        instabilitycheck(current, next);
        current = next;
        if(i % this->skips){
            result.push_back(current);
        }
    }
    return result;
}

//ODE2 function
template<>
std::vector<float> FiniteDifferenceMethod<float>::ODE2Solve(std::function<float(float,float,float)> ode){
    try{
        if (this->initialconditions.size() < 2){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    float past = this->initialconditions.front();
    this->initialconditions.pop_front();
    float current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> range = this->Rangelist[0];
    float step = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((range[1] - range[0]) / step);

    std::vector<float> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        float next = 0;
        next = ode(past,current,step);
        instabilitycheck(current, next);
        past = current;
        current = next;
        if(i % this->skips){
            result.push_back(current);
        }
    }
    return result;
}

//ODE3 function
template<>
std::vector<float> FiniteDifferenceMethod<float>::ODE3Solve(std::function<float(float,float,float,float)> ode){
    try{
        if (this->initialconditions.size() < 3){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    float past1 = this->initialconditions.front(); // past 2 is more recent than past 1
    this->initialconditions.pop_front();
    float past2 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> range = this->Rangelist[0];
    float step = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((range[1] - range[0]) / step);

    std::vector<float> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        float next = 0;
        next = ode(past1,past2,current,step);
        instabilitycheck(current, next);
        past1 = past2;
        past2 = current;
        current = next;
        if(i % this->skips){
            result.push_back(current);
        }
    }
    return result;
}

//ODE4 function
template<>
std::vector<float> FiniteDifferenceMethod<float>::ODE4Solve(std::function<float(float,float,float,float,float)> ode){
    try{
        if (this->initialconditions.size() < 4){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    float past1 = this->initialconditions.front(); // past 2 is more recent than past 1
    this->initialconditions.pop_front();
    float past2 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float past3 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> range = this->Rangelist[0];
    float step = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((range[1] - range[0]) / step);

    std::vector<float> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        float next = 0;
        next = ode(past1,past2,past3,current,step);
        instabilitycheck(current, next);
        past1 = past2;
        past2 = past3;
        past3 = current;
        current = next;
        if(i % this->skips){
            result.push_back(current);
        }
    }
    return result;
}

//ODE5 function
template<>
std::vector<float> FiniteDifferenceMethod<float>::ODE5Solve(std::function<float(float,float,float,float,float,float)> ode){
    try{
        if (this->initialconditions.size() < 5){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    float past1 = this->initialconditions.front(); // past 2 is more recent than past 1
    this->initialconditions.pop_front();
    float past2 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float past3 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float past4 = this->initialconditions.front();
    this->initialconditions.pop_front();
    float current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> range = this->Rangelist[0];
    float step = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((range[1] - range[0]) / step);

    std::vector<float> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        float next = 0;
        next = ode(past1,past2,past3,past4,current,step);
        instabilitycheck(current, next);
        past1 = past2;
        past2 = past3;
        past3 = past4;
        past4 = current;
        current = next;
        if(i % this->skips){
            result.push_back(current);
        }
    }
    return result;
}

//PDE1 functions

template<typename T>
void FiniteDifferenceMethod<T>::BoundaryCondition(std::function<void(T&)> boundarycondition){
    this->boundaries.emplace(this->boundaries.end(),boundarycondition);
}

template<>                                                                                                                                      //orders of the dimension
std::vector<std::vector<float>> FiniteDifferenceMethod<std::vector<float>>::PDE1Solve(std::function<void(float,float,long long,std::vector<std::vector<float>>&)> pde, int dim1, int dim2){
    std::ofstream* fileresult = nullptr;
    if(_WRITETOFILE){
        fileresult = new std::ofstream(filename,std::ios::app);
    }
    try{
        if (this->initialconditions.empty()){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    std::vector<float> current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> rangea = this->Rangelist[0];
    std::array<float,2> rangeb = this->Rangelist[1];
    float stepa = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepb = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((rangea[1] - rangea[0]) / stepa);
    long long grids = round((rangeb[1] - rangeb[0]) / stepb);

    std::vector<std::vector<float>> result;
    result.push_back(current);

    std::vector<std::vector<float>> gridlist;
    for(int i = 0; i < dim1-1; i++){
        std::vector<float> temp = current;
        gridlist.push_back(temp);
    }
    gridlist.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        showProgressBar(i,num_iterations,i * stepa);
        std::vector<float> next = gridlist[gridlist.size()-1];
        //apply boundary conditions
        for(const auto& boundarycondition : this->boundaries){
            boundarycondition(next);
        }
        gridlist.push_back(next);

        #pragma omp parallel for
        for(long long j = floor(dim2 / 2); j < grids - static_cast<long long>(ceil(dim2 / 2)); j++){
            pde(stepa,stepb,j,gridlist);
        }

        for(int j = 0; j < gridlist.size()-1; j++){
            gridlist[j] = gridlist[j+1];
        }
        gridlist.pop_back(); // removes "next"

        if(i % this->skips){
            if (_WRITETOFILE)
            {
                for(const auto& yes : gridlist[gridlist.size()-1]){
                    *fileresult << yes;
                    *fileresult << ",";
                }
                *fileresult << "\n";
            }
            else{
                result.push_back(gridlist[gridlist.size()-1]);
            }
        }
    }
    if(_WRITETOFILE && fileresult){
        delete fileresult;
        fileresult = nullptr;
    }
    return result;
}


template<>                                                                                                                                      //orders of the dimension
std::vector<std::vector<std::vector<float>>> FiniteDifferenceMethod<std::vector<std::vector<float>>>::PDE2Solve(std::function<void(float,float,float,long long,long long,std::vector<std::vector<std::vector<float>>>&)> pde, int dim1, int dim2, int dim3){
    std::ofstream* fileresult = nullptr;
    if(_WRITETOFILE){
        fileresult = new std::ofstream(filename,std::ios::app);
        std::cout << "Created file.\n";
    }
    try{
        if (this->initialconditions.empty()){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    std::cout << "Entered function\n";
    std::vector<std::vector<float>> current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> rangea = this->Rangelist[0];
    std::array<float,2> rangeb = this->Rangelist[1];
    std::array<float,2> rangec = this->Rangelist[2];
    float stepa = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepb = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepc = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((rangea[1] - rangea[0]) / stepa);
    long long gridsa = round((rangeb[1] - rangeb[0]) / stepb);
    long long gridsb = round((rangec[1] - rangec[0]) / stepc);

    std::vector<std::vector<std::vector<float>>> gridlist;
    for(int i = 0; i < dim1-1; i++){
        std::vector<std::vector<float>> temp = current;
        gridlist.push_back(temp);
    }

    gridlist.push_back(current);
    
    std::vector<std::vector<std::vector<float>>> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        showProgressBar(i,num_iterations,i * stepa);
        std::vector<std::vector<float>> next = gridlist[gridlist.size()-1];
        //apply boundary conditions
        for(const auto& boundarycondition : this->boundaries){
            boundarycondition(next);
        } 
        gridlist.push_back(next);

        #pragma omp parallel for collapse(2)
        for(long long j = floor(dim2 / 2); j < gridsa - static_cast<long long>(ceil(dim2 / 2)); j++){
            for(long long k = floor(dim3 / 2); k < gridsb - static_cast<long long>(ceil(dim3 / 2)); k++){
                pde(stepa,stepb,stepc,j,k,gridlist);
            }
        }

        for(int j = 0; j < gridlist.size()-1; j++){
            gridlist[j] = gridlist[j+1];
        }
        gridlist.pop_back(); // removes "next"

        if(i % this->skips){
            if (_WRITETOFILE)
            {
                for(const auto& yes : gridlist[gridlist.size()-1]){
                    for(const auto& no : yes){
                        *fileresult << no;
                        *fileresult << ",";
                    }
                }
                *fileresult << "\n";
            }
            else{
                result.push_back(gridlist[gridlist.size()-1]);
            }
        }
    }
    if(_WRITETOFILE && fileresult){
        delete fileresult;
        fileresult = nullptr;
    }
    return result;
}

template<>
std::vector<std::vector<std::vector<std::vector<float>>>> FiniteDifferenceMethod<std::vector<std::vector<std::vector<float>>>>::PDE3Solve(std::function<void(float,float,float,float,long long, 
long long,long long,std::vector<std::vector<std::vector<std::vector<float>>>>&)> pde,
                                            int dim1, int dim2, int dim3,int dim4){
    std::ofstream* fileresult = nullptr;
    if(_WRITETOFILE){
        fileresult = new std::ofstream(filename,std::ios::app);
        std::cout << "Created file.\n";
    }
    try{
        if (this->initialconditions.empty()){
            throw std::runtime_error("Not enough initial conditions");
        }
        if (this->Rangelist.front()[1] == 0.0f){
            throw std::runtime_error("Invalid range");
        }
        if (_IGNORE_WARNINGS){
        if (this->Rangelist.size() > 1){
            std::cout << "WARNING: More than 1 range is given for an Ordinary Differential Equation. ASBsolver will take the first range.";
        }
        }
    }
    catch(const std::runtime_error& e){
        std::cerr << "Error thrown" << e.what() << std::endl;
    }
    std::cout << "Entered function\n";
    std::vector<std::vector<std::vector<float>>> current = this->initialconditions.front();
    this->initialconditions.pop_front();
    std::array<float,2> rangea = this->Rangelist[0];
    std::array<float,2> rangeb = this->Rangelist[1];
    std::array<float,2> rangec = this->Rangelist[2];
    std::array<float,2> ranged = this->Rangelist[3];
    float stepa = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepb = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepc = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    float stepd = this->stepsize.front();
    this->stepsize.erase(this->stepsize.begin());
    long long num_iterations = round((rangea[1] - rangea[0]) / stepa);
    long long gridsa = round((rangeb[1] - rangeb[0]) / stepb);
    long long gridsb = round((rangec[1] - rangec[0]) / stepc);
    long long gridsc = round((ranged[1] - ranged[0]) / stepd);

    std::vector<std::vector<std::vector<std::vector<float>>>> gridlist;
    for(int i = 0; i < dim1-1; i++){
        std::vector<std::vector<std::vector<float>>> temp = current;
        gridlist.push_back(temp);
    }

    gridlist.push_back(current);
    
    std::vector<std::vector<std::vector<std::vector<float>>>> result;
    result.push_back(current);
    for(long long i = 0; i < num_iterations; i++){
        showProgressBar(i,num_iterations,i * stepa);
        std::vector<std::vector<std::vector<float>>> next = gridlist[gridlist.size()-1];
        //apply boundary conditions
        for(const auto& boundarycondition : this->boundaries){
            boundarycondition(next);
        } 
        gridlist.push_back(next);

        #pragma omp parallel for collapse(3)
        for(long long j = floor(dim2 / 2); j < gridsa - static_cast<long long>(ceil(dim2 / 2)); j++){
            for(long long k = floor(dim3 / 2); k < gridsb - static_cast<long long>(ceil(dim3 / 2)); k++){
                for(long long l = floor(dim4 / 2); l < gridsc - static_cast<long long>(ceil(dim4 / 2)); l++){
                    pde(stepa,stepb,stepc,stepd,j,k,l,gridlist);
                }
            }
        }

        for(int j = 0; j < gridlist.size()-1; j++){
            gridlist[j] = gridlist[j+1];
        }
        gridlist.pop_back(); // removes "next"

        if(i % this->skips){
            if (_WRITETOFILE)
            {
                for(const auto& yes : gridlist[gridlist.size()-1]){
                    for(const auto& no : yes){
                        for(const auto& maybe : no){
                            *fileresult << maybe;
                            *fileresult << ",";
                        }
                    }
                }
                *fileresult << "\n";
            }
            else{
                result.push_back(gridlist[gridlist.size()-1]);
            }
        }
    }
    if(_WRITETOFILE && fileresult){
        delete fileresult;
        fileresult = nullptr;
    }
    return result;
}

//explicit template instantitation

template class FiniteDifferenceMethod<float>;
template class FiniteDifferenceMethod<std::vector<float>>;
template class FiniteDifferenceMethod<std::vector<std::vector<float>>>;
template class FiniteDifferenceMethod<std::vector<std::vector<std::vector<float>>>>;