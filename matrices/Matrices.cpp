#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include "Matrices.hpp"


using namespace ASBsolver;

Matrix Matrix::multiply(float num){
    auto result = this->raw;
    for(int i = 0; i < this->raw.size(); i++){
        for(int j = 0; j < this->raw[0].size(); j++){
            result[i][j] = this->raw[i][j] * num;
        }
    }
    return result;
}

Matrix::Matrix(std::vector<std::vector<float>> matrix){
    this->raw = matrix;
    this->rows = matrix.size();
    this->columns = matrix[0].size();
}

std::vector<std::vector<float>>& Matrix::get(){
    return this->raw;
}

float& Matrix::get(int i, int j){
    return this->raw[i][j];
}


Matrix Matrix::multiply(Matrix matrix){
    std::vector<std::vector<float>> result;
    int resultrows = this->rows;
    int resultcolumns = matrix.columns;
    for(int i = 0; i < resultrows; i++){
        std::vector<float> temp;
        for(int j = 0; j < resultcolumns; j++){
            temp.push_back(0);
        }
        result.push_back(temp);
    }
    Matrix resultmatrix(result);
    for(int i = 0; i < this->raw.size(); i++){
        for(int j = 0; j < this->raw[0].size(); j++){
            resultmatrix.get(i,j) += this->get(i,j) * matrix.get(j,i);
        }
    }
    return resultmatrix;
}

Matrix Matrix::add(Matrix matrix){
    std::vector<std::vector<float>> result = this->get();
    Matrix resultmatrix(result);
    for(int i = 0; i < this->rows; i++){
        for(int j = 0; j < this->columns; j++){
            resultmatrix.get(i,j) += matrix.get(i,j);
        }
    }
    return resultmatrix;
}

Matrix Matrix::subtract(Matrix matrix){
    std::vector<std::vector<float>> result = this->get();
    Matrix resultmatrix(result);
    for(int i = 0; i < this->rows; i++){
        for(int j = 0; j < this->columns; j++){
            resultmatrix.get(i,j) -= matrix.get(i,j);
        }
    }
    return resultmatrix;
}

Matrix Matrix::cross(Matrix m2){
    std::vector<std::vector<float>> result;
    std::vector<float> temp;
    for(int i = 0; i < 3; i++){
        temp.push_back(0);
    }
    result.push_back(temp);
    std::vector<std::vector<float>> M1 = this->get();
    std::vector<std::vector<float>> M2 = m2.get();
    if(M1[0].size() > 2 && M2[0].size() > 2){
        result[0][0] = M1[0][1] * M2[0][2] - M1[0][2] * M2[0][1];
        result[0][1] = M1[0][2] * M2[0][0] - M1[0][0] * M2[0][2];
        result[0][2] = M1[0][0] * M2[0][1] - M1[0][1] * M2[0][0];
        Matrix resultmatrix(result);
        return resultmatrix;
    }
    else{
        std::cout << "Cross product is defined specifically for 1x3 matrices\n";
        Matrix resultmatrix(result);
        return resultmatrix;
    }
}

float Matrix::magnitude(){
    std::vector<std::vector<float>> M1 = this->get();
    if(M1[0].size() > 2){
        return sqrt(pow(M1[0][0],2) + pow(M1[0][1],2) + pow(M1[0][2],2));
    }
    else{
        return sqrt(pow(M1[0][0],2) + pow(M1[0][1],2));
    }
}



float Matrix::dot(Matrix m2){
    std::vector<std::vector<float>> M1 = this->get();
    std::vector<std::vector<float>> M2 = m2.get();
    float result = 0;
    for (size_t i = 0; i < M1[0].size(); ++i) {
        result += M1[0][i] * M2[0][i];
    }
    return result;
}

void Matrix::clear(){
    for(int i = 0; i < this->rows; i++){
        for(int j = 0; j < this->columns; j++){
            this->get(i,j) = 0;
        }
    }
}

void Matrix::identitym(){
    for(int i = 0; i < this->rows; i++){
        for(int j = 0; j < this->columns; j++){
            if(i == j){
                this->get(i,j) = 1;
            }
            else{
                this->get(i,j) = 0;
            }
        }
    }
}

float Matrix::subtendedangle(Matrix m2) {
    float thing = this->dot(m2);
    return acos((thing / this->magnitude()) / m2.magnitude());
}

Matrix Matrix::functoarray(std::function<float(float)> func, float start, float end, float step) {
    std::vector<float> functoarray;
    for (float i = start; i < end; i += step) {
        functoarray.push_back(func(i));
    }
    std::vector<std::vector<float>> m1;
    m1.emplace_back(functoarray);
    Matrix M1(m1);
    return M1;
}

Matrix functoarray(std::function<float(float, float)> func, float start, float start2, float end, float end2, float step, float step2) {
    std::vector<std::vector<float>> functoarray;
    for (float i = start; i < end; i += step) {
        std::vector<float> temporary;
        for (float j = start2; j < end2; j += step2) {
            temporary.push_back(func(i, j));
        }
        functoarray.push_back(temporary);
    }
    Matrix M1(functoarray);
    return M1;
}