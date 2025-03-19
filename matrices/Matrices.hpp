#pragma once
#include <vector>
#include <array>

namespace ASBsolver {
    class Matrix{
        private:
            std::vector<std::vector<float>> raw;
            
        public:
            int columns;

            int rows;

            Matrix(std::vector<std::vector<float>> matrix);
            
            Matrix multiply(float num);

            Matrix multiply(Matrix matrix);

            Matrix operator*(Matrix matrix){
                return multiply(matrix);
            }

            Matrix operator*(float num){
                return multiply(num);
            }

            void operator=(Matrix m2){
                raw = m2.get();
            }

            Matrix add(Matrix matrix);

            Matrix operator+(Matrix matrix){
                return add(matrix);
            }

            Matrix subtract(Matrix matrix);

            Matrix operator-(Matrix matrix){
                return subtract(matrix);
            }

            std::vector<std::vector<float>>& get();

            float& get(int i, int j);

            Matrix cross(Matrix m2);

            float magnitude();

            float dot(Matrix m2);

            void clear();

            void identitym();
    };

}