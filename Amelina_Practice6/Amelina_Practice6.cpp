#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

// Функція для обчислення норми матриці ||A||нескінченність
double matrixNormInfinity(const std::vector<std::vector<double>>& matrix) {
    double norm = 0.0;

#pragma omp parallel for reduction(max: norm)
    for (size_t i = 0; i < matrix.size(); i++) {
        double rowSum = 0.0;

        for (size_t j = 0; j < matrix[i].size(); j++) {
            rowSum += std::abs(matrix[i][j]);
        }

        norm = std::max(norm, rowSum);
    }

    return norm;
}

// Функція для обчислення норми матриці ||A||F
double matrixNormFrobenius(const std::vector<std::vector<double>>& matrix) {
    double norm = 0.0;

#pragma omp parallel for reduction(+: norm)
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            norm += std::pow(matrix[i][j], 2);
        }
    }

    norm = std::sqrt(norm);

    return norm;
}

int main() {
    std::vector<std::vector<double>> matrix = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };

    double normInfinity = matrixNormInfinity(matrix);
    double normFrobenius = matrixNormFrobenius(matrix);

    std::cout << "Matrix norm (Infinity): " << normInfinity << std::endl;
    std::cout << "Matrix norm (Frobenius): " << normFrobenius << std::endl;

    return 0;
}