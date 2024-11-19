#include <iostream>
#include "Eigen/Dense"

int main() {
    // Define a 1x3 row vector
    Eigen::RowVector3d rowVector;
    rowVector << 1.0, 2.0, 3.0;

    // Define a 3x3 matrix
    Eigen::Matrix3d matrix;
    matrix << 1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0;

    // Multiply the vector with the matrix
    Eigen::RowVector3d result = rowVector * matrix;

    std::cout << "Result of vector-matrix multiplication:" << std::endl;
    std::cout << result << std::endl;

    return 0;
}