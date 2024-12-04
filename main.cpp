#include <iostream>
#include "Eigen/Dense"
#include "node/node.hpp"


int main() {
    // Define a 1x3 row vector
    int problem_dimension=1;
    int element_order=1;
    int domain_size=1;
    int partition=10;
    int lambda=100;
    int mu=10;
    double d=0.5;
    int steps=0.5;
    std::string deformation_type= "EXT";
    int max_iteration=10;
    double tol=1e-10;
    std::string boundary_condition="DBC";
    std::string gauss_points_values="On";

    switch (problem_dimension) {
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
    }
    // Initialize vectors with appropriate sizes and values



    return 0;
}