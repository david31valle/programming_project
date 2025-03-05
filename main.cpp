#include <iostream>
#include "Eigen/Dense"
#include "preprocess/mesh.hpp"
#include "initialization/Initialize.hpp"


int main() {
    // Input parameters
    int problem_dimension = 3;              // 1D, 2D, or 3D mesh
    int element_order = 1;                  // Linear (1), Quadratic (2), Cubic (3)
    int domain_size = 1;                    // Physical size of the domain
    int partition = 10;                      // Number of partitions
    int lambda = 100;                       // Material property (Lame's First Parameter)
    int mu = 10;                            // Material property (Shear modulus)
    double d = 0.5;                         // Deformation parameter
    double steps = 0.5;                     // Deformation steps
    std::string deformation_type = "EXT";   // Deformation type ("EXT", "COMP", etc.)
    int max_iteration = 10;                 // Maximum Newton-Raphson iterations
    double tol = 1e-10;                     // Tolerance for solver convergence
    std::string boundary_condition = "DBC"; // Boundary condition ("DBC", "PBC")
    std::string gauss_points_values = "On"; // Output Gauss point values ("On"/"Off")


    auto [node_list, element_list]=generate_mesh(domain_size, partition, element_order, problem_dimension);
    auto test=initialize(problem_dimension, node_list, element_list, domain_size, element_order, lambda, mu);
    return 0;
}