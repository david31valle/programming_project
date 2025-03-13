#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "preprocess/mesh.hpp"
#include "initialization/Initialize.hpp"
#include "Problem/problem.hpp"
#include "element/element.hpp"

int main() {
    // Input parameters
    int problem_dimension = 3;              // 1D, 2D, or 3D mesh
    int element_order = 1;                  // Linear (1), Quadratic (2), Cubic (3)
    int domain_size = 1;                    // Physical size of the domain
    int partition = 3;                      // Number of partitions
    double lambda = 100;                       // Material property (Lame's First Parameter)
    double mu = 10;                            // Material property (Shear modulus)
    double d = 0.5;                         // Deformation parameter
    double steps = 0.5;                     // Deformation steps
    std::string deformation_type = "EXT";   // Deformation type ("EXT", "COMP", etc.)
    int max_iteration = 10;                 // Maximum Newton-Raphson iterations
    double tol = 1e-10;                     // Tolerance for solver convergence
    std::string boundary_condition = "DBC"; // Boundary condition ("DBC", "PBC")
    std::string gauss_points_values = "Off"; // Output Gauss point values ("On"/"Off")

    auto start = std::chrono::high_resolution_clock::now();

    auto [node_list, element_list]=generate_mesh(domain_size, partition, element_order, problem_dimension);
    element_list.array() -= 1;
    auto [Node_List, Element_List]=initialize(problem_dimension, node_list, element_list, domain_size, element_order, lambda, mu);

    //Element_List[0].printElementData();

    problem fem_problem(problem_dimension, Node_List, Element_List, domain_size,
                        boundary_condition, deformation_type, element_order, d,
                        steps, max_iteration, tol, gauss_points_values);
    return 0;
}