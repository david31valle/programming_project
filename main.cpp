#include <iostream>
#include "Eigen/Dense"
#include "node/node.hpp"
#include "preprocess/mesh_1D.hpp"
#include "preprocess/mesh_2D.hpp"
#include "preprocess/mesh_3D.hpp"

int main() {
    // Input parameters
    int problem_dimension = 3;              // 1D, 2D, or 3D mesh
    int element_order = 1;                  // Linear (1), Quadratic (2), Cubic (3)
    int domain_size = 1;                    // Physical size of the domain
    int partition = 2;                      // Number of partitions
    int lambda = 100;                       // Material property (Lame's First Parameter)
    int mu = 10;                            // Material property (Shear modulus)
    double d = 0.5;                         // Deformation parameter
    double steps = 0.5;                     // Deformation steps
    std::string deformation_type = "EXT";   // Deformation type ("EXT", "COMP", etc.)
    int max_iteration = 10;                 // Maximum Newton-Raphson iterations
    double tol = 1e-10;                     // Tolerance for solver convergence
    std::string boundary_condition = "DBC"; // Boundary condition ("DBC", "PBC")
    std::string gauss_points_values = "On"; // Output Gauss point values ("On"/"Off")

    try {
        switch (problem_dimension) {
            case 1: {
                    //------------ 1D Mesh ------------
                std::cout << "Generating 1D Mesh...\n";
                std::vector<int> element_orders = {element_order};

                // Create and generate the 1D mesh
                Mesh_1D mesh_1D(domain_size, partition, element_orders);
                mesh_1D.generateMesh();
                mesh_1D.printMesh();
                break;
            }

            case 2: {
                    //------------ 2D Mesh ------------
                std::cout << "Generating 2D Mesh...\n";
                std::vector<int> element_orders = {element_order, element_order};

                // Create and generate the 2D mesh (Assuming Mesh2D class exists)
                Mesh_2D mesh_2D(domain_size, partition, element_orders);
                mesh_2D.generateMesh();
                mesh_2D.printMesh();
                break;
            }

            case 3: {
                std::cout << "Generating 3D Mesh...\n";
                std::vector<int> element_orders = {element_order, element_order, element_order};

                // Create and generate the 3D mesh (Assuming Mesh3D class exists)
                Mesh_3D mesh_3D(domain_size, partition, element_orders);
                mesh_3D.generateMesh();
                mesh_3D.printMesh();
                break;
            }

            default:
                throw std::invalid_argument("Invalid problem dimension! Use 1 (1D), 2 (2D), or 3 (3D).");
        }

        std::cout << "Mesh generation complete.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }

    return 0;
}