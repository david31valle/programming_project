#include <iostream>
#include "Eigen/Dense"
#include "node/node.hpp"
#include "preprocess/mesh_1D.hpp"

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

    //---------------- Mesh ----------------
    mesh_1D mesh(problem_dimension, domain_size, partition, element_order);
    mesh.generate_mesh();

    // Retrieve nodes and elements
    //const Eigen::VectorXd nodes = mesh.getNodes();
    //const Eigen::MatrixXd elements = mesh.getElements();
    std::vector<double> nodes = mesh.getNodeList();
    std::vector<std::pair<int, int>> elements = mesh.getElementList();

    // Print Nodes and Elements (for testing purposes)
    //std::cout << "Nodes:\n" << nodes << std::endl;
    //std::cout << "Elements:\n" << elements << std::endl;
    std::cout << "Nodes:\n";
    for (const auto& node : nodes) {
        std::cout << node << " ";
    }
    std::cout << "\nElements:\n";
    for (const auto& element : elements)
    {
        std::cout << '[' << element.first << ", " << element.second << "] "<< std::endl;
    }
}