#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include <vector>

// Mesh_1D class declaration
class mesh_1D {
    // Initializing the variables
    int problem_dimension = 0, partition = 0, element_order = 0;
    double domain_size = 0.0;
    std::vector<double> node_list;
    std::vector<std::pair<int, int>> element_list;

    void generateIndividualMesh();

public:
    // Constructor for Mesh_1D
    mesh_1D(int problem_dimension, double domain_size, int partition, int element_order);

    // Generation of mesh
    void generate_mesh();

    // Getters for nodes and elements
    std::vector<double>& getNodeList() ;                 // Access node coordinates
    std::vector<std::pair<int, int>>& getElementList(); // Access elements
};

#endif //MESH_1D_HPP