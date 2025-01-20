//
// Created by trueno on 10/13/24.
//

#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include <vector>

// Mesh_1D class declaration
class Mesh_1D {
    // Initializing the variables
    int problem_dimension = 0, partition = 0, element_order = 0;
    double domain_size = 0.0;
    std::vector<double> node_list;
    std::vector<std::pair<int, int>> element_list;

    void generateIndividualMesh();
    const std::vector<double>& getNodeList() const;

public:
    // Constructor for Mesh_1D
    Mesh_1D(int problem_dimension, double domain_size, int partition, int element_order);

    // Generation of mesh
    void generate_mesh();

    // Getters for nodes and elements
    [[nodiscard]] static const std::vector<double>& getNodeList() ;                 // Access node coordinates
    [[nodiscard]] const std::vector<std::pair<int, int>>& getElementList() const; // Access elements
};

#endif //MESH_1D_HPP