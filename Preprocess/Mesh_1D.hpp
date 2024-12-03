//
// Created by trueno on 10/13/24.
//

#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include <vector>

// Mesh_1D class declaration
class Mesh_1D {
private:
    // Initializing the variables
    int PD = 0, partition = 0, element_order = 0;
    double domain_size = 0.0;
    std::vector<double> NL, EL_1, EL_2, EL_3;

    std::vector<double> nodes;
    std::vector<std::vector<double>> elements;

    void generateIndividualMesh();

public:
    // Constructor for Mesh_1D
    Mesh_1D(int PD, double domain_size, int partition, int element_order);

    // Generation of mesh
    void generate_mesh();

    // Getters for nodes and elements
    std::vector<double> getNodes();
    std::vector<std::vector<double>> getElements();
};

#endif //MESH_1D_HPP