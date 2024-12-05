//
// Created by trueno on 10/13/24.
//

#include "Mesh_1D.hpp"
#include <vector>

// Constructor definition
Mesh_1D::Mesh_1D(const int problem_dimension, const double domain_size, const int partition, const int element_order)
    : problem_dimension(problem_dimension), partition(partition), element_order(element_order), domain_size(domain_size) {}

void Mesh_1D::generate_mesh() {
    generateIndividualMesh();
}

void Mesh_1D::generateIndividualMesh() {

    const int degree = element_order;
    const int number_of_nodes = degree * partition + 1;
    const int number_of_elements = partition;
    const int nodes_per_element = degree + 1;

    const double dx = domain_size / (degree * partition);

    // Generate Nodes
    nodes.resize(number_of_nodes);
    for (int i = 0; i < number_of_nodes; ++i) {
        nodes[i] = i * dx;
    }

    // Generate Elements
    elements.resize(number_of_elements, std::vector<double>(nodes_per_element));
    for (int i = 0; i < number_of_elements; ++i) {
        for (int j = 0; j < nodes_per_element; ++j) {
            if (i == 0 && j == 0) {
                elements[i][j] = 1;
            }
            else if (j == 0) {
                elements[i][j] = elements[i - 1][nodes_per_element - 1];
            }
            else {
                elements[i][j] = elements[i][j - 1] + 1;
            }
        }
    }
}

std::vector<double> Mesh_1D::getNodes() {
    return nodes;
}

std::vector<std::vector<double>> Mesh_1D::getElements() {
    return elements;
}