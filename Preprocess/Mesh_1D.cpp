//
// Created by trueno on 10/13/24.
//

#include "Mesh_1D.hpp"
#include <vector>

// Constructor definition
Mesh_1D::Mesh_1D(const int PD, const double domain_size, const int partition, const int element_order)
    : PD(PD), domain_size(domain_size), partition(partition), element_order(element_order) {}

void Mesh_1D::generate_mesh() {
    generateIndividualMesh();
}

void Mesh_1D::generateIndividualMesh() {

    const int degree = element_order;
    const int NoN = degree * partition + 1;
    const int NoE = partition;
    const int NPE = degree + 1;

    const double dx = domain_size / (degree * partition);

    // Generate Nodes
    nodes.resize(NoN);
    for (int i = 0; i < NoN; ++i) {
        nodes[i] = i * dx;
    }

    // Generate Elements
    elements.resize(NoE, std::vector<double>(NPE));
    for (int i = 0; i < NoE; ++i) {
        for (int j = 0; j < NPE; ++j) {
            if (i == 0 && j == 0) {
                elements[i][j] = 1;
            }
            else if (j == 0) {
                elements[i][j] = elements[i - 1][NPE - 1];
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