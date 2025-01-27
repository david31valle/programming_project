//
// Created by Maitreya Limkar on 24-10-2024.
//

#include "Mesh_2D.hpp"
#include <vector>

// Constructor definition
Mesh_2D::Mesh_2D(const int PD, const double domain_size, const int partition, const int element_order)
    : PD(PD), domain_size(domain_size), partition(partition), element_order(element_order) {}

void Mesh_2D::generate_mesh() {
    generateIndividualMesh();
}

void Mesh_2D::generateIndividualMesh() {

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

std::vector<double> Mesh_2D::getNodes() {
    return nodes;
}

std::vector<std::vector<double>> Mesh_2D::getElements() {
    return elements;
}