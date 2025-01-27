#include "mesh_1D.hpp"
#include <vector>

// Constructor definition
mesh_1D::mesh_1D(const int problem_dimension, const double domain_size, const int partition, const int element_order)
    : problem_dimension(problem_dimension), domain_size(domain_size), partition(partition), element_order(element_order) {}

void mesh_1D::generate_mesh() {
    generateIndividualMesh();
}

void mesh_1D::generateIndividualMesh() {

    const int degree = element_order;
    const int number_of_nodes = degree * partition + 1;
    const int number_of_elements = partition;
    const int nodes_per_element = degree + 1;

    const double dx = domain_size / (degree * partition);

    // Generate Nodes
    node_list.resize(number_of_nodes);
    for (int i = 0; i < number_of_nodes; ++i) {
        node_list[i] = i * dx;
    }

    // Generate Elements
    element_list.clear();
    for (int i = 0; i < partition; ++i) {
        for (int j = 0; j < degree; ++j) {
            int start_node = i * degree + j;
            element_list.emplace_back(start_node, start_node + 1); // Connectivity
        }
    }
}

std::vector<double>& mesh_1D::getNodeList() {
    return node_list;
}

std::vector<std::pair<int, int>>& mesh_1D::getElementList() {
    return element_list;
}