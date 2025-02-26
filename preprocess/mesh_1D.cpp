//
// Created by Maitreya Limkar on 17-02-2025.
//

#include "Mesh_1D.hpp"
#include <iostream>
#include <algorithm>
#include <unordered_map>

// Constructor
Mesh_1D::Mesh_1D(const double domain_size, const int partition, const std::vector<int>& element_orders)
    : domain_size(domain_size), partition(partition), element_orders(element_orders) {}

// Main function to generate the mesh
void Mesh_1D::generateMesh() {
    for (const auto& degree : element_orders) {
        NodeList_1D nl;
        ElementList el;
        generateIndividualMesh(degree, nl, el);
        node_lists.push_back(nl);
        element_lists.push_back(el);
    }

    if (element_orders.size() > 1) {
        mergeNodeLists();
        updateElementLists();
    } else {
        merged_node_list = node_lists[0];
    }
}

// Generate nodes and elements for a given degree
void Mesh_1D::generateIndividualMesh(const int degree, NodeList_1D& nl, ElementList& el) const {
    const int num_nodes = degree * partition + 1;
    const int num_elements = partition;
    const int nodes_per_element = degree + 1;

    double dx = domain_size / (degree * partition);
    nl = NodeList_1D::LinSpaced(num_nodes, 0.0, domain_size);

    el.reserve(num_elements);
    for (int i = 0; i < num_elements; ++i) {
        std::vector<int> element(nodes_per_element);
        if (i == 0) {
            for (int j = 0; j < nodes_per_element; ++j) {
                element[j] = j + 1;
            }
        } else {
            element[0] = el[i - 1].back();
            for (int j = 1; j < nodes_per_element; ++j) {
                element[j] = element[j - 1] + 1;
            }
        }
        el.push_back(element);
    }
}

// Combine node lists without duplicates
void Mesh_1D::mergeNodeLists() {
    std::vector<double> combined_nodes;
    for (const auto& nl : node_lists) {
        for (int i = 0; i < nl.size(); ++i) {
            combined_nodes.push_back(nl(i));
        }
    }

    std::sort(combined_nodes.begin(), combined_nodes.end());
    combined_nodes.erase(std::unique(combined_nodes.begin(), combined_nodes.end()), combined_nodes.end());
    merged_node_list = Eigen::Map<NodeList_1D>(combined_nodes.data(), combined_nodes.size());
}

// Update element indices based on merged node list
void Mesh_1D::updateElementLists() {
    std::unordered_map<double, int> node_map;
    for (int i = 0; i < merged_node_list.size(); ++i) {
        node_map[merged_node_list(i)] = i + 1;  // 1-based indexing
    }

    for (auto& el : element_lists) {
        for (auto& element : el) {
            for (auto& node : element) {
                node = node_map[node_lists[&el - &element_lists[0]](node - 1)];
            }
        }
    }
}

// Print the mesh
void Mesh_1D::printMesh() const {
    std::cout << "Nodes (x):\n" << merged_node_list.transpose() << "\n\n";

    for (size_t i = 0; i < element_lists.size(); ++i) {
        std::cout << "Element Connectivity (Order " << element_orders[i] << "):\n";
        for (const auto& element : element_lists[i]) {
            for (const auto& node : element) {
                std::cout << node << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

// Accessor Methods
NodeList_1D Mesh_1D::getNodeList() const {
    return merged_node_list;
}

std::vector<ElementList> Mesh_1D::getElementLists() const {
    return element_lists;
}