//
// Created by Maitreya Limkar on 17-02-2025.
//

#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include "../Eigen/Dense"
#include <vector>

// Type aliases for clarity
using NodeList_1D = Eigen::VectorXd;
using ElementList = std::vector<std::vector<int>>;

// Mesh_1D class declaration
class Mesh_1D {
public:
    Mesh_1D(double domain_size, int partition, const std::vector<int>& element_orders);
    void generateMesh();
    void printMesh() const;

    [[nodiscard]] NodeList_1D getNodeList() const;
    [[nodiscard]] std::vector<ElementList> getElementLists() const;

private:
    double domain_size;
    int partition;
    std::vector<int> element_orders;

    std::vector<NodeList_1D> node_lists;       // Nodes for each element order
    std::vector<ElementList> element_lists; // Element connectivity
    NodeList_1D merged_node_list;              // Unified node list

    void generateIndividualMesh(int degree, NodeList_1D& nl, ElementList& el) const;
    void mergeNodeLists();
    void updateElementLists();
};
#endif //MESH_1D_HPP
