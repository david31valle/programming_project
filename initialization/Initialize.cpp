//
// Created by David Valle on 04-Dec-24.
//

#include "Initialize.hpp"


std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
initialize(int problem_dimension, Eigen::MatrixXd &node_list, Eigen::MatrixXd &element_list, int domain_size,
           int element_order, int lamda, int mu) {
    double tol = 1e-6;
    int NGP = 0;
    int number_of_nodes=node_list.rows();
    int number_of_elements=element_list.rows();

    switch (problem_dimension) {
        case 1:
            switch (element_order) {
                case 1: NGP = 2; break;
                case 2: NGP = 3; break;
                case 3: NGP = 4; break;
                case 4: NGP = 5; break;
                default: std::cerr << "Invalid element_order for PD=1\n"; break;
            }
            break;
        case 2:
            switch (element_order) {
                case 1: NGP = 4; break;
                case 2: NGP = 9; break;
                case 3: NGP = 16; break;
                case 4: NGP = 25; break;
                default: std::cerr << "Invalid element_order for PD=2\n"; break;
            }
            break;
        case 3:
            switch (element_order) {
                case 1: NGP = 8; break;
                case 2: NGP = 27; break;
                case 3: NGP = 64; break;
                case 4: NGP = 125; break;
                default: std::cerr << "Invalid element_order for PD=3\n"; break;
            }
            break;
        default:
            std::cerr << "Invalid PD value\n";
            break;
    }

    int rows=element_list.rows();
    int cols=element_list.cols();
    std::cout<<element_list;
    std::vector<node> Node_List;

    for (int i = 0; i < number_of_nodes; ++i) {
        // Find all occurrences of 'i+1' in 'element_list'
        std::vector<int> element_indices;
        for (int row = 0; row < rows ; ++row) {
            for (int col = 0; col < cols; ++col) {
                if (element_list(row, col) == i+1) {
                    element_indices.push_back(col*rows+row);  // Store element index
                }
            }
        }

        Eigen::VectorXd X = node_list.row(i);
        Node_List.push_back(node(i, problem_dimension, X, X, element_indices));
    }


    return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>();
}
