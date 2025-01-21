//
// Created by David Valle on 09-Dec-24.
//

#include "element.hpp"

element::element(int element_number, int problem_dimension, std::vector<double> node_list,
                 Eigen::VectorXd spatial_coordinate, Eigen::VectorXd material_coordinate, int number_gauss_point,
                 int element_order, double lambda, double mu) {
    this->element_number=element_number;
    this->problem_dimension=problem_dimension;
    this->node_list=node_list;



}

Eigen::VectorXd element::Residual() {
    return Eigen::VectorXd();
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> element::residual_K() {
    return std::pair<Eigen::VectorXd, Eigen::MatrixXd>();
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> element::residual_gauss_K() {
    return std::pair<Eigen::VectorXd, Eigen::MatrixXd>();
}
