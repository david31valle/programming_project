//
// Created by David Valle on 09-Dec-24.
//

#ifndef PROGRAMMING_PROJECT_ELEMENT_HPP
#define PROGRAMMING_PROJECT_ELEMENT_HPP
#include <vector>
#include "../Eigen/Dense"
#include "../utils/utils.hpp"
class element {
public:
    int element_number;
    int degree;
    std::vector<double>  node_list;
    Eigen::MatrixXd node_per_element;
    Eigen::MatrixXd spatial_coordinate;
    Eigen::MatrixXd material_coordinate;
    double lambda;
    double mu;
    int number_GP;
    std::vector<std::vector<double>> gauss_points;
    int problem_dimension;
    std::vector<std::vector<double>> shape_functions_N;
    std::vector<std::vector<double>> gradient_N_xi;
    element(int element_number, int problem_dimension, std::vector<double>  node_list,    Eigen::MatrixXd spatial_coordinate,
            Eigen::MatrixXd material_coordinate, int number_gauss_point, int element_order, double lambda,
            double mu);

    //methods

    Eigen::MatrixXd compute_J(const Eigen::MatrixXd& material_coordinate, int number_GP, int problem_dimension, const std::vector<std::vector<double>>& gradient_N_xi);
    Eigen::MatrixXd compute_GradN(const Eigen::MatrixXd& J_gp, int number_GP, int problem_dimension, const std::vector<std::vector<double>>& gradient_N_xi);
    Eigen::VectorXd Residual();
    std::pair <Eigen::VectorXd, Eigen::MatrixXd> residual_K();
    std::pair <Eigen::VectorXd, Eigen::MatrixXd> residual_gauss_K();

    //private
    Eigen::MatrixXd shape_function;
    Eigen::MatrixXd gradient_shape_functions;
    Eigen::MatrixXd Jacobian;
    Eigen::MatrixXd gradient_N;

    //Eigen::MatrixXd compute_GP(int number_GP, int problem_dimension);
    //std::vector<std::vector<double>> compute_gp(int number_gauss_point, int problem_dimension);

};


#endif //PROGRAMMING_PROJECT_ELEMENT_HPP
