//
// Created by David Valle on 09-Dec-24.
//

#ifndef PROGRAMMING_PROJECT_ELEMENT_HPP
#define PROGRAMMING_PROJECT_ELEMENT_HPP
#include <vector>
#include "../Eigen/Dense"

class element {
public:
    int element_number;
    int degree;
    std::vector<double>  node_list;
    Eigen::MatrixXd node_per_element;

    Eigen::VectorXd spatial_coordinate;
    Eigen::VectorXd material_coordinate;
    double lambda;
    double mu;
    int number_gauss_point;
    Eigen::MatrixXd gauss_point;
    int problem_dimension;
    element(int element_number, int problem_dimension, std::vector<double>  node_list,    Eigen::VectorXd spatial_coordinate,
            Eigen::VectorXd material_coordinate, int number_gauss_point, int element_order, double lambda,
            double mu);

    //methods

    Eigen::VectorXd Residual();
    std::pair <Eigen::VectorXd, Eigen::MatrixXd> residual_K();
    std::pair <Eigen::VectorXd, Eigen::MatrixXd> residual_gauss_K();




private:
    Eigen::MatrixXd shape_function;
    Eigen::MatrixXd gradient_shape_functions;
    Eigen::MatrixXd Jacobian;

    Eigen::MatrixXd compute_GP(int number_GP, int problem_dimension);


};


#endif //PROGRAMMING_PROJECT_ELEMENT_HPP
