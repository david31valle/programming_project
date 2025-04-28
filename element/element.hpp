

#ifndef PROGRAMMING_PROJECT_ELEMENT_HPP
#define PROGRAMMING_PROJECT_ELEMENT_HPP
#include <vector>
#include "../Eigen/Dense"
#include "../utils/utils.hpp"

#include <iostream>
class element {
public:
    int element_number;
    int degree;
    Eigen::MatrixXd  node_list;
    int node_per_element;
    Eigen::MatrixXd spatial_coordinate;
    Eigen::MatrixXd material_coordinate;
    double lambda;
    double mu;
    int number_GP;
    int problem_dimension;
    element(int element_number,
            int problem_dimension,
            Eigen::MatrixXd  node_list,
            Eigen::MatrixXd spatial_coordinate,
            Eigen::MatrixXd material_coordinate,
            int number_gauss_point,
            int element_order,
            double lambda,
            double mu);

    //methods

    Eigen::MatrixXd compute_J(const Eigen::MatrixXd &material_coordinate, int number_GP, int problem_dimension);
    Eigen::MatrixXd compute_GradN(const Eigen::MatrixXd& J_gp, int number_GP, int problem_dimension, const std::vector<std::vector<double>>& gradient_N_xi);
    static Eigen::MatrixXd convertToEigenMatrix(const std::vector<std::vector<double>>& vec);
    Eigen::VectorXd Residual(double dt) ;
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> RK(double dt);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> RK_GP(int NGP_val) ;

    void compute_R();

    //private
    Eigen::MatrixXd shape_function;
    Eigen::MatrixXd gradient_shape_functions;
    Eigen::MatrixXd Jacobian;
    Eigen::MatrixXd gradient_N;

    static std::vector<std::vector<double>> gauss_points_vector;
    static std::vector<std::vector<double>> shape_functions_N_vector;
    static std::vector<std::vector<double>> gradient_N_xi_vector;

    static Eigen::MatrixXd gauss_points;
    static Eigen::MatrixXd shape_functions_N;
    static Eigen::MatrixXd gradient_N_xi;
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> compute_RK();
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_RK_GP(int NGP_val);

    static void initialize_static_data(int node_per_element, int number_gauss_point, int problem_dimension, int element_order);
    void printElementData() const;
    //Eigen::MatrixXd compute_GP(int number_GP, int problem_dimension);
    //std::vector<std::vector<double>> compute_gp(int number_gauss_point, int problem_dimension);

};


#endif //PROGRAMMING_PROJECT_ELEMENT_HPP
