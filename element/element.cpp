

#include "element.hpp"
#include "../utils/utils.hpp"

element::element(int element_number, int problem_dimension, std::vector<double> node_list,
                 Eigen::MatrixXd spatial_coordinate, Eigen::MatrixXd material_coordinate, int number_gauss_point,
                 int element_order, double lambda, double mu) {
    this->element_number=element_number;
    this->problem_dimension=problem_dimension;
    this->node_list=node_list;
    node_per_element.resize(node_list.size(),2);
    this->spatial_coordinate.resizeLike(spatial_coordinate);
    this->spatial_coordinate=spatial_coordinate;
    this->material_coordinate.resizeLike(material_coordinate);
    this->material_coordinate=material_coordinate;
    this->lambda=lambda;
    this->mu=mu;
    number_GP=number_gauss_point;
    gauss_points = compute_gp(number_GP, problem_dimension);
    auto shape_function= compute_N_xi_gp(degree, gauss_points, problem_dimension);
    shape_functions_N=shape_function.first;
    gradient_N_xi=shape_function.second;
    Jacobian= compute_J(material_coordinate, number_GP, problem_dimension, gradient_N_xi);
    gradient_N= compute_GradN(Jacobian, number_gauss_point, problem_dimension, gradient_N_xi);

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

Eigen::MatrixXd element::compute_J(const Eigen::MatrixXd &material_coordinate, int number_GP, int problem_dimension,
                                   const std::vector<std::vector<double>> &gradient_N_xi) {
    // Determine matrix dimensions
    int rows = material_coordinate.rows();
    int cols = number_GP * problem_dimension;

    // Initialize J_gp as a zero matrix
    Eigen::MatrixXd J_gp = Eigen::MatrixXd::Zero(rows, cols);

    for (int gp = 0; gp < number_GP; ++gp) {
        // Extract the corresponding gradient block from the std::vector representation
        Eigen::MatrixXd GradN_xi(gradient_N_xi.size(), problem_dimension);

        for (size_t i = 0; i < gradient_N_xi.size(); ++i) {
            for (int j = 0; j < problem_dimension; ++j) {
                GradN_xi(i, j) = gradient_N_xi[i][gp * problem_dimension + j];
            }
        }

        // Compute Jacobian
        Eigen::MatrixXd J = material_coordinate * GradN_xi;

        // Store J into J_gp at the correct position
        J_gp.block(0, gp * problem_dimension, J.rows(), J.cols()) = J;
    }

    return J_gp; // Return the final Jacobian matrix
}

Eigen::MatrixXd element::compute_GradN(const Eigen::MatrixXd &J_gp, int number_GP, int problem_dimension,
                                       const std::vector<std::vector<double>> &gradient_N_xi) {
    // Determine matrix dimensions
    int numNodes = gradient_N_xi.size(); // Number of nodes (rows in GradN_xi)
    int cols = number_GP * problem_dimension; // Total columns in GradN_X_gp

    // Initialize GradN_X_gp as a zero matrix
    Eigen::MatrixXd GradN_X_gp = Eigen::MatrixXd::Zero(numNodes, cols);

    for (int gp = 0; gp < number_GP; ++gp) {
        // Extract the corresponding gradient block
        Eigen::MatrixXd GradN_xi(numNodes, problem_dimension);
        for (int i = 0; i < numNodes; ++i) {
            for (int j = 0; j < problem_dimension; ++j) {
                GradN_xi(i, j) = gradient_N_xi[i][gp * problem_dimension + j];
            }
        }

        // Extract the Jacobian block
        Eigen::MatrixXd J = J_gp.block(0, gp * problem_dimension, problem_dimension, problem_dimension);

        // Compute the inverse transpose of J
        Eigen::MatrixXd JinvT = J.inverse().transpose();

        // Compute GradN_X_gp at this Gauss point
        Eigen::MatrixXd GradN_X = (JinvT * GradN_xi.transpose()).transpose();

        // Store result in GradN_X_gp
        GradN_X_gp.block(0, gp * problem_dimension, GradN_X.rows(), GradN_X.cols()) = GradN_X;
    }

    return GradN_X_gp;
}
/*
Eigen::MatrixXd element::compute_GP(int number_GP, int problem_dimension) {
Eigen::MatrixXd gauss_point;


    return Eigen::MatrixXd();
}
*/


