//
// Created by David Valle on 28-Jan-25.
//


#include <iostream>

#include "utils/utils.hpp"
#include "element/element.hpp"
int main(){
    std::cout<<"hello"<<std::endl;

    //std::vector<std::vector<double>>gp=compute_gp(2,1);
    //auto N_grad=compute_N_xi_gp(1,gp,1);

    int element_number = 1;
    int problem_dimension = 2;
    std::vector<double> node_list = {0.0, 1.0, 2.0}; // Example node list
    int number_gauss_point = 2;
    int element_order = 1;
    double lambda = 1.0;
    double mu = 0.5;

    // Example spatial and material coordinate matrices
    Eigen::MatrixXd spatial_coordinate(3, 2);
    spatial_coordinate << 0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0;

    Eigen::MatrixXd material_coordinate = spatial_coordinate; // Assume the same for testing

    // Initialize the element object
    element test_element(element_number, problem_dimension, node_list, spatial_coordinate, material_coordinate,
                         number_gauss_point, element_order, lambda, mu);
    return 0;


}