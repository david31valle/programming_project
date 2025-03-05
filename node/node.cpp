//
// Created by David Valle on 01-Dec-24.
//

#include "node.hpp"

node::node(int node_number, int problem_dimension,  Eigen::VectorXd X_material_position,
           Eigen::VectorXd x_spatial_position, std::vector<int> element_list) {
    this->node_number = node_number;
    this->problem_dimension = problem_dimension;

    // Resize member vectors and assign passed values
    this->X_material_position.resize(problem_dimension, problem_dimension);
    this->X_material_position = X_material_position;

    this->x_spatial_position.resize(problem_dimension, problem_dimension);
    this->x_spatial_position = x_spatial_position;

    this->element_list=element_list;

    initialization(problem_dimension);

}

void node::initialization(int PD) {
    boundary_condition.resize(PD);
    boundary_condition.setOnes();

    boundary_condition_value.resize(PD);
    boundary_condition_value.setZero();

    DOF.resize(PD);
    DOF.setZero();

    degree_constrained.resize(PD);
    degree_constrained.setZero();

    global_index.resize(PD);
    global_index.setZero();

    size_t size = PD * PD + PD * PD;

    gauss_point_BC.resize(size);
    gauss_point_BC.setOnes();

    gauss_point_DOF.resize(size);
    gauss_point_DOF.setZero();

    gauss_point_values.resize(size);
    gauss_point_values.setZero();
}

