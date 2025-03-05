
#ifndef PROGRAMMING_PROJECT_INITIALIZE_HPP
#define PROGRAMMING_PROJECT_INITIALIZE_HPP


#include <iostream>
#include "../Eigen/Dense"
#include "../node/node.hpp"

std::pair<Eigen::MatrixXd,Eigen::MatrixXd> initialize(int problem_dimension,
                                                      Eigen::MatrixXd& node_list,
                                                      Eigen::MatrixXd& element_list,
                                                      int domain_size,
                                                      int element_order,
                                                      int lamda,
                                                      int mu
);



#endif //PROGRAMMING_PROJECT_INITIALIZE_HPP
