#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include "../Eigen/Dense"
#include "../Eigen/Sparse"
#include "../Eigen/IterativeLinearSolvers"
#include "../node/node.hpp"
#include "../element/element.hpp"
#include "iostream"
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>


using namespace Eigen;
class problem {
public:
    // Member variables
    int problem_dimension;
    std::vector<node> Node_List;
    std::vector<element> Element_List;
    std::string boundary_condition;
    int domain_size;
    std::string deformation_type;
    int element_order;
    double d;
    int steps;
    int max_iter;
    double tol;
    std::string gauss_points_values;
    Eigen::MatrixXd F;
    std::string filename;
    double load_factor = 0;
    Eigen::MatrixXd m;

    // DOF variables
    int DOFs=0;
    int DOCs=0;
    std::vector<int> CNL; // Constrained node list
    std::vector<int> PNL; // Prescribed node list

    Eigen::VectorXd Rtot;
    Eigen::SparseMatrix<double> Ktot;
    //Eigen::SparseMatrix<double> Kpu;
    //Eigen::SparseMatrix<double> Kpp;

    Eigen::MatrixXd Kuu;
    Eigen::MatrixXd Kpu;
    Eigen::MatrixXd Kpp;

    // Constructor
    problem(int problem_dimension,
            const std::vector<node>& Node_List,
            const std::vector<element>& Element_List,
            int domain_size,
            const std::string& boundary_condition,
            const std::string& deformation_type,
            int element_order,
            double d,
            int steps,
            int max_iter,
            double tol,
            const std::string& gauss_points_values);


private:
    void Assign_BC() ;
    void Assign_DOF_DBC();
    void initialize_F();
    void Assign_GP_DOFs();
    void problem_info();
    void prescribe();
    void assemble();
    void update(const Eigen::VectorXd& dx);
    void Residual(double dt);

};

#endif
