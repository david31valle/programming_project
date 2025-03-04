#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include "../Eigen/Dense"
#include "../Eigen/src/SparseCore/SparseMatrix.h"
#include "../node/node.hpp"
#include "../Element/Element.cpp"

using namespace Eigen;
class Problem {
public:
    Problem(int PD, std::vector<Node>& NL, std::vector<Element>& EL,
            double domain_size, const Eigen::MatrixXd& BC, const std::string& DEF,
            const std::vector<int>& element_order, double d, int steps,
            int max_iter, double tol, const std::string& GP_vals);
            Eigen::MatrixXd m;

private:
    Eigen::MatrixXd AssignBC(std::vector<Node>& NL, double domain_size, const Eigen::MatrixXd& BC, const Eigen::MatrixXd& F);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Assemble(int PD, std::vector<Node>& NL, std::vector<Element>& EL, int DOFs);
    Eigen::MatrixXd Residual(int PD, std::vector<Node>& NL, std::vector<Element>& EL, int DOFs);
    void update(std::vector<Node>& NL, std::vector<Element>& EL, const Eigen::VectorXd& dx);
    void prescribe(std::vector<Node>& NL, std::vector<Element>& EL, double LF);
    void PostProcess(int PD, const std::vector<Node>& NL, const std::vector<Element>& EL, int step);
};

#endif
