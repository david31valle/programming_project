#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include "../Eigen/Dense"
#include "../Eigen/Sparse"
#include "../node/node.hpp"
#include "../element/element.cpp"

using namespace Eigen;
class Problem {
public:
    Problem(int PD, std::vector<node>& NL, std::vector<element>& EL,
            double domain_size, const Eigen::MatrixXd& BC, const std::string& DEF,
            const std::vector<int>& element_order, double d, int steps,
            int max_iter, double tol, const std::string& GP_vals);
            Eigen::MatrixXd m;

private:
    Eigen::MatrixXd AssignBC(std::vector<node>& NL, double domain_size, const Eigen::MatrixXd& BC, const Eigen::MatrixXd& F);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Assemble(int PD, std::vector<node>& NL, std::vector<element>& EL, int DOFs);
    Eigen::MatrixXd Residual(int PD, std::vector<node>& NL, std::vector<element>& EL, int DOFs);
    void update(std::vector<node>& NL, std::vector<element>& EL, const Eigen::VectorXd& dx);
    void prescribe(std::vector<node>& NL, std::vector<element>& EL, double LF);
    void PostProcess(int PD, const std::vector<node>& NL, const std::vector<element>& EL, int step);
};

#endif
