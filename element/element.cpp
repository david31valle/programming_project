

#include "element.hpp"
#include "../utils/utils.hpp"

// Define static members
Eigen::MatrixXd element::gauss_points;
Eigen::MatrixXd element::shape_functions_N;
Eigen::MatrixXd element::gradient_N_xi;
std::vector<std::vector<double>> element::gauss_points_vector;
std::vector<std::vector<double>> element::shape_functions_N_vector;
std::vector<std::vector<double>> element::gradient_N_xi_vector;


element::element(int element_number, int problem_dimension,
                 Eigen::MatrixXd  node_list,
                 Eigen::MatrixXd spatial_coordinate,
                 Eigen::MatrixXd material_coordinate,
                 int number_gauss_point, int element_order,
                 double lambda, double mu) {
    this->element_number = element_number;
    this->problem_dimension = problem_dimension;
    this->node_list=node_list;
    this->node_per_element=node_list.cols();
    this->spatial_coordinate.resizeLike(spatial_coordinate);
    this->spatial_coordinate = spatial_coordinate;
    this->material_coordinate.resizeLike(material_coordinate);
    this->material_coordinate = material_coordinate;
    this->lambda = lambda;
    this->mu = mu;
    number_GP = number_gauss_point;

    // Initialize shared data (computed only once)
    initialize_static_data(node_per_element, number_gauss_point, problem_dimension, element_order);

    // Use shared static members
    this->shape_functions_N = shape_functions_N;
    this->gradient_N_xi = gradient_N_xi;

    Jacobian = compute_J(material_coordinate, number_GP, problem_dimension);
    gradient_N = compute_GradN(Jacobian, number_gauss_point, problem_dimension, gradient_N_xi_vector);
}

Eigen::MatrixXd element::compute_J(const Eigen::MatrixXd &material_coordinate, int number_GP, int problem_dimension) {
    // Determine matrix dimensions
    int rows = material_coordinate.rows();
    int cols = number_GP * problem_dimension;

    // Initialize J_gp as a zero matrix
    Eigen::MatrixXd J_gp = Eigen::MatrixXd::Zero(rows, cols);

    for (int gp = 0; gp < number_GP; ++gp) {
        // Extract the corresponding gradient block directly from element::gradient_N_xi
        // Adjusting for MATLAB's 1-based indexing (gp+1 → gp in C++)
        //auto test1=gp * problem_dimension;
        //auto  test2= gradient_N_xi.rows();
        Eigen::MatrixXd GradN_xi = gradient_N_xi.block(0, gp * problem_dimension, gradient_N_xi.rows(), problem_dimension);

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

void element::initialize_static_data(int node_per_element, int number_gauss_point, int problem_dimension, int element_order) {
    if (gauss_points.size() == 0) {

        // Compute only once
        gauss_points_vector = compute_gp(number_gauss_point, problem_dimension);
        gauss_points = convertToEigenMatrix(gauss_points_vector);
        //auto testing=build(node_per_element, problem_dimension, number_gauss_point, gauss_points);
        //std::cout<< " number of GP" << number_gauss_point<<std::endl;
        //std::cout <<" number of dp" << problem_dimension<<std::endl;
        //std::cout <<" GP" << gauss_points<<std::endl;

        auto shape_function = compute_N_xi_gp(element_order, gauss_points_vector, problem_dimension);
        shape_functions_N_vector = shape_function.first;
        gradient_N_xi_vector = shape_function.second;

        std::cout << "shape_functions_N_vector (" << shape_functions_N_vector.size() << " x "
                 << (shape_functions_N_vector.empty() ? 0 : shape_functions_N_vector[0].size()) << "):" << std::endl;

        /*for (const auto& row : shape_functions_N_vector) {
            for (double value : row) {
                std::cout << value << "\t";  // Use tab for better spacing
            }
            std::cout << std::endl;  // New line for each row
        }*/


         gauss_points = convertToEigenMatrix(gauss_points_vector);
         shape_functions_N = convertToEigenMatrix(shape_functions_N_vector);
         gradient_N_xi = convertToEigenMatrix(gradient_N_xi_vector);
        //std::cout <<"in the creation"<<std::endl;
         //std::cout<< shape_functions_N<< std::endl;



    }
}

Eigen::MatrixXd element::convertToEigenMatrix(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return Eigen::MatrixXd();  // Return an empty matrix if input is empty
    size_t rows = vec.size();
    size_t cols = vec[0].size();
    Eigen::MatrixXd mat(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mat(i, j) = vec[i][j];  // Assign each element individually
        }
    }

    return mat;
}


void element::compute_R() {
    Eigen::VectorXd x_gp;
    Eigen::MatrixXd Gradx_gp;
    x_gp = spatial_coordinate * number_GP;
    Gradx_gp = spatial_coordinate * gradient_N_xi;

    Eigen::MatrixXd II = Eigen::MatrixXd::Identity(problem_dimension, problem_dimension);
    Eigen::VectorXd R = Eigen::VectorXd::Zero(node_per_element* problem_dimension);

    int NGP = gauss_points.cols();
    Eigen::VectorXd wp = gauss_points.row(gauss_points.rows() - 1);

    for (int gp = 0; gp < NGP; ++gp) {
        Eigen::VectorXd N = shape_functions_N.col(gp);   // Access shape function
        Eigen::MatrixXd GradN = gradient_N_xi.block(0, gp * problem_dimension, gradient_N_xi.rows(), problem_dimension);
        double JxW = Jacobian.block(0, gp * problem_dimension, Jacobian.rows(), problem_dimension).determinant() * wp(gp);

        Eigen::VectorXd x = x_gp.segment(gp, 1);
        Eigen::MatrixXd Gradx = Gradx_gp.block(0, gp * problem_dimension, Gradx_gp.rows(), problem_dimension);

        Eigen::MatrixXd F = Gradx;
        Eigen::MatrixXd Finv = F.inverse();
        Eigen::MatrixXd FinvT = Finv.transpose();
        double J = F.determinant();
        Eigen::MatrixXd P = mu * (F - FinvT) + 0.5 * lambda * (J * J - 1.0) * FinvT;

        Eigen::VectorXd R1 = Eigen::VectorXd::Zero(problem_dimension);
        Eigen::MatrixXd R2 = P;

        for (int I = 0; I < node_per_element; ++I) {
            Eigen::VectorXd R_tmp = (R1 * N(I) + R2 * GradN.row(I).transpose()) * JxW;
            R.segment(I * problem_dimension, problem_dimension) += R_tmp;
        }
    }
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> element::compute_RK() {
    Eigen::VectorXd R = Eigen::VectorXd::Zero(node_per_element * problem_dimension);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(node_per_element * problem_dimension,
                                              node_per_element * problem_dimension);

    Eigen::MatrixXd x_gp = spatial_coordinate * shape_functions_N;
    Eigen::MatrixXd Gradx_gp = spatial_coordinate * gradient_N;
    Eigen::MatrixXd II = Eigen::MatrixXd::Identity(problem_dimension, problem_dimension);

    int NGP = gauss_points.cols();
    Eigen::VectorXd wp = gauss_points.row(gauss_points.rows() - 1);

    for (int gp = 0; gp < NGP; ++gp) {
        Eigen::VectorXd N = shape_functions_N.col(gp);
        Eigen::MatrixXd GradN = gradient_N.block(0, gp * problem_dimension, gradient_N.rows(), problem_dimension);
        double JxW =
                Jacobian.block(0, gp * problem_dimension, Jacobian.rows(), problem_dimension).determinant() * wp(gp);

        Eigen::MatrixXd x = x_gp.col(gp);
        Eigen::MatrixXd Gradx = Gradx_gp.block(0, gp * problem_dimension, Gradx_gp.rows(), problem_dimension);

        Eigen::MatrixXd F = Gradx;
        Eigen::MatrixXd Finv = F.inverse();
        Eigen::MatrixXd FinvT = Finv.transpose();
        double J = F.determinant();
        Eigen::MatrixXd P = mu * (F - FinvT) + 0.5 * lambda * (J * J - 1.0) * FinvT;

        Eigen::VectorXd R1 = Eigen::VectorXd::Zero(problem_dimension);
        Eigen::MatrixXd dR1_dx = Eigen::MatrixXd::Zero(problem_dimension, problem_dimension);
        std::vector<Eigen::MatrixXd> dR1_dGradx(problem_dimension,
                                                Eigen::MatrixXd::Zero(problem_dimension, problem_dimension));

        Eigen::MatrixXd R2 = P;
        std::vector<Eigen::MatrixXd> dR2_dx(problem_dimension,
                                            Eigen::MatrixXd::Zero(problem_dimension, problem_dimension));
        Eigen::MatrixXd dR2_dGradx = Eigen::MatrixXd::Zero(problem_dimension * problem_dimension,
                                                           problem_dimension * problem_dimension);

        Eigen::MatrixXd dP_dF = Eigen::MatrixXd::Zero(problem_dimension * problem_dimension,
                                                      problem_dimension * problem_dimension);

        for (int i = 0; i < problem_dimension; ++i) {
            for (int j = 0; j < problem_dimension; ++j) {
                for (int k = 0; k < problem_dimension; ++k) {
                    for (int l = 0; l < problem_dimension; ++l) {
                        dP_dF(i * problem_dimension + j, k * problem_dimension + l) =
                                lambda * J * J * FinvT(i, j) * FinvT(k, l)
                                - 0.5 * lambda * (J * J - 1) * FinvT(i, l) * FinvT(k, j)
                                + mu * (II(i, k) * II(j, l) + FinvT(i, l) * FinvT(k, j));
                    }
                }
            }
        }

        dR2_dGradx = dP_dF;

        for (int I = 0; I < node_per_element; ++I) {
            R.segment(I * problem_dimension, problem_dimension) +=
                    (R1 * N(I) + R2 * GradN.row(I).transpose()) * JxW;
        }

        for (int I = 0; I < node_per_element; ++I) {
            for (int J = 0; J < node_per_element; ++J) {
                Eigen::MatrixXd K1 = dR1_dx * N(I) * N(J) * JxW;

                Eigen::MatrixXd K2 = Eigen::MatrixXd::Zero(problem_dimension, problem_dimension);
                for (int i = 0; i < problem_dimension; ++i) {
                    for (int j = 0; j < problem_dimension; ++j) {
                        for (int l = 0; l < problem_dimension; ++l) {
                            K2(i, j) += dR1_dGradx[l](i, j) * N(I) * GradN(J, l) * JxW;
                        }
                    }
                }

                Eigen::MatrixXd K3 = Eigen::MatrixXd::Zero(problem_dimension, problem_dimension);
                for (int i = 0; i < problem_dimension; ++i) {
                    for (int k = 0; k < problem_dimension; ++k) {
                        for (int j = 0; j < problem_dimension; ++j) {
                            K3(i, k) += dR2_dx[j](i, k) * GradN(I, j) * N(J) * JxW;
                        }
                    }
                }

                Eigen::MatrixXd K4 = Eigen::MatrixXd::Zero(problem_dimension, problem_dimension);
                for (int i = 0; i < problem_dimension; ++i) {
                    for (int k = 0; k < problem_dimension; ++k) {
                        for (int j = 0; j < problem_dimension; ++j) {
                            for (int m = 0; m < problem_dimension; ++m) {
                                K4(i, k) +=
                                        dR2_dGradx(i * problem_dimension + j, k * problem_dimension + m) * GradN(I, j) *
                                        GradN(J, m) * JxW;
                            }
                        }
                    }
                }

                Eigen::MatrixXd K_tmp = K1 + K2 + K3 + K4;
                K.block(I * problem_dimension, J * problem_dimension, problem_dimension, problem_dimension) += K_tmp;
            }
        }

        /*if (element_number == 0) {
            std::cout << "R" << std::endl;
            std::cout << R << std::endl;
            std::cout << "K" << std::endl;
            std::cout << K << std::endl;
        }*/
    }
    return {R, K};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> element::compute_RK_GP(int NGP_val) {
    Eigen::MatrixXd R_GP = Eigen::MatrixXd::Zero(node_per_element, NGP_val);
    Eigen::MatrixXd K_GP = Eigen::MatrixXd::Zero(node_per_element, node_per_element);

    Eigen::MatrixXd II = Eigen::MatrixXd::Identity(problem_dimension, problem_dimension);
    int NGP = gauss_points.cols();
    Eigen::VectorXd wp = gauss_points.row(gauss_points.rows() - 1);

    Eigen::VectorXd x_gp = spatial_coordinate * number_GP;
    Eigen::MatrixXd Gradx_gp = spatial_coordinate * gradient_N_xi;

    for (int gp = 0; gp < NGP; ++gp) {
        Eigen::VectorXd N = shape_functions_N.col(gp);
        Eigen::MatrixXd GradN = gradient_N_xi.block(0, gp * problem_dimension, gradient_N_xi.rows(), problem_dimension);
        double JxW = Jacobian.block(0, gp * problem_dimension, Jacobian.rows(), problem_dimension).determinant() * wp(gp);

        Eigen::VectorXd x = x_gp.segment(gp, 1);
        Eigen::MatrixXd Gradx = Gradx_gp.block(0, gp * problem_dimension, Gradx_gp.rows(), problem_dimension);

        Eigen::MatrixXd F = Gradx;
        Eigen::MatrixXd Finv = F.inverse();
        Eigen::MatrixXd FinvT = Finv.transpose();
        double J = F.determinant();
        Eigen::MatrixXd P = mu * (F - FinvT) + 0.5 * lambda * (J * J - 1.0) * FinvT;

        Eigen::VectorXd R1 = Eigen::Map<Eigen::VectorXd>(F.transpose().data(), F.size());
        Eigen::VectorXd R2 = Eigen::Map<Eigen::VectorXd>(P.transpose().data(), P.size());
        Eigen::VectorXd R = Eigen::VectorXd(R1.size() + R2.size());
        R << R1, R2;

        for (int s = 0; s < R.size(); ++s) {
            for (int I = 0; I < node_per_element; ++I) {
                double R_tmp = R(s) * N(I) * JxW;
                R_GP(I, s) += R_tmp;
            }
        }

        for (int I = 0; I < node_per_element; ++I) {
            for (int J = 0; J < node_per_element; ++J) {
                double K_tmp = N(I) * N(J) * JxW;
                K_GP(I, J) += K_tmp;
            }
        }
    }
    return {R_GP, K_GP};
}


Eigen::VectorXd element::Residual(double dt) {
    return compute_RK().first;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> element::RK(double dt) {
    return compute_RK();
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> element::RK_GP(int NGP_val) {
    return compute_RK_GP(NGP_val);
}


void element::printElementData() const {
    std::cout << "Element Number: " << element_number << std::endl;
    std::cout << "Problem Dimension: " << problem_dimension << std::endl;

    std::cout << "Node List:" << std::endl;
    std::cout << node_list << std::endl;

    std::cout << "Spatial Coordinates: " << std::endl << spatial_coordinate << std::endl;
    std::cout << "Material Coordinates: " << std::endl << material_coordinate << std::endl;

    std::cout << "Lambda: " << lambda << std::endl;
    std::cout << "Mu: " << mu << std::endl;

    std::cout << "Number of Gauss Points: " << number_GP << std::endl;

    std::cout << "Gauss Points:" << std::endl << gauss_points << std::endl;
    std::cout << "Shape Functions (N): " << std::endl << shape_functions_N << std::endl;
    std::cout << "Gradient N xi: " << std::endl << gradient_N_xi << std::endl;
    std::cout << "Jacobian: " << std::endl << Jacobian << std::endl;
    std::cout << "gradient_N " << std::endl << gradient_N<< std::endl;
}


