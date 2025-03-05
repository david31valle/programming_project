
#include "../Eigen/Dense"
#include "../Eigen/Sparse"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "problem.hpp"
#include "../node/node.hpp"
#include "../element/element.cpp"


//Auxiliary functions
std::pair<std::vector<node>, std::vector<int>> Assign_BC(
        std::vector<node>& NL,
        double domain_size,
        std::string BC_type,
        const Eigen::MatrixXd& F
);

std::pair<std::vector<node>, std::vector<int>> Assign_DOF_DBC(std::vector<node>& NL);

std::pair<std::vector<node>, std::vector<element>> prescribe(std::vector<node>& NL, std::vector<element>& EL, double LF);

void problemInfo(int problem_dimension,
                 const std::vector<node>& NL,
                 const std::vector<element>& EL,
                 const std::vector<int>& DOFs,  // Changed DOFs to vector<int>
                 int element_order,
                 const std::string& filename);

std::pair<std::vector<node>, int> Assign_GP_DOFs(std::vector<node>& NL);

std::pair<Eigen::MatrixXd, Eigen::SparseMatrix<double>> Assemble_GP(
        int problem_dimension,
        const std::vector<node>& NL,
        const std::vector<element>& EL,
        int GP_DOFs);

std::vector<double> Residual(
        int problem_dimension,
        const std::vector<node>& NL,
        const std::vector<element>& EL,
        int DOFs,
        double dt);


void Problem(
        int problem_dimension,
        std::vector<node>& NL,
        const std::vector<element>& EL,
        const int domain_size,
        const std::string BC,
        const std::string DEF,
        int element_order,
        double d,
        int steps,
        int max_iter,
        double tol,
        const std::string GP_vals
) {

    //Output file Name
    std::string filename = std::to_string(problem_dimension) + "D_Normal_" + std::to_string(EL.size()) +
                           "_EL=_[" + std::to_string(element_order) + "]_.txt";

    //---------------------------------------------------------
    //Create BC conditions
    // Create identity matrix I (problem_dimension x problem_dimension)
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(problem_dimension, problem_dimension);
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(problem_dimension, problem_dimension);
    if (DEF == "EXT") {
        // Case 'EXT': F = I; F(1,1) = 1 + d;
        F = I;
        F(0, 0) = 1 + d;
    }
    else if (DEF == "EXP") {
        // Case 'EXP': F = (1+d) * I;
        F = (1 + d) * I;
    }
    else if (DEF == "SHR") {
        // Case 'SHR': F = I; F(1,2) = d;
        F = I;
        F(0, 1) = d;
    };

    std::vector<int> DOFs;
    std::tie(NL, DOFs) = Assign_BC (NL, domain_size, BC, F);


    // Assigne Gauss Points Degrees of freedoms
    int GP_DOFs;
    if (GP_vals == "On") {
        std::tie(NL, GP_DOFs) = Assign_GP_DOFs(NL);
    }

    // Problem information

    problemInfo(problem_dimension, NL, EL, DOFs, element_order, filename);



}

//-----------------------------------Assign BC--------------------------------------------------------
std::pair<std::vector<node>, std::vector<int>> Assign_BC(
        std::vector<node>& NL,
        double domain_size,
        std::string BC_type,
        const Eigen::MatrixXd& F
) {
    int NoNs = NL.size();  // Number of nodes
    int problem_dimension = NL[0].problem_dimension;     // Problem dimension (assuming all nodes share the same problem_dimension)
    double tol = 1e-6;

    double W = domain_size;  // X-direction (Width)
    double H = domain_size;  // Y-direction (Height)
    double D = domain_size;  // Z-direction (Depth)

    for (int i = 0; i < NoNs; ++i) {
        Eigen::MatrixXd& X = NL[i].X_material_position;  // node position

        // Check based on problem dimension problem_dimension
        if (problem_dimension == 1) {
            if (std::fabs(X(0) - 0) < tol || std::fabs(X[0] - W) < tol) {
                NL[i].boundary_condition = std::vector<int>(problem_dimension, 0);  // Assign BC as zeros vector (problem_dimension x 1)

                // Eigen multiplication F * X - X
                Eigen::VectorXd X_eigen(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) X_eigen(j) = X[j];
                Eigen::VectorXd BCval_eigen = F * X_eigen - X_eigen;

                NL[i].boundary_condition_value.resize(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) NL[i].boundary_condition_value[j] = BCval_eigen(j);
            }
        }
        else if (problem_dimension == 2) {
            if (std::fabs(X[0] - 0) < tol || std::fabs(X[0] - W) < tol) {
                NL[i].boundary_condition = std::vector<int>(problem_dimension, 0);  // Assign BC as zeros vector (problem_dimension x 1)

                // Eigen multiplication F * X - X
                Eigen::VectorXd X_eigen(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) X_eigen(j) = X[j];
                Eigen::VectorXd BCval_eigen = F * X_eigen - X_eigen;

                NL[i].boundary_condition_value.resize(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) NL[i].boundary_condition_value[j] = BCval_eigen(j);
            }
        }
        else if (problem_dimension == 3) {
            if (std::fabs(X[0] - 0) < tol || std::fabs(X[0] - W) < tol) {
                NL[i].boundary_condition = std::vector<int>(problem_dimension, 0);  // Assign BC as zeros vector (problem_dimension x 1)

                // Eigen multiplication F * X - X
                Eigen::VectorXd X_eigen(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) X_eigen(j) = X[j];
                Eigen::VectorXd BCval_eigen = F * X_eigen - X_eigen;

                NL[i].boundary_condition_value.resize(problem_dimension);
                for (int j = 0; j < problem_dimension; ++j) NL[i].boundary_condition_value[j] = BCval_eigen(j);
            }
        }
    }

    // Call Assign_DOF_DBC function to process DOFs (degrees of freedom)
    std::vector<int> DOFs;
    std::tie(NL, DOFs) = Assign_DOF_DBC(NL);

    return std::make_pair(NL, DOFs);
}


//-------------------------------Assign DOF DBC---------------------------------------------
std::pair<std::vector<node>, std::vector<int>> Assign_DOF_DBC(std::vector<node>& NL) {
    int problem_dimension = NL[0].problem_dimension;  // Problem dimension
    int NoN = NL.size();  // Number of nodes

    std::vector<int> DOFs;  // List to store degrees of freedom

    for (int i = 0; i < NoN; ++i) {
        int DOF_counter = 0;

        for (int j = 0; j < problem_dimension; ++j) {
            if (NL[i].boundary_condition[j] == 0) {  // If BC (boundary condition) is free
                DOF_counter++;
                NL[i].DOF[j] = DOFs.size() + DOF_counter;  // Assign DOF number
            } else {
                NL[i].DOF[j] = 0;  // Constrained, so DOF is set to zero
            }
        }

        // Add non-zero DOF values to the DOFs vector
        for (int j = 0; j < problem_dimension; ++j) {
            if (NL[i].DOF[j] > 0) {
                DOFs.push_back(NL[i].DOF[j]);
            }
        }
    }

    return std::make_pair(NL, DOFs);  // Return updated NL and DOFs list
}

//----------------------------------------Assign GP Degrees of Freedoms-----------------------------------------------
std::pair<std::vector<node>, int> Assign_GP_DOFs(std::vector<node>& NL) {
    int NoNs = NL.size();  // Number of nodes
    int GP_DOFs = 0;  // Initialize GP_DOFs counter

    // Iterate through all nodes
    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd& BC = NL[i].gauss_point_BC;  // Boundary condition for Gauss points
        std::vector<int>& DOF = NL[i].GP_DOF;  // Degrees of freedom for Gauss points

        // If any element of BC is 1, increment GP_DOFs and assign to DOF
        if (std::find(BC.begin(), BC.end(), 1) != BC.end()) {
            GP_DOFs++;
            std::fill(DOF.begin(), DOF.end(), GP_DOFs);  // Assign the current GP_DOFs to all DOF elements
        }
    }

    return std::make_pair(NL, GP_DOFs);
}

//----------------------------Problem info-------------------------------

void problemInfo(int problem_dimension,
                 const std::vector<node>& NL,
                 const std::vector<element>& EL,
                 const std::vector<int>& DOFs,  // Changed DOFs to vector<int>
                 int element_order,
                 const std::string& filename)
{
    // Calculate the number of elements and nodes
    int NoEs = EL.size();  // Number of elements
    int NoNs = NL.size();  // Number of nodes

    // Display problem information on the console
    std::cout << "======================================================" << std::endl;
    std::cout << "================  Problem information  ===============" << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Problem dimension                   : " << problem_dimension << std::endl;
    std::cout << "Number of nodes                     : " << NoNs << std::endl;
    std::cout << "Number of bulk elements             : " << NoEs << std::endl;
    std::cout << "Number of DOFs                      : ";
    for (int dof : DOFs) {
        std::cout << dof << " ";
    }
    std::cout << std::endl;
    std::cout << "element order                       : [ "<< element_order << "\n";
    std::cout << "======================================================" << std::endl;

    // Open the file and write the problem information
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "======================================================\n";
        file << "================  Problem information  ===============\n";
        file << "======================================================\n";
        file << "Problem dimension                   : " << problem_dimension << "\n";
        file << "Number of nodes                     : " << NoNs << "\n";
        file << "Number of bulk elements             : " << NoEs << "\n";
        file << "Number of DOFs                      : ";
        for (int dof : DOFs) {
            file << dof << " ";
        }
        file << "\n";
        std::cout << "element order                       : [ "<< element_order << "\n";
        file << "======================================================\n\n\n";
        file.close();  // Close the file
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}


//---------------------------------Prescribe----------------------------------
std::pair<std::vector<node>, std::vector<element>> prescribe(std::vector<node>& NL, std::vector<element>& EL, double LF) {
    int NoEs = EL.size();
    int NoNs = NL.size();
    int NPE = EL[0].NPE;
    int problem_dimension = NL[0].problem_dimension;

    // Loop through all nodes to prescribe displacements
    for (int i = 0; i < NoNs; ++i) {
        std::vector<int>& BC = NL[i].boundary_condition;
        std::vector<double>& BCval = NL[i].boundary_condition_value;
        std::vector<double>& X = NL[i].X_material_position;
        std::vector<double>& x = NL[i].x;

        for (int p = 0; p < problem_dimension; ++p) {
            if (BC[p] == 0) {
                x[p] = X[p] + LF * BCval[p];
            }
        }

        NL[i].x = x;
    }

    // Loop through all elements to update spatial coordinates
    for (int e = 0; e < NoEs; ++e) {
        Eigen::MatrixXd x(problem_dimension, NPE);
        const std::vector<int>& NdL = EL[e].NdL;

        for (int i = 0; i < NPE; ++i) {
            for (int p = 0; p < problem_dimension; ++p) {
                x(p, i) = NL[NdL[i]].x[p];
            }
        }

        EL[e].x = std::vector<double>(x.data(), x.data() + x.size());
    }

    return std::make_pair(NL, EL);
}

//---------------------------------Assemble-------------------------------------
std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>> Assemble(
        int problem_dimension,
        const std::vector<node>& NL,
        const std::vector<element>& EL,
        int DOFs)
{
    int NoEs = EL.size();
    int NoNs = NL.size();
    int NPE = EL[0].NPE;

    Eigen::VectorXd Rtot = Eigen::VectorXd::Zero(DOFs);
    std::vector<Eigen::Triplet<double>> tripletList;
    int sprC = 0; // sparse counter

    for (int e = 0; e < NoEs; ++e) {
        // Retrieve residual vector and stiffness matrix for current element
        auto [R, K] = EL[e].RK();
        const std::vector<int>& NdL = EL[e].NdL;

        int dim = problem_dimension;
        int first = 0;
        int last = dim - 1;

        // Update Rtot
        for (int i = 0; i < NPE; ++i) {
            const Eigen::VectorXd & BC = NL[NdL[i]].boundary_condition;
            const std::vector<int>& DOF = NL[NdL[i]].DOF;

            for (int p = first; p <= last; ++p) {
                if (BC[p] == 1) {
                    Rtot(DOF[p]) += R(i * dim + p - first);
                }
            }
        }

        // Update Ktot using triplets for sparse storage
        int dim_i = problem_dimension;
        int first_i = 0;
        int last_i = dim_i - 1;

        for (int i = 0; i < NPE; ++i) {
            const Eigen::VectorXd & BC_i = NL[NdL[i]].boundary_condition;
            const Eigen::VectorXd& DOF_i = NL[NdL[i]].DOF;

            for (int p = first_i; p <= last_i; ++p) {
                if (BC_i[p] == 1) {
                    int dim_j = problem_dimension;
                    int first_j = 0;
                    int last_j = dim_j - 1;

                    for (int j = 0; j < NPE; ++j) {
                        const std::vector<int>& BC_j = NL[NdL[j]].boundary_condition;
                        const std::vector<int>& DOF_j = NL[NdL[j]].DOF;

                        for (int q = first_j; q <= last_j; ++q) {
                            if (BC_j[q] == 1) {
                                double value = K(i * dim_i + p - first_i, j * dim_j + q - first_j);
                                tripletList.emplace_back(DOF_i[p], DOF_j[q], value);
                                sprC++;
                            }
                        }
                    }
                }
            }
        }
    }

    // Assemble the sparse stiffness matrix from triplets
    Eigen::SparseMatrix<double> Ktot(DOFs, DOFs);
    Ktot.setFromTriplets(tripletList.begin(), tripletList.end());

    return {Rtot, Ktot};
}

//----------------------------Assemble Gauss Points----------------------------------------
// Function definition for Assemble_GP
std::pair<Eigen::MatrixXd, Eigen::SparseMatrix<double>> Assemble_GP(
        int problem_dimension,
        const std::vector<node>& NL,
        const std::vector<element>& EL,
        int GP_DOFs)
{
    int NoEs = EL.size();
    int NoNs = NL.size();
    int NPE = EL[0].NPE;
    int NGP_val = NL[0].gauss_point_BC.size(); // Number of Gauss point values

    Eigen::MatrixXd Rtot_GP = Eigen::MatrixXd::Zero(GP_DOFs, NGP_val);
    std::vector<Eigen::Triplet<double>> tripletList;
    int sprC = 0; // sparse counter

    for (int e = 0; e < NoEs; ++e) {
        // Get residual and stiffness for current element at Gauss points
        auto [R, K] = EL[e].RK_gp(NGP_val);  // Assuming RK_GP returns a pair of std::vector types
        const std::vector<int>& NdL = EL[e].NdL;

        // Assemble Rtot_GP
        for (int s = 0; s < NGP_val; ++s) {
            for (int i = 0; i < NPE; ++i) {
                const std::vector<int>& BC = NL[NdL[i]].gauss_point_BC;
                const std::vector<int>& DOF = NL[NdL[i]].GP_DOF;
                int dim = 1;  // dim is 1 for Gauss point quantities

                for (int p = 0; p < dim; ++p) {
                    if (BC[p] == 1) {
                        Rtot_GP(DOF[p], s) += R[i * dim + p];  // Access R as a vector
                    }
                }
            }
        }

        // Assemble Ktot_GP using triplets for sparse storage
        for (int i = 0; i < NPE; ++i) {
            const std::vector<int>& BC_i = NL[NdL[i]].gauss_point_BC;
            const std::vector<int>& DOF_i = NL[NdL[i]].GP_DOF;
            int dim_i = 1;

            for (int p = 0; p < dim_i; ++p) {
                if (BC_i[p] == 1) {
                    for (int j = 0; j < NPE; ++j) {
                        const Eigen::VectorXdstd::vector<int>& BC_j = NL[NdL[j]].gauss_point_BC;
                        const Eigen::VectorXd& DOF_j = NL[NdL[j]].GP_DOF;
                        int dim_j = 1;

                        for (int q = 0; q < dim_j; ++q) {
                            if (BC_j[q] == 1) {
                                double value = K[i * dim_i + p][j * dim_j + q];  // Access K as a 2D vector
                                tripletList.emplace_back(DOF_i[p], DOF_j[q], value);
                                sprC++;
                            }
                        }
                    }
                }
            }
        }
    }

    // Assemble the sparse stiffness matrix from triplets
    Eigen::SparseMatrix<double> Ktot_GP(GP_DOFs, GP_DOFs);
    Ktot_GP.setFromTriplets(tripletList.begin(), tripletList.end());

    return {Rtot_GP, Ktot_GP};
}

// ----------------------------------Residual-----------------------------------
std::vector<double> Residual(
        int problem_dimension,
        const std::vector<node>& NL,
        const std::vector<element>& EL,
        int DOFs,
        double dt)
{
    int NoEs = EL.size();
    int NoNs = NL.size();
    int NPE = EL[0].NPE;

    // Initialize Rtot as a standard vector of zeros
    std::vector<double> Rtot(DOFs, 0.0);

    for (int e = 0; e < NoEs; ++e) {
        // Obtain the residual vector for the current element
        std::vector<double>  R = EL[e].Residual();
        const std::vector<int>& NdL = EL[e].NdL;

        int dim = problem_dimension; // Dimension based on the problem
        int first = 1;
        int last = dim;

        for (int i = 0; i < NPE; ++i) {
            const std::vector<int>& BC = NL[NdL[i]].boundary_condition;
            const std::vector<int>& DOF = NL[NdL[i]].DOF;

            for (int p = first; p <= last; ++p) {
                if (BC[p - 1] == 1) {
                    Rtot[DOF[p - 1]] += R[(i * dim) + (p - first)];
                }
            }
        }
    }

    return Rtot;
}


//---------------------------------Update----------------------------------------------------------------
void update(std::vector<node>& NL, std::vector<element>& EL, const std::vector<double>& dx) {
    int NoEs = EL.size();
    int NoNs = NL.size();
    int problem_dimension = NL[0].problem_dimension;
    int NPE = EL[1].NPE;
    // Update the nodes
    for (int i = 0; i < NoNs; i++) {
        node& node = NL[i];
        for (size_t p = 0; p < node.boundary_condition.size(); p++) {
            if (node.boundary_condition[p] == 1) {
                node.x[p] += dx[node.DOF[p]];
            }
        }
    }


    // Update the elements
    for (int e = 0; e < NoEs; e++) {
        element& elem = EL[e];
        elem.x.resize(problem_dimension, std::vector<double>(NPE));
        for (int i = 0; i < NPE; i++) {
            int nodeIndex = elem.NdL[i];
            elem.x[i] = NL[nodeIndex].x;
        }
    }
}
