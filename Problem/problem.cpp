//
// Created by David Valle on 05-Mar-25.
//

#include "problem.hpp"

problem::problem(int problem_dimension, const std::vector<node>& Node_List, const std::vector<element>& Element_List,
                 int domain_size, const std::string& boundary_condition, const std::string& deformation_type,
                 int element_order, double d, int steps, int max_iter, double tol, const std::string& gauss_points_values)
        : problem_dimension(problem_dimension), Node_List(Node_List), Element_List(Element_List),
          boundary_condition(boundary_condition), domain_size(domain_size), deformation_type(deformation_type),
          element_order(element_order), d(d), steps(steps), max_iter(max_iter), tol(tol),
          gauss_points_values(gauss_points_values) {

    filename = std::to_string(problem_dimension) + "D_Normal_" + std::to_string(Element_List[0].node_per_element.cols()) +  // Equivalent to `size(EL,2)` in MATLAB
                           "_EL=[" + std::to_string(element_order) + "]_.txt";

    initialize_F();
    Assign_BC();
    if(gauss_points_values=="On"){
        Assign_GP_DOFs();
    }
    problem_info();
    //Newton-Raphson method

    int counter = 1;
    double load_step=1/steps;
    assemble();
    while (load_factor <= 1.0 + 1e-8) { // Load step iteration
        std::cout << "\nLoad factor: " << load_factor << std::endl;

        // Apply load step
        prescribe();

        std::ofstream file(filename, std::ios::app);
        int error_counter = 1;
        bool isNotAccurate = true;

        while (isNotAccurate) { // Newton-Raphson loop
            assemble(); // Assemble global stiffness matrix and residual

            // Compute initial residual norm
            if (error_counter == 1) {
                double Norm_R0 = Rtot.norm();
                std::cout << "Residual Norm at Predictor: " << Norm_R0 << " , normalized: 1\n";
                if (file.is_open()) {
                    file << "Residual Norm at Predictor: " << Norm_R0 << " , normalized: 1\n";
                }
            }

            // Compute displacement increment and update nodal values
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(Ktot);
            Eigen::VectorXd dx = -solver.solve(Rtot);
            Eigen::VectorXd f_reaction = Kpu * dx;
            update(dx);

            // Recompute residual with updated values
            Residual(1.0);

            std::cout << "Residual Norm @ Increment " << error_counter << " at Corrector: " << Rtot.norm()
                      << " , normalized: " << Rtot.norm() << std::endl;
            if (file.is_open()) {
                file << "Residual Norm @ Increment " << error_counter << " at Corrector: " << Rtot.norm()
                     << " , normalized: " << Rtot.norm()  << "\n";
            }

            // Check for convergence
            if (Rtot.norm() < tol) {
                isNotAccurate = false;
            } else if (error_counter > max_iter || Rtot.norm() > 1e6) {
                isNotAccurate = false;
                std::cout << "Convergence is not obtained!" << std::endl;
            }

            error_counter++;
        } // End Newton-Raphson loop

        file.close();
        load_factor += load_step;
        counter++;

        /*if (gauss_points_values == "On") {
            //assemble_GP(); // Assemble at Gauss points
            Eigen::VectorXd dx_GP = Ktot_GP.colPivHouseholderQr().solve(Rtot_GP);
            //update_GP(dx_GP);
        }*/
    } // End


    std::cout<<"instance created"<<std::endl;


}

void problem::Assign_BC() {
    int NoNs = Node_List.size();  // Number of node
    double W = domain_size;  // X-direction ... width

    // Iterate over nodes
    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd X = Node_List[i].X_material_position;  // Node position

        if (std::abs(X(0) - 0) < tol || std::abs(X(0) - W) < tol) {
            Node_List[i].boundary_condition = Eigen::VectorXd::Zero(problem_dimension);
            Node_List[i].boundary_condition_value = F * X - X;
        }
    }

    // Call Assign_DOF_DBC (implemented in problem class)
    Assign_DOF_DBC();
}

void problem::Assign_DOF_DBC() {
    int NoNs = Node_List.size();  // Number of nodes
    int index = 0;
    CNL.clear();
    PNL.clear();

    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd& BC = Node_List[i].boundary_condition;  // Reference to avoid copying
        Eigen::VectorXd& DOF = Node_List[i].DOF;
        Eigen::VectorXd& DOC = Node_List[i].degree_constrained;
        Eigen::VectorXd& Numbering = Node_List[i].global_index;

        // Resize vectors to match problem dimension if needed
        if (DOF.size() != BC.size()) DOF = Eigen::VectorXd::Zero(BC.size());
        if (DOC.size() != BC.size()) DOC = Eigen::VectorXd::Zero(BC.size());
        if (Numbering.size() != BC.size()) Numbering = Eigen::VectorXd::Zero(BC.size());

        for (int p = 0; p < BC.size(); ++p) {
            index++;
            Numbering(p) = index;  // Assign global index

            if (BC(p) == 1) {
                DOFs++;
                DOF(p) = DOFs;  // Assign DOF number
            }
            else if (BC(p) == 0 && Node_List[i].boundary_condition_value.isZero()) {
                DOCs++;
                DOC(p) = DOCs;  // Assign DOC number
                CNL.push_back(index);
            }
            else {
                PNL.push_back(index);
            }
        }
    }
}

void problem::initialize_F() {
    // Identity matrix of size PD × PD
    F = Eigen::MatrixXd::Identity(problem_dimension, problem_dimension);

    if (deformation_type == "EXT") {
        F(0, 0) = 1 + d;  // Modify the first diagonal element
    }
    else if (deformation_type == "EXP") {
        F *= (1 + d);  // Scale the whole matrix
    }
    else if (deformation_type == "SHR") {
        if (problem_dimension > 1) {
            F(0, 1) = d;  // Modify shear component
        }
    }

    ;
}

void problem::Assign_GP_DOFs() {
    int NoNs = Node_List.size();  // Number of nodes
    int GP_DOFs = 0;

    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd& BC = Node_List[i].gauss_point_BC;  // Reference to avoid copying
        Eigen::VectorXd& DOF = Node_List[i].gauss_point_DOF;

        // Resize DOF if needed
        if (DOF.size() != BC.size()) {
            DOF = Eigen::VectorXd::Zero(BC.size());
        }

        // Check if any value in BC is 1 (equivalent to MATLAB's `ismember(1, BC)`)
        if ((BC.array() == 1.0).any()) {
            GP_DOFs++;
            DOF.array() = GP_DOFs;  // Assign the same GP_DOFs to all elements
        }
    }
}

void problem::problem_info() {
    int NoEs = Element_List.size();  // Number of elements
    int NoNs = Node_List.size();     // Number of nodes
      // Get filename

    // Print to console
    std::cout << "======================================================\n";
    std::cout << "================  Problem information  ===============\n";
    std::cout << "======================================================\n";
    std::cout << "Problem dimension                   : " << problem_dimension << "\n";
    std::cout << "Number of nodes                     : " << NoNs << "\n";
    std::cout << "Number of bulk elements             : " << NoEs << "\n";
    std::cout << "Number of DOFs                      : " << DOFs << "\n";
    std::cout << "Element order                       : [" << element_order << "]\n";
    std::cout << "======================================================\n";

    // Write to file
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "======================================================\n";
        file << "================  Problem information  ===============\n";
        file << "======================================================\n";
        file << "Problem dimension                   : " << problem_dimension << "\n";
        file << "Number of nodes                     : " << NoNs << "\n";
        file << "Number of bulk elements             : " << NoEs << "\n";
        file << "Number of DOFs                      : " << DOFs << "\n";
        file << "Element order                       : [" << element_order << "]\n";
        file << "======================================================\n\n\n";
        file.close();
        std::cout << "Problem information saved to: " << filename << "\n";
    } else {
        std::cerr << "Error: Unable to open file for writing: " << filename << "\n";
    }

}

void problem::prescribe() {
    int NoEs = Element_List.size();  // Number of elements
    int NoNs = Node_List.size();     // Number of nodes
    int PD = problem_dimension;      // Problem dimension

    // Update the new position of nodes based on BC and BCval
    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd& BC = Node_List[i].boundary_condition;
        Eigen::VectorXd& BCval = Node_List[i].boundary_condition_value;
        Eigen::MatrixXd X = Node_List[i].X_material_position;
        Eigen::MatrixXd x = Node_List[i].x_spatial_position;

        for (int p = 0; p < PD; ++p) {
            if (BC(p) == 0) {
                x(p) = X(p) + load_factor * BCval(p);
            }
        }

        Node_List[i].x_spatial_position = x;  // Update node position
    }

    // Copy updated node positions to the element spatial configuration
    for (int e = 0; e < NoEs; ++e) {
        int NPE = Element_List[e].node_per_element.cols();  // Nodes per element
        Eigen::MatrixXd x(PD, NPE);  // Temporary matrix for element positions

        Eigen::VectorXd NdL = Element_List[e].node_per_element.row(0).transpose();  // Extract node indices

        for (int i = 0; i < NPE; ++i) {
            x.col(i) = Node_List[NdL(i)].x_spatial_position;  // Assign node positions
        }

        Element_List[e].material_coordinate = x;  // Update element's spatial position
    }
}


void problem::assemble() {
    int NoEs = Element_List.size();  // Number of elements
    int NoNs = Node_List.size();     // Number of nodes
    int PD = problem_dimension;      // Problem dimension

    int DOFs_total = DOFs;           // Total degrees of freedom

    // Resize global force vector
    Rtot = Eigen::VectorXd::Zero(DOFs_total);

    // Sparse matrix assembly lists
    std::vector<Eigen::Triplet<double>> tripletList;

    // Generate all node indices
    std::vector<int> allIndices(NoNs * PD);
    std::iota(allIndices.begin(), allIndices.end(), 0);  // Fill with 0 to NoNs*PD-1

    // Identify unknown indices (excluding constrained and prescribed nodes)
    std::vector<int> unknown_indices;
    std::set_difference(allIndices.begin(), allIndices.end(), CNL.begin(), CNL.end(),
                        std::back_inserter(unknown_indices));
    std::set_difference(unknown_indices.begin(), unknown_indices.end(), PNL.begin(), PNL.end(),
                        std::back_inserter(unknown_indices));

    for (int e = 0; e < NoEs; ++e) {
        // Retrieve element stiffness matrix and residual vector
        Eigen::VectorXd R;
        Eigen::MatrixXd K;
        std::tie(R, K) = Element_List[e].compute_RK();  // Calls residual and stiffness matrix

        // Get element node list
        Eigen::VectorXi NdL = Element_List[e].node_per_element.row(0).transpose().cast<int>();

        // === ASSEMBLE GLOBAL RESIDUAL FORCE VECTOR ===
        for (int i = 0; i < NdL.size(); ++i) {
            Eigen::VectorXd& BC = Node_List[NdL(i)].boundary_condition;
            Eigen::VectorXd& DOF = Node_List[NdL(i)].DOF;

            for (int p = 0; p < PD; ++p) {
                if (BC(p) == 1) {
                    int dofIndex = DOF(p);
                    Rtot(dofIndex) += R((i * PD) + p);
                }
            }
        }

        // === ASSEMBLE GLOBAL STIFFNESS MATRIX ===
        for (int i = 0; i < NdL.size(); ++i) {
            Eigen::VectorXd& DOF_i = Node_List[NdL(i)].global_index;
            for (int p = 0; p < PD; ++p) {
                for (int j = 0; j < NdL.size(); ++j) {
                    Eigen::VectorXd& DOF_j = Node_List[NdL(j)].global_index;
                    for (int q = 0; q < PD; ++q) {
                        tripletList.emplace_back(DOF_i(p), DOF_j(q), K((i * PD) + p, (j * PD) + q));
                    }
                }
            }
        }
    }

// Build the sparse stiffness matrix using triplets
    Eigen::SparseMatrix<double> Total_stiffness_matrix(DOFs_total, DOFs_total);
    Total_stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    // Extract submatrices correctly using block operations
    Ktot = Total_stiffness_matrix.block(0, 0, unknown_indices.size(), unknown_indices.size());
    Kpu = Total_stiffness_matrix.block(PNL[0], unknown_indices[0], PNL.size(), unknown_indices.size());
    Kpp = Total_stiffness_matrix.block(PNL[0], PNL[0], PNL.size(), PNL.size());
}


void problem::update(const Eigen::VectorXd& dx) {
    int NoNs = Node_List.size();  // Number of nodes
    int NoEs = Element_List.size(); // Number of elements
    int PD = problem_dimension;   // Problem dimension

    // === Update nodal positions ===
    for (int i = 0; i < NoNs; ++i) {
        Eigen::VectorXd& BC = Node_List[i].boundary_condition;
        Eigen::VectorXd& DOF = Node_List[i].DOF;
        Eigen::VectorXd x = Node_List[i].x_spatial_position;

        for (int p = 0; p < BC.size(); ++p) {
            if (BC(p) == 1) {
                x(p) += dx(static_cast<int>(DOF(p)));
            }
        }

        Node_List[i].x_spatial_position = x;  // Store updated position
    }

    // === Update element spatial coordinates ===
    for (int e = 0; e < NoEs; ++e) {
        int NPE = Element_List[e].node_per_element.cols();  // Nodes per element
        Eigen::MatrixXd x(PD, NPE);  // Temporary matrix for element positions

        Eigen::VectorXi NdL = Element_List[e].node_per_element.row(0).transpose().cast<int>();  // Extract node indices

        for (int i = 0; i < NPE; ++i) {
            x.col(i) = Node_List[NdL(i)].x_spatial_position;  // Assign updated node positions
        }

        Element_List[e].spatial_coordinate = x;  // Store updated spatial positions
    }
}

void problem::Residual(double dt) {
    int NoEs = Element_List.size();  // Number of elements
    int NoNs = Node_List.size();     // Number of nodes
    int PD = problem_dimension;      // Problem dimension
    int NPE = Element_List[0].node_per_element.cols();  // Nodes per element

    // Initialize global residual vector
    Rtot = Eigen::VectorXd::Zero(DOFs);

    for (int e = 0; e < NoEs; ++e) {
        // Get element residual vector
        Eigen::VectorXd R_e = Element_List[e].Residual(dt);

        // Get the nodes of the element
        Eigen::VectorXi NdL = Element_List[e].node_per_element.row(0).transpose().cast<int>();

        // Assemble global residual vector
        for (int i = 0; i < NPE; ++i) {
            Eigen::VectorXd& BC = Node_List[NdL(i)].boundary_condition;
            Eigen::VectorXd& DOF = Node_List[NdL(i)].DOF;

            for (int p = 0; p < PD; ++p) {
                if (BC(p) == 1) {
                    int dofIndex = static_cast<int>(DOF(p));
                    Rtot(dofIndex) += R_e((i * PD) + p);
                }
            }
        }
    }
}