#include "mesh.hpp"

std::pair<Eigen::MatrixXd , Eigen::MatrixXd> generate_mesh(int domain_size, int partition, int element_order, int problem_dimension){
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> generated_mesh;
    try {

        switch (problem_dimension) {
            case 1: {
                //------------ 1D Mesh ------------
                std::cout << "Generating 1D Mesh...\n";
                std::vector<int> element_orders = {element_order};

                // Create and generate the 1D mesh
                Mesh_1D mesh(domain_size, partition, element_orders);
                mesh.generateMesh();
                mesh.printMesh();
                auto elements=mesh.getElementLists()[0];


                int rows = elements.size();
                int cols = elements[0].size();  // Assume all inner vectors have the same size

                Eigen::MatrixXd element_list(rows, cols);

                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        element_list(i, j) = elements[i][j];
                    }
                }

                generated_mesh.first=mesh.getNodeList();
                generated_mesh.second=element_list;
                break;
            }

            case 2: {
                //------------ 2D Mesh ------------
                std::cout << "Generating 2D Mesh...\n";
                std::vector<int> element_orders = {element_order, element_order};

                // Create and generate the 2D mesh (Assuming Mesh2D class exists)
                Mesh_2D mesh(domain_size, partition, element_orders);
                mesh.generateMesh();
                //mesh.printMesh();
                auto elements=mesh.getElementLists()[0];
                auto y=mesh.getNodeList();

                int rows = elements.size();
                int cols = elements[0].size();  // Assume all inner vectors have the same size

                Eigen::MatrixXd element_list(rows, cols);

                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        element_list(i, j) = elements[i][j];
                    }
                }
                //std::cout << "Eigen Matrix:\n" << element_list << std::endl;
                generated_mesh.first=mesh.getNodeList();
                generated_mesh.second=element_list;
                break;
            }

            case 3: {
                std::cout << "Generating 3D Mesh...\n";
                std::vector<int> element_orders = {element_order, element_order, element_order};

                // Create and generate the 3D mesh (Assuming Mesh3D class exists)
                Mesh_3D mesh(domain_size, partition, element_orders);
                mesh.generateMesh();
                //mesh.printMesh();

                auto elements=mesh.getElementLists()[0];

                int rows = elements.size();
                int cols = elements[0].size();  // Assume all inner vectors have the same size

                Eigen::MatrixXd element_list(rows, cols);
                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        element_list(i, j) = elements[i][j];
                    }
                }
                //std::cout << "Eigen Matrix:\n" << element_list << std::endl;
                generated_mesh.first=mesh.getNodeList();
                generated_mesh.second=element_list;
                break;
            }

            default:
                throw std::invalid_argument("Invalid problem dimension! Use 1 (1D), 2 (2D), or 3 (3D).");
        }

        std::cout << "Mesh generation complete.\n";
        return generated_mesh;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return generated_mesh;
    }

}

