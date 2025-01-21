cmake_minimum_required(VERSION 3.26)
project(programming_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(programming_project
        main.cpp
        Eigen
        Preprocess/Mesh_1D.hpp
        Preprocess/Mesh_1D.cpp
        node/node.hpp
        node/node.cpp
        initialization/Initialize.hpp
        element/element.hpp
        element/element.cpp)