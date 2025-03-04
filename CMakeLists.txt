cmake_minimum_required(VERSION 3.26)
project(programming_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(programming_project
        main.cpp
        Eigen
        preprocess/mesh_1D.hpp
        preprocess/mesh_1D.cpp
        preprocess/mesh_2D.hpp
        preprocess/mesh_2D.cpp
        node/node.hpp
        node/node.cpp
        initialization/Initialize.hpp
        element/element.hpp
        element/element.cpp
        utils/utils.hpp
        utils/utils.cpp
        utils/utils.cpp
        Problem/Problem.cpp
        Problem/Problem.cpp
        Problem/problem.hpp)

add_executable(testing
        tests.cpp
        utils/utils.hpp
        utils/utils.cpp
        element/element.hpp
        element/element.cpp
        Problem/Problem.cpp
        Problem/problem.hpp)