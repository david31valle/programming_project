cmake_minimum_required(VERSION 3.16)
project(programming_project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

add_executable(programming_project
        main.cpp
        Eigen

        preprocess/mesh_1D.hpp
        preprocess/mesh_1D.cpp
        preprocess/mesh_2D.hpp
        preprocess/mesh_3D.hpp
        preprocess/mesh_2D.cpp
        preprocess/mesh_3D.cpp
        preprocess/mesh.cpp
        preprocess/mesh.hpp
        initialization/Initialize.cpp
        initialization/Initialize.hpp
        node/node.hpp
        node/node.cpp
        element/element.hpp
        element/element.cpp
        utils/utils.hpp
        utils/utils.cpp
        Problem/problem.hpp
        Problem/problem.cpp


)

add_executable(
        tests.cpp
        utils/utils.hpp
        utils/utils.cpp
        tests.cpp
        utils/NaturalElement.h

)