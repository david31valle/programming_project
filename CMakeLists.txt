project(programming_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(programming_project
        main.cpp
        Eigen
        initialization/Initialize.hpp
        preprocess/mesh_1D.hpp
        preprocess/mesh_1D.cpp
        preprocess/mesh_2D.hpp
        preprocess/mesh_3D.hpp
        preprocess/mesh_2D.cpp
        preprocess/mesh_3D.cpp
        preprocess/mesh.cpp
        preprocess/mesh.hpp
        initialization/Initialize.cpp
        node/node.cpp

)

add_executable(testing
        tests.cpp
        utils/utils.hpp
        utils/utils.cpp
        element/element.hpp
        element/element.cpp
        Problem/Problem.cpp
        Problem/problem.hpp
)
