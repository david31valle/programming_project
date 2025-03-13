//
// Created by David Valle on 28-Jan-25.
//


#include <iostream>

#include "utils/utils.hpp"
#include "utils/NaturalElement.h"
int main(){
    std::cout<<"hello"<<std::endl;

    std::vector<std::vector<double>>gp=compute_gp(8,3);
    auto N_grad=compute_N_xi_gp(1,gp,3);

    return 0;


}