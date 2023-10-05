#include "econ_functions.hpp"

#include <cmath>
#include <iostream>

double utility(double c, double c_lb, double gamma){
    if (c>=(c_lb/2)){
        return pow(c, 1-gamma)/(1-gamma);
    }
    else{
        std::cout<< "Utility is not defined for c<c_lb" <<std::endl;
        return -1/c_lb - 100;
    }
}


