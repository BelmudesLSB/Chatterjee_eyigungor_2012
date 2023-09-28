#include <iostream>
#include <cmath>
#include "aux_functions.hpp"

void displayV(double *v, int M){

    std::cout << "............" << std::endl;
    for (int i=0;i<M;i++){
        std::cout << v[i] << "    ";
    }
    std::cout << "............" << std::endl;
}

void displayQ(double *q, int ny, int nb){
    for (int i=0;i<ny;i++){
        std::cout << "state i = " << i << ": ";
        for (int j=0;j<nb;j++){
            std::cout << q[i*nb+j] << " ";
        }
        std::cout << "|||||||||||||| " << std::endl;
    }
}

void copy_values(double *f_to_fill, double *g_get_values, int R){
    // This function copies the values on vector g on vector f.
    // Note, R is the length of the vector.
    for (int i=0;i<R;i++){
        f_to_fill[i] = g_get_values[i];
    }
}

void clear_values(double *f_to_clear, int R){
    for (int i=0;i<R;i++){
        f_to_clear[i] = 0.00;
    }
}