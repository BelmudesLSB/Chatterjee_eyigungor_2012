#include <iostream>
#include <cmath>
#include <string>
#include "economy.hpp"
#include "math_functions.hpp"
#include "aux_functions.hpp"
#include "ggq_algorithm.hpp"
#include "econ_functions.hpp"

int main(){

    /*
    Section 1: Parameters for this economy using the baseline calibration in Chatterjee and Eyigungor (2019):
    */

    // Preference parameters:
    double beta = 0.95402;      // discount factor.
    double gamma = 2;           // risk aversion.

    // Income process parameters:
    double rho = 0.948503;      // persistence of income.
    double sigma_e = 0.027092;  // standard deviation of income shocks.
    double sigma_m = 0.003;     // standard deviation of income shocks.
    double t = 3;               // number of standard deviations for Tauchen (1986).    
    double d_0 = -0.18819;      // Output loss at default parameter.
    double d_1 = 0.24558;       // Output loss at default parameter. 

    // Grid Parameters:
    int ny = 21;                // Y grid points. (?)
    int nb = 152;               // B grid points. (?)
    double b_min = -0.8;        // Minimum bond holdings.
    double b_max = 0.0;         // Maximum bond holdings.
    double m_bar = 0.006;       // Maximum value for m.
    
    // Convergence parameters:
    int max_iter = 50000;        // Maximum number of iterations.
    double tol = 1e-6;          // Tolerance for convergence.

    // Term structure parameters:
    double r = 0.01;            // Risk free interest rate.
    double lambda = 0.05;       // Reciprocal of average maturity.
    double z = 0.03;            // Coupon payments.
    double xi = 0.035;          // Probability of reentry.
    double eta_q = 0.995;       // weight on the old bond price.
    double eta_w = 0.9;         // weight on the new continuation value.
    double eta_vd = 0;          // weight on the new value of default.

    /*
    Section 2: Allocate memory space to store the results of the model:
    */

    // Pointers to store the grids and transition matrices:
    double* ygridPtr = new double[ny];      // Memory for the income grid.
    double* bgridPtr = new double[nb];      // Memory for the bond grid.
    double* pgridPtr = new double[ny*ny];   // Memory for the transition matrix.

    // Pointers to store the value functions:

    double* Vd_Ptr = new double[ny];          // Memory for the value function of default.    
    double* W_ptr = new double[ny*nb];        // Memory for the expected continuation value.

    // Pointers to store the policy functions and prices:

    double* q_Ptr = new double[ny*nb];      // Memory for the bond price.

    /* 
    Section 3: Create the instance of the class and initialize the economy:
    */

    CE_economy c_economy(beta, gamma, rho, sigma_e, sigma_m, d_0, d_1, ny, nb, b_min, b_max, m_bar, t, max_iter, tol, r, lambda, z, xi, eta_q, eta_w, eta_vd, ygridPtr, bgridPtr, pgridPtr, Vd_Ptr, W_ptr, q_Ptr);
    c_economy.initialize();

    /*
    Section 4: Solve the model:
    */  

    double* Vd_0 = new double[ny];
    double* E_mV = new double[nb*ny]; // Store the results of the ggq_algorithm.
    double* q_1 = new double[nb*ny];  // Store the results of the ggq_algorithm regarding the price, using q_0 as argument.
    double* q_0 = new double[nb*ny]; 
    double* w_0 = new double[nb*ny];
    double* C_vector = new double[nb]; // Create the vector for consumption:
    double* W_vector = new double[nb]; // Create the vector for continuation values:
    double dis_q = 1.00;
    double dis_w = 1.00;

    int iter = 0;

    while (dis_q > tol || dis_w > tol){ 

            copy_values(Vd_0, c_economy.Vd, ny);
            double aux = 0.00;
            for (int i=0; i<ny; i++){
                aux = 0.00;
                for (int j=0; j<ny; j++){
                    aux += c_economy.P[i*ny+j] * Vd_0[j];
                }
                c_economy.Vd[i] =  (utility(c_economy.Ygrid[i] - c_economy.phi(i) - m_bar, tol, gamma) + beta * (1-xi) * aux  + xi * c_economy.W[i*nb+(nb-1)]); 
            } 
            
            clear_values(q_1, nb*ny);         // Clear the values of the vector.
            clear_values(E_mV, nb*ny);        // Clear the values of the vector.

            for (int i=0; i<ny; i++){
                for (int j=0; j<nb; j++){
                    clear_values(C_vector, nb); // Clear the values of the vector.
                    clear_values(W_vector, nb); // Clear the values of the vector.

                    for (int x=0; x<nb; x++){
                        C_vector[x] = c_economy.consumption_x(x, i, j, c_economy.Q);
                        W_vector[x] = c_economy.W[i*nb + x];
                    }               
                    E_mV[i*nb+j] = ggq_topdown(c_economy.Vd[i], ny, nb, C_vector, W_vector, lambda, z, r, 1-gamma, tol, -m_bar, m_bar, 0, sigma_m, tol, max_iter, i, j, c_economy.Q, q_1, c_economy.P);         
                }
            }

            copy_values(q_0, c_economy.Q, ny*nb); // Store the old bond price.
            copy_values(w_0, c_economy.W, ny*nb); // Store the old continuation value.

            clear_values(c_economy.W, ny*nb);       // Clear the values of the vector.
            clear_values(c_economy.Q, ny*nb);       // Clear the values of the vector.

            for (int i=0; i<ny; i++){
                for (int j=0; j<nb; j++){
                    for (int i_prime=0; i_prime<ny; i_prime++){
                        c_economy.W[i*nb + j] += beta * c_economy.P[i*ny + i_prime] * E_mV[i_prime*nb + j];
                        c_economy.Q[i*nb + j] += c_economy.P[i*ny + i_prime] * q_1[i_prime*nb + j];
                    }
                }
            }

            dis_q = 0;
            dis_w = 0;
            double aux_q = 0;
            double aux_w = 0;

            for (int i=0; i<ny; i++){
                for (int j=0; j<nb; j++){
                    aux_q = std::abs(c_economy.Q[i*nb + j] - q_0[i*nb + j]);
                    aux_w = std::abs(c_economy.W[i*nb + j] - w_0[i*nb + j]);
                }
                if (aux_q > dis_q){
                    dis_q = aux_q;
                }
                if (aux_w > dis_w){
                    dis_w = aux_w;
                }
            }

            for (int i=0; i<ny; i++){
                for (int j=0; j<nb; j++){
                    c_economy.W[i*nb + j] = eta_w * w_0[i*nb + j] + (1-eta_w) * c_economy.W[i*nb + j];
                    c_economy.Q[i*nb + j] = eta_q * q_0[i*nb + j] + (1-eta_q) * c_economy.Q[i*nb + j];
                    
                }
            }

            iter += 1;

            if (iter % 250 == 0){
                std::cout << "Iteration: " << iter << std::endl;
                std::cout << "Distances| Q:" << dis_q << " and W:" << dis_w << std::endl;
            }

            if (iter>max_iter){
                std::cout << "The algorithm did not converge" << std::endl;
                break;
            }

            if (dis_q < tol && dis_w < tol){
                std::cout<< "The algorithm converged" << std::endl;
                std::cout << "The algorithm converged in " << iter << " iterations" << std::endl;
                std::cout << "Distances: " << dis_q << " and " << dis_w << std::endl;
                std::cout << "The bond price is: " << std::endl;
                displayQ(c_economy.Q, ny, nb);
                std::cout << "The continuation value is: " << std::endl;
                displayQ(c_economy.W, ny, nb);
                std::cout << "The value of default is: " << std::endl;
                displayV(c_economy.Vd, ny);
                break;
            }

        }
   

    delete [] C_vector;
    delete [] W_vector;
    delete [] E_mV;
    delete [] Vd_0;
    delete [] q_1;
    delete [] q_0;
    delete [] w_0;

    
    //displayQ(c_economy.W, ny, nb);
    //displayV(c_economy.Vd, ny);
    
}

