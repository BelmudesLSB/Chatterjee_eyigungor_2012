#include <iostream>
#include <cmath>
#include <string>
#include <mex.h>

#include "economy.hpp"
#include "math_functions.hpp"
#include "aux_functions.hpp"
#include "ggq_algorithm.hpp"
#include "econ_functions.hpp"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    
    // Check for the proper number of arguments
    if (nrhs != 1 || nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "One input required.");
    }

    // Read the input parameters from MATLAB
    const mxArray* paramsStruct = prhs[0];

    // Preference parameters:
    //int nthreads = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "nthreads")));
    double beta = mxGetScalar(mxGetField(paramsStruct, 0, "beta"));
    double gamma =  mxGetScalar(mxGetField(paramsStruct, 0, "gamma"));

    // Income process parameters:
    double rho = mxGetScalar(mxGetField(paramsStruct, 0, "rho"));
    double sigma_e = mxGetScalar(mxGetField(paramsStruct, 0, "sigma_e"));
    double sigma_m = mxGetScalar(mxGetField(paramsStruct, 0, "sigma_m"));
    double t = mxGetScalar(mxGetField(paramsStruct, 0, "t"));
    double d_0 = mxGetScalar(mxGetField(paramsStruct, 0, "d_0"));
    double d_1 = mxGetScalar(mxGetField(paramsStruct, 0, "d_1"));

    // Grid Parameters:
    int ny = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "ny")));
    int nb = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "nb")));
    double b_min = mxGetScalar(mxGetField(paramsStruct, 0, "b_min"));
    double b_max = mxGetScalar(mxGetField(paramsStruct, 0, "b_max"));
    double m_bar = mxGetScalar(mxGetField(paramsStruct, 0, "m_bar"));


    // Convergence parameters:
    int max_iter = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "max_iter")));
    double tol = mxGetScalar(mxGetField(paramsStruct, 0, "tol"));

    // Term structure parameters:

    double r = mxGetScalar(mxGetField(paramsStruct, 0, "r"));
    double lambda = mxGetScalar(mxGetField(paramsStruct, 0, "lambda"));
    double z = mxGetScalar(mxGetField(paramsStruct, 0, "z"));
    double xi = mxGetScalar(mxGetField(paramsStruct, 0, "xi"));
    double eta_q = mxGetScalar(mxGetField(paramsStruct, 0, "eta_q"));
    double eta_w = mxGetScalar(mxGetField(paramsStruct, 0, "eta_w"));
    double eta_vd = mxGetScalar(mxGetField(paramsStruct, 0, "eta_vd"));

    //New! set the number of threads:
    //omp_set_num_threads(nthreads);

    /*
    Section 2: Allocate memory space to store the results of the model:
    */

   // Create vectors to store the output of the model:

    // Pointers to store the grids and transition matrices:
    mxArray* ygrid = mxCreateDoubleMatrix(ny, 1, mxREAL);
    mxArray* bgrid = mxCreateDoubleMatrix(nb, 1, mxREAL);
    mxArray* pgrid = mxCreateDoubleMatrix(ny*ny, 1, mxREAL);

    // Pointers to store the value functions:

    mxArray* Vd = mxCreateDoubleMatrix(ny, 1, mxREAL);
    mxArray* W = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);

    // Pointers to store the policy functions and prices:

    mxArray* Q = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);

    double* ygridPtr = mxGetPr(ygrid);
    double* bgridPtr = mxGetPr(bgrid);
    double* pgridPtr = mxGetPr(pgrid);
    double* Vd_Ptr = mxGetPr(Vd);
    double* W_ptr = mxGetPr(W);
    double* q_Ptr = mxGetPr(Q);


    /* 
    Section 3: Create the instance of the class and initialize the economy:
    */
    CE_economy c_economy(beta, gamma, rho, sigma_e, sigma_m, d_0, d_1, ny, nb, b_min, b_max, m_bar, t, max_iter, tol, r, lambda, z, xi, eta_q, eta_w, eta_vd, ygridPtr, bgridPtr, pgridPtr, Vd_Ptr, W_ptr, q_Ptr);
    mexPrintf("Initializing economy... \n");
    c_economy.initialize();

    /*
    Section 4: Solve the model:
    */  
    mexPrintf("Solving the model... \n");
    c_economy.solve();

    /*
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

            if (iter % 1000 == 0){
                std::cout << "Iteration: " << iter << std::endl;
                std::cout << "Distances| Q:" << dis_q << " and W:" << dis_w << std::endl;
                mexPrintf("Iteration: %d\n", iter);
                mexPrintf("Current Distance: Q %f and W %f\n", dis_q, dis_w);
            }

            if (iter>max_iter){
                std::cout << "The algorithm did not converge" << std::endl;
                mexPrintf("The algorithm did not converge\n");
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
                mexPrintf("Solution found\n");
                mexPrintf("The algorithm converged in %d iterations\n", iter);
                mexPrintf("Distances: Q %f and W %f\n", dis_q, dis_w);
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

    /*
    
    /*
    Section 5: Copy and Export to MATLAB
    */

    copy_values(c_economy.Ygrid, ygridPtr, ny);
    copy_values(c_economy.Bgrid, bgridPtr, nb);
    copy_values(c_economy.P, pgridPtr, ny*ny);
    copy_values(c_economy.Vd, Vd_Ptr, ny);
    copy_values(c_economy.W, W_ptr, ny*nb);
    copy_values(c_economy.Q, q_Ptr, ny*nb);

    const char* fieldNames[6]= {"ygrid", "bgrid", "pgrid", "Vd", "W", "Q"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fieldNames);
    mxSetField(plhs[0], 0, "ygrid", ygrid);
    mxSetField(plhs[0], 0, "bgrid", bgrid);
    mxSetField(plhs[0], 0, "pgrid", pgrid);
    mxSetField(plhs[0], 0, "Vd", Vd);
    mxSetField(plhs[0], 0, "W", W);
    mxSetField(plhs[0], 0, "Q", Q);
}
    

