#include <iostream>
#include <cmath>
#include <string>
#include <mex.h>
#include <omp.h>

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

    // Parallelization parameters:
    
    int nthreads = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "nthreads")));
    omp_set_num_threads(nthreads);
    
    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            mexPrintf("Number of threads: %d\n", omp_get_num_threads());
        }
    }

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
    

