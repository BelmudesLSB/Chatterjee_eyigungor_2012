#ifndef economy
#define economy

#include <cmath>
#include <iostream>

// This header file, defines all the elements that the class is going to have as well as all the methods that we can apply to the class.

class CE_economy{

public:

    // Parameters for this economy:

    double Beta;    // discount factor.
    double Gamma;   // risk aversion.
    double Rho;     // persistence of income.
    double Sigma_e; // standard deviation of income shocks.
    double Sigma_m; // standard deviation of income shocks.
    double D_0;     // Output loss at default parameter.
    double D_1;     // Output loss at default parameter.
    int Ny;         // number of income states.
    int Nb;         // number of bond states.
    double Bmin;    // minimum bond holdings.
    double Bmax;    // maximum bond holdings.
    double M_bar;   // Maximum value for m.
    double T;       // number of standard deviations for Tauchen (1986).
    int Max_iter;   // Maximum number of iterations.
    double Tol;     // Tolerance for convergence.
    double R;       // gross interest rate.
    double Lambda;  // Reciprocal of average maturity.
    double Z;       // Coupon payments.
    double Xi;      // Probability of reentry.
    double Eta_q;     // Relaxation parameter for bond prices.
    double Eta_w;     // Relaxation parameter for continuation value.
    double Eta_vd;    // Relaxation parameter for value of default.
    
    // Pointers to operate on:

    double *Ygrid;  // pointer for grid of income.
    double *Bgrid;  // pointer for grid of bonds.
    double *P;      // pointer for transition matrix.
    double *Vd;     // pointer for value function of default.
    double *W;      // pointer for the expected continuation value.
    double *Q;      // pointer for bond price.

    // Constructor:
    CE_economy(double beta, double gamma, double rho, double sigma_e, double sigma_m, double d_0, double d_1, int ny, int nb, double b_min, double b_max, double m_bar, double t, int max_iter, double tol, double r, double lambda, double z, double xi, double eta_q, double eta_w, double eta_vd, double* ygridPtr, double* bgridPtr, double* pgridPtr, double* Vd_Ptr, double* W_ptr, double* q_Ptr);
    
    //// Methods associated with this class:

    // Grids and value for default:
    int bond_grid();
    int income_grid();
    int prob_grid();
    int guess_Vd();
    int guess_W();
    int guess_Q();
    int initialize();
    void solve();

    // Value of default:
    double phi(int i);

    // Consumption: 

    double consumption_x(int x, int i, int j, double* q);



    
    
};

#endif
