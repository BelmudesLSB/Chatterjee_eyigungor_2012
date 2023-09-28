#include <cmath>
#include <iostream>
#include "economy.hpp"
#include "math_functions.hpp"
#include "aux_functions.hpp"
#include "econ_functions.hpp"


CE_economy::CE_economy(double beta, double gamma, double rho, double sigma_e, double sigma_m, double d_0, double d_1, int ny, int nb, double b_min, double b_max, double m_bar, double t, int max_iter, double tol, double r, double lambda, double z, double xi, double eta_q, double eta_w, double eta_vd, double* ygridPtr, double* bgridPtr, double* pgridPtr, double* Vd_Ptr, double* W_ptr, double* q_Ptr){
    
    //// Set the parameters:

    Beta = beta;
    Gamma = gamma;
    Rho = rho;
    Sigma_e = sigma_e;
    Sigma_m = sigma_m;
    D_0 = d_0;
    D_1 = d_1;
    Ny = ny;
    Nb = nb;
    Bmin = b_min;
    Bmax = b_max;
    M_bar = m_bar;
    T = t;
    Max_iter = max_iter;
    Tol = tol;
    R = r;
    Lambda = lambda;
    Z = z;
    Xi = xi;
    Eta_q = eta_q;
    Eta_w = eta_w;
    Eta_vd = eta_vd;


    //// Set the pointers:

    Ygrid = ygridPtr;
    Bgrid = bgridPtr;
    P = pgridPtr;
    Vd = Vd_Ptr;
    W = W_ptr;
    Q = q_Ptr;   
}

int CE_economy::bond_grid(){
    //// This method creates the bond grid that is store in the object of the class Bgrid.

    double b_step = (Bmax - Bmin)/(Nb - 1);
    for (int i = 0; i < Nb; i++){
        Bgrid[i] = Bmin + i*b_step;
    }
    return EXIT_SUCCESS;
}

int CE_economy::income_grid(){
    //// This method creates the income grid that is store in the object of the class Ygrid.

    double sigma_y = sqrt(pow(Sigma_e,2)/(1-pow(Rho,2)));
    Ygrid[Ny-1] = T*sigma_y;
    Ygrid[0] = -T*sigma_y;
    double omega = (Ygrid[Ny-1]-Ygrid[0])/(Ny-1);

    for (int i=1; i<Ny-1; i++){
        Ygrid[i] = Ygrid[i-1] + omega;
    }
    return EXIT_SUCCESS;
}

int CE_economy::prob_grid(){
    //// This method creates the transition matrix for the income process that is store in the object of the class P.

    double sigma_y = sqrt(pow(Sigma_e,2)/(1-pow(Rho,2)));
    double omega = (2*T*sigma_y)/(Ny-1);

    for (int i=0; i<Ny; i++){
        for (int j=0; j<Ny; j++){
            // We need to treat endpoints separately:
            if (j==0 || j==Ny-1){
                if (j==0){
                    P[i*Ny+j] = normcdf((Ygrid[0]-Rho*Ygrid[i]+omega/2)/Sigma_e, 1);
                }
                else {
                    P[i*Ny+j] = 1-normcdf((Ygrid[Ny-1]-Rho*Ygrid[i]-omega/2)/Sigma_e,1);
                }
            // Probability distribution inside the y_grid:
            } else {
                P[i*Ny+j] = normcdf((Ygrid[j]-Rho*Ygrid[i]+omega/2)/Sigma_e, 1)-normcdf((Ygrid[j]-Rho*Ygrid[i]-omega/2)/Sigma_e, 1);
            }
        }
    }
    return EXIT_SUCCESS;
}

int CE_economy::guess_Vd(){

    // Guess that the default function as the value with no reentry:
    double* Vd_0 = new double[Ny];
    double* Vd_1 = new double[Ny];
    bool convergence = false;

    for (int i=0; i<Ny; i++){
        if (Ygrid[i] - phi(i) - M_bar <=0){
            //std::cout<< "Check calibration!!" << std::endl;
        }
        Vd_0[i] = utility(Ygrid[i] - phi(i) - M_bar, Tol, Gamma);
        //std::cout << "Vd_0[" << i << "] = " << Vd_0[i] << std::endl;
    }
    
    double distance = 1;
    int iter = 0;

    while (distance > Tol && iter < Max_iter){
        distance = 0;
        for (int i=0; i<Ny; i++){
            double aux = 0;
            for (int j=0; j<Ny; j++){
                aux += P[i*Ny+j] * Vd_0[j];
            }
            Vd_1[i] = utility(Ygrid[i] - phi(i) - M_bar, Tol, Gamma) + Beta * aux;
        }
        for(int i=0; i<Ny; i++){
            double temp = std::abs(Vd_1[i]-Vd_0[i]);
            if (temp > distance){
                distance = temp;
            }
        }

        copy_values(Vd_0, Vd_1, Ny);
        iter += 1;
    }

    copy_values(Vd, Vd_0, Ny);

    delete [] Vd_0;
    delete [] Vd_1;

    if (distance > Tol){
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }          
}

int CE_economy::guess_Q(){

    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nb; j++){
            Q[i*Nb+j] = (1/(1+R)) * ((double)j/((double)Nb-1));
        }
    }
    return EXIT_SUCCESS;
}

int CE_economy::guess_W(){
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nb; j++){
            W[i*Nb+j] = Beta * Vd[i];
        }
    }
    return EXIT_SUCCESS;
}

int CE_economy::initialize(){
    // This functions initializes grids, and verifies results.
    CE_economy::bond_grid();
    CE_economy::income_grid();
    CE_economy::prob_grid();
    CE_economy::guess_Q();

    for (int i=0; i<Ny; i++){
        Ygrid[i] = exp(Ygrid[i]);
    }

    // Check results:

    if (Bgrid[0] != Bmin){
        std::cout << "ERROR !!!!! :: Bmin is not the first element of the bond grid." << std::endl;
        return EXIT_FAILURE;
    }

    if (Bgrid[Nb-1] != Bmax){
        std::cout << "ERROR !!!!! :: Bmax is not the last element of the bond grid." << std::endl;
        return EXIT_FAILURE;
    }

    if (CE_economy::guess_Vd()==EXIT_FAILURE){
        std::cout << "ERROR !!!!! :: V_d did not converge. Check initialization." << std::endl;
        return EXIT_FAILURE;
    }

    CE_economy::guess_W();

    // Check wether the probability grid sums to one:

    double sum = 0;
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Ny; j++){
            sum += P[i*Ny+j];
        }
    }

    if (sum < Ny - Tol || sum > Ny + Tol){
        std::cout << "ERROR !!!!! :: Check probability grid" << std::endl;
        //std::cout << "The sum of the elements is: " << sum << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

double CE_economy::phi(int i){
    // Value of default:
    /*
    This function returns the value of phi given the income index.
    The argument is:
    i: income index.
    The output is:
    phi: value of output under default.
    */
    double aux = D_0 * Ygrid[i] + D_1 * pow(Ygrid[i],2);
    if (aux>0){
        return aux;
    } else {
        return 0.00;
    }
}

double CE_economy::consumption_x(int x, int i, int j, double* q){
    // Consumption as a function of income today, asset holdings today, asset holdings tomorrow and bond price:
    /* 
    This method computes consumption given the bond choice for tomorrow.
    The arguments are:
    x: choice of debt for the next period as an index at the grid.
    i: income as an index at the grid.
    j: current bond holdings as an index at the grid.
    q: bond price.   
    The output is:
    c: consumption value.
    */
    double c = Ygrid[i] + (Lambda + (1-Lambda)*Z) * Bgrid[j] - q[i*Nb+x] * (Bgrid[x] - (1-Lambda) * Bgrid[j]);
    return c;
}

