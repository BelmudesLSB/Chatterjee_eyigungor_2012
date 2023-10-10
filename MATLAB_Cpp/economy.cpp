#include <cmath>
#include <iostream>
#include <omp.h>
#include <mex.h>
#include "economy.hpp"
#include "math_functions.hpp"
#include "aux_functions.hpp"
#include "econ_functions.hpp"
#include "ggq_algorithm.hpp"


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

void CE_economy::solve(){

    double* Vd_0 = new double[Ny];
    double* E_mV = new double[Nb*Ny]; // Store the results of the ggq_algorithm.
    double* q_1 = new double[Nb*Ny];  // Store the results of the ggq_algorithm regarding the price, using q_0 as argument.
    double* q_0 = new double[Nb*Ny]; 
    double* w_0 = new double[Nb*Ny];
    //double* C_vector = new double[Nb]; // Create the vector for consumption:
    //double* W_vector = new double[Nb]; // Create the vector for continuation values:
    double dis_q = 1.00;
    double dis_w = 1.00;

    int iter = 0;

    while (dis_q > Tol || dis_w > Tol){ 

            copy_values(Vd_0, Vd, Ny);
            double aux = 0.00;
            int j = 0;

            #pragma omp parallel for private(aux, j)
            for (int i=0; i<Ny; i++){
                aux = 0.00;
                j = 0;
                for (j; j<Ny; j++){
                    aux += P[i*Ny+j] * Vd_0[j];
                }
                Vd[i] =  (utility(Ygrid[i] - phi(i) - M_bar, Tol, Gamma) + Beta * (1-Xi) * aux  + Xi * W[i*Nb+(Nb-1)]); 
            }
                

            clear_values(q_1, Nb*Ny);           // Clear the values of the vector.
            clear_values(E_mV, Nb*Ny);          // Clear the values of the vector.

            int x;
            #pragma omp parallel for private(j, x) collapse(2) schedule(dynamic) shared(E_mV, P, Q, Vd)
            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    double* C_vector = new double[Nb]; // Create the vector for consumption:
                    double* W_vector = new double[Nb]; 
                    clear_values(C_vector, Nb);  // Clear the values of the vector.
                    clear_values(W_vector, Nb);  // Clear the values of the vector.
                    for (int x=0; x<Nb; x++){
                        C_vector[x] = consumption_x(x, i, j, Q);
                        W_vector[x] = W[i*Nb + x];
                    }
                    E_mV[i*Nb+j] = ggq_topdown(Vd[i], Ny, Nb, C_vector, W_vector, Lambda, Z, R, 1-Gamma, Tol, -M_bar, M_bar, 0, Sigma_m, Tol, Max_iter, i, j, Q, P).EMV;
                    q_1[i*Nb+j] = ggq_topdown(Vd[i], Ny, Nb, C_vector, W_vector, Lambda, Z, R, 1-Gamma, Tol, -M_bar, M_bar, 0, Sigma_m, Tol, Max_iter, i, j, Q, P).Q1;
                    delete [] C_vector;
                    delete [] W_vector;         
                }
            }

            copy_values(q_0, Q, Ny*Nb); // Store the old bond price.
            copy_values(w_0, W, Ny*Nb); // Store the old continuation value.

            clear_values(W, Ny*Nb);       // Clear the values of the vector.
            clear_values(Q, Ny*Nb);       // Clear the values of the vector.

            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    for (int i_prime=0; i_prime<Ny; i_prime++){
                        W[i*Nb + j] += Beta * P[i*Ny + i_prime] * E_mV[i_prime*Nb + j];
                        Q[i*Nb + j] += P[i*Ny + i_prime] * q_1[i_prime*Nb + j];
                    }
                }
            }

            dis_q = 0;
            dis_w = 0;
            double aux_q = 0;
            double aux_w = 0;

            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    aux_q = std::abs(Q[i*Nb + j] - q_0[i*Nb + j]);
                    aux_w = std::abs(W[i*Nb + j] - w_0[i*Nb + j]);
                }
                if (aux_q > dis_q){
                    dis_q = aux_q;
                }
                if (aux_w > dis_w){
                    dis_w = aux_w;
                }
            }

            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    W[i*Nb + j] = Eta_w * w_0[i*Nb + j] + (1-Eta_w) * W[i*Nb + j];
                    Q[i*Nb + j] = Eta_q * q_0[i*Nb + j] + (1-Eta_q) * Q[i*Nb + j];
                    
                }
            }

            iter += 1;

            if (iter % 500 == 0){
                Eta_q = Eta_q * 0.1 + 0.9 * 0.995;
                Eta_w = Eta_w * 0.1 + 0.9 * 0.9;
            }

            if (iter % 5000 == 0){
                mexPrintf("Iteration: %d\n", iter);
                mexPrintf("Current Distance: Q %f and W %f\n", dis_q, dis_w);
            }

            if (iter> Max_iter){
                std::cout << "The algorithm did not converge" << std::endl;
                mexPrintf("The algorithm did not converge\n");
                break;
            }

            if (dis_q < Tol && dis_w < Tol){
                mexPrintf("Solution found\n");
                mexPrintf("The algorithm converged in %d iterations\n", iter);
                mexPrintf("Distances: Q %f and W %f\n", dis_q, dis_w);
                break;
            }
    }
}

