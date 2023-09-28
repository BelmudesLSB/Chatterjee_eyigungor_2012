#ifndef GGQ_H
#define GGQ_H

double ggq_topdown(double VD, int ny, int Nlist, double *clist, double *Wlist, double Lambda, double Z, double R, double ONE_MINUS_RRA, double C_LB, double GGQ_MLB, double GGQ_MUB, double GGQ_MMEAN, double GGQ_MSTD,
    double ERRTOL_BISECT, int MAXITERS_BISECT, int i, int j, double* Q_0, double * Q_1, double* P);
    
double bisect_zero(double a, double b, double t, int maxiter, double c1, double c2, double W1_minus_W2, double ONE_MINUS_RRA);
    
double m_root_fun(double m, double U1, double U2, double W1_minus_W2, double ONE_MINUS_RRA);

double CENDupdatefun(double m, double U, double W, double GGQ_MMEAN, double GGQ_MSTD, double ONE_MINUS_RRA);
    
double gauss_legendre_CENDupdate(double m1, double m2, double U, double W, double GGQ_MMEAN, double GGQ_MSTD, double ONE_MINUS_RRA);

#endif