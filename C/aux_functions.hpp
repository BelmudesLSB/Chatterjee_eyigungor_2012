#ifndef aux_functions
#define aux_functions

void displayV(double *v, int M);

void displayQ(double *q, int ny, int nb);

bool convergence_check(double* f, double* g, int R, double error_tol);

void copy_values(double *f_to_fill, double *g_get_values, int R);

void clear_values(double *f_to_clear, int R);


#endif