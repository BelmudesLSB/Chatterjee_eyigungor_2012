#include "math_functions.hpp"
#include <cmath>

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif 

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

double normpdf(double x, double sigma){
	return exp(-0.5*x*x/(sigma*sigma))*M_2_SQRTPI/(2*sigma*M_SQRT2);
}

double normcdf(double x, double sigma){
	return 0.5*(1+erf(x/(sigma*M_SQRT2)));
}