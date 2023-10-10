#include <cmath>
#include <cstdlib>
#include <iostream>

#include "math_functions.hpp"
#include "aux_functions.hpp"
#include "econ_functions.hpp"
#include "ggq_algorithm.hpp"

Result ggq_topdown(double VD, int ny, int Nlist, double *clist, double *wlist, double Lambda, double Z, double R, double ONE_MINUS_RRA, double C_LB, double GGQ_MLB, double GGQ_MUB, double GGQ_MMEAN, double GGQ_MSTD,
	double ERRTOL_BISECT, int MAXITERS_BISECT, int y, int b, double* Q_0, double* P){

	/* The arguments are:
		VD: value of default
		Nlist: number of points in the list.
		clist: list of consumption values.
		wlist: list of continuation values.
		ONE_MINUS_RRA: 1-RRA
		C_LB: lower bound on consumption
		GGQ_MLB: lower bound on m
		GGQ_MUB: upper bound on m
		GGQ_MMEAN: mean of m
		GGQ_MSTD: standard deviation of m
		ERRTOL_BISECT: error tolerance for bisection
		MAXITERS_BISECT: maximum number of iterations for bisection
        The output is:
        CENDupdate: the CEND update
        y : the index of the state in the ygrid (This will help update the price of bonds)
        b : the index of the state in the bgrid (This will help update the price of bonds)
        q_0 : the price of bonds in the previous iteration.
        q_1 : the price of bonds in the current iteration.
	*/

	// initialize at mub

	double Cnow, Wnow, Vnow, Ccandidate, Wcandidate, Vcandidate;
	int i, ichoice;
	double mnow = GGQ_MUB;	// m_1
	bool bContinue = false;
	double Q_1 = 0.00;

	// initiate Vnow at a really small number
	Vnow = wlist[0];
	for (i = 1; i < Nlist; i++)
	{
		if (wlist[i]<Vnow)
		{
			Vnow = wlist[i];
		}
	}
	Vnow += pow(C_LB,ONE_MINUS_RRA)/ONE_MINUS_RRA - 1000.0;

	// initialize at m_1
	for (i = 0; i < Nlist; i++){
		Ccandidate = clist[i];
		// Check feasibility at m_1:
		if (Ccandidate + mnow >= C_LB)
		// NOTE: q'>=qmin should only apply if b'-(1-m)*b>0 
		{
            Wcandidate = wlist[i];
			Vcandidate = pow(Ccandidate + mnow, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wcandidate; //mnow=GGQ_MLB
			// U_SCALING = 1/(1-RRA)
			if (Vcandidate > Vnow)
			{
				Cnow = Ccandidate;
				Wnow = Wcandidate;
				Vnow = Vcandidate;
				// record choices
				ichoice = i;
				bContinue = true;
			}
			else if (Vcandidate == Vnow)
			{
				if (Ccandidate < Cnow)//choose lowest c
				{
					Cnow = Ccandidate;
					Wnow = Wcandidate;
					Vnow = Vcandidate;
					// record choices
					ichoice = i;
					bContinue = true;
				}
			}
		}
	} // initialize

	if (!bContinue)
	{
		// no feasible choice; always default
		// CENDupdate = VD;
		// We dont need to integrate anything E_m[v(y,b,m)]=VD(y)
		Result result;
		result.EMV = VD;
		result.Q1 = 0;
		return result;
	}
	else
	{	// There exists at least one feasible choice at m_1

		int ichoice_prev = ichoice;	// x_1(y,b,m_1)
		double M1, Vnow_M1, mnext, Cnext, Wnext;
		double CENDupdate = 0.0;
		bContinue = true;
		double temp;

		while (bContinue){

			M1 = std::max(GGQ_MLB, C_LB - Cnow); // feasible m for current Cnow
			Vnow_M1 = pow(Cnow + M1, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wnow;
			mnext = M1;
			Cnext = Cnow; // initialize at Cnow for tie-breaks
			bContinue = false;

			// go through choices
			for (i = 0; i < Nlist; i++)
			{	
				if (i != ichoice_prev)
				{
					Ccandidate = clist[i];
					
					if (Ccandidate>Cnow) // only larger c is possibly a better choice, feasibility is automatically satisfied!
					{
						Wcandidate = wlist[i];
						// Evaluate at the lowest possible point! Using strict concavity.
						Vcandidate = pow(Ccandidate + M1, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wcandidate;
						if (Vcandidate>Vnow_M1)
						{
							// find bisection point, save a register-- store in Vcandidate
							Vcandidate = bisect_zero(M1, mnow, ERRTOL_BISECT, MAXITERS_BISECT, Cnow, Ccandidate, Wnow - Wcandidate , ONE_MINUS_RRA);
							if (Vcandidate > mnext)
							{
								// NOTE: if the difference is too small, the root from the bisection may be equal to MUB. This will produce an error.
								ichoice = i;
								mnext = Vcandidate;
								Cnext = Ccandidate;
								Wnext = Wcandidate;
								bContinue = true;
							}
							else if (Vcandidate == mnext && Ccandidate<Cnext) // choose lowest consumption
							{
								ichoice = i;
								mnext = Vcandidate;
								Cnext = Ccandidate;
								Wnext = Wcandidate;
								bContinue = true;
							}
						} // is a better choice
					} // possible better choice
				}
			} // note, it's possible for this loop to have no feasible choice (e.g. q=0)
			if (VD>=pow(Cnow + mnow, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wnow)
			{
				// always default
				temp = VD*(normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
				CENDupdate += temp;
				bContinue = false;
			}
			else if (VD<=pow(Cnow + mnext, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wnow)
			{
				// never default on the interval we are analysing.
				if (bContinue)
				// there was something better than c_now at some point! (bContinue is true)
				{   
					temp = gauss_legendre_CENDupdate(mnext, mnow, Cnow, Wnow, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA);
                    double prob_m = (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(mnext - GGQ_MMEAN, GGQ_MSTD))/(normcdf(GGQ_MUB - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
                    Q_1 += prob_m * (Lambda + (1-Lambda)*(Z + Q_0[y*Nlist+ichoice_prev])) * (1/(1+R));
					//std::cout << "2. prob_m = " << prob_m << std::endl;
					CENDupdate += temp;
					// update and continue algorithm
					
                    ichoice_prev = ichoice;
					mnow = mnext;
					Cnow = Cnext;
					Wnow = Wnext;
				}
				else
				// There is nothing better than cnow in the interval (bContinue is false !!) 
				{
					if (C_LB - Cnow > GGQ_MLB) // (c_lb > m_lb + c_i) And there is nothing better than cnow for the rest of the interval
					{
						// (Cnow,Wnow) applies until M1=C_LB-Cnow>MLB.
						temp = gauss_legendre_CENDupdate(M1, mnow, Cnow, Wnow, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA);
                        double prob_m = (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf(M1 - GGQ_MMEAN, GGQ_MSTD))/(normcdf(GGQ_MUB - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
                        Q_1 += prob_m * (Lambda + (1-Lambda)*(Z + Q_0[y*Nlist+ichoice_prev])) * (1/(1+R));
						////std::cout << "3. prob_m = " << prob_m << std::endl;
						CENDupdate += temp;

						// reinitialize algorithm at the maximal value at M1. why?
						
						bContinue = false; // Stop algorithm after this
						mnow = M1;
						Vnow = Vnow_M1;
						for (i = 0; i < Nlist; i++)
						{
							Ccandidate = clist[i];
							if (i != ichoice && Ccandidate + M1 > C_LB) // ichoice is the previous choice
							{
								Wcandidate = wlist[i];
								Vcandidate = pow(Ccandidate + M1, ONE_MINUS_RRA)/ONE_MINUS_RRA + Wcandidate;
								if (Vcandidate > Vnow)
								{
									Cnow = Ccandidate;
									Wnow = Wcandidate;
									Vnow = Vcandidate;
									// record choices
									ichoice_prev = i;
									bContinue = true;
								}
								else if (Vcandidate == Vnow)
								{
									if (Ccandidate < Cnow)
									{
										Cnow = Ccandidate;
										Wnow = Wcandidate;
										Vnow = Vcandidate;
										// record choices
										ichoice_prev = i;
										bContinue = true;
									}
								}
							}
						}
						
						if (!bContinue)
						{
							// all choices are not feasible, VD applies for m<=M1
							temp = VD*(normcdf(M1 - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
							CENDupdate += temp;
						}
					}
					else
					{ 
						// (Cnow,Wnow) applies until MLB. The policy is well defined over all the interval.
						temp = gauss_legendre_CENDupdate(GGQ_MLB, mnow, Cnow, Wnow, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA);
                        double prob_m = (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf( GGQ_MLB - GGQ_MMEAN, GGQ_MSTD))/(normcdf(GGQ_MUB - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
						Q_1 += prob_m * (Lambda + (1-Lambda)*(Z + Q_0[y*Nlist+ichoice_prev])) * (1/(1+R));
						//std::cout << "4. prob_m = " << prob_m << std::endl;
						CENDupdate += temp;
					}
				}
			}
			else
			{
				// intermediate default
				Vcandidate = pow(ONE_MINUS_RRA*(VD - Wnow),1/ONE_MINUS_RRA) - Cnow; // def threshold
				temp = gauss_legendre_CENDupdate(Vcandidate, mnow, Cnow, Wnow, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
					+ VD*(normcdf(Vcandidate - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
                double prob_m = (normcdf(mnow - GGQ_MMEAN, GGQ_MSTD) - normcdf( Vcandidate - GGQ_MMEAN, GGQ_MSTD))/(normcdf(GGQ_MUB - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));
                Q_1 += prob_m * (Lambda + (1-Lambda)*(Z + Q_0[y*Nlist+ichoice_prev])) * (1/(1+R));
				CENDupdate += temp;
				//std::cout << "1. prob_m = " << prob_m << std::endl;
				// Lucas:
				bContinue = false;
			}

		} // while loop

		CENDupdate = CENDupdate/(normcdf(GGQ_MUB - GGQ_MMEAN, GGQ_MSTD) - normcdf(GGQ_MLB - GGQ_MMEAN, GGQ_MSTD));

		Result result;
		result.EMV = CENDupdate;
		result.Q1 = Q_1;

		return result;
	}
}

double bisect_zero(double a, double b, double t, int maxiter, double c1, double c2, double W1_minus_W2,	double ONE_MINUS_RRA)
{
	double xlb = a, xub = b, xmid, fmid;
	double flb = m_root_fun(xlb, c1, c2, W1_minus_W2, ONE_MINUS_RRA);
	for (int i = 0; i < maxiter; i++)
	{
		xmid = 0.5*(xlb + xub);
		fmid = m_root_fun(xmid, c1, c2, W1_minus_W2, ONE_MINUS_RRA);
		if (std::signbit(flb) == std::signbit(fmid))
		{
			// same sign
			xlb = xmid;
		}
		else
		{
			// different signs
			xub = xmid;
		}
		if (xub - xlb <= t)
		{
			break;
		}
	}
	return 0.5*(xlb + xub);
}

double m_root_fun(double m, double U1, double U2, double W1_minus_W2, double ONE_MINUS_RRA)
{
	return pow(U1 + m, ONE_MINUS_RRA)/ONE_MINUS_RRA - pow(U2 + m, ONE_MINUS_RRA)/ONE_MINUS_RRA + W1_minus_W2;
}

double CENDupdatefun(double m, double U, double W, double GGQ_MMEAN, double GGQ_MSTD, double ONE_MINUS_RRA){
	return normpdf(m - GGQ_MMEAN, GGQ_MSTD)*(pow(U + m, ONE_MINUS_RRA)/ONE_MINUS_RRA+W);
}

double gauss_legendre_CENDupdate(double m1, double m2, double U, double W, double GGQ_MMEAN, double GGQ_MSTD, double ONE_MINUS_RRA)
{
	double alpha = 0.5*(m1 + m2);
	double beta = 0.5*(m2 - m1);

	return 	beta*(0.1460811336496904*CENDupdatefun(alpha + beta*0.0000000000000000, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1445244039899700*CENDupdatefun(alpha + beta*-0.1455618541608951, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1445244039899700*CENDupdatefun(alpha + beta*	0.1455618541608951, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1398873947910731*CENDupdatefun(alpha + beta* -0.2880213168024011, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1398873947910731*CENDupdatefun(alpha + beta*	0.2880213168024011, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1322689386333375*CENDupdatefun(alpha + beta* -0.4243421202074388, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1322689386333375*CENDupdatefun(alpha + beta*	0.4243421202074388, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1218314160537285*CENDupdatefun(alpha + beta* -0.5516188358872198, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1218314160537285*CENDupdatefun(alpha + beta*	0.5516188358872198, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1087972991671484*CENDupdatefun(alpha + beta* -0.6671388041974123, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.1087972991671484*CENDupdatefun(alpha + beta*	0.6671388041974123, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0934444234560339*CENDupdatefun(alpha + beta* -0.7684399634756779, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0934444234560339*CENDupdatefun(alpha + beta*	0.7684399634756779, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0761001136283793*CENDupdatefun(alpha + beta* -0.8533633645833173, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0761001136283793*CENDupdatefun(alpha + beta*	0.8533633645833173, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0571344254268572*CENDupdatefun(alpha + beta* -0.9200993341504008, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0571344254268572*CENDupdatefun(alpha + beta*	0.9200993341504008, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0369537897708525*CENDupdatefun(alpha + beta* -0.9672268385663063, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0369537897708525*CENDupdatefun(alpha + beta*	0.9672268385663063, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0160172282577743*CENDupdatefun(alpha + beta* -0.9937521706203895, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA)
		+	0.0160172282577743*CENDupdatefun(alpha + beta*	0.9937521706203895, U, W, GGQ_MMEAN, GGQ_MSTD, ONE_MINUS_RRA));
}