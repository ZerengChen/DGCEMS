#include "eqstate.h"

/* The following EOS by Jackett and McDougall.
* Details can be found in paper "Algorithms for Density, 
* Potential Temperature, Conservative Temperature, 
* and the Freezing Temperature of Seawater".
* We note that the pressure effect is not included at present.
*/
void EosByFeistel(double *dest, double th, double s) {

	double anum, aden;

	double th2 = th*th, sqrts = sqrt(s);

	anum = 9.9984085444849347*pow(10,2.0) + \
		th*(7.3471625860981584*pow(10,0.0) + \
			th*(-5.3211231792841769*pow(10,-2.0) + \
				th*3.6492439109814549*pow(10,-4.0))) + \
		s*(2.5880571023991390*pow(10,0.0) - \
			th*6.7168282786692355*pow(10,-3.0) + \
			s*1.9203202055760151*pow(10,-3.0));

	aden = 1.0000000000000000*pow(10,0.0) + \
		th*(7.2815210113327091*pow(10,-3.0) + \
			th*(-4.4787265461983921*pow(10,-5.0) + \
				th*(3.3851002965802430*pow(10,-7.0) + \
					th*1.3651202389758572*pow(10,-10.0)))) + \
		s*(1.7632126669040377*pow(10,-3.0) - \
			th*(8.8066583251206474*pow(10,-6.0) + \
				th2*1.8832689434804897*pow(10,-10.0)) + \
			sqrts*(5.7463776745432097*pow(10,-6.0) + \
				th2*1.4716275472242334*pow(10,-9.0)));

	(*dest) = anum / aden;
}

/* The following EOS by UNESCO(1983).
* We note that the pressure effect is not included at present.
*/
void EosByUNESCO(double *dest, double T, double S) {
	double T2 = T*T;
	double T3 = T*T2;
	double T4 = T2*T2;
	double T5 = T*T4;
	double S15 = pow(S,1.5);
	double S2 = S*S;
	double S3 = S*S2;
	double x;
	x = 999.842594 + 6.793952*pow(10,-2.0)*T - 9.09529*pow(10,-3.0)*T2 + 1.001685*pow(10,-4.0)*T3;
	x = x - 1.120083*pow(10,-6.0)*T4 + 6.536332*pow(10,-9.0)*T5;
	x = x + S*(0.824493 - 4.0899*pow(10,-3.0)*T + 7.6438*pow(10,-5.0)*T2 - 8.2467*pow(10,-7.0)*T3);
	x = x + S*5.3875*pow(10,-9.0)*T4;
	x = x + sqrt(S3)*(-5.72466*pow(10,-3.0) + 1.0227*pow(10,-4.0)*T - 1.6546*pow(10,-6.0)*T2);
	x = x + 4.8314*pow(10,-4.0)*S2;
	(*dest) = x;
}

// Calculate the density with the linear equation of state, see Tuomas(GMD,2017) for details.
void EosByLinear(double *dest, double T, double S, double rho0, double T0, \
	double S0, double alphaT, double betaS) 
{
	(*dest) = rho0 - alphaT*(T - T0) + betaS * (S - S0);
}