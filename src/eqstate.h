#pragma once
#ifndef _EQSTATE_H
#define _EQSTATE_H
#include <cmath>

/* The following EOS by Jackett and McDougall.
* Details can be found in paper "Algorithms for Density, 
* Potential Temperature, Conservative Temperature, 
* and the Freezing Temperature of Seawater".
* We note that the pressure effect is not included at present.
*/
void EosByFeistel(double *dest, double th, double s);

/* The following EOS by UNESCO(1983).
* We note that the pressure effect is not included at present.
*/
void EosByUNESCO(double *dest, double T, double S);

// Calculate the density with the linear equation of state, see Tuomas(GMD,2017) for details.
void EosByLinear(double *dest, double T, double S, double rho0, double T0, \
	double S0, double alphaT, double betaS);

#endif