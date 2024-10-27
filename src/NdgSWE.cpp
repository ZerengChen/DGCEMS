#include "NdgSWE.h"

/*This function is used to impose the boundary condition for the pure hydrulic problem*/
void ImposeBoundaryCondition(double *gra, NdgEdgeType type, double *nx, double *ny, double *fm, double *fp, \
	double *zM, double *zP, double *fext, int Nfp, int Nvar, int Ne){
	/*All other fields except H, Hu, Hv are set to be the inner value at the boundaries unless they are prescribed out of the program */
	// assign the local node values
	double *huM = fm, *hvM = fm + Nfp*Ne, *hM = fm + 2 * Nfp*Ne;
	double *huP = fp, *hvP = fp + Nfp*Ne, *hP = fp + 2 * Nfp*Ne;
	double *huE = fext, *hvE = fext + Nfp*Ne, *hE = fext + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		zP[i] = zM[i];
	}

	// get next node values
	if (type == NdgEdgeInner) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huM[i];
			hvP[i] = hvM[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeZeroGrad) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huM[i];
			hvP[i] = hvM[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeClamped) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hE[i];
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fext[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeClampedDepth) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huM[i] / hE[i];
			const double vE = hvM[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];   // outward normal flux
			const double uvE = -uE * ny[i] + vE * nx[i];  // tangential flux
			const double un = unE - sqrt(-*gra / zP[i]) * (hE[i] - hM[i]);
			const double uv = uvE;
			hP[i] = hE[i];
			huP[i] = (un * nx[i] - uv * ny[i]) * hM[i];
			hvP[i] = (un * ny[i] + uv * nx[i]) * hM[i];
			for (int n = 3; n < Nvar; n++) {
				fp[n*Nfp*Ne + i] = fext[n*Nfp*Ne + i];
		    }
		}
	}
	else if (type == NdgEdgeClampedVel) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fext[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeSlipWall) {
		for (int i = 0; i < Nfp; i++){
			const double qxM = huM[i];
			const double qyM = hvM[i];
			double qnM = qxM * nx[i] + qyM * ny[i];   // outward normal flux
			double qvM = -qxM * ny[i] + qyM * nx[i];  // outward tangential flux
			// adjacent value
			hP[i] = hM[i];
			huP[i] = (-qnM) * nx[i] - qvM * ny[i];
			hvP[i] = (-qnM) * ny[i] + qvM * nx[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonSlipWall) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = -huM[i];
			hvP[i] = -hvM[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
#ifdef _BAROCLINIC
			//for Plume--------
			if (Nvar = 5) {
				fp[3 * Nfp*Ne + i] = hP[i] * 20.0;
			}
			//-----------------
#endif
		}
	}
	else if (type == NdgEdgeFlather) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];   // outward normal flux
			const double uvE = -uE * ny[i] + vE * nx[i];  // tangential flux
			const double un = unE - sqrt(-*gra / zP[i]) * (hE[i] - hM[i]);
			const double uv = uvE;
			hP[i] = hM[i];
			huP[i] = (un * nx[i] - uv * ny[i]) * hM[i];
			hvP[i] = (un * ny[i] + uv * nx[i]) * hM[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonLinearFlather) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];  // outward normal flux
			const double un = huM[i] / hM[i] * nx[i] + hvM[i] / hM[i] * ny[i];
			const double temp = 0.5 * (un - unE) + sqrt(*gra * hE[i]);
			hP[i] = temp * temp / *gra;
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonLinearFlatherFlow) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];   // outward normal flux
			const double uvE = -uE * ny[i] + vE * nx[i];  // tangential flux
			const double un = unE - 2 * sqrt(*gra * hE[i]) + 2 * sqrt(*gra * hM[i]);
			const double uv = uvE;
			hP[i] = hM[i];
			huP[i] = (un * nx[i] - uv * ny[i]) * hM[i];
			hvP[i] = (un * ny[i] + uv * nx[i]) * hM[i];
			//   } else if (type == NdgEdgeNonReflectingFlux) {
			//     const double unM = surf->huM / surf->hM * nx + surf->hvM / surf->hM *
			//     ny; const double RLP = unM + 2 * sqrt(gra * surf->hM); const double
			//     hs = fext->h[idM]; surf->hP = hs; const double RRM = RLP - 4 *
			//     sqrt(gra * hs); const double un = 0.5 * (RRM + RLP); const double uvM
			//     = -surf->huM / surf->hM * ny + surf->hvM / surf->hM * nx; surf->huP =
			//     (un * nx - uvM * ny) * surf->hM; surf->hvP = (un * ny + uvM * nx) *
			//     surf->hM;
			for (int n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fext[n*Nfp*Ne + i];
			}
		}
	}
	else {
		printf("Matlab:%s:Unknown boundary type: %d\n", __FILE__, type);
	}
	return;
}

void EvaluateHydroStaticReconstructValue(double hmin, double *fm, double *fp, double *zM, double *zP, int Nfp, int Nvar, int Ne)
{
	/*Passive transport substances, such as temperature and salt are not reconstructed*/
	double um, vm, up, vp, etaM, etaP, zstar;
	double *huM = fm, *huP = fp;
	double *hvM = fm + Nfp*Ne, *hvP = fp + Nfp*Ne;
	double *hM = fm + 2 * Nfp*Ne, *hP = fp + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		zstar = fmax(zM[i], zP[i]);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, huM + i, &um);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, hvM + i, &vm);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, huP + i, &up);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, hvP + i, &vp);
		etaM = hM[i] + zM[i];
		etaP = hP[i] + zP[i];
		zstar = fmin(etaM, zstar);
		huM[i] = hM[i] * um;
		hvM[i] = hM[i] * vm;
		huP[i] = hP[i] * up;
		hvP[i] = hP[i] * vp;
		for (int n = 3; n < Nvar; n++){
			double variable;
			/*Passive transport substances, such as temperature and salt are reconstructed following the same strategy for hu and hv*/
			EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, fm + n*Nfp*Ne + i, &variable);
			fm[n*Nfp*Ne + i] = hM[i] * variable;
			EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, fp + n*Nfp*Ne + i, &variable);
			fp[n*Nfp*Ne + i] = hP[i] * variable;
		}
		hM[i] = etaM - zstar;
		hP[i] = fmax(0, etaP - zstar) - fmax(0, zP[i] - zstar);
		zM[i] = zstar;
		zP[i] = zstar;
	}
}
/*
void EvaluateFlowRateByDeptheThreshold(double hmin, double *h, double *hu, double *hv, double *um, double *vm)
{
	if (*h > hmin) {
		//     const double sqrt2 = 1.414213562373095;
		//     double h4 = pow(h, 4);
		//     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
		//     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
		*um = *hu / *h;
		*vm = *hv / *h;
	}
	else {
		*um = 0.0;
		*vm = 0.0;
	}
}
*/

/*This function is used to calcualte the numerical flux in the primitive continuity equation(PCE) */
void GetPCENumericalFluxTerm_HLLC_LAI(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Ne){
	/*HU ranks first, HV second and H third*/
	/*The HLLC numerical flux presented in LAI(2012, Modeling one- and two-dimensional shallow water flows with discontinuous Galerkin method) is adopted*/
	double HR, HL, UR, UL;
	double SL, SM, SR;
	double USTAR, HSTAR;
	double *hum = fm, *hvm = fm + Nfp*Ne, *hm = fm + 2 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *hp = fp + 2 * Nfp*Ne;
	/*H, HUN, HVT*/
	double *QL = (double *)malloc(3*sizeof(double)), *QR = (double *)malloc(3*sizeof(double));
	/*H, UN, VT*/
	double *VARIABLEL = (double *)malloc(3 * sizeof(double)), *VARIABLER = (double *)malloc(3 * sizeof(double));
	double FL , FR ;
	double FSTARL, FSTARR ;

	double QSTARL, QSTARR;
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);

		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < 3; n++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		/*Here, we use the critical value to categorize different condition, I wonder whether this value can be
		directly set to zero, if this value is set to zero, then they calculation in the further part can be looked
		to be very smooth*/
		if (!(HL < Hcrit && HR < Hcrit))
		{

			/*Compute FL AND FR*/
			FL = HL*UL;

			FR = HR*UR;

			if ((HL>Hcrit) & (HR > Hcrit)){
				USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				HSTAR = 1.0 / (*gra)*pow(0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR), 2);
				SL = fmin(UL - sqrt(*gra * HL), USTAR - sqrt(*gra*HSTAR));
				SR = fmax(UR + sqrt(*gra * HR), USTAR + sqrt(*gra*HSTAR));
				SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
				/*Compute QSTARL AND QSTARR*/
				QSTARL = HL*((SL - UL) / (SL - SM)) * 1;
				QSTARR = HR*((SR - UR) / (SR - SM)) * 1;

				FSTARL = FL + SL * (QSTARL - QL[0]);
				FSTARR = FR + SR * (QSTARR - QR[0]);

				/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
				if (SL >= 0){
					dest[i] = FL;
					SPY = 1;
				}
				else if (SM >= 0 && SL < 0){
					dest[i] = FSTARL;
					SPY = 1;
				}
				else if (SM < 0 && SR > 0){
					dest[i] = FSTARR;
					SPY = 1;
				}
				else if (0 >= SR){
					dest[i] = FR;
					SPY = 1;
				}

			}
			else if (HL>Hcrit && HR <= Hcrit){
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
				//是否将SM的计算移到这里来
				//				SM = SR;
				if (SL >= 0){
					dest[i] = FL;
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					dest[i] = (SR * FL - SL * FR + SL * SR * (HR - HL)) / (SR - SL);
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = FR;
					SPY = 1;
				}
			}
			else if (HL <= Hcrit && HR > Hcrit){
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
				//				SM = SL;
				if (SL >= 0){
					dest[i] = FL;
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					dest[i] = (SR * FL - SL * FR + SL * SR * (HR - HL)) / (SR - SL);
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = FR;
					SPY = 1;
				}
			}
			//SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
			if (SPY == 0){
				printf("Error in PCE HLLC flux!");
			}
		}
	}
	free(QL), free(QR);
	free(VARIABLEL), free(VARIABLER);
}


void GetPCENumericalFluxTerm_HLL(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Ne){
	/*HU ranks first, HV second and H third*/
	/*The HLLC numerical flux presented in LU(2020, Computer Methods in Applied Mechanics and Engineering) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SR;
	//double HSTAR, USTAR;
	double *hum = fm, *hvm = fm + Nfp*Ne, *hm = fm + 2 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *hp = fp + 2 * Nfp*Ne;

	double FL, FR;
	double *QL = (double *)malloc(3 * sizeof(double)), *QR = (double *)malloc(3 * sizeof(double));
	double *VARIABLEL = (double *)malloc(3 * sizeof(double)), *VARIABLER = (double *)malloc(3 * sizeof(double));
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);


		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < 3; n++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];

		if (!(HL < Hcrit && HR < Hcrit))
		{
			if (HL > Hcrit && HR>Hcrit){
				double USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				double CSTAR = 0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR);
				SL = fmin(UL - sqrt(*gra*HL), USTAR - CSTAR);
				SR = fmax(UR + sqrt(*gra*HR), USTAR + CSTAR);
			}
			else if (HL>Hcrit && HR <= Hcrit){
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
			}
			else if (HL <= Hcrit && HR > Hcrit){
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
			}
			else{
				SL = 0.0;
				SR = 0.0;
			}


			/*Compute FL AND FR*/
			FL = HL*UL;

			FR = HR*UR;

			/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
			if (SL >= 0 && SR > 0){
				dest[i] = FL;
				SPY = 1;
			}
			else if (SL < 0 && SR > 0){
				dest[i] = (SR * FL - SL * FR + SL * SR * (HR - HL)) / (SR - SL);
				SPY = 1;
			}
			else if (SL < 0 && SR <= 0){
				dest[i] = FR;
				SPY = 1;
			}
			else if ((fabs(SL) < EPS) & (fabs(SR) < EPS)){
				dest[i] = FL;
				SPY = 1;
			}
			if (SPY == 0){
				printf("Error occured when calculating PCE HLL flux! check please!");
			}
		}
	}
	free(QL), free(QR);
	free(VARIABLEL), free(VARIABLER);
}
/*
* Purpose: This function is used to evaluate the local Riemann problem at the element surface to impose the boundary condition for slope limiter
*
* Input:
*      double[Nfp x Ne x (Nvar+1)] fm value at the local surface, water depth is included in this variable
* 	   double[Nfp x Ne x (Nvar+1)] fp value at the adjacent surface, water depth is included in this variable
*      double[Nfp x Ne]  nx the direction vector in x direction
* 	   double[Nfp x Ne]  ny the direction vector in y direction
* 	   double  gra the accelaration due to gravity
*      double  Hcrit the critical water depth used to identify whether the studied point should be treated as a wet or dry point
*      double[2 x Ne] BEFToE the topological relation at the boundary
*      int Ne number of the boundary face
*      int Nfp number of nodes on the face
*      int Nvar number of variable
*      int Nh approximation order in horizontal direction
*      int Nz approximation order in vertical direction
* Output:
* 		double[4 x Ne x Nvar] dest the Riemann variable value at the interface, it contains momentum components in normal direction and 
*                             tangential direction and contains concentrantion (multiplied by water depth) of passive tranport substance.
*                             The returned value is of size 4 in the first dimension, because only 4 vertex are included i
*/

void EvaluateVerticalFaceRiemannProblem(double *dest, double *Tempfm, double *Tempfp, double *Tempnx, double *Tempny, \
	double *gra, double Hcrit, int Ne, int Nfp, int Nvar, int Nh, int Nz){
	/*The HLLC numerical flux presented in LAI(2012, Modeling one- and two-dimensional shallow water flows with discontinuous Galerkin method) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SM, SR;
	double USTAR, HSTAR;
	/*This part is used to extract the vertex value*/
	double *fm = (double *)malloc(4 * (Nvar + 1)*sizeof(double)), *fp = (double *)malloc(4 * (Nvar + 1)*sizeof(double));
	double *hum = fm, *hvm = fm + 4, *hup = fp, *hvp = fp + 4, *hm = fm + 2 * 4, *hp = fp + 2 * 4;
	double *nx = (double *)malloc(4 * sizeof(double)), *ny = (double *)malloc(4 * sizeof(double));
	for (int i = 0; i < 2; i++){
		nx[i] = Tempnx[i*Nh], nx[i + 2] = Tempnx[Nz*(Nh + 1) + i*Nh];
		ny[i] = Tempny[i*Nh], ny[i + 2] = Tempny[Nz*(Nh + 1) + i*Nh];
		for (int n = 0; n < Nvar + 1; n++){
			fm[n * 4 + i] = Tempfm[n*Nfp*Ne + i*Nh];
			fp[n * 4 + i] = Tempfp[n*Nfp*Ne + i*Nh];
			fm[n * 4 + i + 2] = Tempfm[n*Nfp*Ne + Nz*(Nh + 1) + i*Nh];
			fp[n * 4 + i + 2] = Tempfp[n*Nfp*Ne + Nz*(Nh + 1) + i*Nh];
		}
	}
	double *QL = (double *)malloc((Nvar + 1)*sizeof(double)), *QR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *VARIABLEL = (double *)malloc((Nvar + 1)*sizeof(double)), *VARIABLER = (double *)malloc((Nvar + 1)*sizeof(double));
	double *QSTARL = (double *)malloc((Nvar + 1)*sizeof(double)), *QSTARR = (double *)malloc((Nvar + 1)*sizeof(double));
	for (int i = 0; i < 4; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar + 1; n++)
		{
			QL[n] = *(fm + 4*n + i);
			QR[n] = *(fp + 4*n + i);
		}

		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < Nvar + 1; n++){
			/*In this part, we calculate the original variable. If the water depth is smaller than the critical value,
			we directly set the value of the original variable to be zero, and the original value is further used to calculate
			the flux term in left and right side. The question here is whether the water depth should be set to zero if it is
			smaller than the critical value, since it is used when calcualte the flux term. In the current version, we leave them
			unchanged*/
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];
		/*Here, we use the critical value to categorize different condition, I wonder whether this value can be
		directly set to zero, if this value is set to zero, then they calculation in the further part can be looked
		to be very smooth*/
		if (!(HL < Hcrit && HR < Hcrit))
		{
			if ((HL>Hcrit) & (HR > Hcrit)){
				USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				HSTAR = 1.0 / (*gra)*pow(0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR), 2);
				SL = fmin(UL - sqrt(*gra * HL), USTAR - sqrt(*gra*HSTAR));
				SR = fmax(UR + sqrt(*gra * HR), USTAR + sqrt(*gra*HSTAR));
				SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
				/*Compute QSTARL AND QSTARR*/
				QSTARL[0] = HL*((SL - UL) / (SL - SM)) * 1;
				QSTARL[1] = HL*((SL - UL) / (SL - SM)) * SM;
				QSTARR[0] = HR*((SR - UR) / (SR - SM)) * 1;
				QSTARR[1] = HR*((SR - UR) / (SR - SM)) * SM;
				for (int i = 2; i<Nvar + 1; i++){
					QSTARL[i] = HL*((SL - UL) / (SL - SM)) * VARIABLEL[i];
					QSTARR[i] = HR*((SR - UR) / (SR - SM)) * VARIABLER[i];
				}
				/*AND FINALLY THE HLLC FLUX (AFTER ROTATION)*/
				if (SL >= 0){
					//*Fhn = FL[0];
					dest[i] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4*Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SM >= 0 && SL < 0){
					dest[i] = (QSTARL[1] * nx[i] - QSTARL[2] * ny[i]);
					dest[i + 4*Ne] = (QSTARL[1] * ny[i] + QSTARL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QSTARL[n];
					}
					SPY = 1;
				}
				else if (SM < 0 && SR > 0){
					dest[i] = (QSTARR[1] * nx[i] - QSTARR[2] * ny[i]);
					dest[i + 4*Ne] = (QSTARR[1] * ny[i] + QSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QSTARR[n];
					}
					SPY = 1;
				}
				else if (0 >= SR){
					dest[i] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4*Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QR[n];
					}
					SPY = 1;
				}

			}
			else if (HL>Hcrit && HR <= Hcrit){
				/*For this situation, the three wave structure has degenerated to the two wave structure.
				*For this special case, i.e. HR <= Hcrit, SM = SR.
				*/
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
				if (SL >= 0){
					dest[i] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4*Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					QSTARL[0] = HL *(SL - UL) / (SL - SR);
					QSTARL[1] = HL *(SL - UL) / (SL - SR)*SR;
					QSTARL[2] = HL *(SL - UL) / (SL - SR)*VL;
					for (int i = 3; i<Nvar + 1; i++){
						QSTARL[i] = HL*((SL - UL) / (SL - SR)) * VARIABLEL[i];
					}
					dest[i] = QSTARL[1] * nx[i] - QSTARL[2] * ny[i];
					dest[i + 4*Ne] = QSTARL[1] * ny[i] + QSTARL[2] * nx[i];
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QSTARL[n];
					}
					SPY = 1;
				}
				else if (SR <= 0){
					/*Actually, I wonder whether this situation would occur*/
					dest[i] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4*Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QR[n];
					}
					SPY = 1;
				}
			}
			else if (HL <= Hcrit && HR > Hcrit){
				/*
				 * For this special case, i.e. HL <= Hcrit, SM = SL.
				 */
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
				if (SL >= 0){
					dest[i] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4*Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4*Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					QSTARR[0] = HR *(SR - UR) / (SR - SL);
					QSTARR[1] = HR *(SR - UR) / (SR - SL)*SL;
					QSTARR[2] = HR *(SR - UR) / (SR - SL)*VR;
					for (int i = 3; i<Nvar + 1; i++){
						QSTARR[i] = HR*((SR - UR) / (SR - SL)) * VARIABLER[i];
					}
					dest[i] = (QSTARR[1] * nx[i] - QSTARR[2] * ny[i]);
					dest[i + 4 * Ne] = (QSTARR[1] * ny[i] + QSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QSTARR[n];
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4*Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * 4 * Ne] = QR[n];
					}
					SPY = 1;
				}
			}
			//SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
			if (SPY == 0){
				printf("Error occured when calculating the HLLC Riemann problem! check please!");
			}
		}
	}
	free(QL), free(QR);
	free(QSTARL), free(QSTARR);
	free(VARIABLEL), free(VARIABLER);
	free(fm), free(fp);
	free(nx), free(ny);
}


void EvaluatePhysicalVariableByDepthThreshold(double hmin, double *h, double *variable, double *outPut){
	if (*h > hmin){
		*outPut = *variable / *h;
	}
	else{
		*outPut = 0;
	}
}

/** Rotate flux to outward normal direction */
void RotateFluxToNormal2d(double *hu, ///< flux at x component
	double *hv, ///< flux at y component
	double *nx, ///< outward normal vector
	double *ny, ///< outward normal vector
	double *qn,      ///< normal flux
	double *qv       ///< tangent flux
	) {
	*qn = +*hu * *nx + *hv * *ny;
	*qv = -*hu * *ny + *hv * *nx;
	return;
}

/*HLLC numerical flux is adopted in vertical face, */
void EvaluateVerticalFaceNumFlux_HLLC_LAI(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	/*The HLLC numerical flux presented in LAI(2012, Modeling one- and two-dimensional shallow water flows with discontinuous Galerkin method) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SM, SR;
	double USTAR, HSTAR;
	double *hum = fm, *hvm = fm + Nfp*Ne, *hm = fm + 2 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *hp = fp + 2 * Nfp*Ne;

	double *QL = (double *)malloc((Nvar + 1)*sizeof(double)), *QR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *FL = (double *)malloc((Nvar + 1)*sizeof(double)), *FR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *FSTARL = (double *)malloc((Nvar + 1)*sizeof(double)), *FSTARR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *VARIABLEL = (double *)malloc((Nvar + 1)*sizeof(double)), *VARIABLER = (double *)malloc((Nvar + 1)*sizeof(double));
	double *QSTARL = (double *)malloc((Nvar + 1)*sizeof(double)), *QSTARR = (double *)malloc((Nvar + 1)*sizeof(double));
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar + 1; n++)
		{
			QL[n] = *(fm + n*Nfp*Ne + i);
			QR[n] = *(fp + n*Nfp*Ne + i);
		}

		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < Nvar + 1; n++){
			/*In this part, we calculate the original variable. If the water depth is smaller than the critical value,
			we directly set the value of the original variable to be zero, and the original value is further used to calculate
			the flux term in left and right side. The question here is whether the water depth should be set to zero if it is
			smaller than the critical value, since it is used when calcualte the flux term. In the current version, we leave them
			unchanged*/
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];
		/*Here, we use the critical value to categorize different condition, I wonder whether this value can be
		directly set to zero, if this value is set to zero, then they calculation in the further part can be looked
		to be very smooth*/
		if (!(HL < Hcrit && HR < Hcrit))//等号能不能取？
		{

			/*Compute FL AND FR*/
			FL[0] = HL*UL;
			FL[1] = HL*pow(UL, 2) + 0.5*(*gra)*pow(HL, 2);
			FL[2] = HL*UL*VL;
			for (int n = 3; n < Nvar + 1; n++){
				FL[n] = FL[0] * VARIABLEL[n];
			}

			FR[0] = HR*UR;
			FR[1] = HR*pow(UR, 2) + 0.5*(*gra)*pow(HR, 2);
			FR[2] = HR*UR*VR;
			for (int n = 3; n < Nvar + 1; n++){
				FR[n] = FR[0] * VARIABLER[n];
			}

			if ((HL>Hcrit) & (HR > Hcrit)){
				USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				HSTAR = 1.0 / (*gra)*pow(0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR), 2);
				SL = fmin(UL - sqrt(*gra * HL), USTAR - sqrt(*gra*HSTAR));
				SR = fmax(UR + sqrt(*gra * HR), USTAR + sqrt(*gra*HSTAR));
				SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
				/*Compute QSTARL AND QSTARR*/
				QSTARL[0] = HL*((SL - UL) / (SL - SM)) * 1;
				QSTARL[1] = HL*((SL - UL) / (SL - SM)) * SM;
				QSTARR[0] = HR*((SR - UR) / (SR - SM)) * 1;
				QSTARR[1] = HR*((SR - UR) / (SR - SM)) * SM;
				for (int i = 2; i<Nvar + 1; i++){
					QSTARL[i] = HL*((SL - UL) / (SL - SM)) * VARIABLEL[i];
					QSTARR[i] = HR*((SR - UR) / (SR - SM)) * VARIABLER[i];
				}
				/*Compute FSTARL and FSTARR*/
				for (int n = 0; n < Nvar + 1; n++){
					FSTARL[n] = FL[n] + SL * (QSTARL[n] - QL[n]);
					FSTARR[n] = FR[n] + SR * (QSTARR[n] - QR[n]);
				}
				/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
				if (SL >= 0){
					//*Fhn = FL[0];
					dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SM >= 0 && SL < 0){
					dest[i] = (FSTARL[1] * nx[i] - FSTARL[2] * ny[i]);
					dest[i + Nfp*Ne] = (FSTARL[1] * ny[i] + FSTARL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FSTARL[n];
					}
					SPY = 1;
				}
				else if (SM < 0 && SR > 0){
					dest[i] = (FSTARR[1] * nx[i] - FSTARR[2] * ny[i]);
					dest[i + Nfp*Ne] = (FSTARR[1] * ny[i] + FSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FSTARR[n];
					}
					SPY = 1;
				}
				else if (0 >= SR){
					dest[i] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}

			}
			else if (HL>Hcrit && HR <= Hcrit){
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
				//是否将SM的计算移到这里来
				//				SM = SR;
				if (SL >= 0){
					dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					double Fqxn = (SR * FL[1] - SL * FR[1] + SL * SR * (HR*UR - HL*UL)) / (SR - SL);
					double Fqyn = (SR * FL[2] - SL * FR[2] + SL * SR * (HR*VR - HL*VL)) / (SR - SL);
					dest[i] = (Fqxn * nx[i] - Fqyn * ny[i]);
					dest[i + Nfp*Ne] = (Fqxn * ny[i] + Fqyn * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1)*Nfp*Ne] = (SR * FL[n] - SL * FR[n] + SL * SR * (QR[n] - QL[n])) / (SR - SL);
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}
			}
			else if (HL <= Hcrit && HR > Hcrit){
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
				//				SM = SL;
				if (SL >= 0){
					dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					double Fqxn = (SR * FL[1] - SL * FR[1] + SL * SR * (HR*UR - HL*UL)) / (SR - SL);
					double Fqyn = (SR * FL[2] - SL * FR[2] + SL * SR * (HR*VR - HL*VL)) / (SR - SL);
					dest[i] = (Fqxn * nx[i] - Fqyn * ny[i]);
					dest[i + Nfp*Ne] = (Fqxn * ny[i] + Fqyn * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1)*Nfp*Ne] = (SR * FL[n] - SL * FR[n] + SL * SR * (QR[n] - QL[n])) / (SR - SL);
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar + 1; n++){
						dest[i + (n - 1) * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}
			}
			//SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
			if (SPY == 0){
				printf("Error in ADV HLLC flux!");
			}
		}
	}
	free(QL), free(QR);
	free(QSTARL), free(QSTARR);
	free(FL), free(FR);
	free(FSTARL), free(FSTARR);
	free(VARIABLEL), free(VARIABLER);
}

/*HLLC numerical flux is adopted in vertical face, */
void EvaluateVerticalFaceNumFlux_HLL(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	/*The HLL numerical flux presented in Toro(2020, Shock-capturing methods for free-surface shallow flows) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SR;
	//double HSTAR, USTAR;
	double *hum = fm, *hvm = fm + Nfp*Ne, *hm = fm + 2 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *hp = fp + 2 * Nfp*Ne;

	double *QL = (double *)malloc((Nvar + 1)*sizeof(double)), *QR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *FL = (double *)malloc((Nvar + 1)*sizeof(double)), *FR = (double *)malloc((Nvar + 1)*sizeof(double));
	double *VARIABLEL = (double *)malloc((Nvar + 1)*sizeof(double)), *VARIABLER = (double *)malloc((Nvar + 1)*sizeof(double));
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar + 1; n++)
		{
			QL[n] = *(fm + n*Nfp*Ne + i);
			QR[n] = *(fp + n*Nfp*Ne + i);
		}

		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < Nvar + 1; n++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];

		if (!(HL < Hcrit && HR < Hcrit))
		{

			if (HL > Hcrit && HR>Hcrit){
				double USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				double CSTAR = 0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR);
				SL = fmin(UL - sqrt(*gra*HL), USTAR - CSTAR);
				SR = fmax(UR + sqrt(*gra*HR), USTAR + CSTAR);
			}
			else if (HL>Hcrit && HR <= Hcrit){
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
			}
			else if (HL <= Hcrit && HR > Hcrit){
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
			}
			else{
				SL = 0.0;
				SR = 0.0;
			}
			/*Compute FL AND FR*/
			FL[0] = HL*UL;
			FL[1] = HL*pow(UL, 2) + 0.5*(*gra)*pow(HL, 2);
			FL[2] = HL*UL*VL;
			/*
			for (int n = 3; n < Nvar + 1; n++){
			FL[n] = FL[0] * VARIABLEL[n];
			}
			*/

			FR[0] = HR*UR;
			FR[1] = HR*pow(UR, 2) + 0.5*(*gra)*pow(HR, 2);
			FR[2] = HR*UR*VR;
			/*
			for (int n = 3; n < Nvar + 1; n++){
			FR[n] = FR[0] * VARIABLER[n];
			}
			*/

			/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
			if (SL >= 0 && SR > 0){
				//*Fhn = FL[0];
				dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
				dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
				/*
				for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FL[n];
				}
				*/
				SPY = 1;
			}
			else if (SL < 0 && SR > 0){
				double Fqxn = (SR * FL[1] - SL * FR[1] + SL * SR * (HR*UR - HL*UL)) / (SR - SL);
				double Fqyn = (SR * FL[2] - SL * FR[2] + SL * SR * (HR*VR - HL*VL)) / (SR - SL);
				dest[i] = (Fqxn * nx[i] - Fqyn * ny[i]);
				dest[i + Nfp*Ne] = (Fqxn * ny[i] + Fqyn * nx[i]);
				/*
				for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FSTARL[n];
				}
				*/
				SPY = 1;
			}
			else if (SL < 0 && SR <= 0){
				dest[i] = (FR[1] * nx[i] - FR[2] * ny[i]);
				dest[i + Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
				/*
				for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FSTARR[n];
				}
				*/
				SPY = 1;
			}
			else if ((fabs(SL) < EPS) & (fabs(SR) < EPS)){
				dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
				dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
				/*
				for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FR[n];
				}
				*/
				SPY = 1;
			}
			if (SPY == 0){
				printf("Error occured when calculating the HLL flux! check please!");
			}
		}
	}
	free(QL), free(QR);
	free(FL), free(FR);
	free(VARIABLEL), free(VARIABLER);
}
