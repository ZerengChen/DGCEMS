#include "NdgSWE3D.h"
#include "NdgSWE.h"
#include "NdgMath.h"

void FetchFacialValue(double *dest, double *source, int size, int *fpIndex){
	for (int i = 0; i < size; i++){
		dest[i] = source[fpIndex[i]];
	}
}

void EvaluateVerticalFaceSurfFlux(double *dest, double *fm, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	double u, v, theta;
	double *hu = fm, *hv = fm + Nfp*Ne, *h = fm + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hu + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hv + i, &v);
		dest[i] = (hu[i] * u + 0.5 * (*gra)*pow(h[i],2))*nx[i] + hu[i] * v*ny[i];
		dest[i + Ne*Nfp] = hu[i] * v*nx[i] + (hv[i] * v + 0.5 * (*gra) * pow(h[i], 2))*ny[i];
		for (int field = 3; field < Nvar+1; field++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, fm + field*Ne*Nfp + i, &theta);
			dest[i + (field-1)*Ne*Nfp] = hu[i] * theta * nx[i] + hv[i] * theta * ny[i];
		}
	}
}

void EvaluateHorizontalFaceSurfFlux(double *flux, double *fm, double *nz, double Hcrit, int Nfp, int Nvar, int Ne){
	double *hum = fm, *hvm = fm + Nfp*Ne, *omega = fm + 2 * Nfp*Ne, *H = fm + 3 * Nfp*Ne;
	double u, v, variable, theta;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, hum + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, hvm + i, &v);
		*(flux + i) = u*omega[i]*nz[i];
		*(flux + Ne*Nfp + i) = v*omega[i]*nz[i];
		for (int n = 2; n < Nvar; n++){
			/*Here 2*Ne*Nfp stands for the memory occupied by h and omega*/
			variable = *(fm + 2 * Ne*Nfp + n*Ne*Nfp + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, &variable, &theta);
			*(flux + n*Ne*Nfp + i) = theta*omega[i]*nz[i];
		}
	}
}

void EvaluateHorizontalFaceNumFlux(double *FluxS, double *fm, double *fp, double *nz, double Hcrit, int Nfp, int Nvar, int Ne){
	double *hum = fm, *hvm = fm + Nfp*Ne, *omegam = fm + 2 * Nfp*Ne, *Hm = fm + 3 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *omegap = fp + 2 * Nfp*Ne, *Hp = fp + 3 * Nfp*Ne;
	double um, vm, up, vp, variablem, variablep, thetam, thetap;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hum + i, &um);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hvm + i, &vm);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hup + i, &up);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hvp + i, &vp);
		*(FluxS + i) = 0.5*(um*omegam[i] + up*omegap[i]) * nz[i];
		*(FluxS + Ne*Nfp + i) = 0.5*(vm*omegam[i] + vp*omegap[i]) * nz[i];
		for (int n = 2; n < Nvar; n++){
			/*Here 2*Ne*Nfp stands for the memory occupied by h and omega*/
			variablem = *(fm + 2 * Ne*Nfp + n*Ne*Nfp + i);
			variablep = *(fp + 2 * Ne*Nfp + n*Ne*Nfp + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, &variablem, &thetam);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, &variablep, &thetap);
			*(FluxS + n*Ne*Nfp + i) = 0.5*(thetam*omegam[i] + thetap*omegap[i]) * nz[i];
		}
	}

}

void EvaluatePrebalanceVolumeTerm(double *Edest, double *Gdest, double *Hdest, double *fphys, \
	int *varIndex, int Nvar, double *gra, int Np, int K, double Hcrit)
{
	double *hu = fphys, *hv = fphys + Np*K, *omega = fphys + 2 * Np*K;
	double *h = fphys + 3 * Np*K, *z = fphys + 5 * Np*K;
	double u, v, variable, theta;
	for (int i = 0; i < Np; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hu + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hv + i, &v);
		/*NOTE THAT: THE WET AND DRY NEED TO BE TACKLED HERE*/
		Edest[i] = hu[i] * u + 0.5 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
		Gdest[i] = hu[i] * v;
		Hdest[i] = u*omega[i];
		Edest[i + Np*K] = hu[i] * v;
		Gdest[i + Np*K] = hv[i] * v + 0.5 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
		Hdest[i + Np*K] = v*omega[i];
		/*If temperature, salt, sediment or other passive transport material is included, the following part is used to calculate
		the volume term corresponding to these terms*/
		for (int n = 2; n < Nvar; n++){
			variable = *(fphys + ((int)varIndex[n] - 1)*Np*K + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, &variable, &theta);
			Edest[i + n*Np*K] = hu[i] * theta;
			Gdest[i + n*Np*K] = hv[i] * theta;
			Hdest[i + n*Np*K] = omega[i] * theta;
		}
	}
}
//For radiation and roller
void EvaluatePrebalanceVolumeSourceTerm(double *Edest, double *Gdest, double *fxx, double *fxy, double *fyy, int Np, int K)
{
	for (int i = 0; i < Np; i++) {
		Edest[i] = Edest[i] - fxx[i];
		Gdest[i] = Gdest[i] - fxy[i];
		Edest[i + Np * K] = Edest[i + Np * K] - fxy[i];
		Gdest[i + Np * K] = Gdest[i + Np * K] - fyy[i];
		
	}
}

void GetModCoefficient(double *dest, double *source, double *InvV, int Np, int NLayer){
	int RowOPA = Np;
	int ColOPB = NLayer;
	int ColOPA = Np;
	double alpha = 1.0;
	int LDA = Np;
	int LDB = Np;
	double Beta = 0.0;
	int LDC = Np;
	double *B = source;
	MatrixMultiply(InvV, B, dest, RowOPA, ColOPB, ColOPA, alpha);
	//dgemm("N", "N", &RowOPA, &ColOPB, &ColOPA, &alpha, InvV, &LDA, B, &LDB, &Beta, dest, &LDC);
}

void GetIntegralValue(double *dest, int Np2d, double *V, double *source){
	int RowOPA = Np2d;
	int ColOPA = Np2d;
	int ColOPB = 1;
	double alpha = 1.0;
	int LDA = Np2d;
	double *B = source;
	int LDB = Np2d;
	double Beta = 0.0;
	int LDC = Np2d;
	MatrixMultiply(V, B, dest, RowOPA, ColOPB, ColOPA, alpha);
	//dgemm("N", "N", &RowOPA, &ColOPB, &ColOPA, &alpha, V, &LDA, B, &LDB, &Beta, dest, &LDC);
	MultiplyByConstant(dest, dest, sqrt(2), Np2d);
}

/*
* Purpose: This function is used to get the depth-averaged value
*
* Input:
*      double * dest the position to store the depth averaged value
* 	   int[1] Np2d the number of the interpolation point of the two-dimensional master cell
*      double * V2d the start position of the dimensional Vandermode matrix
* 	   double * Tempdest the tempora position to store the result
* 	   double * Tempfield3d the three dimensional field with the Jacobian considered
*      double * field3d the three dimensional field
* 	   double * Jz the three dimensional Jacobian coefficient
*      double * fmod the space to store the modal coefficient
*      double * InvV3d the inverse of the three dimensional Vandermonde matrix
*      int[1] Np3d the number of the interpolation point of the three-dimensional master cell
*      int[1] NLayer the number of vertical layer
* Output:
* 		double * dest the position to store the depth averaged value
*/

void VerticalColumnIntegralField3d(double *dest, int Np2d, double *V2d, double *Tempdest, \
	double *Tempfield3d, double *field3d, double * Jz, \
	double *fmod, double *InvV3d, int Np3d, int NLayer){

	DotProduct(Tempfield3d, field3d, Jz, Np3d*NLayer);
	GetModCoefficient(fmod, Tempfield3d, InvV3d, Np3d, NLayer);
	for (int L = 0; L < NLayer; L++){
		Add(Tempdest, Tempdest, fmod + L*Np3d, Np2d);
	}
	GetIntegralValue(dest, Np2d, V2d, Tempdest);
}


/*
* Purpose: This function is used to integrate the facial physical field from bottom to surface
*
* Input:
*      double[LNfp2d] dest the integral field, with LNfp2d the number of the intepolation point of the one-dimensional line in horizontal direction
* 	   double[Nfp2d x Nlayer] source the physical field to be integrated, with Nfp2d the number of the intepolation point of the 2d face, and Nlayer the number of layers
*      double[Nfp2d] fmod space used to store the Vandemonde coefficient
* 	   double[Nfp2d x Nfp2d] InvV2d the inverse two dimensional Vandermode matrix
* 	   int[1] Nfp2d the number of the intepolation point
*      double[Nfp2d x Nlayer] Jz the facial Jacobian of the vertical face
* 	   int[1] Nlayer the number of layers
*      double[LNfp2d x LNfp2d] V1d the one dimensional Vandermode matrix
*      int[1] LNfp2d the number of the intepolation point of the one-dimensional line in horizontal direction
*      int[1] FToF the face to face topological relation
* Output:
* 		double[LNfp] dest the integral field, with LNfp2d the number of the intepolation point of the one-dimensional line in horizontal direction
*/

void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF){
	double *TempSource = (double*)malloc((int)Nfp2d*sizeof(double));
	double *Tempdest = (double*)malloc((int)LNfp2d*sizeof(double));
	memset(TempSource, 0, (int)Nfp2d * sizeof(double));
	memset(Tempdest, 0, (int)LNfp2d*sizeof(double));
	/*ptrdiff_t one = 1;
	double beta = 0, alpha = 1.0;*/
	int one = 1;
	for (int i = 0; i < Nlayer; i++){
		DotProduct(TempSource, source + i*(int)Nfp2d, Jz + i*(int)Nfp2d, (int)Nfp2d);
		//dgemm("N", "N", &RowOPA, &ColOPB, &ColOPA, &alpha, V, &LDA, B, &LDB, &Beta, dest, &LDC);
		MatrixMultiply(InvV2d, TempSource, fmod, Nfp2d, one, Nfp2d, 1.0);
		//dgemm("N", "N", &Nfp2d, &one, &Nfp2d, &alpha, InvV2d, &Nfp2d, TempSource, &Nfp2d, &beta, fmod, &Nfp2d);
		for (int j = 0; j < LNfp2d; j++){
			Tempdest[j] += fmod[j];
		}
	}
	MatrixMultiply(V1d, Tempdest, dest, LNfp2d, one, LNfp2d, 1.0);
	//dgemm("N", "N", &LNfp2d, &one, &LNfp2d, &alpha, V1d, &LNfp2d, Tempdest, &LNfp2d, &beta, dest, &LNfp2d);
	DotDivideByConstant(dest, dest, sqrt(2.0) / 2, LNfp2d);
	if (FToF >= 3) Flip(dest, LNfp2d);
	
	free(TempSource);
	free(Tempdest);
}

/*
* Purpose: This function is used to integrate the three dimensional physical field from bottom to a given elevation
*
* Input:
*      double[Np x Nlayer] dest the pointer to the start address of the integral field, from bottom to any given elevation
* 	   double[Np x Nlayer] source three dimensional field to be integrated
*      double[Np x Nlayer] Jz the jacobian coefficient in vertical direction
* 	   double[Np] fmod the space used to store the mode coefficient
* 	   int[1] NLayer number of layers in vertical direction
*      ptrdiff_t[1] Np the number of interpolation points for the 3d master cell
* 	   double[Np x Np] InvV3d inverse Vandermonde matrix 
*      int[1] Np2d the number of interpolation points for the 2d master cell
*      int[1] Npz number of interpolation points in vertical direction
*      double[Np x Np] Vint the integral matrix used to integrate the physical field
* Output:
* 		double[Np x Nlayer] the pointer to the start address of the integral field, from bottom to any given elevation
*/
void VerticalIntegralFromBottom(double *dest, double *source, double *Jz, double *fmod, int NLayer, int Np, double *InvV3d, int Np2d, int Npz, double *Vint){
	
	memset(dest, 0, NLayer * Np * 1 * sizeof(double));
	double *TempSource = (double*)malloc((int)Np*sizeof(double));
	int *UpEid = (int*)malloc(Np2d*sizeof(int));
	for (int i = 0; i < Np2d; i++){
		UpEid[i] = (int)Np - Np2d + i;
	}
	double *FacialData = (double*)malloc(Np2d*sizeof(double));
	double *TempAverageData = (double*)malloc(Np*sizeof(double));
	int one = 1;
	double beta = 0, alpha = 1.0;

	DotProduct(TempSource, source + (NLayer - 1)*(int)Np, Jz + (NLayer - 1)*(int)Np, (int)Np);
	MatrixMultiply(InvV3d, TempSource, fmod, Np, one, Np, alpha);
	MatrixMultiply(Vint, fmod, dest + (NLayer - 1)*Np, Np, one, Np, alpha);
	/*dgemm("N", "N", &Np, &one, &Np, &alpha, InvV3d, &Np, TempSource, &Np, &beta, fmod, &Np);
	dgemm("N", "N", &Np, &one, &Np, &alpha, Vint, &Np, fmod, &Np, &beta, dest + (NLayer - 1)*(int)Np, &Np);*/

	FetchFacialValue(FacialData, dest + (NLayer - 1)*(int)Np, Np2d, UpEid);

	RepmatValue(TempAverageData, FacialData, Np2d, Npz);

	for (int L = 1; L < NLayer; L++){
		DotProduct(TempSource, source + (NLayer - L - 1)*(int)Np, Jz + (NLayer - L - 1)*(int)Np, (int)Np);
		MatrixMultiply(InvV3d, TempSource, fmod, Np, one, Np, alpha);
		MatrixMultiply(Vint, fmod, dest + (NLayer - L - 1)*Np, Np, one, Np, alpha);
		/*dgemm("N", "N", &Np, &one, &Np, &alpha, InvV3d, &Np, TempSource, &Np, &beta, fmod, &Np);
		dgemm("N", "N", &Np, &one, &Np, &alpha, Vint, &Np, fmod, &Np, &beta, dest + (NLayer - L - 1)*(int)Np, &Np);*/
		Add(dest + (NLayer - L - 1)*(int)Np, dest + (NLayer - L - 1)*(int)Np, TempAverageData, Np);
		FetchFacialValue(FacialData, dest + (NLayer - L - 1)*(int)Np, Np2d, UpEid);
		RepmatValue(TempAverageData, FacialData, Np2d, Npz);
	}

	free(TempSource);
	free(UpEid);
	free(FacialData);
	free(TempAverageData);
}

/*
* Purpose: This function is used to repmat the two-dimensional edge value in vertical direction
*
* Input:
*      double[Nfp x Nlayer] dest the pointer to the start address of the vertical extended value
* 	   double[LNfp] source the two-dimensional edge value
* 	   int[1] NLayer number of layers in vertical direction
*      int[1] Nfp the number of interpolation points for the 2d master cell
*      int[1] LNfp the number of interpolation points for the 1d master cell
*      int[1] Npz number of interpolation points in vertical direction
*      int[1] FToF the face index
* Output:
* 		double[Nfp x Nlayer] dest the pointer to the start address of the vertical extended value
*/

void VerticalIntegralFromSurface(double *dest, double *source, double *Jz, double *fmod, int NLayer, int Np, double *InvV3d, int Np2d, int Npz, double *Vint) {

	memset(dest, 0, NLayer * Np * sizeof(double));
	double *TempSource = (double*)malloc((int)Np * sizeof(double));
	int *BotEid = (int*)malloc(Np2d * sizeof(int));
	for (int i = 0; i < Np2d; i++) {
		BotEid[i] = i;
	}
	double *FacialData = (double*)malloc(Np2d * sizeof(double));
	double *TempAverageData = (double*)malloc(Np * sizeof(double));
	int one = 1;
	double beta = 0, alpha = 1.0;
	// For the surface most element
	DotProduct(TempSource, source + (1 - 1)*(int)Np, Jz + (1 - 1)*(int)Np, (int)Np);
	//dgemm("N", "N", &Np, &one, &Np, &alpha, InvV3d, &Np, TempSource, &Np, &beta, fmod, &Np);
	//dgemm("N", "N", &Np, &one, &Np, &alpha, Vint, &Np, fmod, &Np, &beta, dest + (1 - 1)*(int)Np, &Np);
	MatrixMultiply(InvV3d, TempSource, fmod, Np, one, Np, alpha);
	MatrixMultiply(Vint, fmod, dest + (1 - 1)*(int)Np, Np, one, Np, alpha);

	FetchFacialValue(FacialData, dest + (1 - 1)*(int)Np, Np2d, BotEid);

	RepmatValue(TempAverageData, FacialData, Np2d, Npz);

	for (int L = 1; L < NLayer; L++) {
		DotProduct(TempSource, source + L * (int)Np, Jz + L * (int)Np, (int)Np);
		//dgemm("N", "N", &Np, &one, &Np, &alpha, InvV3d, &Np, TempSource, &Np, &beta, fmod, &Np);
		//dgemm("N", "N", &Np, &one, &Np, &alpha, Vint, &Np, fmod, &Np, &beta, dest + L * (int)Np, &Np);
		MatrixMultiply(InvV3d, TempSource, fmod, Np, one, Np, alpha);
		MatrixMultiply(Vint, fmod, dest + L * (int)Np, Np, one, Np, alpha);

		Add(dest + L * (int)Np, dest + L * (int)Np, TempAverageData, Np);
		FetchFacialValue(FacialData, dest + L * (int)Np, Np2d, BotEid);
		RepmatValue(TempAverageData, FacialData, Np2d, Npz);
	}

	free(TempSource);
	free(BotEid);
	free(FacialData);
	free(TempAverageData);
}

void VerticalRepmatFacialValue_(double *dest, double *source, int NLayer, int Nfp, int LNfp, int Npz, int FToF)
{
	double *FlipValue = (double*)malloc(LNfp*sizeof(double));
	double *PerFaceValue = (double*)malloc(Nfp*sizeof(double));
	if (FToF >= 3){
		memcpy(FlipValue, source, LNfp*sizeof(double));
		Flip(FlipValue, LNfp);
		RepmatValue(PerFaceValue, FlipValue, LNfp, Npz);
		RepmatValue(dest, PerFaceValue, Nfp, NLayer);
	}
	else{
		RepmatValue(PerFaceValue, source, LNfp, Npz);
		RepmatValue(dest, PerFaceValue, Nfp, NLayer);
	}
	free(FlipValue);
	free(PerFaceValue);
}

//Purpose: These functions are used to calculate the roller and radiation
void EvaluateVerticalFaceSurfFluxRandS(double *dest, double *xxfm, double *xyfm, double *yyfm, double *nx, double *ny, int Nfp, int Ne) {
	for (int i = 0; i < Nfp; i++) {
		dest[i] = xxfm[i] *nx[i] + xyfm[i] * ny[i];
		dest[i + Ne * Nfp] = xyfm[i] *nx[i] + yyfm[i] *ny[i];
	}
}
void EvaluateHorizontalFaceSurfFluxRandS(double *dest, double *xfm, double *yfm, double *nz, int Nfp, int Ne) {
	for (int i = 0; i < Nfp; i++) {
		dest[i] = xfm[i] * nz[i];
		dest[i + Ne * Nfp] = yfm[i] * nz[i];
	}
}
void EvaluateHorizontalFaceCentralNumFlux(double *FluxS, double *xxfm, double *xyfm, double *yyfm, double *xxfp, double *xyfp, double *yyfp, double *nx, double *ny, int Nfp, int Ne) {
	for (int i = 0; i < Nfp; i++) {		
		//*(FluxS + i) = 0.5*(xxfm[i] + xxfp[i] + xyfm[i] + xyfp[i]) * nx[i];
		//*(FluxS + Ne * Nfp + i) = 0.5*(xyfm[i] + xyfp[i] + yyfm[i] + yyfp[i]) * ny[i];
		*(FluxS + i) = 0.5*(xxfm[i] + xxfp[i]) * nx[i] + 0.5 * (xyfm[i] + xyfp[i]) * ny[i];
		*(FluxS + Ne * Nfp + i) = 0.5*(xyfm[i] + xyfp[i]) * nx[i] + 0.5 * (yyfm[i] + yyfp[i]) * ny[i];
	}
}
void EvaluateVerticalFaceCentralNumFlux(double *FluxS, double *fm, double *fp, double *nz, int Nfp, int Ne) {
	for (int i = 0; i < Nfp; i++) {
		*(FluxS + i) = 0.5*(fm[i] + fp[i]) * nz[i];
	}
}
/** Evaluate flux term in surface integration */
/*
void EvaluateFluxTerm2d( double hmin, ///< water threshold
	double *gra,  ///< gravity acceleration
	double *h,    ///< water depth
	double *hu,   ///< flux variable
	double *hv,   ///< flux variable
	double *E        ///< surface integral flux term
	) {
	double u, v;
	if (h > hmin) {
		u = hu / h;
		v = hv / h;
	}
	else {
		u = 0.0;
		v = 0.0;
	}
	double huv = h * u * v;
    double h2 = h * h;
	E[0] = hu;
	E[1] = h * u * u + 0.5 * gra * h2;
	return;
}
*/