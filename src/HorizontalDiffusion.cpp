#include "HorizontalDiffusion.h"

void GetFacialFluxTerm(double *dest, double *source, double *n, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = source[i] * n[i];
}
void GetCentralNumFlux(double *dest, double *fm, double *fp, double *n, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * (fm[i] + fp[i]) / 2;
}
void GetIPNumFlux(double *dest, double *fm, double *fp, double *n, double *vfm, \
	double *vfp, double *Jumpn, int Nfp, double *Tau, double Coefficient){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * (fm[i] + fp[i]) / 2 - n[i] * Tau[i] / Coefficient*(Jumpn[i] * (vfm[i] - vfp[i]));
}

void GetIPBoundaryNumFlux(double *dest, double *fm, double *n, double *vfm, \
	double *vfp, double *Jumpn, int Nfp, double *Tau, double Coefficient){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * fm[i] - n[i] * Tau[i] / Coefficient*(Jumpn[i] * (vfm[i] - vfp[i]));
}

/*
* Purpose: This function is used to calculate the inner edge contribution to the auxialary variable
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       int[ 1 ]  face index of the studied inner face
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       int[ 1 ]  field index of studied physical field, i.e. index corresponds to $u,v$ and $\theta$
*       int[Nfp x 1] fm the local physical value on the studied inner edge
*       int[Nfp x 1] fp the adjacent physical value on the studied inner edge
*       double[2 x 1] IEFToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] IEFToN1 index of the local facial interpolation points
*       double[Nfp x 1] IEFToN2 index of the adjacent facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
*       double[Nfp x 1] FluxM the local flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxP the adjacent flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n the direction vector of the studied face
*       double[Nfp x Nfp] Mb the mass matrix of the standard cell
*       double[Nfp x 1]  Js the faical Jacobian of the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetIEContributionToAuxialaryVariable(double *dest, int face, int Ne, int Nfp, int field, double *fm, double *fp, double *IEFToE, double *IEFToF,\
	double *IEFToN1, double *IEFToN2, int Np, int K, double *FluxM, double *FluxP, double *FluxS, double *n, double *Mb, double *Js){
	/*This */
	GetFacialFluxTerm(FluxM + field*Ne*Nfp + face*Nfp, fm + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	GetFacialFluxTerm(FluxP + field*Ne*Nfp + face*Nfp, fp + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	GetCentralNumFlux(FluxS + field*Ne*Nfp + face*Nfp, fm + field*Ne*Nfp + face*Nfp, fp + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, Nfp, IEFToN1, IEFToN2, FluxM + field*Ne*Nfp, FluxP + field*Ne*Nfp, FluxS + field*Ne*Nfp, Js, Mb, dest);
}
/*
* Purpose: This function is used to calculate the boundary edge contribution to the auxialary variable
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       int[ 1 ]  face index of the studied boundary  face
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       int[ 1 ]  field index of studied physical field, i.e. index corresponds to $u,v$ and $\theta$
*       int[Nfp x 1] fm the local physical value on the studied inner edge
*       int[Nfp x 1] fp the boundary condition on the studied boundary edge
*       double[2 x 1] BEFToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] BEFToN1 index of the local facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
*       double[Nfp x 1] FluxM the local flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n the direction vector of the studied face
*       double[Nfp x Nfp] Mb the mass matrix of the standard cell
*       double[Nfp x 1]  Js the faical Jacobian of the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetBEContributionToAuxialaryVariable(double *dest, int face, int Ne, int Nfp, int field, double *fm, double *fp, double *BEFToE, double *BEFToF, \
	double *BEFToN1, int Np, int K, double *FluxM, double *FluxS, double *n, double *Mb, double *Js){
	GetFacialFluxTerm(FluxM + field*Ne*Nfp + face*Nfp, fm + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	GetFacialFluxTerm(FluxS + field*Ne*Nfp + face*Nfp, fp + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, Nfp, BEFToN1, FluxM + field*Ne*Nfp, FluxS + field*Ne*Nfp, Js, Mb, dest);
}

/*
* Purpose: This function is used to calculate the inner edge contribution to the right hand side
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       double[Nfp x 1]  LPDTfm local value of the local partial diffrential term on the studied inner face
*       double[Nfp x 1]  LPDTfp adjacent value of the local partial diffrential term on the studied inner face
* 		int[ 1 ]  face index of the studied inner face
*       double[Np x K]  LPDiffTerm the local partial differential term
*       double[2 x IENe] FToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] FToN1 index of the local facial interpolation points
*       double[Nfp x 1] FToN2 index of the adjacent facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n direction vector used to multiply the numerical flux
*       double[Nfp x 1] fm the local physical value on the studied inner edge, used to calculate the jump term
*       double[Nfp x 1] fp the boundary condition on the studied boundary edge, used to calculate the jump term
*       double[Nfp x 1] Jumpn direction vector used to calculate the jump term
*       double[Nfp x 1] Tau the penalty parameter used in interior penalty flux
*       double[ 1 ] Coefficient the parameter stands for prantl number or 1.0
*       double[Nfp x 1] AVfm the local value of the auxialary variable
*       double[Nfp x 1] AVfp the adjacent value of the auxialary variable
*       double[Np x K] AV the auxialary variable
*       double[Nfp x 1] FluxM the local flux term relates to AVfm
*       double[Nfp x 1] FluxP the adjacent flux term relates to AVfp
*       double[Nfp x 1] Js the facial Jacobian parameter
*       double[Nfp x Nfp] M the mass matrix of the master cell corresponds to the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetIEContributionToRHS(double *dest, double *LPDTfm, double *LPDTfp, int face, double *LPDiffTerm, double *FToE, double *FToF, double *FToN1, \
	double *FToN2, int Np, int K, int Ne, int Nfp, int field, double *FluxS, double *n, double *fm, double *fp, double *Jumpn, double *Tau, double Coefficient, \
	double *AVfm, double *AVfp, double *AV, double *FluxM, double *FluxP, double *Js, double *M){
	/*FToE is the start stress of the property FToE contained in inner edge*/
	FetchInnerEdgeFacialValue(LPDTfm + Ne*Nfp*field + face*Nfp, LPDTfp + Ne*Nfp*field + face*Nfp, LPDiffTerm, FToE + 2 * face, FToN1 + face*Nfp, FToN2 + face*Nfp, Np, Nfp);
	GetIPNumFlux(FluxS + field*Ne*Nfp + face*Nfp, LPDTfm + Ne*Nfp*field + face*Nfp, LPDTfp + Ne*Nfp*field + face*Nfp, n + face*Nfp, fm + face*Nfp, fp + face*Nfp, Jumpn + face*Nfp, Nfp, Tau + face*Nfp, Coefficient);
	FetchInnerEdgeFacialValue(AVfm + field*Ne*Nfp + face*Nfp, AVfp + field*Ne*Nfp + face*Nfp, AV, FToE + 2 * face, FToN1 + face*Nfp, FToN2 + face*Nfp, Np, Nfp);
	GetFacialFluxTerm(FluxM + field*Ne*Nfp + face*Nfp, AVfm + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	GetFacialFluxTerm(FluxP + field*Ne*Nfp + face*Nfp, AVfp + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	StrongFormInnerEdgeRHS(face, FToE, FToF, Np, K, Nfp, FToN1, FToN2, FluxM + field*Ne*Nfp, FluxP + field*Ne*Nfp, FluxS + field*Ne*Nfp, Js, M, dest);
}

/*
* Purpose: This function is used to calculate the boundary edge contribution to the right hand side
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       double[Nfp x 1]  LPDTfm local value of the local partial diffrential term on the studied inner face
* 		int[ 1 ]  face index of the studied inner face
*       double[Np x K]  LPDiffTerm the local partial differential term
*       double[2 x IENe] FToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] FToN1 index of the local facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n direction vector used to multiply the numerical flux
*       double[Nfp x 1] fm the local physical value on the studied inner edge, used to calculate the jump term
*       double[Nfp x 1] fp the boundary condition on the studied boundary edge, used to calculate the jump term
*       double[Nfp x 1] Jumpn direction vector used to calculate the jump term
*       double[Nfp x 1] Tau the penalty parameter used in interior penalty flux
*       double[ 1 ] Coefficient the parameter stands for prantl number or 1.0
*       double[Nfp x 1] AVfm the local value of the auxialary variable
*       double[Np x K] AV the auxialary variable
*       double[Nfp x 1] FluxM the local flux term relates to AVfm
*       double[Nfp x 1] Js the facial Jacobian parameter
*       double[Nfp x Nfp] M the mass matrix of the master cell corresponds to the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/

void GetBEContributionToRHS(double *dest, double *LPDTfm, int face, double *LPDiffTerm, double *FToE, double *FToF, double *FToN1, \
	int Np, int K, int Ne,  int Nfp, int field, double *FluxS, double *n, double *fm, double *fp, double *Jumpn, double *Tau, double Coefficient, \
	double *AVfm, double *AV, double *FluxM, double *Js, double *M){
	/*FToE is the start adress of the property FToE contained in inner edge*/
	FetchBoundaryEdgeFacialValue(LPDTfm + field*Ne*Nfp + face*Nfp, LPDiffTerm, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
	GetIPBoundaryNumFlux(FluxS + field*Ne*Nfp + face*Nfp, LPDTfm + field*Ne*Nfp + face*Nfp, n + face*Nfp, fm + face*Nfp, fp + face*Nfp, Jumpn + face*Nfp, Nfp, Tau + face*Nfp, Coefficient);
	FetchBoundaryEdgeFacialValue(AVfm + field*Ne*Nfp + face*Nfp, AV, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
	GetFacialFluxTerm(FluxM + field*Ne*Nfp + face*Nfp, AVfm + field*Ne*Nfp + face*Nfp, n + face*Nfp, Nfp);
	StrongFormBoundaryEdgeRHS(face, FToE, FToF, Np, K, Nfp, FToN1, FluxM + field*Ne*Nfp, FluxS + field*Ne*Nfp, Js, M, dest);
}