/*******************************************************************************
 * File:        Roe_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _ROE_FLUXES_H
#define	_ROE_FLUXES_H

#include "Vector3D.h"

// Roe Scheme Functions
void Roe_Init(void);
void Roe_Finalize(void);
void Roe_Reset(void);

void Roe_Init_New(void);
void Roe_Finalize_New(void);
void Roe_Reset_New(void);


void Compute_Roe_Transformed_Precondition_Matrix_None(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_CecileVoizat(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_BTW(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_Eriksson(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_Merkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_Turkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix_WeissSmith(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Roe_Transformed_Precondition_Matrix(int nodeID, int reqInv, double **PrecondMatrix);

void Compute_Transformed_Residual_Roe(void);
void Compute_Steady_Residual_Roe_Precondition_CecileVoizat(void);
void Compute_Steady_Residual_Roe_Precondition_BTW(void);
void Compute_Steady_Residual_Roe_Precondition_Eriksson(void);
void Compute_Steady_Residual_Roe_Precondition_Merkel(void);
void Compute_Steady_Residual_Roe_Precondition_Turkel(void);
void Compute_Steady_Residual_Roe_Precondition_WeissSmith(void);

void Compute_Flux_Roe_OneSided(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_LD(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Original(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_CecileVoizat(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_BTW(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_Eriksson(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_WeissSmith(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);

void Compute_Dissipation_Matrix_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Original(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_CecileVoizat(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_BTW(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_Eriksson(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_WeissSmith(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);

void Compute_Residual_Roe_New(int AddTime);


void Compute_RoeVariables(double *Q_L, double *Q_R, double *Q_Roe);
void Compute_RoeAJacobian(double *Q_L, double *Q_R, Vector3D areavec, double **AJacobian_Roe);
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_RoeFlux_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Residual_Roe(int AddTime);
void Compute_Jacobian_Exact_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Approximate_Roe(int AddTime, int Iteration);
void Compute_Jacobian_FiniteDifference_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Roe(int AddTime, int Iteration);

#endif	/* _ROE_FLUXES_H */

