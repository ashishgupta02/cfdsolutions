/*******************************************************************************
 * File:        Roe_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _ROE_FLUXES_H
#define	_ROE_FLUXES_H

#include "Vector3D.h"

// Roe Scheme Functions
void Roe_Init(void);
void Roe_Finalize(void);
void Roe_Reset(void);

void Compute_RoeAverage_Q(double *Q_L, double *Q_R, double *Q_Roe);

void Compute_Transformed_Preconditioner_Matrix_Roe_None(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_Roe_Merkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_Roe_Turkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_Roe(int nodeID, int reqInv, double **PrecondMatrix);

void Compute_Transformed_Residual_Roe(void);
void Compute_Steady_Residual_Roe_Precondition_Merkel(void);
void Compute_Steady_Residual_Roe_Precondition_Turkel(void);

void Compute_Flux_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Thornber(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Original(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Flux_Roe(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);

void Compute_Dissipation_Matrix_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Thornber(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Original(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);
void Compute_Dissipation_Matrix_Roe(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe);

void Compute_Residual_Roe(int AddTime);

#endif	/* _ROE_FLUXES_H */

