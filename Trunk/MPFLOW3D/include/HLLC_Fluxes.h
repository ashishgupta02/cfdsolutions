/*******************************************************************************
 * File:        HLLC_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _HLLC_FLUXES_H
#define	_HLLC_FLUXES_H

#include "Vector3D.h"

// HLLC Scheme Functions
void HLLC_Init(void);
void HLLC_Finalize(void);
void HLLC_Reset(void);

void Compute_Transformed_Preconditioner_Matrix_HLLC_None(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_HLLC_Merkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_HLLC_Turkel(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Transformed_Preconditioner_Matrix_HLLC(int nodeID, int reqInv, double **PrecondMatrix);

void Compute_Transformed_Residual_HLLC(void);
void Compute_Steady_Residual_HLLC_Precondition_Merkel(void);
void Compute_Steady_Residual_HLLC_Precondition_Turkel(void);

void Compute_Flux_HLLC_Original(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_Flux_HLLC_Thornber(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_Flux_HLLC_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_Flux_HLLC_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_Flux_HLLC(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);

void Compute_Residual_HLLC(int AddTime);

#endif	/* _HLLC_FLUXES_H */

