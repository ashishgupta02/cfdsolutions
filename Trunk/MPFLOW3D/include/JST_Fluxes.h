/*******************************************************************************
 * File:        JST_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _JST_FLUXES_H
#define	_JST_FLUXES_H

#include "Vector3D.h"

// JST Scheme Functions
void JST_Init(void);
void JST_Finalize(void);
void JST_Reset(void);
void Compute_Flux_JST(int node_L, int node_R, Vector3D areavec, double *Flux_JST_Conv, double *Flux_JST_Diss, int AddTime); 
void Compute_Residual_JST(int AddTime);

#endif	/* _JST_FLUXES_H */

