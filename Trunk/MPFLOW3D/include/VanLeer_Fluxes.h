/*******************************************************************************
 * File:        VanLeer_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _VANLEER_FLUXES_H
#define	_VANLEER_FLUXES_H

#include "Vector3D.h"

// Van Leer Scheme Functions
void VanLeer_Init(void);
void VanLeer_Finalize(void);
void VanLeer_Reset(void);
void Compute_VanLeerFlux(int node_L, int node_R, Vector3D areavec, double *Flux_VanLeer, int AddTime);
void Compute_Residual_VanLeer(int AddTime);

#endif	/* _VANLEER_FLUXES_H */

