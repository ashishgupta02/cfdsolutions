/*******************************************************************************
 * File:        LDFSS_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _LDFSS_FLUXES_H
#define	_LDFSS_FLUXES_H

#include "Vector3D.h"

// LDFSS Scheme Functions
void LDFSS_Init(void);
void LDFSS_Finalize(void);
void LDFSS_Reset(void);
void Compute_LDFSSFlux(int node_L, int node_R, Vector3D areavec, double *Flux_LDFSS, int AddTime);
void Compute_Residual_LDFSS(int AddTime);

#endif	/* _LDFSS_FLUXES_H */

