/*******************************************************************************
 * File:        Osher_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _OSHER_FLUXES_H
#define	_OSHER_FLUXES_H

#include "Vector3D.h"

// Osher Scheme Functions
void Osher_Init(void);
void Osher_Finalize(void);
void Osher_Reset(void);
void Compute_OsherFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Osher, int AddTime);
void Compute_Residual_Osher(int AddTime);

#endif	/* _OSHER_FLUXES_H */

