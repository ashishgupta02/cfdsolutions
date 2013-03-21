/*******************************************************************************
 * File:        JST_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _JST_FLUXES_H
#define	_JST_FLUXES_H

#include "Vector3D.h"

// JST Scheme Functions
void JST_Init(void);
void JST_Finalize(void);
void JST_Reset(void);
void Compute_JSTFlux(int node_L, int node_R, Vector3D areavec, double *Flux_JST, int AddTime);
void Compute_Residual_JST(int AddTime);

#endif	/* _JST_FLUXES_H */

