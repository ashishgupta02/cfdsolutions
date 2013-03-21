/*******************************************************************************
 * File:        CompressibleUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _COMPRESSIBLEUTILS_H
#define	_COMPRESSIBLEUTILS_H

//------------------------------------------------------------------------------
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void Compute_Flux_Euler_Convective(double *Q, Vector3D areavec, double *Flux);
void Compute_Flux_Euler_Convective(int node_ID, Vector3D areavec, double *Flux);

void Compute_Flux_Jacobian_Euler_Convective(double *Q, Vector3D areavec, double **Jacobian_Conv);
void Compute_Flux_Jacobian_Euler_Convective(int node_ID, Vector3D areavec, double **Jacobian_Conv);

void ConservativeEulerFluxJacobian(double *Qc, Vector3D areavec, double **Jacobian, double gamma);

#endif	/* _COMPRESSIBLEUTILS_H */

