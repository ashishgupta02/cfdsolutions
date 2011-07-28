/*******************************************************************************
 * File:        CompressibleUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _COMPRESSIBLEUTILS_H
#define	_COMPRESSIBLEUTILS_H

double ConservativeGetPressure(double *Qc, double gamma);
double ConservativeGetTemperature(double *Qc, double gamma);

void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma);
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma);
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma);

void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma);
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma);
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma);

void ConservativeEulerFlux(double *Qc, Vector3D areavec, double *Flux);
void ConservativeEulerFluxJacobian(double *Qc, Vector3D areavec, double **Jacobian, double gamma);

#endif	/* _COMPRESSIBLEUTILS_H */

