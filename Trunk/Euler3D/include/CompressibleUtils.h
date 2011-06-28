/*******************************************************************************
 * File:        CompressibleUtils.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef COMPRESSIBLEUTILS_H
#define	COMPRESSIBLEUTILS_H

double ConservativeGetPressure(double *Qc, double gamma);
double ConservativeGetTemperature(double *Qc, double gamma);

void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma);
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma);
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma);

void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma);
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma);
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma);

#endif	/* COMPRESSIBLEUTILS_H */

