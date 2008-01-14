/*[
 * Copyright 2006   Ashish Gupta
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
 * with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

/*
 * File		PPScalarFlowVariables_2D.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/


#ifndef __PPScalarFlowVariables_2D_H__
#define __PPScalarFlowVariables_2D_H__

int CalculateDensity_2D (void);
int CalculatePressure_2D (void);
int CalculateTemperature_2D (void);
int CalculateSpeedSound_2D (void);
int CalculateMachNumber_2D (void);
int CalculateStagnationDensity_2D (void);
int CalculateStagnationPressure_2D (void);
int CalculateStagnationTemperature_2D (void);
int CalculatePressureCoefficient_2D (void);
int CalculateStagnationPressureCoef_2D (void);
int CalculatePitotPressure_2D (void);
int CalculatePitotPressureRatio_2D (void);
int CalculateDynamicPressure_2D (void);
int CalculateEnthalpy_2D (void);
int CalculateStagnationEnthalpy_2D (void);
int CalculateInternalEnergy_2D (void);
int CalculateStagnationEnergy_2D (void);
int CalculateStagnationEnergyDensity_2D (void);
int CalculateKineticEnergy_2D (void);
int CalculateVelocityX_2D (void);
int CalculateVelocityY_2D (void);
int CalculateVelocityMagnitude_2D (void);
int CalculateEquivalentPotentialVelocityRatio_2D (void);
int CalculateMomentumX_2D (void);
int CalculateMomentumY_2D (void);
int CalculateMomentumMagnitude_2D (void);
int CalculateEntropy_2D (void);
int CalculateEntropyMeasureS1_2D (void);
int CalculateVorticityX_2D (void);
int CalculateVorticityY_2D (void);
int CalculateVorticityMagnitude_2D (void);
int CalculateDensityGradientMagnitude_2D (void);
int CalculatePressureGradientMagnitude_2D (void);
int CalculateTemperatureGradientMagnitude_2D (void);
int CalculateDivergenceVelocity_2D (void);
int CalulateIsentropicDensityRatio_2D (void);
int CalculateIsentropicPressureRatio_2D (void);
int CalculateIsentropicTemperatureRatio_2D (void);
int CalculateShock_2D (void);
int CalculateFilteredShock_2D (void);
int CalculateSwirl_2D (void);
int CalculateHelicity_2D (void);
int CalculateRelativeHelicity_2D (void);
int CalculateFilteredRelativeHelicity_2D (void);
int CalculateSutherlandsLaw_2D (void);

int ScalarFlowVariableOptions_2D (void);

#endif

