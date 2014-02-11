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
 * File		PPInitializeSolution_2D.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */


#ifndef __PPInitializeSolution_2D_H__
#define __PPInitializeSolution_2D_H__

#include "CGNSIO.h"

/* User Defined */

typedef struct PP_2D_Data_Un {
    char Name[50];
    int Size;
    double *Data;
} Data_2D_Un;

typedef struct PP_2D_SolutionCGNS_Un {
    unsigned int Size;
    Data_2D_Un *SolCGNS;
} SolutionCGNS_2D_Un;

typedef struct PP_2D_Solution_Un {
    int Location;
    int Size;
    Data_2D_Un **Sols;
} Solution_2D_Un;

/* This Variable is used to store the data present in CGNS file */
extern SolutionCGNS_2D_Un SolutionCGNS2D;

/* This Variable is used to store the data pointers of both the solution
        1) CGNS File Data 
        2) Computed Data
 */
extern Solution_2D_Un Solution2D;

/* Solution Data Declearations */
/* This are the data which are supported by PostCFD */
extern Data_2D_Un Density2D;
extern Data_2D_Un Pressure2D;
extern Data_2D_Un Temperature2D;
extern Data_2D_Un EnergyInternal2D;

extern Data_2D_Un Enthalpy2D;
extern Data_2D_Un Entropy2D;

extern Data_2D_Un DensityStagnation2D;
extern Data_2D_Un PressureStagnation2D;
extern Data_2D_Un TemperatureStagnation2D;
extern Data_2D_Un EnergyStagnation2D;
extern Data_2D_Un EnthalpyStagnation2D;
extern Data_2D_Un EnergyStagnationDensity2D;

extern Data_2D_Un VelocityX2D;
extern Data_2D_Un VelocityY2D;
extern Data_2D_Un VelocityZ2D;
extern Data_2D_Un VelocityMagnitude2D;

extern Data_2D_Un VelocitySound2D;
extern Data_2D_Un VelocitySoundStagnation2D;

extern Data_2D_Un MomentumX2D;
extern Data_2D_Un MomentumY2D;
extern Data_2D_Un MomentumZ2D;
extern Data_2D_Un MomentumMagnitude2D;

extern Data_2D_Un EnergyKinetic2D;
extern Data_2D_Un PressureDynamic2D;

extern Data_2D_Un VorticityX2D;
extern Data_2D_Un VorticityY2D;
extern Data_2D_Un VorticityZ2D;
extern Data_2D_Un VorticityMagnitude2D;

extern Data_2D_Un VelocityAngleX2D;
extern Data_2D_Un VelocityAngleY2D;
extern Data_2D_Un VelocityAngleZ2D;

extern Data_2D_Un VelocityUnitVectorX2D;
extern Data_2D_Un VelocityUnitVectorY2D;
extern Data_2D_Un VelocityUnitVectorZ2D;

extern Data_2D_Un SutherlandLawConstant2D;

extern Data_2D_Un Mach2D;

extern Data_2D_Un CoefPressure2D;

/* Specific to PostCFD */
extern Data_2D_Un DensityGradientX2D;
extern Data_2D_Un DensityGradientY2D;
extern Data_2D_Un DensityGradientZ2D;
extern Data_2D_Un DensityGradientMagnitude2D;

extern Data_2D_Un PressureGradientX2D;
extern Data_2D_Un PressureGradientY2D;
extern Data_2D_Un PressureGradientZ2D;
extern Data_2D_Un PressureGradientMagnitude2D;

extern Data_2D_Un PitotPressure2D;
extern Data_2D_Un PitotPressureRatio2D;

extern Data_2D_Un TemperatureGradientX2D;
extern Data_2D_Un TemperatureGradientY2D;
extern Data_2D_Un TemperatureGradientZ2D;
extern Data_2D_Un TemperatureGradientMagnitude2D;

extern Data_2D_Un PerturbationVelocityX2D;
extern Data_2D_Un PerturbationVelocityY2D;
extern Data_2D_Un PerturbationVelocityZ2D;

extern Data_2D_Un VelocityDivergence2D;

/* Function Declearattions */
void DummySolutions_2D(void);
void SolutionMinMaxVector_2D(int, int, double *, double *);
void SolutionMinMax_2D(int, double *, double *);
int SolutionCC2NodeSelect_2D(void);
int InitializeSolutionCGNS_2D(ZONE *);
int InitializeSolution_2D(void);
int GetSolutionList_2D(void);
int CheckSolutionFieldExistance_2D(const char *);
int GetSolutionFieldLocation_2D(const char *);
int UpdateSolution_2D(Data_2D_Un *);

#endif
