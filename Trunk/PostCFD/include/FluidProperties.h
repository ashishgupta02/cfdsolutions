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
 * File		FluidProperties.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/


#ifndef __FluidProperties_H__
#define __FluidProperties_H__

#include "PPInitializeSolution_2D.h"

#define	InCompressible	0
#define Compressible	1

/* FluidFlag = 0 => Incompressible, FluidFlag = 1 => Compressible */
extern int FluidFlag;

/* NOTE: Use variables
When Incompressible:
	Always Constant: Specific Density
	Variable or Constant: Specific Heat, Viscosity, Conductivity
When Compressible:
	Variable or Constant: Gas Constant, SpecificHeatRatio, Viscosity, Conductivity
*/

/* If Properties are constant */
extern double GasConstant_Const, SpecificDensity_Const, Viscosity_Const;
extern double SpecificHeatRatio_Const, SpecificHeat_Const, Conductivity_Const;

/* Properties Flag */
extern int GasConstant_Flag, Viscosity_Flag;
extern int SpecificHeatRatio_Flag, SpecificHeat_Flag, Conductivity_Flag;

/* If Properties is variable */
extern Data_2D_Un *GasConstant_Var;
extern Data_2D_Un *Viscosity_Var;
extern Data_2D_Un *SpecificHeatRatio_Var;
extern Data_2D_Un *SpecificHeat_Var;
extern Data_2D_Un *Conductivity_Var;

/* Functions Declearation */
int GetFluidFlag (void);
int GetSpecificDensityConst (void);
int GetSpecificHeat (void);
int GetViscosity (void);
int GetConductivity (void);
int GetGasConstant (void);
int GetSpecificHeatRatio (void);
int GetIncompressibleFluidProperties (void);
int GetCompressibleFluidProperties (void);
int GetFluidProperties (void);

#endif

