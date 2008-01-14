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
 * File		ReferenceQuantities.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Algebra.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPInitializeSolution_2D.h"
#include "FluidProperties.h"
#include "ReferenceQuantities.h"

/* Externs from ReferenceQuantities.h */
int RefFlag;
double MachNo_Ref = 0.0;
double AoA_Ref = 0.0;
double VelocityX_Ref = 0.0;
double VelocityY_Ref = 0.0;
double VelocityZ_Ref = 0.0;
double Density_Ref = 0.0;
double Pressure_Ref = 0.0;
double SpeedSound_Ref = 0.0;
double Temperature_Ref = 0.0;
double VelocityMagnitude_Ref = 0.0;

/*--------------------------------------------------------------*/
static void ERROR1 (const char *errmsg) {
	char cmd[129];

	sprintf (cmd, "error_exit {%s}", errmsg);
	printf ("%s\n", cmd);
	exit (1);
}

/*--------------------------------------------------------------*/
static int CalculateAoA_Ref (void) {
	if ((VelocityX_Ref >= 0.0) && (VelocityY_Ref >= 0.0)) {
		 if ((VelocityX_Ref == 0.0) && (VelocityY_Ref == 0.0))
		 	ERROR1("Invalid Reference Velocities");
		
		if (VelocityX_Ref == 0.0)
			AoA_Ref = 90.0;
		else
			AoA_Ref = (180/PI) * atan(VelocityY_Ref/VelocityX_Ref);
	}
	
	if (VelocityX_Ref < 0.0)
		AoA_Ref = 180.0 + ((180/PI) * atan(VelocityY_Ref/VelocityX_Ref));
	
	if ((VelocityX_Ref >= 0.0) && (VelocityY_Ref < 0.0))
		AoA_Ref = 360.0 - ((180/PI) * atan(VelocityY_Ref/VelocityX_Ref));
	
	return 0;
}

/*--------------------------------------------------------------*/
int CalculateReferenceQuantities (int flag) {
	switch (flag) {
		case 111:
			/* Available:
				Mach No, AoA
				Density, Temperature */
			
			/* Calculate Pressure_Ref */
			Pressure_Ref = GasConstant_Const * (Density_Ref * Temperature_Ref);
			
			/* Calulate Speed Sound */
			SpeedSound_Ref = sqrt((SpecificHeatRatio_Const * (GasConstant_Const * Temperature_Ref)));
			
			/* Calculate VelocityX and VelocityY */
			VelocityX_Ref = (MachNo_Ref * SpeedSound_Ref) * cos((PI/180) * AoA_Ref);
			VelocityY_Ref = (MachNo_Ref * SpeedSound_Ref) * sin((PI/180) * AoA_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			break;
		case 112:
			/* Available:
				Mach No, AoA
				Density, SpeedSound */
			
			/* Calculate Temperature */
			Temperature_Ref = pow(SpeedSound_Ref, 2)/(SpecificHeatRatio_Const * GasConstant_Const);
			
			/* Calculate Pressure */
			Pressure_Ref = GasConstant_Const * (Density_Ref * Temperature_Ref);
			
			/* Calculate VelocityX and VelocityY */
			VelocityX_Ref = (MachNo_Ref * SpeedSound_Ref) * cos((PI/180) * AoA_Ref);
			VelocityY_Ref = (MachNo_Ref * SpeedSound_Ref) * sin((PI/180) * AoA_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
						
			break;
		case 121:
			/* Available:
				Mach No, AoA
				Pressure, Temperature */
			
			/* Calculate Density */
			Density_Ref = Pressure_Ref/(GasConstant_Const * Temperature_Ref);
			
			/* Calculate Speed Sound */
			SpeedSound_Ref = sqrt((SpecificHeatRatio_Const * (GasConstant_Const * Temperature_Ref)));
			
			/* Calculate VelocityX and VelocityY */
			VelocityX_Ref = (MachNo_Ref * SpeedSound_Ref) * cos((PI/180) * AoA_Ref);
			VelocityY_Ref = (MachNo_Ref * SpeedSound_Ref) * sin((PI/180) * AoA_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			break;
		case 122:
			/* Available:
				Mach No, AoA
				Pressure, SpeedSound */
			
			/* Calculate Temperature */
			Temperature_Ref = pow(SpeedSound_Ref, 2)/(SpecificHeatRatio_Const * GasConstant_Const);
			
			/* Calculate Density */
			Density_Ref = Pressure_Ref/(GasConstant_Const * Temperature_Ref);

			/* Calculate VelocityX and VelocityY */
			VelocityX_Ref = (MachNo_Ref * SpeedSound_Ref) * cos((PI/180) * AoA_Ref);
			VelocityY_Ref = (MachNo_Ref * SpeedSound_Ref) * sin((PI/180) * AoA_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			break;
		case 211:
			/* Available:
				VelocityX, VelocityY
				Density, Temperature */
			
			/* Calculate Pressure */
			Pressure_Ref = GasConstant_Const * (Density_Ref * Temperature_Ref);
			
			/* Calculate Speed Sound */
			SpeedSound_Ref = sqrt((SpecificHeatRatio_Const * (GasConstant_Const * Temperature_Ref)));
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			/* Calculate Mach Number */
			MachNo_Ref = VelocityMagnitude_Ref/SpeedSound_Ref;
			
			/* Calculate Angle of Attack */
			CalculateAoA_Ref();
			break;
		case 212:
			/* Available:
				VelocityX, VelocityY
				Density, SpeedSound */
			
			/* Calculate Temperature */
			Temperature_Ref = pow(SpeedSound_Ref, 2)/(SpecificHeatRatio_Const * GasConstant_Const);
			
			/* Calculate Pressure */
			Pressure_Ref = GasConstant_Const * (Density_Ref * Temperature_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			/* Calculate Mach No */
			MachNo_Ref = VelocityMagnitude_Ref/SpeedSound_Ref;
			
			/* Calculate Angle of Attack */
			CalculateAoA_Ref();
			break;
		case 221:
			/* Available:
				VelocityX, VelocityY
				Pressure, Temperature */
			
			/* Calculate Density */
			Density_Ref = Pressure_Ref/(GasConstant_Const * Temperature_Ref);
			
			/* Calculate Speed Sound */
			SpeedSound_Ref = sqrt((SpecificHeatRatio_Const * (GasConstant_Const * Temperature_Ref)));

			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			/* Calculate Mach No */
			MachNo_Ref = VelocityMagnitude_Ref/SpeedSound_Ref;
			
			/* Calculate Angle of Attack */
			CalculateAoA_Ref();
			break;
		case 222:
			/* Available:
				VelocityX, VelocityY
				Pressure, SpeedSound */
			
			/* Calculate Temperature */
			Temperature_Ref = pow(SpeedSound_Ref, 2)/(SpecificHeatRatio_Const * GasConstant_Const);
			
			/* Calculate Density */
			Density_Ref = Pressure_Ref/(GasConstant_Const * Temperature_Ref);
			
			/* Calculate VelocityMagnitude */
			VelocityMagnitude_Ref = sqrt((pow(VelocityX_Ref, 2) + pow(VelocityY_Ref, 2)));
			
			/* Calculate Mach No */
			MachNo_Ref = VelocityMagnitude_Ref/SpeedSound_Ref;
			
			/* Calculate Angle of Attack */
			CalculateAoA_Ref();
			break;
		default:
			ERROR1("Invalid Option");
	}
	return 0;
}

/*--------------------------------------------------------------*/
int GetReferenceQuantities (void) {
	int flag, flag1;
	
	printf("Reference Quantities:\n");
	printf("Give either Ref:\n\t1) Mach No / Angle of Attack\n\t2)Velocities\n");
	printf("Select Option (1/2) = ");
	scanf("%d", &flag);
	
	switch (flag) {
		case 1:
			/* Mach No and AoA */
			printf("\tMach No Reference = ");
			scanf("%lf", &MachNo_Ref);
			
			printf("\tAngle of Attack Reference = ");
			scanf("%lf", &AoA_Ref);
			
			break;
		case 2:
			/* Velocities */
			printf("\tVelocityX Reference = ");
			scanf("%lf", &VelocityX_Ref);
			
			printf("\tVelocityY Reference = ");
			scanf("%lf", &VelocityY_Ref);
			
			break;
		default:
			ERROR1("Invalid Option");
	}
	
	flag1 = flag;
	
	printf("Give either Ref: 1) Density or 2) Pressure : (1/2): ");
	scanf("%d", &flag);
	
	switch (flag) {
		case 1:
			/* Density Reference */
			printf("\tDensity Reference = ");
			scanf("%lf", &Density_Ref);
			break;
		case 2:
			/* Pressure Reference */
			printf("\tPressure Reference = ");
			scanf("%lf", &Pressure_Ref);
			break;
		default:
			ERROR1("Invalid Option");
	}
	
	flag1 = (flag1 * 10) + flag ;
	
	printf("Give either Ref: 1) Temperature or 2) Speed of Sound : (1/2): ");
	scanf("%d", &flag);
	
	switch (flag) {
		case 1:
			/* Temperature */
			printf("\tTemperature Reference = ");
			scanf("%lf", &Temperature_Ref);
			break;
		case 2:
			/* Speed of Sound */
			printf("\tSpeed of Sound Reference = ");
			scanf("%lf", &SpeedSound_Ref);
			break;
		default:
			ERROR1("Invalid Option");
	}
	
	flag1 = (flag1 * 10) + flag;
	
	CalculateReferenceQuantities(flag1);
	
	RefFlag = 1;
	return 0;
}

/*--------------------------------------------------------------*/
/* CheckRefFlag :
	0 = Not initialized
	1 = Initialized
*/
int CheckRefFlag (void) {
	return (RefFlag);
}
