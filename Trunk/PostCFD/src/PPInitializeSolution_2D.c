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
 * File		PPInitializeSolution_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PrePostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Interpolation_2D.h"
#include "PostAnalysis_2D.h"

/* For CellCenter to Node Convertion */
static int CC2NodeFlag = 0;
static int SolutionCGNSFlag = 0;
static int SolutionFlag = 0;

SolutionCGNS_2D_Un SolutionCGNS2D;
Solution_2D_Un Solution2D;

/* Defination of externs of PPInitializeSolution_2D.h*/
/* Solution Data Definations */
/* This are the data which are supported by PostCFD */
Data_2D_Un Density2D;
Data_2D_Un Pressure2D;
Data_2D_Un Temperature2D;
Data_2D_Un EnergyInternal2D;

Data_2D_Un Enthalpy2D;
Data_2D_Un Entropy2D;

Data_2D_Un DensityStagnation2D;
Data_2D_Un PressureStagnation2D;
Data_2D_Un TemperatureStagnation2D;
Data_2D_Un EnergyStagnation2D;
Data_2D_Un EnthalpyStagnation2D;
Data_2D_Un EnergyStagnationDensity2D;

Data_2D_Un VelocityX2D;
Data_2D_Un VelocityY2D;
Data_2D_Un VelocityZ2D;
Data_2D_Un VelocityMagnitude2D;

Data_2D_Un VelocitySound2D;
Data_2D_Un VelocitySoundStagnation2D;

Data_2D_Un MomentumX2D;
Data_2D_Un MomentumY2D;
Data_2D_Un MomentumZ2D;
Data_2D_Un MomentumMagnitude2D;

Data_2D_Un EnergyKinetic2D;
Data_2D_Un PressureDynamic2D;

Data_2D_Un VorticityX2D;
Data_2D_Un VorticityY2D;
Data_2D_Un VorticityZ2D;
Data_2D_Un VorticityMagnitude2D;

Data_2D_Un VelocityAngleX2D;
Data_2D_Un VelocityAngleY2D;
Data_2D_Un VelocityAngleZ2D;

Data_2D_Un VelocityUnitVectorX2D;
Data_2D_Un VelocityUnitVectorY2D;
Data_2D_Un VelocityUnitVectorZ2D;

Data_2D_Un SutherlandLawConstant2D;

Data_2D_Un Mach2D;

Data_2D_Un CoefPressure2D;

Data_2D_Un DensityGradientX2D;
Data_2D_Un DensityGradientY2D;
Data_2D_Un DensityGradientZ2D;
Data_2D_Un DensityGradientMagnitude2D;

Data_2D_Un PressureGradientX2D;
Data_2D_Un PressureGradientY2D;
Data_2D_Un PressureGradientZ2D;
Data_2D_Un PressureGradientMagnitude2D;

Data_2D_Un PitotPressure2D;
Data_2D_Un PitotPressureRatio2D;

Data_2D_Un TemperatureGradientX2D;
Data_2D_Un TemperatureGradientY2D;
Data_2D_Un TemperatureGradientZ2D;
Data_2D_Un TemperatureGradientMagnitude2D;

Data_2D_Un PerturbationVelocityX2D;
Data_2D_Un PerturbationVelocityY2D;
Data_2D_Un PerturbationVelocityZ2D;

Data_2D_Un VelocityDivergence2D;

/* Dummy Solutions */
/*--------------------------------------------------------------*/
void DummySolutions_2D(void) {
	int i, select;
	
	GetSolutionList_2D();
	printf("Select Solution = ");
	scanf("%d", &select);
	switch (Solution2D.Location) {
		case Vertex:
			for (i = 0; i < NoNodes2D; i++) {
				Solution2D.Sols[select]->Data[i] = 0.0192 * Node2D[i].Coordinate[0] + 3.84 * pow(Node2D[i].Coordinate[1] - 0.5 ,2) + 0.2;
			}
			break;
		case CellCenter:
			for (i = 0; i < NoCells2D; i++) {
				Solution2D.Sols[select]->Data[i] = (27 * (1 - (4 * pow(Cell2D[i].Centroid[1] - 0.5, 2))) + 3)*(1 - Cell2D[i].Centroid[0]/90);
			}
			break;
	}
}

/*--------------------------------------------------------------*/
void SolutionMinMaxVector_2D (int SolLocation1, int SolLocation2, double *MinValueSol, double *MaxValueSol) {
	int i;
	double Sol;
	
	switch (Solution2D.Location) {
		case Vertex:
			*MinValueSol = sqrt((Solution2D.Sols[SolLocation1]->Data[0] * Solution2D.Sols[SolLocation1]->Data[0])
						+ (Solution2D.Sols[SolLocation2]->Data[0] * Solution2D.Sols[SolLocation2]->Data[0]));
			*MaxValueSol = *MinValueSol;
			for (i = 1; i < NoNodes2D; i++) {
				Sol = sqrt((Solution2D.Sols[SolLocation1]->Data[i] * Solution2D.Sols[SolLocation1]->Data[i])
						+ (Solution2D.Sols[SolLocation2]->Data[i] * Solution2D.Sols[SolLocation2]->Data[i]));
				if (*MinValueSol > Sol)
					*MinValueSol = Sol;
				
				if (*MaxValueSol < Sol)
					*MaxValueSol = Sol;
			}
			break;
		case CellCenter:
			*MinValueSol = sqrt((Solution2D.Sols[SolLocation1]->Data[0] * Solution2D.Sols[SolLocation1]->Data[0])
						+ (Solution2D.Sols[SolLocation2]->Data[0] * Solution2D.Sols[SolLocation2]->Data[0]));
			*MaxValueSol = *MinValueSol;
			for (i = 1; i < NoCells2D; i++) {
				Sol = sqrt((Solution2D.Sols[SolLocation1]->Data[i] * Solution2D.Sols[SolLocation1]->Data[i])
						+ (Solution2D.Sols[SolLocation2]->Data[i] * Solution2D.Sols[SolLocation2]->Data[i]));
				if (*MinValueSol > Sol)
					*MinValueSol = Sol;
				
				if (*MaxValueSol < Sol)
					*MaxValueSol = Sol;			
			}
			break;
	}
}

/*--------------------------------------------------------------*/
void SolutionMinMax_2D (int SolLocation, double *MinValueSol, double *MaxValueSol) {
	int i;
	
	switch (Solution2D.Location) {
		case Vertex:
			*MinValueSol = Solution2D.Sols[SolLocation]->Data[0];
			*MaxValueSol = Solution2D.Sols[SolLocation]->Data[0];
			for (i = 1; i < NoNodes2D; i++) {
				if (*MinValueSol > Solution2D.Sols[SolLocation]->Data[i])
					*MinValueSol = Solution2D.Sols[SolLocation]->Data[i];
				
				if (*MaxValueSol < Solution2D.Sols[SolLocation]->Data[i])
					*MaxValueSol = Solution2D.Sols[SolLocation]->Data[i];
			}
			break;
		case CellCenter:
			*MinValueSol = Solution2D.Sols[SolLocation]->Data[0];
			*MaxValueSol = Solution2D.Sols[SolLocation]->Data[0];
			for (i = 1; i < NoCells2D; i++) {
				if (*MinValueSol > Solution2D.Sols[SolLocation]->Data[i])
					*MinValueSol = Solution2D.Sols[SolLocation]->Data[i];
				
				if (*MaxValueSol < Solution2D.Sols[SolLocation]->Data[i])
					*MaxValueSol = Solution2D.Sols[SolLocation]->Data[i];
			}
			break;
	}
}

/*---------------------------------------------------------------*/
int SolutionCC2NodeSelect_2D(void) {
	printf("Solution Location : Cell-Center\n");
	printf("\t Convert Solution to Vertex (0/1):");
	scanf("%d", &CC2NodeFlag);
	
	switch (CC2NodeFlag) {
		case 0:
			/* No Interpolation */
			break;
		case 1:
			/* Interpolation */
			SetInterpolationOptionCC2Node_2D ();
			break;
		default:
			/* Error */
			CC2NodeFlag = 0;
			MSG("Invalid Option: Continuing with no interpolation");
	}
	return 0;
}

/*---------------------------------------------------------------*/
int InitializeSolutionCGNS_2D (ZONE *P) {
	int i, ns, nf;
	SOLUTION *sol;
	FIELD *fld;
	Data_2D_Un tmp;
	
	for (sol = P->sols, ns = 1; ns <= P->nsols; sol++, ns++) {
		if (ns >= 2) {
			MSG("InitializeSolutionCGNS_2D: Only Single Solution Node Supported");
			return 1;
		}
		
		if (SolutionCGNSFlag == 0) {
			/* Ask to convert solution from CellCenter to Vertex */
			CC2NodeFlag = 0;
			if (sol->location == CellCenter) {
				SolutionCC2NodeSelect_2D();
				
				if (CC2NodeFlag == 0)
					Solution2D.Location = CellCenter;
				
				if (CC2NodeFlag == 1)
					Solution2D.Location = Vertex;
			}
			else if (sol->location == Vertex)
				Solution2D.Location = Vertex;
				
			SolutionCGNS2D.Size = sol->nflds;
			
			/* No of Initial Solution Pointers */
			SolutionCGNS2D.SolCGNS = (Data_2D_Un *)malloc(sol->nflds * sizeof(Data_2D_Un));
			if (SolutionCGNS2D.SolCGNS == NULL) {
				MSG("InitializeSolutionCGNS_2D: Memory Failed");
				return 1;
			}
			
			switch (sol->location) {
				case Vertex:
					/* Solution Field Location Vertex */				
					for (fld = sol->flds, nf = 1; nf <= sol->nflds; fld++, nf++){
						/* Copy the solution field name */
						strcpy(SolutionCGNS2D.SolCGNS[nf-1].Name, fld->name);
						
						/* Specify the size of data */
						SolutionCGNS2D.SolCGNS[nf-1].Size = NoNodes2D;
						
						/* Allocate the memory */
						SolutionCGNS2D.SolCGNS[nf-1].Data = (double *) malloc (NoNodes2D * sizeof(double));
						if (SolutionCGNS2D.SolCGNS[nf-1].Data == NULL) {
							MSG("InitializeSolutionCGNS_2D: Memory Failed");
							return 1;
						}
						
						/* Storing the solution field data */
						for (i = 0; i < NoNodes2D; i++)
							SolutionCGNS2D.SolCGNS[nf-1].Data[i] = fld->data[i];
					}
					break;
				case CellCenter:
					/* Solution Field Location Cell-Center */
					for (fld = sol->flds, nf = 1; nf <= sol->nflds; fld++, nf++) {					
						/* Copy the solution field name */
						strcpy(SolutionCGNS2D.SolCGNS[nf-1].Name, fld->name);
						
						if (CC2NodeFlag == 0) {
							/* For Data at Cell-Center */
							/* Specify the size of data */
							SolutionCGNS2D.SolCGNS[nf-1].Size = NoCells2D;
							/* Allocate the memory */
							SolutionCGNS2D.SolCGNS[nf-1].Data = (double *) malloc (NoCells2D * sizeof(double));
							if (SolutionCGNS2D.SolCGNS[nf-1].Data == NULL) {
								MSG("InitializeSolutionCGNS_2D: Memory Failed");
								return 1;
							}
							
							/* Storing the solution field data */
							for (i = 0; i < NoCells2D; i++)
								SolutionCGNS2D.SolCGNS[nf-1].Data[i] = fld->data[i];
								
						}
						else if (CC2NodeFlag == 1) {
							/* For Interpolated Cell-Center Data to Vertex */
							/* Specify the size of data */
							SolutionCGNS2D.SolCGNS[nf-1].Size = NoNodes2D;
							
							/* Allocate the memory */
							SolutionCGNS2D.SolCGNS[nf-1].Data = (double *) malloc (NoNodes2D * sizeof(double));
							if (SolutionCGNS2D.SolCGNS[nf-1].Data == NULL) {
								MSG("InitializeSolutionCGNS_2D: Memory Failed");
								return 1;
							}
							
							/* Create a temporary data */
							tmp.Data = (double *) malloc (NoCells2D * sizeof(double));
							if (tmp.Data == NULL) {
								MSG("InitializSolutionCGNS_2D: Memory Failed");
								return 1;
							}
							
							/* Initializing the tmp variable */
							for (i = 0; i < NoCells2D; i++)
								tmp.Data[i] = fld->data[i];
							
							/* Interpolation of Data */
							CellCenter2Node_2D (&tmp, &SolutionCGNS2D.SolCGNS[nf-1]);
							
							/* Free tmp memory */
							free(tmp.Data);
						}
					}
					break;
				default:
					/* Unsupported Solution Location */
					MSG("Unsupported Solution Location");
					return 1;
			}
			/* To avoid Reinitialization of Solution */
			SolutionCGNSFlag = 1;
		}
	}
	
	return 0;
}

/*---------------------------------------------------------------*/
int InitializeSolution_2D (void) {
	int i;
	
	if (SolutionFlag == 0) {

		Solution2D.Size = SolutionCGNS2D.Size;
		
		Solution2D.Sols = (Data_2D_Un **) malloc (Solution2D.Size * sizeof(Data_2D_Un *));
		if (Solution2D.Sols == NULL) {
			MSG("IntializeSolution_2D: Memory Failed");
			return 1;
		}
		
		for (i = 0; i < Solution2D.Size; i++)
			Solution2D.Sols[i] = &SolutionCGNS2D.SolCGNS[i];
		
		SolutionFlag = 1;
	}
	
	return 0;
}

/*---------------------------------------------------------------*/
int GetSolutionList_2D (void) {
	int i;
	
	printf("Solution Available:\n");
	
	for(i = 0; i < Solution2D.Size; i++)
		printf("\t%d.%s\n", i, Solution2D.Sols[i]->Name);
	
	return 0;
}

/*---------------------------------------------------------------*/
/* Function: CheckSolutionFieldExistance_2D
	Input: name string
	Output: 1 if present else 0 if not present 
*/
int CheckSolutionFieldExistance_2D (const char *str) {
	int i, result;
	
	result = 0;
	
	for (i = 0; i < Solution2D.Size; i++) {
		if (!strcmp(Solution2D.Sols[i]->Name, str)) { 
			result = 1;
			return (result);
		}
	}		
	
	return (result);
}

/*---------------------------------------------------------------*/
/* Function: GetSolutionFieldLocation_2D
	Input: name string
	Output: Solution2D array index if present 
		else -1 if not present
*/
int GetSolutionFieldLocation_2D (const char *str) {
	int i;
	
	for (i = 0; i < Solution2D.Size; i++) {
		if (!strcmp(Solution2D.Sols[i]->Name, str))
			return (i);
	}
	
	return (-1);
}

/*---------------------------------------------------------------*/
/* Function: UpdateSolution_2D
	Input: pointer to new solution
*/
int UpdateSolution_2D (Data_2D_Un *newsol) {
	int i, size;
	Data_2D_Un **tmp;
	
	size = Solution2D.Size;
	
	/* Creat a tmp pointer to store pointers */
	tmp = (Data_2D_Un **) malloc (size * sizeof(Data_2D_Un *));
	if (tmp == NULL) {
		MSG("UpdateSolution_2D: Memory Failed: 1");
		return 1;
	}
	
	/* Storing pointers location */
	for (i = 0; i < size; i++)
		tmp[i] = Solution2D.Sols[i];
	
	/* Free Solution Sols Array */
	free (Solution2D.Sols);
	
	/* Increase Size by One */
	Solution2D.Size = Solution2D.Size + 1;
	
	/* ReAllocation of Solution2D.Sols */
	Solution2D.Sols = (Data_2D_Un **) malloc (Solution2D.Size * sizeof(Data_2D_Un *));
	if (Solution2D.Sols == NULL) {
		MSG("UpdateSolution_2D: Memory Failed: 2");
		return 1;
	}
	
	/* Recopying the pointer location */
	for (i = 0; i < size; i++)
		Solution2D.Sols[i] = tmp[i];
	
	/* Copying Location of new solution */
	Solution2D.Sols[size] = newsol;
	
	/* Freeing Tmp variable */
	free(tmp);
	
	return 0;
}
