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
 * File		PPVectorFlowVariables_2D.c
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
#include "PPInitializeSolution_2D.h"
#include "Gradient_2D.h"
#include "PPVectorFlowVariables_2D.h"
#include "PPSolutionGradients_2D.h"

/* Vector Flow Variable Calculation Option  */
static char *VectorOption[] = {
    "usage  : Vector Flow Variable [options]",
    "options:",
    "	1	= Velocity",
    "	2	= Vorticity",
    "	3	= Momentum",
    "	4	= Perturbation Velocity",
    "	5	= Solution Gradients",
    NULL
};

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static int usage(char **usgmsg) {
    int n;

    for (n = 0; NULL != usgmsg[n]; n++)
        fprintf(stderr, "%s\n", usgmsg[n]);

    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVelocity_2D(void) {

    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePerturbationVelocity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateMomentum_2D(Solution_2D_Un *sol1) {
    sol1 = NULL; /* remove this */
    //	int i, ns;
    //	SOLUTION *sol;

    /* Initialize Velocity Vector */
    //	InitializeSolutionField_2D(P, 1);

    /* Initialize Density */
    //	InitializeSolutionField_2D(P, 2);


    //	for (sol = P->sols, ns = 1; ns <= P->nsols; sol++, ns++) {
    //		if (ns >= 2) {
    //			Warn("Momentum_2D: Only Single Solution Node Supported");
    //			return 1;
    //		}

    //		switch(sol->location) {
    //			case Vertex:
    /* To Calculate Momentum at Node */
    //				MomentumNode = (Data_Vector_2D *) malloc (NoNodes2D * sizeof (Data_Vector_2D));
    //				if (MomentumNode == NULL) {
    //					Warn("Momentum_2D: Memory Failure : 1");
    //					return 1;
    //				}

    //				for (i = 0; i < NoNodes2D; i++) {
    //					MomentumNode[i].Data[0] = VelocityNode[i].Data[0] * DensityNode[i].Data;
    //					MomentumNode[i].Data[1] = VelocityNode[i].Data[1] * DensityNode[i].Data;
    //				} 

    //				break;
    //			case CellCenter:
    //				/* To Calculate Momentun at Cell-Center */
    //				MomentumCC = (Data_Vector_2D *) malloc (NoCells2D * sizeof (Data_Vector_2D));
    //				if (MomentumCC == NULL) {
    //					Warn("Momentum_2D: Memory Failure : 2");
    //					return 1;
    //				}

    //				for (i = 0; i < NoCells2D; i++) {
    //					MomentumCC[i].Data[0] = VelocityCC[i].Data[0] * DensityCC[i].Data;
    //					MomentumCC[i].Data[1] = VelocityCC[i].Data[1] * DensityCC[i].Data;
    //				}

    //				break;
    //			default:
    /* Unsupported Solution Location */
    //				Warn("Momentum_2D: Unsupported Solution Location");
    //				return 1;
    //		}
    //		printf("Momentum Calculation Completed\n");
    //	}

    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVorticity_2D(ZONE *P) {
    P = NULL; /* remove this */
    //	int i, ns;
    //	Data_2D *velocityx;
    //	Data_2D *velocityy;
    //	Data_Vector_2D *gx;
    //	Data_Vector_2D *gy;
    //	SOLUTION *sol;

    /* Initialize Velocity Vector */
    //	InitializeSolutionField_2D(P, 1);

    //	for (sol = P->sols, ns = 1; ns <= P->nsols; sol++, ns++) {
    //		if (ns >= 2) {
    //			Warn("Vorticity_2D: Only Single Solution Node Supported");
    //			return 1;
    //		}

    //		switch(sol->location) {
    //			case Vertex:
    /* To Calculate Vorticity at Node */
    //				velocityx = (Data_2D *) malloc (NoNodes2D * sizeof (Data_2D));
    //				if (velocityx == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 1");
    //					return 1;
    //				}

    //				velocityy = (Data_2D *) malloc (NoNodes2D * sizeof (Data_2D));
    //				if (velocityy == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 2");
    //					return 1;
    //				}

    //				for (i = 0; i < NoNodes2D; i++) {
    //					velocityx[i].Data = VelocityNode[i].Data[0];
    //					velocityy[i].Data = VelocityNode[i].Data[1];
    //				}
    //				
    //				gx = (Data_Vector_2D *) malloc (NoNodes2D * sizeof (Data_Vector_2D));
    //				if (gx == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 3");
    //					return 1;
    //				}

    //				gy = (Data_Vector_2D *) malloc (NoNodes2D * sizeof (Data_Vector_2D));
    //				if (gy == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 4");
    //					return 1;
    //				}

    //				NodeGradientLSM_2D (velocityx, gx);
    //				NodeGradientLSM_2D (velocityy, gy);

    //				VorticityNode = (Data_2D *) malloc (NoNodes2D * sizeof (Data_2D));
    //				if (VorticityNode == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 5");
    //					return 1;
    //				}

    //				for (i = 0; i < NoNodes2D; i++)
    //					VorticityNode[i].Data = (gy[i].Data[0] - gx[i].Data[1]);

    //				break;
    //			case CellCenter:
    /* To Calculate Vorticity at Cell-Center */
    //				velocityx = (Data_2D *) malloc (NoCells2D * sizeof (Data_2D));
    //				if (velocityx == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 6");
    //					return 1;
    //				}

    //				velocityy = (Data_2D *) malloc (NoCells2D * sizeof (Data_2D));
    //				if (velocityy == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 7");
    //					return 1;
    //				}

    //				for (i = 0; i < NoCells2D; i++) {
    //					velocityx[i].Data = VelocityCC[i].Data[0];
    //					velocityy[i].Data = VelocityCC[i].Data[1];
    //				}

    //				gx = (Data_Vector_2D *) malloc (NoCells2D * sizeof (Data_Vector_2D));
    //				if (gx == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 8");
    //					return 1;
    //				}

    //				gy = (Data_Vector_2D *) malloc (NoCells2D * sizeof (Data_Vector_2D));
    //				if (gy == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 9");
    //					return 1;
    //				}

    //				CellGradientLSM_2D (velocityx, gx);
    //				CellGradientLSM_2D (velocityy, gy);

    //				VorticityCC = (Data_2D *) malloc (NoCells2D * sizeof (Data_2D));
    //				if (VorticityCC == NULL) {
    //					Warn("Vorticity_2D: Memory Failure : 10");
    //					return 1;
    //				}

    //				for (i = 0; i < NoCells2D; i++)
    //					VorticityCC[i].Data = (gy[i].Data[0] - gx[i].Data[1]);

    //				break;
    //			default:
    /* Unsupported Solution Location */
    //				Warn("Vorticity_2D: Unsupported Solution Location");
    //				return 1;
    //		}
    //		printf("Vorticity Calculation Completed\n");
    //	}

    return 0;
}

/*---------------------------------------------------------------*/
int VectorFlowVariableOptions_2D(void) {
    int Option;

    printf("Vector Flow Variable Module:\n");
    usage(VectorOption);
    printf("Select Option: ");
    scanf("%d", &Option);

    switch (Option) {
        case 1:
            /* Velocity */
            CalculateVelocity_2D();
            break;
        case 2:
            /* Vorticity */
            //			CalculateVorticity_2D();
            break;
        case 3:
            /* Momentum */
            //			CalculateMomentum_2D();
            break;
        case 4:
            /* Perturbation Velocity */
            CalculatePerturbationVelocity_2D();
            break;
        case 5:
            /* Solution Gradient: Density, Pressure, Temperature */
            SolutionGradient_2D();
            break;
        default:
            /* Invalid Option */
            Warn("Invalid Option");
    }

    return 0;
}
