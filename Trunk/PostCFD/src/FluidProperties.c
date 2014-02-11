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
 * File		FluidProperties.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPInitializeSolution_2D.h"
#include "FluidProperties.h"

/* Externs from FluidProperties.h */
int FluidFlag = -1;

/* If Properties are constant */
double GasConstant_Const, SpecificDensity_Const, Viscosity_Const;
double SpecificHeatRatio_Const, SpecificHeat_Const, Conductivity_Const;

/* Properties Flag */
int GasConstant_Flag, Viscosity_Flag;
int SpecificHeatRatio_Flag, SpecificHeat_Flag, Conductivity_Flag;

/* If Properties is variable */
Data_2D_Un *GasConstant_Var;
Data_2D_Un *Viscosity_Var;
Data_2D_Un *SpecificHeatRatio_Var;
Data_2D_Un *SpecificHeat_Var;
Data_2D_Un *Conductivity_Var;

/*--------------------------------------------------------------*/
static void ERROR(char *errmsg) {
    char cmd[129];

    sprintf(cmd, "error_exit {%s}", errmsg);
    printf("%s\n", cmd);
    exit(1);
}

/*---------------------------------------------------------------*/
int GetFluidFlag(void) {
    printf("Select a Fluid Property Type\n");
    printf("Select\n\tIncompressible\t0\n\tCompressible\t1\n");
    printf("Fluid Type (0/1):");
    scanf("%d", &FluidFlag);

    return 0;
}

/*---------------------------------------------------------------*/
int GetSpecificDensityConst(void) {

    printf("Specific Density\n");
    printf("\tSpecific Density = ");
    scanf("%lf", &SpecificDensity_Const);

    return 0;
}

/*---------------------------------------------------------------*/
int GetSpecificHeat(void) {
    int flag, select;

    printf("Specific Heat\n");
    printf("\tConstant or Variable (0/1):");
    scanf("%d", &flag);

    switch (flag) {
        case 0:
            /* Constant */
            printf("Specific Heat = ");
            scanf("%lf", &SpecificHeat_Const);
            SpecificHeat_Flag = 0;
            break;
        case 1:
            /* Variable */
            SpecificHeat_Flag = 1;
            printf("Constant Specific Heat Reference Value:\n");
            printf("Specific Heat (free-stream) = ");
            scanf("%lf", &SpecificHeat_Const);

            printf("Variable Specific Heat:\n");
            GetSolutionList_2D();
            printf("Select the Variable = ");
            scanf("%d", &select);

            SpecificHeat_Var = Solution2D.Sols[select];
            if (SpecificHeat_Var == NULL)
                ERROR("GetSpecificHeat: Solution Memory Location Failed");

            break;
        default:
            printf("Error: Invalid Option\n");
            exit(1);
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetViscosity(void) {
    int flag, select;

    printf("Viscosity\n");
    printf("\tConstant or Variable (0/1):");
    scanf("%d", &flag);

    switch (flag) {
        case 0:
            /* Constant */
            Viscosity_Flag = 0;
            printf("Viscosity = ");
            scanf("%lf", &Viscosity_Const);
            break;
        case 1:
            /* Variable */
            Viscosity_Flag = 1;
            printf("Constant Viscosity Reference Value:\n");
            printf("Viscosity (free-stream) = ");
            scanf("%lf", &Viscosity_Const);

            printf("Variable Viscosity:\n");
            GetSolutionList_2D();
            printf("Select the Variable = ");
            scanf("%d", &select);

            Viscosity_Var = Solution2D.Sols[select];
            if (Viscosity_Var == NULL)
                ERROR("GetViscosity: Solution Memory Location Failed");

            break;
        default:
            printf("Error: Invalid Option\n");
            exit(1);
    }


    return 0;
}

/*---------------------------------------------------------------*/
int GetConductivity(void) {
    int flag, select;

    printf("Conductivity\n");
    printf("\tConstant or Variable (0/1):");
    scanf("%d", &flag);

    switch (flag) {
        case 0:
            /* Constant */
            Conductivity_Flag = 0;
            printf("Conductivity = ");
            scanf("%lf", &Conductivity_Const);
            break;
        case 1:
            /* Variable */
            Conductivity_Flag = 1;
            printf("Constant Conductivity Reference Value:\n");
            printf("Conductivity (free-stream) = ");
            scanf("%lf", &Conductivity_Const);

            printf("Variable Conductivity:\n");
            GetSolutionList_2D();
            printf("Select the Variable = ");
            scanf("%d", &select);

            Conductivity_Var = Solution2D.Sols[select];
            if (Conductivity_Var == NULL)
                ERROR("GetConductivity: Solution Memory Location Failed");

            break;
        default:
            printf("Error: Invalid Option\n");
            exit(1);
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetGasConstant(void) {
    int flag, select;

    printf("Gas Constant\n");
    printf("\tConstant or Variable (0/1):");
    scanf("%d", &flag);

    switch (flag) {
        case 0:
            /* Constant */
            GasConstant_Flag = 0;
            printf("Gas Constant = ");
            scanf("%lf", &GasConstant_Const);
            break;
        case 1:
            /* Variable */
            GasConstant_Flag = 1;
            printf("Constant Gas Constant Reference Value:\n");
            printf("Gas Constant (free-stream) = ");
            scanf("%lf", &GasConstant_Const);

            printf("Variable GasConstant:\n");
            GetSolutionList_2D();
            printf("Select the Variable = ");
            scanf("%d", &select);

            GasConstant_Var = Solution2D.Sols[select];
            if (GasConstant_Var == NULL)
                ERROR("GetGasConstant: Solution Memory Location Failed");

            break;
        default:
            printf("Error: Invalid Option\n");
            exit(1);
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetSpecificHeatRatio(void) {
    int flag, select;

    printf("Specific Heat Ratio\n");
    printf("\tConstant or Variable (0/1):");
    scanf("%d", &flag);

    switch (flag) {
        case 0:
            /* Constant */
            SpecificHeatRatio_Flag = 0;
            printf("Specific Heat Ratio = ");
            scanf("%lf", &SpecificHeatRatio_Const);
            break;
        case 1:
            /* Variable */
            SpecificHeatRatio_Flag = 1;
            printf("Constant Specific Heat Ratio Reference Value:\n");
            printf("Specific Heat Ratio (free-stream) = ");
            scanf("%lf", &SpecificHeatRatio_Const);

            printf("Variable Specific Heat Ratio:\n");
            GetSolutionList_2D();
            printf("Select the Variable = ");
            scanf("%d", &select);

            SpecificHeatRatio_Var = Solution2D.Sols[select];
            if (SpecificHeatRatio_Var == NULL)
                ERROR("GetSpecificHeatRatio: Solution Memory Location Failed");

            break;
        default:
            printf("Error: Invalid Option\n");
            exit(1);
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetIncompressibleFluidProperties(void) {
    GetSpecificDensityConst();
    GetSpecificHeat();
    GetViscosity();
    GetConductivity();

    return 0;
}

/*---------------------------------------------------------------*/
int GetCompressibleFluidProperties(void) {
    GetGasConstant();
    GetSpecificHeatRatio();
    GetViscosity();
    GetConductivity();

    return 0;
}

/*---------------------------------------------------------------*/
int GetFluidProperties(void) {
    GetFluidFlag();

    switch (FluidFlag) {
        case InCompressible:
            /* Incompressible */
            GetIncompressibleFluidProperties();
            break;
        case Compressible:
            /* Compressible */
            GetCompressibleFluidProperties();
            break;
        default:
            /* Error */
            printf("Invalid Fluid Type: Neither InCompressible nor Compressible\n");
            exit(1);
    }

    return 0;
}
