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
 * File		PPScalarFlowVariables_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Gradient_2D.h"
#include "FluidProperties.h"
#include "ReferenceQuantities.h"
#include "PPScalarFlowVariables_2D.h"

/* Scalar Flow Variable Calculation Option  */
static char *ScalarOption[] = {
    "usage  : Scalar Flow Variable [options]",
    "options:",
    "	1	= Density",
    "	2	= Pressure",
    "	3	= Temperature",
    "	4	= Speed of Sound",
    "	5	= Mach Number",
    "	6	= Stagnation Density",
    "	7	= Stagnation Pressure",
    "	8	= Stagnation Temperature",
    "	9	= Pressure Coefficient",
    "	10	= Stagnation Pressure Coefficient",
    "	11	= Pitot Pressure",
    "	12	= Pitot Pressure Ratio",
    "	13	= Dynamic Pressure",
    "	14	= Enthalpy",
    "	15	= Stagnation Enthalpy",
    "	16	= Internal Energy",
    "	17	= Stagnation Energy",
    "	18	= Stagnation Energy Density",
    "	19	= Kinetic Energy",
    "	20	= VelocityX",
    "	21	= VelocityY",
    "	22	= Velocity Magnitude",
    "	23	= Equivalent Potential Velocity Ratio",
    "	24	= MomentumX",
    "	25	= MomentumY",
    "	26	= Momentum Magnitude",
    "	27	= Entropy",
    "	28	= Entropy Measure S1",
    "	29	= VorticityX",
    "	30	= VorticityY",
    "	31	= Vorticity Magnitude",
    "	32	= Density Gradient Magnitude",
    "	33	= Pressure Gradient Magnitude",
    "	34	= Temperature Gradient Magnitude",
    "	35	= Divergence of Velocity",
    "	36	= Isentropic Density Ratio",
    "	37	= Isentropic Pressure Ratio",
    "	38	= Isentropic Temperature Ratio",
    "	39	= Shock",
    "	40	= Filtered Shock",
    "	41	= Swirl",
    "	42	= Helicity",
    "	43	= Relative Helicity",
    "	44	= Filtered Relative Helicity",
    "	45	= Sutherlands Law",
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

/*--------------------------------------------------------------*/
static void MSG(char *msg) {
    char cmd[129];

    sprintf(cmd, "Warning: {%s}", msg);
    printf("%s\n", cmd);
}

/*---------------------------------------------------------------*/
int CalculateDensity_2D(void) {
    int i, flag1, flag2, flag3;

    switch (FluidFlag) {
        case -1:
            /* Set Fluid Properties */
            printf("Set Fluid Properties\n");
            return 1;
            break;
        case 0:
            /* Incompressible Fluid */
            printf("Incompressible Flow: Constant Density\n");
            return 1;
            break;
        case 1:
            /* Compressible Fluid */
            flag1 = CheckSolutionFieldExistance_2D("Density");
            /* Density Missing */
            if (flag1 == 0) {
                /* Check Pressure */
                flag2 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag2 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag2 = GetSolutionFieldLocation_2D("Pressure");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Name Identifier */
                strcpy(Density2D.Name, "Density");

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Specify Size */
                        Density2D.Size = NoNodes2D;
                        /* Allocating Memory */
                        Density2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                        if (Density2D.Data == NULL) {
                            MSG("Memory Failed: 1");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Density */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Density2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Density2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Density2D);
                        break;
                    case CellCenter:
                        /* Specify Size */
                        Density2D.Size = NoCells2D;
                        /* Allocating Memory */
                        Density2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                        if (Density2D.Data == NULL) {
                            MSG("Memory Failed: 2");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Density */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Density2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Density2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Density2D);
                        break;
                }
            }
            /* If Density Exist */
            if (flag1 == 1) {
                /* Get Density Location */
                flag1 = GetSolutionFieldLocation_2D("Density");
                if (flag1 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Get replace flag */
                MSG("Density Exists");
                printf("Replace (0/1):");
                scanf("%d", &flag2);

                if (flag2 == 0) {
                    printf("Calculation Skipped\n");
                    return 0;
                }

                flag2 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag2 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag2 = GetSolutionFieldLocation_2D("Pressure");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Calculate Density */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                    case CellCenter:
                        /* Calculate Density */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                }
            }
            break;
    }
    printf("Density Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePressure_2D(void) {
    int i, flag1, flag2, flag3;

    switch (FluidFlag) {
        case -1:
            /* Set Fluid Properties */
            printf("Set Fluid Properties\n");
            return 1;
            break;
        case 0:
            /* Incompressible Fluid */
            flag1 = CheckSolutionFieldExistance_2D("Pressure");
            /* Pressure Missing */
            if (flag1 == 0) {
                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Name Identifier */
                strcpy(Pressure2D.Name, "Pressure");

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Specify Size */
                        Pressure2D.Size = NoNodes2D;
                        /* Allocating Memory */
                        Pressure2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                        if (Pressure2D.Data == NULL) {
                            MSG("Memory Failed: 1");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Pressure2D.Data[i] = SpecificDensity_Const * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Pressure2D.Data[i] = SpecificDensity_Const * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Pressure2D);
                        break;
                    case CellCenter:
                        /* Specify Size */
                        Pressure2D.Size = NoCells2D;
                        /* Allocating Memory */
                        Pressure2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                        if (Pressure2D.Data == NULL) {
                            MSG("Memory Failed: 2");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Pressure2D.Data[i] = SpecificDensity_Const * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Pressure2D.Data[i] = SpecificDensity_Const * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Pressure2D);
                        break;
                }
            }
            /* If Pressure Exist */
            if (flag1 == 1) {
                /* Get Pressure Location */
                flag1 = GetSolutionFieldLocation_2D("Pressure");
                if (flag1 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Get replace flag */
                MSG("Pressure Exists");
                printf("Replace (0/1):");
                scanf("%d", &flag2);

                if (flag2 == 0) {
                    printf("Calculation Skipped\n");
                    return 0;
                }

                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = SpecificDensity_Const * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = SpecificDensity_Const * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                    case CellCenter:
                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = SpecificDensity_Const * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = SpecificDensity_Const * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                }
            }
            break;
        case 1:
            /* Compressible Fluid */
            flag1 = CheckSolutionFieldExistance_2D("Pressure");
            /* Pressure Missing */
            if (flag1 == 0) {
                /* Check Density */
                flag2 = CheckSolutionFieldExistance_2D("Density");
                /* Calculate Density */
                if (flag2 == 0)
                    CalculateDensity_2D();

                /* Density Field Location */
                flag2 = GetSolutionFieldLocation_2D("Density");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Name Identifier */
                strcpy(Pressure2D.Name, "Pressure");

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Specify Size */
                        Pressure2D.Size = NoNodes2D;
                        /* Allocating Memory */
                        Pressure2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                        if (Pressure2D.Data == NULL) {
                            MSG("Memory Failed: 1");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Pressure2D.Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Pressure2D.Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Pressure2D);
                        break;
                    case CellCenter:
                        /* Specify Size */
                        Pressure2D.Size = NoCells2D;
                        /* Allocating Memory */
                        Pressure2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                        if (Pressure2D.Data == NULL) {
                            MSG("Memory Failed: 2");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Pressure2D.Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Pressure2D.Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Pressure2D);
                        break;
                }
            }
            /* If Pressure Exist */
            if (flag1 == 1) {
                /* Get Pressure Location */
                flag1 = GetSolutionFieldLocation_2D("Pressure");
                if (flag1 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Get replace flag */
                MSG("Pressure Exists");
                printf("Replace (0/1):");
                scanf("%d", &flag2);

                if (flag2 == 0) {
                    printf("Calculation Skipped\n");
                    return 0;
                }

                flag2 = CheckSolutionFieldExistance_2D("Density");
                /* Calculate Density */
                if (flag2 == 0)
                    CalculateDensity_2D();

                /* Density Field Location */
                flag2 = GetSolutionFieldLocation_2D("Density");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Temperature */
                flag3 = CheckSolutionFieldExistance_2D("Temperature");
                /* Calculate Temperature */
                if (flag3 == 0)
                    CalculateTemperature_2D();

                /* Temperature Field Location */
                flag3 = GetSolutionFieldLocation_2D("Temperature");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                    case CellCenter:
                        /* Calculate Pressure */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] * (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                }
            }
            break;
    }
    printf("Pressure Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateTemperature_2D(void) {
    int i, flag1, flag2, flag3;

    switch (FluidFlag) {
        case -1:
            /* Set Fluid Properties */
            printf("Set Fluid Properties\n");
            return 1;
            break;
        case 0:
            /* Incompressible Fluid */
            flag1 = CheckSolutionFieldExistance_2D("Temperature");
            /* Temperature Missing */
            if (flag1 == 0) {
                /* Check Pressure */
                flag3 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag3 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag3 = GetSolutionFieldLocation_2D("Pressure");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Name Identifier */
                strcpy(Temperature2D.Name, "Temperature");

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Specify Size */
                        Temperature2D.Size = NoNodes2D;
                        /* Allocating Memory */
                        Temperature2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                        if (Temperature2D.Data == NULL) {
                            MSG("Memory Failed: 1");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Const);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Var->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Temperature2D);
                        break;
                    case CellCenter:
                        /* Specify Size */
                        Temperature2D.Size = NoCells2D;
                        /* Allocating Memory */
                        Temperature2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                        if (Temperature2D.Data == NULL) {
                            MSG("Memory Failed: 2");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Const);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Var->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Temperature2D);
                        break;
                }
            }
            /* If Temperature Exist */
            if (flag1 == 1) {
                /* Get Temperature Location */
                flag1 = GetSolutionFieldLocation_2D("Temperature");
                if (flag1 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Get replace flag */
                MSG("Temperature Exists");
                printf("Replace (0/1):");
                scanf("%d", &flag2);

                if (flag2 == 0) {
                    printf("Calculation Skipped\n");
                    return 0;
                }

                /* Check Pressure */
                flag3 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag3 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag3 = GetSolutionFieldLocation_2D("Pressure");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Const);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Var->Data[i]);
                                break;
                        }
                        break;
                    case CellCenter:
                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Const);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag3]->Data[i] / (SpecificDensity_Const * GasConstant_Var->Data[i]);
                                break;
                        }
                        break;
                }
            }
            break;
        case 1:
            /* Compressible Fluid */
            flag1 = CheckSolutionFieldExistance_2D("Temperature");
            /* Temperature Missing */
            if (flag1 == 0) {
                /* Check Pressure */
                flag2 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag2 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag2 = GetSolutionFieldLocation_2D("Pressure");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Density */
                flag3 = CheckSolutionFieldExistance_2D("Density");
                /* Calculate Density */
                if (flag3 == 0)
                    CalculateDensity_2D();

                /* Density Field Location */
                flag3 = GetSolutionFieldLocation_2D("Density");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Name Identifier */
                strcpy(Pressure2D.Name, "Temperature");

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Specify Size */
                        Temperature2D.Size = NoNodes2D;
                        /* Allocating Memory */
                        Temperature2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                        if (Temperature2D.Data == NULL) {
                            MSG("Memory Failed: 1");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Temperature2D);
                        break;
                    case CellCenter:
                        /* Specify Size */
                        Temperature2D.Size = NoCells2D;
                        /* Allocating Memory */
                        Temperature2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                        if (Temperature2D.Data == NULL) {
                            MSG("Memory Failed: 2");
                            MSG("Calculation Unsuccessful");
                            return (-1);
                        }

                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Temperature2D.Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }

                        UpdateSolution_2D(&Temperature2D);
                        break;
                }
            }
            /* If Temperature Exist */
            if (flag1 == 1) {
                /* Get Temperature Location */
                flag1 = GetSolutionFieldLocation_2D("Temperature");
                if (flag1 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Get replace flag */
                MSG("Temperature Exists");
                printf("Replace (0/1):");
                scanf("%d", &flag2);

                if (flag2 == 0) {
                    printf("Calculation Skipped\n");
                    return 0;
                }

                flag2 = CheckSolutionFieldExistance_2D("Pressure");
                /* Calculate Pressure */
                if (flag2 == 0)
                    CalculatePressure_2D();

                /* Pressure Field Location */
                flag2 = GetSolutionFieldLocation_2D("Pressure");
                if (flag2 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Check Density */
                flag3 = CheckSolutionFieldExistance_2D("Density");
                /* Calculate Density */
                if (flag3 == 0)
                    CalculateDensity_2D();

                /* Density Field Location */
                flag3 = GetSolutionFieldLocation_2D("Density");
                if (flag3 == -1) {
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                switch (Solution2D.Location) {
                    case Vertex:
                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoNodes2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                    case CellCenter:
                        /* Calculate Temperature */
                        switch (GasConstant_Flag) {
                            case 0:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Const * Solution2D.Sols[flag3]->Data[i]);
                                break;
                            case 1:
                                for (i = 0; i < NoCells2D; i++)
                                    Solution2D.Sols[flag1]->Data[i] = Solution2D.Sols[flag2]->Data[i] / (GasConstant_Var->Data[i] * Solution2D.Sols[flag3]->Data[i]);
                                break;
                        }
                        break;
                }
            }
            break;
    }
    printf("Temperature Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateSpeedSound_2D(void) {
    printf("Speed of Sound Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateMachNumber_2D(void) {
    printf("Mach Number Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationDensity_2D(void) {
    printf("Stagnation Density Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationPressure_2D(void) {
    printf("Stagnation Pressure Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationTemperature_2D(void) {
    printf("Stagnation Temperature Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePressureCoefficient_2D(void) {
    int i, flag1, flag2;

    flag1 = CheckSolutionFieldExistance_2D("CoefPressure");
    /* If CoefPressure Not Exist */
    if (flag1 == 0) {
        flag2 = CheckSolutionFieldExistance_2D("Pressure");
        /* Calculate Pressure */
        if (flag2 == 0)
            CalculatePressure_2D();

        /* Pressure Field Location */
        flag2 = GetSolutionFieldLocation_2D("Pressure");
        if (flag2 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Name Identifier */
        strcpy(CoefPressure2D.Name, "CoefPressure");

        switch (Solution2D.Location) {
            case Vertex:
                /* Specify Size */
                CoefPressure2D.Size = NoNodes2D;
                /* Allocating Memory */
                CoefPressure2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (CoefPressure2D.Data == NULL) {
                    MSG("Memory Failed: 1");
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Calculate Coefficient of Pressure */
                for (i = 0; i < NoNodes2D; i++)
                    CoefPressure2D.Data[i] = (Solution2D.Sols[flag2]->Data[i] - Pressure_Ref) / ((0.5 * Density_Ref) * pow(VelocityMagnitude_Ref, 2));

                UpdateSolution_2D(&CoefPressure2D);
                break;
            case CellCenter:
                /* Specify Size */
                CoefPressure2D.Size = NoCells2D;
                /* Allocating Memory */
                CoefPressure2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (CoefPressure2D.Data == NULL) {
                    MSG("Memory Failed: 2");
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Calculate Coefficient of Pressure */
                for (i = 0; i < NoCells2D; i++)
                    CoefPressure2D.Data[i] = (Solution2D.Sols[flag2]->Data[i] - Pressure_Ref) / ((0.5 * Density_Ref) * pow(VelocityMagnitude_Ref, 2));

                UpdateSolution_2D(&CoefPressure2D);
                break;
        }
    }
    /* If CoefPressure Exist */
    if (flag1 == 1) {
        /* Get CoefPressure Location */
        flag1 = GetSolutionFieldLocation_2D("CoefPressure");
        if (flag1 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Get replace flag */
        MSG("Coefficient of Pressure Exists");
        printf("Replace (0/1):");
        scanf("%d", &flag2);

        if (flag2 == 0) {
            printf("Calculation Skipped\n");
            return 0;
        }

        flag2 = CheckSolutionFieldExistance_2D("Pressure");
        /* Calculate Pressure */
        if (flag2 == 0)
            CalculatePressure_2D();

        /* Pressure Field Location */
        flag2 = GetSolutionFieldLocation_2D("Pressure");
        if (flag2 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        switch (Solution2D.Location) {
            case Vertex:
                /* Calculate Coefficient of Pressure */
                for (i = 0; i < NoNodes2D; i++)
                    Solution2D.Sols[flag1]->Data[i] = (Solution2D.Sols[flag2]->Data[i] - Pressure_Ref) / ((0.5 * Density_Ref) * pow(VelocityMagnitude_Ref, 2));

                break;
            case CellCenter:
                /* Calculate Coefficient of Pressure */
                for (i = 0; i < NoCells2D; i++)
                    Solution2D.Sols[flag1]->Data[i] = (Solution2D.Sols[flag2]->Data[i] - Pressure_Ref) / ((0.5 * Density_Ref) * pow(VelocityMagnitude_Ref, 2));

                break;
        }
    }
    printf("Coefficient of Pressure Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationPressureCoef_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePitotPressure_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePitotPressureRatio_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateDynamicPressure_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateEnthalpy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationEnthalpy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateInternalEnergy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationEnergy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateStagnationEnergyDensity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateKineticEnergy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVelocityX_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVelocityY_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVelocityMagnitude_2D(void) {
    int i, flag1, flag2, flag3;

    flag1 = CheckSolutionFieldExistance_2D("VelocityMagnitude");
    /* VelocityMagnitude Missing */
    if (flag1 == 0) {
        /* Check VelocityX */
        flag2 = CheckSolutionFieldExistance_2D("VelocityX");
        /* Calculate VelocityX */
        if (flag2 == 0)
            CalculateVelocityX_2D();

        /* VelocityX Field Location */
        flag2 = GetSolutionFieldLocation_2D("VelocityX");
        if (flag2 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Check VelocityY */
        flag3 = CheckSolutionFieldExistance_2D("VelocityY");
        /* Calculate VelocityY */
        if (flag3 == 0)
            CalculateVelocityY_2D();

        /* VelocityY Field Location */
        flag3 = GetSolutionFieldLocation_2D("VelocityY");
        if (flag3 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Name Identifier */
        strcpy(VelocityMagnitude2D.Name, "VelocityMagnitude");

        switch (Solution2D.Location) {
            case Vertex:
                /* Specify Size */
                VelocityMagnitude2D.Size = NoNodes2D;
                /* Allocating Memory */
                VelocityMagnitude2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (VelocityMagnitude2D.Data == NULL) {
                    MSG("Memory Failed: 1");
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Calculate VelocityMagnitude */
                for (i = 0; i < NoNodes2D; i++)
                    VelocityMagnitude2D.Data[i] = sqrt(pow(Solution2D.Sols[flag2]->Data[i], 2) + pow(Solution2D.Sols[flag3]->Data[i], 2));

                UpdateSolution_2D(&VelocityMagnitude2D);
                break;
            case CellCenter:
                /* Specify Size */
                VelocityMagnitude2D.Size = NoCells2D;
                /* Allocating Memory */
                VelocityMagnitude2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (VelocityMagnitude2D.Data == NULL) {
                    MSG("Memory Failed: 2");
                    MSG("Calculation Unsuccessful");
                    return (-1);
                }

                /* Calculate VelocityMagnitude */
                for (i = 0; i < NoCells2D; i++)
                    VelocityMagnitude2D.Data[i] = sqrt(pow(Solution2D.Sols[flag2]->Data[i], 2) + pow(Solution2D.Sols[flag3]->Data[i], 2));

                UpdateSolution_2D(&VelocityMagnitude2D);
                break;
        }
    }
    /* If VelocityMagnitude Exist */
    if (flag1 == 1) {
        /* Get VelocityMagnitude Location */
        flag1 = GetSolutionFieldLocation_2D("VelocityMagnitude");
        if (flag1 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Get replace flag */
        MSG("VelocityMagnitude Exists");
        printf("Replace (0/1):");
        scanf("%d", &flag2);

        if (flag2 == 0) {
            printf("Calculation Skipped\n");
            return 0;
        }

        /* Check VelocityX */
        flag2 = CheckSolutionFieldExistance_2D("VelocityX");
        /* Calculate VelocityX */
        if (flag2 == 0)
            CalculateVelocityX_2D();

        /* VelocityX Field Location */
        flag2 = GetSolutionFieldLocation_2D("VelocityX");
        if (flag2 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        /* Check VelocityY */
        flag3 = CheckSolutionFieldExistance_2D("VelocityY");
        /* Calculate VelocityY */
        if (flag3 == 0)
            CalculateVelocityY_2D();

        /* VelocityY Field Location */
        flag3 = GetSolutionFieldLocation_2D("VelocityY");
        if (flag3 == -1) {
            MSG("Calculation Unsuccessful");
            return (-1);
        }

        switch (Solution2D.Location) {
            case Vertex:
                /* Calculate VelocityMagnitude */
                for (i = 0; i < NoNodes2D; i++)
                    Solution2D.Sols[flag1]->Data[i] = sqrt(pow(Solution2D.Sols[flag2]->Data[i], 2) + pow(Solution2D.Sols[flag3]->Data[i], 2));
                break;
            case CellCenter:
                /* Calculate VelocityMagnitude */
                for (i = 0; i < NoCells2D; i++)
                    Solution2D.Sols[flag1]->Data[i] = sqrt(pow(Solution2D.Sols[flag2]->Data[i], 2) + pow(Solution2D.Sols[flag3]->Data[i], 2));
                break;
        }
    }
    printf("Velocity Magnitude Calculation Successful\n");
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateEquivalentPotentialVelocityRatio_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateMomentumX_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateMomentumY_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateMomentumMagnitude_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateEntropy_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateEntropyMeasureS1_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVorticityX_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVorticityY_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateVorticityMagnitude_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateDensityGradientMagnitude_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePressureGradientMagnitude_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateTemperatureGradientMagnitude_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateDivergenceVelocity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalulateIsentropicDensityRatio_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateIsentropicPressureRatio_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateIsentropicTemperatureRatio_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateShock_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateFilteredShock_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateSwirl_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateHelicity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateRelativeHelicity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateFilteredRelativeHelicity_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int CalculateSutherlandsLaw_2D(void) {
    return 0;
}

/*---------------------------------------------------------------*/
int ScalarFlowVariableOptions_2D(void) {
    int Option;

    printf("Scalar Flow Variable Module:\n");
    usage(ScalarOption);
    printf("Select Option: ");
    scanf("%d", &Option);

    switch (Option) {
        case 1:
            /* Density */
            CalculateDensity_2D();
            break;
        case 2:
            /* Pressure */
            CalculatePressure_2D();
            break;
        case 3:
            /* Temperature */
            CalculateTemperature_2D();
            break;
        case 4:
            /* Speed of Sound */
            CalculateSpeedSound_2D();
            break;
        case 5:
            /* Mach Number */
            CalculateMachNumber_2D();
            break;
        case 6:
            /* Stagnation Density */
            CalculateStagnationDensity_2D();
            break;
        case 7:
            /* Stagnation Pressure */
            CalculateStagnationPressure_2D();
            break;
        case 8:
            /* Stagnation Temperature */
            CalculateStagnationTemperature_2D();
            break;
        case 9:
            /* Pressure Coefficient */
            CalculatePressureCoefficient_2D();
            break;
        case 10:
            /* Stagnation Pressure Coefficient */
            CalculateStagnationPressureCoef_2D();
            break;
        case 11:
            /* Pitot Pressure */
            CalculatePitotPressure_2D();
            break;
        case 12:
            /* Pitot Pressure Ratio */
            CalculatePitotPressureRatio_2D();
            break;
        case 13:
            /* Dynamic Pressure */
            CalculateDynamicPressure_2D();
            break;
        case 14:
            /* Enthalpy */
            CalculateEnthalpy_2D();
            break;
        case 15:
            /* Stagnation Enthalpy */
            CalculateStagnationEnthalpy_2D();
            break;
        case 16:
            /* Internal Energy */
            CalculateInternalEnergy_2D();
            break;
        case 17:
            /* Stagnation Energy */
            CalculateStagnationEnergy_2D();
            break;
        case 18:
            /* Stagnation Energy Density */
            CalculateStagnationEnergyDensity_2D();
            break;
        case 19:
            /* Kinetic Energy */
            CalculateKineticEnergy_2D();
            break;
        case 20:
            /* VelocityX */
            CalculateVelocityX_2D();
            break;
        case 21:
            /* VelocityY */
            CalculateVelocityY_2D();
            break;
        case 22:
            /* Velocity Magnitude */
            CalculateVelocityMagnitude_2D();
            break;
        case 23:
            /* Equivalent Potential Velocity Ratio */
            CalculateEquivalentPotentialVelocityRatio_2D();
            break;
        case 24:
            /* MomentumX */
            CalculateMomentumX_2D();
            break;
        case 25:
            /* MomentumY */
            CalculateMomentumY_2D();
            break;
        case 26:
            /* Momentum Magnitude */
            CalculateMomentumMagnitude_2D();
            break;
        case 27:
            /* Entropy */
            CalculateEntropy_2D();
            break;
        case 28:
            /* Entropy Measure S1 */
            CalculateEntropyMeasureS1_2D();
            break;
        case 29:
            /* VorticityX */
            CalculateVorticityX_2D();
            break;
        case 30:
            /* VorticityY */
            CalculateVorticityY_2D();
            break;
        case 31:
            /* Vorticity Magnitude */
            CalculateVorticityMagnitude_2D();
            break;
        case 32:
            /* Density Gradient Magnitude */
            CalculateDensityGradientMagnitude_2D();
            break;
        case 33:
            /* Pressure Gradient Magnitude */
            CalculatePressureGradientMagnitude_2D();
            break;
        case 34:
            /* Temperature Gradient Magnitude */
            CalculateTemperatureGradientMagnitude_2D();
            break;
        case 35:
            /* Divergence of Velocity */
            CalculateDivergenceVelocity_2D();
            break;
        case 36:
            /* Isentropic Density Ratio */
            CalulateIsentropicDensityRatio_2D();
            break;
        case 37:
            /* Isentropic Pressure Ratio */
            CalculateIsentropicPressureRatio_2D();
            break;
        case 38:
            /* Isentropic Temperature Ratio */
            CalculateIsentropicTemperatureRatio_2D();
            break;
        case 39:
            /* Shock */
            CalculateShock_2D();
            break;
        case 40:
            /* Filtered Shock */
            CalculateFilteredShock_2D();
            break;
        case 41:
            /* Swirl */
            CalculateSwirl_2D();
            break;
        case 42:
            /* Helicity */
            CalculateHelicity_2D();
            break;
        case 43:
            /* Relative Helicity */
            CalculateRelativeHelicity_2D();
            break;
        case 44:
            /* Filtered Relative Helicity */
            CalculateFilteredRelativeHelicity_2D();
            break;
        case 45:
            /* Sutherlands Law */
            CalculateSutherlandsLaw_2D();
            break;
        default:
            /* Error Message */
            MSG("Invalid Option");
    }
    return 0;
}
