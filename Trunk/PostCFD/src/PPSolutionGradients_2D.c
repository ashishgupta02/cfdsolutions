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
 * File		PPSolutionGradients_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "PostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Interpolation_2D.h"
#include "Gradient_2D.h"
#include "PPSolutionGradients_2D.h"

/* Gradient Calculation Option  */
static char *GradientOption[] = {
    "usage  : Gradient [options]",
    "options:",
    "	1	= Line Integral",
    "	2	= Least Square Method",
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
int LSMSolutionFieldGradient_2D(int FieldLocation, Data_2D_Un *OutX, Data_2D_Un *OutY) {

    switch (Solution2D.Location) {
        case Vertex:
            /* To Calculate Node Gradient */
            NodeGradientLSM_2D(Solution2D.Sols[FieldLocation], OutX, OutY);
            break;
        case CellCenter:
            /* To Calculate Cell Gradient */
            CellGradientLSM_2D(Solution2D.Sols[FieldLocation], OutX, OutY);
            break;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int LISolutionFieldGradient_2D(int FieldLocation, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    Data_2D_Un tmp1;
    Data_2D_Un tmp2;

    switch (Solution2D.Location) {
        case Vertex:
            /* To Calculate Node Gradient */
            /* Pre-Calculate Gradient at Cell-Center */
            tmp1.Size = NoCells2D;
            tmp1.Data = (double *) malloc(NoCells2D * sizeof (double));
            if (tmp1.Data == NULL)
                FATAL(NULL, "LISolutionFieldGradient_2D: Memory Failure : 1");

            tmp2.Size = NoCells2D;
            tmp2.Data = (double *) malloc(NoCells2D * sizeof (double));
            if (tmp2.Data == NULL)
                FATAL(NULL, "LISolutionFieldGradient_2D: Memory Failure : 2");

            /* Cell Gradient */
            CellGradientLI_2D(Solution2D.Sols[FieldLocation], &tmp1, &tmp2);

            /* Node Gradient */
            NodeGradientLI_2D(&tmp1, &tmp2, OutX, OutY);

            free(tmp1.Data);
            free(tmp2.Data);
            break;
        case CellCenter:
            /* To Calculate Cell Gradient */
            tmp1.Size = NoNodes2D;
            tmp1.Data = (double *) malloc(NoNodes2D * sizeof (double));
            if (tmp1.Data == NULL)
                FATAL(NULL, "LISolutionFieldGradient_2D: Memory Failure : 3");

            /* Interpolate Data */
            CellCenter2Node_2D(Solution2D.Sols[FieldLocation], &tmp1);

            /* Calculate Gradient at Cell-Center */
            CellGradientLI_2D(&tmp1, OutX, OutY);

            free(tmp1.Data);
            break;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CalculateSolutionFieldGradient_2D(int FieldLocation, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    int Option;

    usage(GradientOption);
    printf("Give Option: ");
    scanf("%d", &Option);

    switch (Option) {
        case 1:
            /* Line Integral */
            LISolutionFieldGradient_2D(FieldLocation, OutX, OutY);
            break;
        case 2:
            /* Least Square Method */
            LSMSolutionFieldGradient_2D(FieldLocation, OutX, OutY);
            break;
        default:
            /* Error */
            MSG("Invalid Option");
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CalculateTemperatureGradient_2D(int TemperatureLocation) {
    int L1, L2;
    int Option1, Option2;

    Option2 = 1;

    /* Check Data Presence (In Solution2D) */
    if (CheckSolutionFieldExistance_2D("TemperatureGradientX")
            || CheckSolutionFieldExistance_2D("TemperatureGradientY")) {
        MSG("TemperatureGradient exist");
        printf("Replace Data (0/1) : ");
        scanf("%d", &Option1);

        switch (Option1) {
            case 0:
                Option2 = 0;
                break;
            case 1:
                Option2 = 2;
                break;
            default:
                Option2 = 0;
                MSG("Invalid Option: Gradient Not Calculated");
        }
    }

    switch (Option2) {
        case 0:
            /* Escape Calculation */
            break;
        case 1:
            /* Calculate Gradient */
            /* Initialize Name */
            strcpy(TemperatureGradientX2D.Name, "TemperatureGradientX");
            strcpy(TemperatureGradientY2D.Name, "TemperatureGradientY");

            /* If Vertex memory allocation */
            if (Solution2D.Location == Vertex) {
                /* Size of array */
                TemperatureGradientX2D.Size = NoNodes2D;
                TemperatureGradientY2D.Size = NoNodes2D;

                /* Memory allocation for X component */
                TemperatureGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (TemperatureGradientX2D.Data == NULL) {
                    MSG("CalculateTemperatureGradient_2D: Memory Failed: 1");
                    break;
                }

                /* Memory allocation for Y component */
                TemperatureGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (TemperatureGradientY2D.Data == NULL) {
                    MSG("CalculateTemperatureGradient_2D: Memory Failed: 2");
                    break;
                }
            }

            /* If CellCenter memory allocation */
            if (Solution2D.Location == CellCenter) {
                /* Size of array */
                TemperatureGradientX2D.Size = NoCells2D;
                TemperatureGradientY2D.Size = NoCells2D;

                /* Memory allocation for X component */
                TemperatureGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (TemperatureGradientX2D.Data == NULL) {
                    MSG("CalculateTemperatureGradient_2D: Memory Failed: 3");
                    break;
                }

                /* Memory allocation for Y component */
                TemperatureGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (TemperatureGradientY2D.Data == NULL) {
                    MSG("CalculateTemperatureGradient_2D: Memory Failed: 4");
                    break;
                }
            }

            /* Calculate Gradient */
            CalculateSolutionFieldGradient_2D(TemperatureLocation, &TemperatureGradientX2D, &TemperatureGradientY2D);

            /* Update Solution2D Struct */
            UpdateSolution_2D(&TemperatureGradientX2D);
            UpdateSolution_2D(&TemperatureGradientY2D);
            break;
        case 2:
            /* Replace Gradient */
            L1 = GetSolutionFieldLocation_2D("TemperatureGradientX");
            L2 = GetSolutionFieldLocation_2D("TemperatureGradientY");

            /* If Both GradientX and GradientY are Present */
            if ((L1 != -1) && (L2 != -1)) {
                CalculateSolutionFieldGradient_2D(TemperatureLocation, Solution2D.Sols[L1], Solution2D.Sols[L2]);
                break;
            }

            /* If Only GradientX is present */
            if (L1 != -1) {
                /* Initialize Name */
                strcpy(TemperatureGradientY2D.Name, "TemperatureGradientY");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    TemperatureGradientY2D.Size = NoNodes2D;

                    /* Memory allocation for Y component */
                    TemperatureGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (TemperatureGradientY2D.Data == NULL) {
                        MSG("CalculateTemperatureGradient_2D: Memory Failed: 5");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    TemperatureGradientY2D.Size = NoCells2D;

                    /* Memory allocation for Y component */
                    TemperatureGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (TemperatureGradientY2D.Data == NULL) {
                        MSG("CalculateTemperatureGradient_2D: Memory Failed: 6");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(TemperatureLocation, Solution2D.Sols[L1], &TemperatureGradientY2D);
                UpdateSolution_2D(&TemperatureGradientY2D);
                break;
            }

            /* If Only GradientY is present */
            if (L2 != -1) {
                /* Initialize Name */
                strcpy(TemperatureGradientX2D.Name, "TemperatureGradientX");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    TemperatureGradientX2D.Size = NoNodes2D;

                    /* Memory allocation for X component */
                    TemperatureGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (TemperatureGradientX2D.Data == NULL) {
                        MSG("CalculateTemperatureGradient_2D: Memory Failed: 7");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    TemperatureGradientX2D.Size = NoCells2D;

                    /* Memory allocation for X component */
                    TemperatureGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (TemperatureGradientX2D.Data == NULL) {
                        MSG("CalculateTemperatureGradient_2D: Memory Failed: 8");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(TemperatureLocation, Solution2D.Sols[L2], &TemperatureGradientX2D);
                UpdateSolution_2D(&TemperatureGradientX2D);
                break;
            }

            break;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CalculatePressureGradient_2D(int PressureLocation) {
    int L1, L2;
    int Option1, Option2;

    Option2 = 1;

    /* Check Data Presence (In Solution2D) */
    if (CheckSolutionFieldExistance_2D("PressureGradientX")
            || CheckSolutionFieldExistance_2D("PressureGradientY")) {
        MSG("PressureGradient exist");
        printf("Replace Data (0/1) : ");
        scanf("%d", &Option1);

        switch (Option1) {
            case 0:
                Option2 = 0;
                break;
            case 1:
                Option2 = 2;
                break;
            default:
                Option2 = 0;
                MSG("Invalid Option: Gradient Not Calculated");
        }
    }

    switch (Option2) {
        case 0:
            /* Escape Calculation */
            break;
        case 1:
            /* Calculate Gradient */
            /* Initialize Name */
            strcpy(PressureGradientX2D.Name, "PressureGradientX");
            strcpy(PressureGradientY2D.Name, "PressureGradientY");

            /* If Vertex memory allocation */
            if (Solution2D.Location == Vertex) {
                /* Size of array */
                PressureGradientX2D.Size = NoNodes2D;
                PressureGradientY2D.Size = NoNodes2D;

                /* Memory allocation for X component */
                PressureGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (PressureGradientX2D.Data == NULL) {
                    MSG("CalculatePressureGradient_2D: Memory Failed: 1");
                    break;
                }

                /* Memory allocation for Y component */
                PressureGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (PressureGradientY2D.Data == NULL) {
                    MSG("CalculatePressureGradient_2D: Memory Failed: 2");
                    break;
                }
            }

            /* If CellCenter memory allocation */
            if (Solution2D.Location == CellCenter) {
                /* Size of array */
                PressureGradientX2D.Size = NoCells2D;
                PressureGradientY2D.Size = NoCells2D;

                /* Memory allocation for X component */
                PressureGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (PressureGradientX2D.Data == NULL) {
                    MSG("CalculatePressureGradient_2D: Memory Failed: 3");
                    break;
                }

                /* Memory allocation for Y component */
                PressureGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (PressureGradientY2D.Data == NULL) {
                    MSG("CalculatePressureGradient_2D: Memory Failed: 4");
                    break;
                }
            }

            /* Calculate Gradient */
            CalculateSolutionFieldGradient_2D(PressureLocation, &PressureGradientX2D, &PressureGradientY2D);

            /* Update Solution2D Struct */
            UpdateSolution_2D(&PressureGradientX2D);
            UpdateSolution_2D(&PressureGradientY2D);
            break;
        case 2:
            /* Replace Gradient */
            L1 = GetSolutionFieldLocation_2D("PressureGradientX");
            L2 = GetSolutionFieldLocation_2D("PressureGradientY");

            /* If Both GradientX and GradientY are Present */
            if ((L1 != -1) && (L2 != -1)) {
                CalculateSolutionFieldGradient_2D(PressureLocation, Solution2D.Sols[L1], Solution2D.Sols[L2]);
                break;
            }

            /* If Only GradientX is present */
            if (L1 != -1) {
                /* Initialize Name */
                strcpy(PressureGradientY2D.Name, "PressureGradientY");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    PressureGradientY2D.Size = NoNodes2D;

                    /* Memory allocation for Y component */
                    PressureGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (PressureGradientY2D.Data == NULL) {
                        MSG("CalculatePressureGradient_2D: Memory Failed: 5");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    PressureGradientY2D.Size = NoCells2D;

                    /* Memory allocation for Y component */
                    PressureGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (PressureGradientY2D.Data == NULL) {
                        MSG("CalculatePressureGradient_2D: Memory Failed: 6");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(PressureLocation, Solution2D.Sols[L1], &PressureGradientY2D);
                UpdateSolution_2D(&PressureGradientY2D);
                break;
            }

            /* If Only GradientY is present */
            if (L2 != -1) {
                /* Initialize Name */
                strcpy(PressureGradientX2D.Name, "PressureGradientX");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    PressureGradientX2D.Size = NoNodes2D;

                    /* Memory allocation for X component */
                    PressureGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (PressureGradientX2D.Data == NULL) {
                        MSG("CalculatePressureGradient_2D: Memory Failed: 7");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    PressureGradientX2D.Size = NoCells2D;

                    /* Memory allocation for X component */
                    PressureGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (PressureGradientX2D.Data == NULL) {
                        MSG("CalculatePressureGradient_2D: Memory Failed: 8");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(PressureLocation, Solution2D.Sols[L2], &PressureGradientX2D);
                UpdateSolution_2D(&PressureGradientX2D);
                break;
            }

            break;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CalculateDensityGradient_2D(int DensityLocation) {
    int L1, L2;
    int Option1, Option2;

    Option2 = 1;

    /* Check Data Presence (In Solution2D) */
    if (CheckSolutionFieldExistance_2D("DensityGradientX")
            || CheckSolutionFieldExistance_2D("DensityGradientY")) {
        MSG("DensityGradient exist");
        printf("Replace Data (0/1) : ");
        scanf("%d", &Option1);

        switch (Option1) {
            case 0:
                Option2 = 0;
                break;
            case 1:
                Option2 = 2;
                break;
            default:
                Option2 = 0;
                MSG("Invalid Option: Gradient Not Calculated");
        }
    }

    switch (Option2) {
        case 0:
            /* Escape Calculation */
            break;
        case 1:
            /* Calculate Gradient */
            /* Initialize Name */
            strcpy(DensityGradientX2D.Name, "DensityGradientX");
            strcpy(DensityGradientY2D.Name, "DensityGradientY");

            /* If Vertex memory allocation */
            if (Solution2D.Location == Vertex) {
                /* Size of array */
                DensityGradientX2D.Size = NoNodes2D;
                DensityGradientY2D.Size = NoNodes2D;

                /* Memory allocation for X component */
                DensityGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (DensityGradientX2D.Data == NULL) {
                    MSG("CalculateDensityGradient_2D: Memory Failed: 1");
                    break;
                }

                /* Memory allocation for Y component */
                DensityGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (DensityGradientY2D.Data == NULL) {
                    MSG("CalculateDensityGradient_2D: Memory Failed: 2");
                    break;
                }
            }

            /* If CellCenter memory allocation */
            if (Solution2D.Location == CellCenter) {
                /* Size of array */
                DensityGradientX2D.Size = NoCells2D;
                DensityGradientY2D.Size = NoCells2D;

                /* Memory allocation for X component */
                DensityGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (DensityGradientX2D.Data == NULL) {
                    MSG("CalculateDensityGradient_2D: Memory Failed: 3");
                    break;
                }

                /* Memory allocation for Y component */
                DensityGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                if (DensityGradientY2D.Data == NULL) {
                    MSG("CalculateDensityGradient_2D: Memory Failed: 4");
                    break;
                }
            }

            /* Calculate Gradient */
            CalculateSolutionFieldGradient_2D(DensityLocation, &DensityGradientX2D, &DensityGradientY2D);

            /* Update Solution2D Struct */
            UpdateSolution_2D(&DensityGradientX2D);
            UpdateSolution_2D(&DensityGradientY2D);
            break;
        case 2:
            /* Replace Gradient */
            L1 = GetSolutionFieldLocation_2D("DensityGradientX");
            L2 = GetSolutionFieldLocation_2D("DensityGradientY");

            /* If Both GradientX and GradientY are Present */
            if ((L1 != -1) && (L2 != -1)) {
                CalculateSolutionFieldGradient_2D(DensityLocation, Solution2D.Sols[L1], Solution2D.Sols[L2]);
                break;
            }

            /* If Only GradientX is present */
            if (L1 != -1) {
                /* Initialize Name */
                strcpy(DensityGradientY2D.Name, "DensityGradientY");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    DensityGradientY2D.Size = NoNodes2D;

                    /* Memory allocation for Y component */
                    DensityGradientY2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (DensityGradientY2D.Data == NULL) {
                        MSG("CalculateDensityGradient_2D: Memory Failed: 5");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    DensityGradientY2D.Size = NoCells2D;

                    /* Memory allocation for Y component */
                    DensityGradientY2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (DensityGradientY2D.Data == NULL) {
                        MSG("CalculateDensityGradient_2D: Memory Failed: 6");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(DensityLocation, Solution2D.Sols[L1], &DensityGradientY2D);
                UpdateSolution_2D(&DensityGradientY2D);
                break;
            }

            /* If Only GradientY is present */
            if (L2 != -1) {
                /* Initialize Name */
                strcpy(DensityGradientX2D.Name, "DensityGradientX");

                /* If Vertex memory allocation */
                if (Solution2D.Location == Vertex) {
                    /* Size of array */
                    DensityGradientX2D.Size = NoNodes2D;

                    /* Memory allocation for X component */
                    DensityGradientX2D.Data = (double *) malloc(NoNodes2D * sizeof (double));
                    if (DensityGradientX2D.Data == NULL) {
                        MSG("CalculateDensityGradient_2D: Memory Failed: 7");
                        break;
                    }
                }

                /* If CellCenter memory allocation */
                if (Solution2D.Location == CellCenter) {
                    /* Size of array */
                    DensityGradientX2D.Size = NoCells2D;

                    /* Memory allocation for X component */
                    DensityGradientX2D.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if (DensityGradientX2D.Data == NULL) {
                        MSG("CalculateDensityGradient_2D: Memory Failed: 8");
                        break;
                    }
                }

                CalculateSolutionFieldGradient_2D(DensityLocation, Solution2D.Sols[L2], &DensityGradientX2D);
                UpdateSolution_2D(&DensityGradientX2D);
                break;
            }

            break;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int SolutionGradient_2D(void) {
    int select, Option;

    printf("Solution Gradient Module:\n");
    printf("Calculate:\n\t1) Density Gradient\n\t2) Pressure Gradient\n\t3) Temperature Gradient\n");
    printf("Select Solution: ");
    scanf("%d", &Option);

    switch (Option) {
        case 1:
            /* Density Gradient */
            printf("Identify Density\n");
            GetSolutionList_2D();
            printf("Density (-1 if not available): ");
            scanf("%d", &select);
            if (select == -1)
                break;

            if ((select >= Solution2D.Size) || (select < -1)) {
                MSG("Invalid Solution Selection");
                break;
            }

            CalculateDensityGradient_2D(select);
            break;
        case 2:
            /* Pressure Gradient */
            printf("Identify Pressure\n");
            GetSolutionList_2D();
            printf("Pressure (-1 if not available): ");
            scanf("%d", &select);
            if (select == -1)
                break;

            if ((select >= Solution2D.Size) || (select < -1)) {
                MSG("Invalid Solution Selection");
                break;
            }

            CalculatePressureGradient_2D(select);
            break;
        case 3:
            /* Temperature Gradient */
            printf("Identify Temperature\n");
            GetSolutionList_2D();
            printf("Temperature (-1 if not available): ");
            scanf("%d", &select);
            if (select == -1)
                break;

            if ((select >= Solution2D.Size) || (select < -1)) {
                MSG("Invalid Solution Selection");
                break;
            }

            CalculateTemperatureGradient_2D(select);
            break;
        default:
            /* Waring Message */
            MSG("Invalid Option");
    }

    return 0;
}
