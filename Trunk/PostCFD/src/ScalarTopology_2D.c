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
 * File		ScalarTopology_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "Vector.h"
#include "PostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "CellPointOperation_2D.h"
#include "Gradient_2D.h"
#include "Interpolation_2D.h"
#include "Stream_2D.h"
#include "Streamline_2D.h"
#include "VectorTopology_2D.h"
#include "ScalarTopology_2D.h"

/* ScalarTopology Argument */
int ScalarTopologyArg = -1;
static int CheckScaTopArg = -1;

static int firstSca = 0;

/* Data Structure to Store Gradients */
Data_2D_Un GraVecX;
Data_2D_Un GraVecY;

/*--------------------------------------------------------------*/
void ScalarTopology_2D(void) {
    int flag = 0;
    Data_2D_Un tmp;

    if (ScalarTopologyArg == -1)
        return;

    /* Should not go for first time */
    if (firstSca) {
        if (ScalarTopologyArg != CheckScaTopArg)
            flag = 1;
    }

    /* If First time */
    if ((!firstSca) || (flag == 1)) {
        if (flag == 1) {
            if (GraVecX.Data != NULL) free(GraVecX.Data);
            if (GraVecY.Data != NULL) free(GraVecY.Data);
        }

        CheckScaTopArg = ScalarTopologyArg;

        switch (Solution2D.Location) {
            case Vertex:
                /* Allocate Gradient Vector Memory */
                GraVecX.Size = NoNodes2D;
                GraVecY.Size = NoNodes2D;

                GraVecX.Data = (double *) malloc(NoNodes2D * sizeof (double));
                GraVecY.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if ((GraVecX.Data == NULL) || (GraVecY.Data == NULL)) {
                    Warn("ScalarTopology_2D: Memory Allocation Failed 1");
                    return;
                }
                /* Get Scalar Gradient */
                NodeGradientLSM_2D(Solution2D.Sols[ScalarTopologyArg], &GraVecX, &GraVecY);

                break;
            case CellCenter:
                /* Initialize tmp variable */
                tmp.Size = NoNodes2D;
                tmp.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (tmp.Data == NULL) {
                    Warn("ScalarTopology_2D: Memory Allocation Failed 2");
                    return;
                }

                /* Interpolate from CellCenter to Vertex */
                ArithmeticCC2Node_2D(Solution2D.Sols[ScalarTopologyArg], &tmp);

                GraVecX.Size = NoNodes2D;
                GraVecY.Size = NoNodes2D;

                GraVecX.Data = (double *) malloc(NoNodes2D * sizeof (double));
                GraVecY.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if ((GraVecX.Data == NULL) || (GraVecY.Data == NULL)) {
                    Warn("ScalarTopology_2D: Memory Allocation Failed 3");
                    return;
                }

                /* Get Scalar Gradient */
                NodeGradientLSM_2D(&tmp, &GraVecX, &GraVecY);

                /* Release tmp memory */
                if (tmp.Data != NULL) free(tmp.Data);
                break;
        }

        firstSca = 1;
        /* Get Critical Point Using Vector Topology */
        VectorTopology_2D();
    }
}
