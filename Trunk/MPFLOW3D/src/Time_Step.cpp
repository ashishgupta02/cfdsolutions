/*******************************************************************************
 * File:        Time_Step.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Compute_DeltaT(int Iteration) {
    int i;
    double denom = 0.0;
    
    // Finally Compute the Local Time Stepping
    if ((CFL_Ramp > 1) && (CFL_MAX > CFL_MIN)) {
        if (Iteration < CFL_Ramp)
            CFL = CFL_MIN + (CFL_MAX - CFL_MIN)*(((double)Iteration)/((double)(CFL_Ramp-1)));
        else
            CFL = CFL_MAX;
    } else
        CFL = CFL_MAX;

    for (i = 0; i < nNode; i++) {
        DeltaT[i] = cVolume[i] * CFL / DeltaT[i];
        MinDeltaT = MIN(MinDeltaT, DeltaT[i]);
        MaxDeltaT = MAX(MaxDeltaT, DeltaT[i]);
    }

    // Check if Global Time Stepping is Required
    if (TimeIntegrationType == TIME_INTEGRATION_TYPE_GLOBAL) {
        denom = DeltaT[0];
        for (i = 0; i < nNode; i++)
            denom = MIN(denom, DeltaT[i]);
        for (i = 0; i < nNode; i++)
            DeltaT[i] = denom;
    }
}

