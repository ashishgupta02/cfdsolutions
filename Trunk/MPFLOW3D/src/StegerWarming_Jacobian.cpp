/*******************************************************************************
 * File:        StegerWarming_Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "CompressibleUtils.h"
#include "Material.h"
#include "Solver.h"

// Static Variable for Speed Up
static int StegerWarmingJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create StegerWarming Jacobian Data Structure
//------------------------------------------------------------------------------
void StegerWarmingJac_Init(void) {
    
    // Check if StegerWarming Data Structure is required
    if (StegerWarmingJac_DB == 0) {
        StegerWarmingJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete StegerWarming Jacobian Data Structure
//------------------------------------------------------------------------------
void StegerWarmingJac_Finalize(void) {
    StegerWarmingJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset StegerWarming Jacobian Data Structure
//------------------------------------------------------------------------------
void StegerWarmingJac_Reset(void) {
    
    if (StegerWarmingJac_DB == 0)
        StegerWarmingJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_StegerWarming_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_StegerWarming_Approximate(int AddTime, int Iteration) {
    
}

