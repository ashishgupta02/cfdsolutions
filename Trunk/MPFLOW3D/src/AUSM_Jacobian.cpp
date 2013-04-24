/*******************************************************************************
 * File:        AUSM_Jacobian.cpp
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
static int AUSMJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create AUSM Jacobian Data Structure
//------------------------------------------------------------------------------
void AUSMJac_Init(void) {
    
    // Check if AUSM Data Structure is required
    if (AUSMJac_DB == 0) {
        AUSMJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete AUSM Jacobian Data Structure
//------------------------------------------------------------------------------
void AUSMJac_Finalize(void) {
    AUSMJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset AUSM Jacobian Data Structure
//------------------------------------------------------------------------------
void AUSMJac_Reset(void) {
    
    if (AUSMJac_DB == 0)
        AUSMJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_AUSM_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_AUSM_Approximate(int AddTime, int Iteration) {
    
}

