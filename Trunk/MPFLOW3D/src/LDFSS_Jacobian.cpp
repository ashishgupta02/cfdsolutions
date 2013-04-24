/*******************************************************************************
 * File:        LDFSS_Jacobian.cpp
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
static int LDFSSJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create LDFSS Jacobian Data Structure
//------------------------------------------------------------------------------
void LDFSSJac_Init(void) {
    
    // Check if LDFSS Data Structure is required
    if (LDFSSJac_DB == 0) {
        LDFSSJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete LDFSS Jacobian Data Structure
//------------------------------------------------------------------------------
void LDFSSJac_Finalize(void) {
    LDFSSJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset LDFSS Jacobian Data Structure
//------------------------------------------------------------------------------
void LDFSSJac_Reset(void) {
    
    if (LDFSSJac_DB == 0)
        LDFSSJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_LDFSS_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_LDFSS_Approximate(int AddTime, int Iteration) {
    
}

