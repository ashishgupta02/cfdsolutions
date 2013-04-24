/*******************************************************************************
 * File:        HLLC_Jacobian.cpp
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
static int HLLCJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create HLLC Jacobian Data Structure
//------------------------------------------------------------------------------
void HLLCJac_Init(void) {
    
    // Check if HLLC Data Structure is required
    if (HLLCJac_DB == 0) {
        HLLCJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete HLLC Jacobian Data Structure
//------------------------------------------------------------------------------
void HLLCJac_Finalize(void) {
    HLLCJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset HLLC Jacobian Data Structure
//------------------------------------------------------------------------------
void HLLCJac_Reset(void) {
    
    if (HLLCJac_DB == 0)
        HLLCJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_HLLC_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_HLLC_Approximate(int AddTime, int Iteration) {
    
}

