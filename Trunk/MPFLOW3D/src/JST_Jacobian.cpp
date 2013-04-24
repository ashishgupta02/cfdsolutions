/*******************************************************************************
 * File:        JST_Jacobian.cpp
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
static int JSTJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create JST Jacobian Data Structure
//------------------------------------------------------------------------------
void JSTJac_Init(void) {
    
    // Check if JST Data Structure is required
    if (JSTJac_DB == 0) {
        JSTJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete JST Jacobian Data Structure
//------------------------------------------------------------------------------
void JSTJac_Finalize(void) {
    JSTJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset JST Jacobian Data Structure
//------------------------------------------------------------------------------
void JSTJac_Reset(void) {
    
    if (JSTJac_DB == 0)
        JSTJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_JST_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_JST_Approximate(int AddTime, int Iteration) {
    
}

