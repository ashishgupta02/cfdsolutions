/*******************************************************************************
 * File:        Osher_Jacobian.cpp
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
static int OsherJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create Osher Jacobian Data Structure
//------------------------------------------------------------------------------
void OsherJac_Init(void) {
    
    // Check if Osher Data Structure is required
    if (OsherJac_DB == 0) {
        OsherJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Osher Jacobian Data Structure
//------------------------------------------------------------------------------
void OsherJac_Finalize(void) {
    OsherJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Osher Jacobian Data Structure
//------------------------------------------------------------------------------
void OsherJac_Reset(void) {
    
    if (OsherJac_DB == 0)
        OsherJac_Init();
    
    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_Osher_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_Osher_Approximate(int AddTime, int Iteration) {
    
}

