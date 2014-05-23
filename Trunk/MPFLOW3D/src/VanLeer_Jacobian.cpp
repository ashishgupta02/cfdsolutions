/*******************************************************************************
 * File:        VanLeer_Jacobian.cpp
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
static int VanLeerJac_DB   = 0;

//------------------------------------------------------------------------------
//! Create VanLeer Jacobian Data Structure
//------------------------------------------------------------------------------
void VanLeerJac_Init(void) {

    // Check if VanLeer Data Structure is required
    if (VanLeerJac_DB == 0) {
        VanLeerJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete VanLeer Jacobian Data Structure
//------------------------------------------------------------------------------
void VanLeerJac_Finalize(void) {
    VanLeerJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset VanLeer Jacobian Data Structure
//------------------------------------------------------------------------------
void VanLeerJac_Reset(void) {

    if (VanLeerJac_DB == 0)
        VanLeerJac_Init();

    // Initialization
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_VanLeer_Exact(int AddTime, int Iteration) {

}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_VanLeer_Approximate(int AddTime, int Iteration) {

}

