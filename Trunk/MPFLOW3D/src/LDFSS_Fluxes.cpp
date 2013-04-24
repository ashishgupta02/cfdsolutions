/*******************************************************************************
 * File:        LDFSS_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 * Reference:   
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"

// Static Variable for Speed Up (if any)


//------------------------------------------------------------------------------
//! Create LDFSS Scheme Data Structure
//------------------------------------------------------------------------------
void LDFSS_Init(void) {
    
}

//------------------------------------------------------------------------------
//! Delete LDFSS Scheme Data Structure
//------------------------------------------------------------------------------
void LDFSS_Finalize(void) {
    
}

//------------------------------------------------------------------------------
//! Reset LDFSS Scheme Data Structure
//------------------------------------------------------------------------------
void LDFSS_Reset(void) {
    
}

//------------------------------------------------------------------------------
//! Compute LDFSS Flux
//------------------------------------------------------------------------------
void Compute_Flux_LDFSS(int node_L, int node_R, Vector3D areavec, double *Flux_LDFSS_Conv, double *Flux_LDFSS_Diss, int AddTime) {
    
}

//------------------------------------------------------------------------------
//! Computes the LDFSS Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_LDFSS(int AddTime) {
    
}
