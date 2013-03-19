/*******************************************************************************
 * File:        AUSM_Fluxes.cpp
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
//! Create AUSM Scheme Data Structure
//------------------------------------------------------------------------------
void AUSM_Init(void) {
    
}

//------------------------------------------------------------------------------
//! Delete AUSM Scheme Data Structure
//------------------------------------------------------------------------------
void AUSM_Finalize(void) {
    
}

//------------------------------------------------------------------------------
//! Reset AUSM Scheme Data Structure
//------------------------------------------------------------------------------
void AUSM_Reset(void) {
    
}

//------------------------------------------------------------------------------
//! Compute AUSM Flux
//------------------------------------------------------------------------------
void Compute_AUSMFlux(int node_L, int node_R, Vector3D areavec, double *Flux_AUSM, int AddTime) {
    
}

//------------------------------------------------------------------------------
//! Computes the AUSM Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_AUSM(int AddTime) {
    
}
