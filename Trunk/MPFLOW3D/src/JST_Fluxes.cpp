/*******************************************************************************
 * File:        JST_Fluxes.cpp
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
//! Create JST Scheme Data Structure
//------------------------------------------------------------------------------
void JST_Init(void) {
    
}

//------------------------------------------------------------------------------
//! Delete JST Scheme Data Structure
//------------------------------------------------------------------------------
void JST_Finalize(void) {
    
}

//------------------------------------------------------------------------------
//! Reset JST Scheme Data Structure
//------------------------------------------------------------------------------
void JST_Reset(void) {
    
}

//------------------------------------------------------------------------------
//! Compute JST Flux
//------------------------------------------------------------------------------
void Compute_Flux_JST(int node_L, int node_R, Vector3D areavec, double *Flux_JST_Conv, double *Flux_JST_Diss, int AddTime) {
    
}

//------------------------------------------------------------------------------
//! Computes the JST Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_JST(int AddTime) {
    
}

