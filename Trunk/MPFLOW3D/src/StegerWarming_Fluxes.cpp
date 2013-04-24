/*******************************************************************************
 * File:        StegerWarming_Fluxes.cpp
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
//! Create StegerWarming Scheme Data Structure
//------------------------------------------------------------------------------
void StegerWarming_Init(void) {
    
}

//------------------------------------------------------------------------------
//! Delete StegerWarming Scheme Data Structure
//------------------------------------------------------------------------------
void StegerWarming_Finalize(void) {
    
}

//------------------------------------------------------------------------------
//! Reset StegerWarming Scheme Data Structure
//------------------------------------------------------------------------------
void StegerWarming_Reset(void) {
    
}

//------------------------------------------------------------------------------
//! Compute StegerWarming Flux
//------------------------------------------------------------------------------
void Compute_Flux_StegerWarming(int node_L, int node_R, Vector3D areavec, double *Flux_StegerWarming_Conv, double *Flux_StegerWarming_Diss, int AddTime){
    
}

//------------------------------------------------------------------------------
//! Computes the StegerWarming Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_StegerWarming(int AddTime) {
    
}
