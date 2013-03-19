/*******************************************************************************
 * File:        Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian(int AddTime, int Iteration) {
    int saveorder;
    
    // Save the Solver Order
    saveorder   = SolverOrder;
    SolverOrder = SOLVER_ORDER_FIRST;
    
    // Check the frequency of Jacobian Update
    if (Iteration%JacobianUpdate != 0)
        return;
    
    // Update the Jacobian
    switch (FluxScheme) {
        case FLUX_SCHEME_ROE: // Roe
            Compute_Jacobian_Roe(AddTime, Iteration);
            break;
        default:
            error("Compute_Jacobian: Invalid Flux Scheme - %d", FluxScheme);
            break;
    }
    
    // For Unsteady Computations
    if (SolverMethod == SOLVER_METHOD_UNSTEADY) {
        
    }
    
    // Switch Back to Original Solver Order
    SolverOrder = saveorder;
}

