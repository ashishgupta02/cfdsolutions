/*******************************************************************************
 * File:        Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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
    int saveorder, savelmroefix;
    
    // Save the Order
    saveorder = Order;
    Order = 1;
    
    // Save the LMRoeFix
    savelmroefix = LMRoeFix;
    LMRoeFix = 0;
    
    // Check the frequency of Jacobian Update
    if (Iteration%JacobianUpdate != 0)
        return;
    
    // Update the Jacobian
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Compute_Jacobian_Roe(AddTime, Iteration);
            break;
        case SOLVER_SCHEME_LMROE: // LMRoe
            Compute_Jacobian_Roe(AddTime, Iteration);
            break;
        default:
            error("Compute_Jacobian: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
    
    // Switch Back to Original Order
    Order    = saveorder;
    LMRoeFix = savelmroefix;
}

