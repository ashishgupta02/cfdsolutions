/*******************************************************************************
 * File:        Residual.cpp
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
//! Computes the Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual(void) {
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Compute_Residual_Roe();
            break;
        case SOLVER_SCHEME_LMROE: // LMRoe
            LMRoeFix = 1;
            Compute_Residual_Roe();
            break;
        default:
            error("Compute_Residual: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
}


