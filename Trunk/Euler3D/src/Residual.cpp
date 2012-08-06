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
void Compute_Residual(int AddTime) {
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Compute_Residual_Roe(AddTime);
            break;
        case SOLVER_SCHEME_HLLC: // HLLC
            Compute_Residual_HLLC(AddTime);
            break;
        case SOLVER_SCHEME_AUSM: // AUSM
            Compute_Residual_AUSM(AddTime);
            break;
        case SOLVER_SCHEME_VANLEER: // Van Leer
            Compute_Residual_VanLeer(AddTime);
            break;
        case SOLVER_SCHEME_LDFSS: // LDFSS
            Compute_Residual_LDFSS(AddTime);
            break;
        case SOLVER_SCHEME_OSHER: // Osher
            Compute_Residual_Osher(AddTime);
            break;
        case SOLVER_SCHEME_STEGERWARMING: // Steger Warming
            Compute_Residual_StegerWarming(AddTime);
            break;
        case SOLVER_SCHEME_JST: // JST
            Compute_Residual_JST(AddTime);
            break;
        default:
            error("Compute_Residual: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
}

