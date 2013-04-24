/*******************************************************************************
 * File:        Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "SolverParameters.h"
#include "Roe_Fluxes.h"
#include "HLLC_Fluxes.h"
#include "VanLeer_Fluxes.h"
#include "AUSM_Fluxes.h"
#include "JST_Fluxes.h"
#include "LDFSS_Fluxes.h"
#include "Osher_Fluxes.h"
#include "StegerWarming_Fluxes.h"

//------------------------------------------------------------------------------
//! Computes the Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Flux(int node_L, int node_R, Vector3D areavec, double *Flux_Conv, double *Flux_Diss, int AddTime) {
    switch (FluxScheme) {
        case FLUX_SCHEME_ROE: // Roe
            Compute_Flux_Roe(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        case FLUX_SCHEME_HLLC: // HLLC
            Compute_Flux_HLLC(node_L, node_R, areavec, Flux_Conv, AddTime);
            for (int i = 0; i < NEQUATIONS; i++)
                Flux_Diss[i] = 0.0;
            break;
        case FLUX_SCHEME_AUSM: // AUSM
            Compute_Flux_AUSM(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        case FLUX_SCHEME_VANLEER: // Van Leer
            Compute_Flux_VanLeer(node_L, node_R, areavec, Flux_Conv, AddTime);
            for (int i = 0; i < NEQUATIONS; i++)
                Flux_Diss[i] = 0.0;
            break;
        case FLUX_SCHEME_LDFSS: // LDFSS
            Compute_Flux_LDFSS(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        case FLUX_SCHEME_OSHER: // Osher
            Compute_Flux_Osher(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        case FLUX_SCHEME_STEGERWARMING: // Steger Warming
            Compute_Flux_StegerWarming(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        case FLUX_SCHEME_JST: // JST
            Compute_Flux_JST(node_L, node_R, areavec, Flux_Conv, Flux_Diss, AddTime);
            break;
        default:
            error("Compute_Flux: Invalid Flux Scheme - %d", FluxScheme);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Transformed_Preconditioner_Matrix(int nodeID, int reqInv, double **PrecondMatrix) {
    switch (FluxScheme) {
        case FLUX_SCHEME_ROE: // Roe
            Compute_Transformed_Preconditioner_Matrix_Roe(nodeID, reqInv, PrecondMatrix);
            break;
        case FLUX_SCHEME_HLLC: // HLLC
            Compute_Transformed_Preconditioner_Matrix_HLLC(nodeID, reqInv, PrecondMatrix);
            break;
        case FLUX_SCHEME_AUSM: // AUSM
            // Nothing to do
            break;
        case FLUX_SCHEME_VANLEER: // Van Leer
            // Nothing to do
            break;
        case FLUX_SCHEME_LDFSS: // LDFSS
            // Nothing to do
            break;
        case FLUX_SCHEME_OSHER: // Osher
            // Nothing to do
            break;
        case FLUX_SCHEME_STEGERWARMING: // Steger Warming
            // Nothing to do
            break;
        case FLUX_SCHEME_JST: // JST
            // Nothing to do
            break;
        default:
            error("Compute_Transformed_Preconditioner_Matrix: Invalid Flux Scheme - %d", FluxScheme);
            break;
    }
}

