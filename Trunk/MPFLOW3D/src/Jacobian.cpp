/*******************************************************************************
 * File:        Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Solver.h"
#include "Roe_Fluxes.h"
#include "HLLC_Fluxes.h"
#include "VanLeer_Fluxes.h"
#include "AUSM_Fluxes.h"
#include "JST_Fluxes.h"
#include "LDFSS_Fluxes.h"
#include "Osher_Fluxes.h"
#include "StegerWarming_Fluxes.h"
#include "Jacobian.h"

//------------------------------------------------------------------------------
//! Initialize the Jacobian Data Structure
//------------------------------------------------------------------------------
void Jacobian_Init(void) {
    switch (JacobianMethod) {
        case JACOBIAN_METHOD_CENTRAL:
            FDJac_Init();
            break;
        case JACOBIAN_METHOD_FORWARD:
            FDJac_Init();
            break;
        case JACOBIAN_METHOD_BACKWARD:
            FDJac_Init();
            break;
        case JACOBIAN_METHOD_ALTERNATE:
            FDJac_Init();
            break;
        case JACOBIAN_METHOD_APPROX:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    RoeJac_Init();
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    RoeJac_Init();
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    AUSMJac_Init();
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    VanLeerJac_Init();
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    LDFSSJac_Init();
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    OsherJac_Init();
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    StegerWarmingJac_Init();
                    break;
                case FLUX_SCHEME_JST: // JST
                    JSTJac_Init();
                    break;
                default:
                    error("Jacobian_Init:1: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        case JACOBIAN_METHOD_EXACT:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    RoeJac_Init();
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    RoeJac_Init();
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    AUSMJac_Init();
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    VanLeerJac_Init();
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    LDFSSJac_Init();
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    OsherJac_Init();
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    StegerWarmingJac_Init();
                    break;
                case FLUX_SCHEME_JST: // JST
                    JSTJac_Init();
                    break;
                default:
                    error("Jacobian_Init:2: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        default:
            error("Jacobian_Init:3: Invalid Jacobian Method - %d", JacobianMethod);
            break;
    }
}

//------------------------------------------------------------------------------
//! Finalize the Jacobian Data Structure
//------------------------------------------------------------------------------
void Jacobian_Finalize(void) {
    switch (JacobianMethod) {
        case JACOBIAN_METHOD_CENTRAL:
            FDJac_Finalize();
            break;
        case JACOBIAN_METHOD_FORWARD:
            FDJac_Finalize();
            break;
        case JACOBIAN_METHOD_BACKWARD:
            FDJac_Finalize();
            break;
        case JACOBIAN_METHOD_ALTERNATE:
            FDJac_Finalize();
            break;
        case JACOBIAN_METHOD_APPROX:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    RoeJac_Finalize();
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    RoeJac_Finalize();
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    AUSMJac_Finalize();
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    VanLeerJac_Finalize();
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    LDFSSJac_Finalize();
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    OsherJac_Finalize();
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    StegerWarmingJac_Finalize();
                    break;
                case FLUX_SCHEME_JST: // JST
                    JSTJac_Finalize();
                    break;
                default:
                    error("Jacobian_Finalize:1: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        case JACOBIAN_METHOD_EXACT:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    RoeJac_Finalize();
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    RoeJac_Finalize();
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    AUSMJac_Finalize();
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    VanLeerJac_Finalize();
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    LDFSSJac_Finalize();
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    OsherJac_Finalize();
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    StegerWarmingJac_Finalize();
                    break;
                case FLUX_SCHEME_JST: // JST
                    JSTJac_Finalize();
                    break;
                default:
                    error("Jacobian_Finalize:2: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        default:
            error("Jacobian_Finalize:3: Invalid Jacobian Method - %d", JacobianMethod);
            break;
    }
}

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
    
    switch (JacobianMethod) {
        case JACOBIAN_METHOD_CENTRAL:
            Compute_Jacobian_FiniteDifference(AddTime, Iteration);
            break;
        case JACOBIAN_METHOD_FORWARD:
            Compute_Jacobian_FiniteDifference(AddTime, Iteration);
            break;
        case JACOBIAN_METHOD_BACKWARD:
            Compute_Jacobian_FiniteDifference(AddTime, Iteration);
            break;
        case JACOBIAN_METHOD_ALTERNATE:
            Compute_Jacobian_FiniteDifference(AddTime, Iteration);
            break;
        case JACOBIAN_METHOD_APPROX:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    Compute_Jacobian_Roe_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    Compute_Jacobian_Roe_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    Compute_Jacobian_AUSM_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    Compute_Jacobian_VanLeer_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    Compute_Jacobian_LDFSS_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    Compute_Jacobian_Osher_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    Compute_Jacobian_StegerWarming_Approximate(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_JST: // JST
                    Compute_Jacobian_JST_Approximate(AddTime, Iteration);
                    break;
                default:
                    error("Compute_Jacobian:1: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        case JACOBIAN_METHOD_EXACT:
            switch (FluxScheme) {
                case FLUX_SCHEME_ROE: // Roe
                    Compute_Jacobian_Roe_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_HLLC: // HLLC
                    Compute_Jacobian_HLLC_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_AUSM: // AUSM
                    Compute_Jacobian_AUSM_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_VANLEER: // Van Leer
                    Compute_Jacobian_VanLeer_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_LDFSS: // LDFSS
                    Compute_Jacobian_LDFSS_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_OSHER: // Osher
                    Compute_Jacobian_Osher_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_STEGERWARMING: // Steger Warming
                    Compute_Jacobian_StegerWarming_Exact(AddTime, Iteration);
                    break;
                case FLUX_SCHEME_JST: // JST
                    Compute_Jacobian_JST_Exact(AddTime, Iteration);
                    break;
                default:
                    error("Compute_Jacobian:2: Invalid Flux Scheme - %d", FluxScheme);
                    break;
            }
            break;
        default:
            error("Compute_Jacobian:3: Invalid Jacobian Method - %d", JacobianMethod);
            break;
    }
    
    // For Unsteady Computations
    if (SolverMethod == SOLVER_METHOD_UNSTEADY) {
        
    }
    
    // Switch Back to Original Solver Order
    SolverOrder = saveorder;
}

