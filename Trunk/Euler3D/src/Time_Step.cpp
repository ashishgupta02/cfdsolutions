/*******************************************************************************
 * File:        Time_Step.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
//void Compute_DeltaT_Old(int Iteration) {
//    int i, nodeL, nodeR;
//    double denom = 0.0;
//    
//    Vector3D areavec;
//    double   area;
//
//    double rho_avg, u_avg, v_avg, w_avg, et_avg, e_avg, p_avg, c_avg, Ubar_avg;
//    double lamda1, lamda4, lamda5;
//    double max_lamda;
//    double eps, alpha, beta, beta_m, beta_p, Ur, mach, sigma;
//    double dtmp;
//    int j, nid;
//    double p_L, p_R, lp, lc, dp_L, dp_R, dp, lmach_L, lmach_R, lmach;
//    
//    // Internal Edges
//    for (i = 0; i < nEdge; i++) {
//        // Get two nodes of edge
//        nodeL = intEdge[i].node[0];
//        nodeR = intEdge[i].node[1];
//
//#ifdef DEBUG
//        assert(nodeR > nodeL);
//#endif
//        
//        // Get area vector
//        areavec = intEdge[i].areav;
//        area  = areavec.magnitude();
//        areavec.normalize();
//
//        // Get the average
//        // Conservative Variable Formulation
//        if (Variable_Type == VARIABLE_CONSERVATIVE) {
//            rho_avg  = 0.5*(Q1[nodeL] + Q1[nodeR]);
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR])/rho_avg;
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR])/rho_avg;
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR])/rho_avg;           
//            et_avg   = 0.5*(Q5[nodeL] + Q5[nodeR])/rho_avg;
//            e_avg    = et_avg - 0.5*(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            p_avg    = rho_avg*e_avg*(Gamma - 1.0);
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        // Primitive Variable Formulation Pressure Velocity Temperature
//        if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//            p_avg    = 0.5*(Q1[nodeL] + Q1[nodeR]) + Gauge_Pressure;
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR]);
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR]);
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR]);
//            rho_avg  = p_avg/(Ref_R*(0.5*(Q5[nodeL] + Q5[nodeR]) + Gauge_Temperature));
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        // Primitive Variable Formulation Density Velocity Pressure
//        if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//            rho_avg  = 0.5*(Q1[nodeL] + Q1[nodeR]);
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR]);
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR]);
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR]);
//            p_avg    = 0.5*(Q5[nodeL] + Q5[nodeR]) + Gauge_Pressure;
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        
//        // Compute Weiss Smith Precondition Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
//            // Compute Precondition Variables
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            eps  = MIN(1.0, MAX(1.0e-10, mach*mach));
//
//            // Compute beta for Ideal Gas
//            beta   = 1.0/(c_avg*c_avg);
//
//            // Compute Ur for Ideal Gas
//            if (dtmp < eps*c_avg)
//                Ur = eps*c_avg;
//            else if ((eps*c_avg < dtmp) && (dtmp < c_avg))
//                Ur = dtmp;
//            else
//                Ur = c_avg;
//
//            // Compute alpha
//            alpha  = 0.5*(1.0 - beta*Ur*Ur);
//            
//            // Compute U'
//            Ubar_avg = Ubar_avg*(1.0 - alpha);
//            
//            // Compute C'
//            c_avg = sqrt(alpha*alpha*Ubar_avg*Ubar_avg + Ur*Ur);
//        }
//        
//        // Compute Cecile Voizat Precondition Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
//            // Conservative Variable Formulation
//            if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                p_L = (Gamma - 1.0)*(Q5[nodeL] - 0.5*((Q2[nodeL]*Q2[nodeL] + Q3[nodeL]*Q3[nodeL] + Q4[nodeL]*Q4[nodeL])/Q1[nodeL]));
//                p_R = (Gamma - 1.0)*(Q5[nodeR] - 0.5*((Q2[nodeR]*Q2[nodeR] + Q3[nodeR]*Q3[nodeR] + Q4[nodeR]*Q4[nodeR])/Q1[nodeR]));
//            }
//            // Primitive Variable Formulation Pressure Velocity Temperature
//            if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                p_L = Q1[nodeL] + Gauge_Pressure;
//                p_R = Q1[nodeR] + Gauge_Pressure;
//            }
//            // Primitive Variable Formulation Density Velocity Pressure
//            if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                p_L = Q5[nodeL] + Gauge_Pressure;
//                p_R = Q5[nodeR] + Gauge_Pressure;
//            }
//            
//            dp_L = dp_R = 0.0;
//            lmach_L = lmach_R = 0.0;
//            for (j = crs_IA_Node2Node[nodeL]; j < crs_IA_Node2Node[nodeL+1]; j++) {
//                nid     = crs_JA_Node2Node[j];
//                // Conservative Variable Formulation
//                if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                    lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                }
//                // Primitive Variable Formulation Pressure Velocity Temperature
//                if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                    lp      = Q1[nid] + Gauge_Pressure;
//                    dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                    lc      = sqrt((Gamma * lp) / dtmp);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Primitive Variable Formulation Density Velocity Pressure
//                if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                    lp      = Q5[nid] + Gauge_Pressure;
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Get the Max Change in Pressure
//                dp_L = MAX(dp_L, fabs(p_L - lp));
//            }
//            dp    = dp_L;
//            lmach = lmach_L;
//            // Check if Right Node is not a Ghost Boundary Node
//            if (nodeR < nNode) {
//                for (j = crs_IA_Node2Node[nodeR]; j < crs_IA_Node2Node[nodeR+1]; j++) {
//                    nid  = crs_JA_Node2Node[j];
//                    // Conservative Variable Formulation
//                    if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                        lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                        lc      = sqrt((Gamma * lp) / Q1[nid]);
//                        lmach_R = MAX(lmach_R, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                    }
//                    // Primitive Variable Formulation Pressure Velocity Temperature
//                    if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                        lp      = Q1[nid] + Gauge_Pressure;
//                        dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                        lc      = sqrt((Gamma * lp) / dtmp);
//                        lmach_R = MAX(lmach_R, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                    }
//                    // Primitive Variable Formulation Density Velocity Pressure
//                    if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                        lp      = Q5[nid] + Gauge_Pressure;
//                        lc      = sqrt((Gamma * lp) / Q1[nid]);
//                        lmach_R = MAX(lmach_R, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                    }
//                    // Get the Max Change in Pressure
//                    dp_R = MAX(dp_R, fabs(p_R - lp));
//                }
//                dp    = 0.5*(dp_L + dp_R);
//                lmach = 0.5*(lmach_L + lmach_R);
//            }
//            
//            // Compute Precondition Variables
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            sigma = 1.0;
//            if (Inf_Mach < 0.5) {
//                sigma = MIN(1.0, MAX(MAX(MAX(1.0e-6, mach), lmach), 2.0*dp/(rho_avg*c_avg*c_avg)));
//               //sigma = MAX(mach, dp/(rho_avg*c_avg*c_avg));
//               //sigma = MAX(1.0e-6, MIN(lmach, dp/(Inf_Mach*Inf_Mach*rho_avg*c_avg*c_avg)));
//                if (sigma > Inf_Mach)
//                    sigma  = Inf_Mach;
//            }
//                       
//            // Compute C'
//            c_avg = 0.5*sqrt(4.0*c_avg*c_avg*sigma*sigma + (1.0 - sigma*sigma)*(1.0 - sigma*sigma)*Ubar_avg*Ubar_avg);
//            
//            // Compute U'
//            Ubar_avg = 0.5*(1.0 + sigma*sigma)*Ubar_avg;
//        }
//        
//        //====================
//        // Compute Briley Taylor Whitfield Pre-Conditioner Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_BTW) {
//            // Conservative Variable Formulation
//            if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                p_L = (Gamma - 1.0)*(Q5[nodeL] - 0.5*((Q2[nodeL]*Q2[nodeL] + Q3[nodeL]*Q3[nodeL] + Q4[nodeL]*Q4[nodeL])/Q1[nodeL]));
//                p_R = (Gamma - 1.0)*(Q5[nodeR] - 0.5*((Q2[nodeR]*Q2[nodeR] + Q3[nodeR]*Q3[nodeR] + Q4[nodeR]*Q4[nodeR])/Q1[nodeR]));
//            }
//            // Primitive Variable Formulation Pressure Velocity Temperature
//            if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                p_L = Q1[nodeL] + Gauge_Pressure;
//                p_R = Q1[nodeR] + Gauge_Pressure;
//            }
//            // Primitive Variable Formulation Density Velocity Pressure
//            if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                p_L = Q5[nodeL] + Gauge_Pressure;
//                p_R = Q5[nodeR] + Gauge_Pressure;
//            }
//            
//            dp_L = dp_R = 0.0;
//            lmach_L = lmach_R = 0.0;
//            for (j = crs_IA_Node2Node[nodeL]; j < crs_IA_Node2Node[nodeL+1]; j++) {
//                nid     = crs_JA_Node2Node[j];
//                // Conservative Variable Formulation
//                if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                    lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                }
//                // Primitive Variable Formulation Pressure Velocity Temperature
//                if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                    lp      = Q1[nid] + Gauge_Pressure;
//                    dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                    lc      = sqrt((Gamma * lp) / dtmp);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Primitive Variable Formulation Density Velocity Pressure
//                if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                    lp      = Q5[nid] + Gauge_Pressure;
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Get the Max Change in Pressure
//                dp_L = MAX(dp_L, fabs(p_L - lp));
//            }
//            dp    = dp_L;
//            lmach = lmach_L;
//            // Check if Right Node is not a Ghost Boundary Node
//            if (nodeR < nNode) {
//                for (j = crs_IA_Node2Node[nodeR]; j < crs_IA_Node2Node[nodeR+1]; j++) {
//                    nid  = crs_JA_Node2Node[j];
//                    // Conservative Variable Formulation
//                    if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                        lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                        lc      = sqrt((Gamma * lp) / Q1[nid]);
//                        lmach_R = MAX(lmach_R, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                    }
//                    // Primitive Variable Formulation Pressure Velocity Temperature
//                    if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                        lp      = Q1[nid] + Gauge_Pressure;
//                        dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                        lc      = sqrt((Gamma * lp) / dtmp);
//                        lmach_R = MAX(lmach_R, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                    }
//                    // Primitive Variable Formulation Density Velocity Pressure
//                    if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                        lp      = Q5[nid] + Gauge_Pressure;
//                        lc      = sqrt((Gamma * lp) / Q1[nid]);
//                        lmach_R = MAX(lmach_R, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                    }
//                    // Get the Max Change in Pressure
//                    dp_R = MAX(dp_R, fabs(p_R - lp));
//                }
//                dp    = 0.5*(dp_L + dp_R);
//                lmach = 0.5*(lmach_L + lmach_R);
//            }
//            
//            // Compute the required parameters
//            beta = 1.0; // Corresponds to no pre-conditioning
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            if (Inf_Mach < 0.5) {
//                beta = MIN(1.0, MAX(MAX(MAX(1.0e-6, mach), lmach), 2.0*dp/(rho_avg*c_avg*c_avg)));
//                if (beta > Inf_Mach)
//                    beta = Inf_Mach;
//            }
//            beta   = beta*beta;
//            beta_m = 0.5*(1.0 - beta);
//            beta_p = 0.5*(1.0 + beta);
//            sigma  = sqrt(Ubar_avg*Ubar_avg*beta_m*beta_m + beta*c_avg*c_avg);
//            
//            // Compute Eigenvalues for BTW
//            lamda1 = fabs(Ubar_avg);
//            lamda4 = fabs(Ubar_avg*beta_p + sigma);
//            lamda5 = fabs(Ubar_avg*beta_p - sigma);
//        }
//        //====================
//        
//        // Compute Eigenvalues
//        if (SolverScheme != SOLVER_SCHEME_ROE_BTW) {
//            lamda1 = fabs(Ubar_avg);
//            lamda4 = fabs(Ubar_avg + c_avg);
//            lamda5 = fabs(Ubar_avg - c_avg);
//        }
//        
//        max_lamda = MAX(lamda1, MAX(lamda4, lamda5));
//        denom = (max_lamda * area);
//
//        // Add to Each Node
//        DeltaT[nodeL] += denom;
//        DeltaT[nodeR] += denom;
//    }
//
//    // Boundary Edges
//    for (i = 0; i < nBEdge; i++) {
//        // Get two nodes of edge
//        nodeL = bndEdge[i].node[0];
//        nodeR = bndEdge[i].node[1];
//
//#ifdef DEBUG
//        assert(nodeR > nodeL);
//#endif
//        
//        // Get area vector
//        areavec = bndEdge[i].areav;
//        area    = areavec.magnitude();
//        areavec.normalize();
//
//        // Note we only need to accumulate for the left node
//        // Conservative Variable Formulation
//        if (Variable_Type == VARIABLE_CONSERVATIVE) {
//            rho_avg  = 0.5*(Q1[nodeL] + Q1[nodeR]);
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR])/rho_avg;
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR])/rho_avg;
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR])/rho_avg;
//            et_avg   = 0.5*(Q5[nodeL] + Q5[nodeR])/rho_avg;
//            e_avg    = et_avg - 0.5*(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            p_avg    = rho_avg*e_avg*(Gamma - 1.0);
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        // Primitive Variable Formulation Pressure Velocity Temperature
//        if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//            p_avg    = 0.5*(Q1[nodeL] + Q1[nodeR]) + Gauge_Pressure;
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR]);
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR]);
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR]);
//            rho_avg  = p_avg/(Ref_R*(0.5*(Q5[nodeL] + Q5[nodeR]) + Gauge_Temperature));
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        // Primitive Variable Formulation Density Velocity Pressure
//        if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//            rho_avg  = 0.5*(Q1[nodeL] + Q1[nodeR]);
//            u_avg    = 0.5*(Q2[nodeL] + Q2[nodeR]);
//            v_avg    = 0.5*(Q3[nodeL] + Q3[nodeR]);
//            w_avg    = 0.5*(Q4[nodeL] + Q4[nodeR]);
//            p_avg    = 0.5*(Q5[nodeL] + Q5[nodeR]) + Gauge_Pressure;
//            c_avg    = sqrt(Gamma*p_avg/rho_avg);
//            Ubar_avg = u_avg*areavec.vec[0] + v_avg*areavec.vec[1] + w_avg*areavec.vec[2];
//        }
//        
//        // Compute Weiss Smith Precondition Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
//            // Compute Precondition Variables
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            eps  = MIN(1.0, MAX(1.0e-10, mach*mach));
//            
//            // Compute beta for Ideal Gas
//            beta   = 1.0/(c_avg*c_avg);
//
//            // Compute Ur for Ideal Gas
//            if (dtmp < eps*c_avg)
//                Ur = eps*c_avg;
//            else if ((eps*c_avg < dtmp) && (dtmp < c_avg))
//                Ur = dtmp;
//            else
//                Ur = c_avg;
//
//            // Compute alpha
//            alpha  = 0.5*(1.0 - beta*Ur*Ur);
//            
//            // Compute U'
//            Ubar_avg = Ubar_avg*(1.0 - alpha);
//            
//            // Compute C'
//            c_avg = sqrt(alpha*alpha*Ubar_avg*Ubar_avg + Ur*Ur);
//        }
//        
//        // Compute Cecile Voizat Precondition Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
//            // Conservative Variable Formulation
//            if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                p_L = (Gamma - 1.0)*(Q5[nodeL] - 0.5*((Q2[nodeL]*Q2[nodeL] + Q3[nodeL]*Q3[nodeL] + Q4[nodeL]*Q4[nodeL])/Q1[nodeL]));
//                p_R = (Gamma - 1.0)*(Q5[nodeR] - 0.5*((Q2[nodeR]*Q2[nodeR] + Q3[nodeR]*Q3[nodeR] + Q4[nodeR]*Q4[nodeR])/Q1[nodeR]));
//            }
//            // Primitive Variable Formulation Pressure Velocity Temperature
//            if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                p_L = Q1[nodeL] + Gauge_Pressure;
//                p_R = Q1[nodeR] + Gauge_Pressure;
//            }
//            // Primitive Variable Formulation Density Velocity Pressure
//            if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                p_L = Q5[nodeL] + Gauge_Pressure;
//                p_R = Q5[nodeR] + Gauge_Pressure;
//            }
//            
//            dp_L = dp_R = 0.0;
//            lmach_L = lmach_R = 0.0;
//            for (j = crs_IA_Node2Node[nodeL]; j < crs_IA_Node2Node[nodeL+1]; j++) {
//                nid     = crs_JA_Node2Node[j];
//                // Conservative Variable Formulation
//                if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                    lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                }
//                // Primitive Variable Formulation Pressure Velocity Temperature
//                if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                    lp      = Q1[nid] + Gauge_Pressure;
//                    dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                    lc      = sqrt((Gamma * lp) / dtmp);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Primitive Variable Formulation Density Velocity Pressure
//                if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                    lp      = Q5[nid] + Gauge_Pressure;
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Get the Max Change in Pressure
//                dp_L = MAX(dp_L, fabs(p_L - lp));
//            }
//            dp    = dp_L;
//            lmach = lmach_L;
//            
//            // Compute Precondition Variables
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            sigma = 1.0;
//            if (Inf_Mach < 0.5) {
//                sigma = MIN(1.0, MAX(MAX(MAX(1.0e-6, mach), lmach), 2.0*dp/(rho_avg*c_avg*c_avg)));
//                //sigma = MAX(mach, dp/(rho_avg*c_avg*c_avg));
//                //sigma = MAX(1.0e-6, MIN(lmach, dp/(Inf_Mach*Inf_Mach*rho_avg*c_avg*c_avg)));
//                if (sigma > Inf_Mach)
//                    sigma  = Inf_Mach;
//            }
//            
//            // Compute C'
//            c_avg = 0.5*sqrt(4.0*c_avg*c_avg*sigma*sigma + (1.0 - sigma*sigma)*(1.0 - sigma*sigma)*Ubar_avg*Ubar_avg);
//            
//            // Compute U'
//            Ubar_avg = 0.5*(1.0 + sigma*sigma)*Ubar_avg;
//        }
//        
//        //====================
//        // Compute Briley Taylor Whitfield Pre-Conditioner Variables
//        if (SolverScheme == SOLVER_SCHEME_ROE_BTW) {
//            // Conservative Variable Formulation
//            if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                p_L = (Gamma - 1.0)*(Q5[nodeL] - 0.5*((Q2[nodeL]*Q2[nodeL] + Q3[nodeL]*Q3[nodeL] + Q4[nodeL]*Q4[nodeL])/Q1[nodeL]));
//                p_R = (Gamma - 1.0)*(Q5[nodeR] - 0.5*((Q2[nodeR]*Q2[nodeR] + Q3[nodeR]*Q3[nodeR] + Q4[nodeR]*Q4[nodeR])/Q1[nodeR]));
//            }
//            // Primitive Variable Formulation Pressure Velocity Temperature
//            if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                p_L = Q1[nodeL] + Gauge_Pressure;
//                p_R = Q1[nodeR] + Gauge_Pressure;
//            }
//            // Primitive Variable Formulation Density Velocity Pressure
//            if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                p_L = Q5[nodeL] + Gauge_Pressure;
//                p_R = Q5[nodeR] + Gauge_Pressure;
//            }
//            
//            dp_L = dp_R = 0.0;
//            lmach_L = lmach_R = 0.0;
//            for (j = crs_IA_Node2Node[nodeL]; j < crs_IA_Node2Node[nodeL+1]; j++) {
//                nid     = crs_JA_Node2Node[j];
//                // Conservative Variable Formulation
//                if (Variable_Type == VARIABLE_CONSERVATIVE) {
//                    lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]);
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
//                }
//                // Primitive Variable Formulation Pressure Velocity Temperature
//                if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
//                    lp      = Q1[nid] + Gauge_Pressure;
//                    dtmp    = lp/(Ref_R*(Q5[nid] + Gauge_Temperature));
//                    lc      = sqrt((Gamma * lp) / dtmp);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Primitive Variable Formulation Density Velocity Pressure
//                if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
//                    lp      = Q5[nid] + Gauge_Pressure;
//                    lc      = sqrt((Gamma * lp) / Q1[nid]);
//                    lmach_L = MAX(lmach_L, sqrt(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/lc);
//                }
//                // Get the Max Change in Pressure
//                dp_L = MAX(dp_L, fabs(p_L - lp));
//            }
//            dp    = dp_L;
//            lmach = lmach_L;
//            
//            // Compute the required parameters
//            beta = 1.0; // Corresponds to no pre-conditioning
//            dtmp = sqrt(u_avg*u_avg + v_avg*v_avg + w_avg*w_avg);
//            mach = dtmp/c_avg;
//            if (Inf_Mach < 0.5) {
//                beta = MIN(1.0, MAX(MAX(MAX(1.0e-6, mach), lmach), 2.0*dp/(rho_avg*c_avg*c_avg)));
//                if (beta > Inf_Mach)
//                    beta = Inf_Mach;
//            }
//            beta   = beta*beta;
//            beta_m = 0.5*(1.0 - beta);
//            beta_p = 0.5*(1.0 + beta);
//            sigma  = sqrt(Ubar_avg*Ubar_avg*beta_m*beta_m + beta*c_avg*c_avg);
//                      
//            // Compute Eigenvalues for BTW
//            lamda1 = fabs(Ubar_avg);
//            lamda4 = fabs(Ubar_avg*beta_p + sigma);
//            lamda5 = fabs(Ubar_avg*beta_p - sigma);
//        }
//        //====================
//        
//        // Compute Eigenvalues
//        if (SolverScheme != SOLVER_SCHEME_ROE_BTW) {
//            lamda1 = fabs(Ubar_avg);
//            lamda4 = fabs(Ubar_avg + c_avg);
//            lamda5 = fabs(Ubar_avg - c_avg);
//        }
//        
//        max_lamda = MAX(lamda1, MAX(lamda4, lamda5));
//        denom = (max_lamda * area);
//        
//        // Add to Left Node Only
//        DeltaT[nodeL] += denom;
//    }
//    
//    // Finally Compute the Local Time Stepping
//    if ((CFL_Ramp > 1) && (CFL_MAX > CFL_MIN)) {
//        if (Iteration < CFL_Ramp)
//            CFL = CFL_MIN + (CFL_MAX - CFL_MIN)*(((double)Iteration)/((double)(CFL_Ramp-1)));
//        else
//            CFL = CFL_MAX;
//    } else
//        CFL = CFL_MAX;
//
//    for (i = 0; i < nNode; i++)
//        DeltaT[i] = cVolume[i] * CFL / DeltaT[i];
//
//    // Check if Global Time Stepping is Required
//    if (TimeAccuracy == 1) {
//        denom = DeltaT[0];
//        for (i = 0; i < nNode; i++)
//            denom = MIN(denom, DeltaT[i]);
//        for (i = 0; i < nNode; i++)
//            DeltaT[i] = denom;
//    }
//}


//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Compute_DeltaT(int Iteration) {
    int i;
    double denom = 0.0;
    
    // Finally Compute the Local Time Stepping
    if ((CFL_Ramp > 1) && (CFL_MAX > CFL_MIN)) {
        if (Iteration < CFL_Ramp)
            CFL = CFL_MIN + (CFL_MAX - CFL_MIN)*(((double)Iteration)/((double)(CFL_Ramp-1)));
        else
            CFL = CFL_MAX;
    } else
        CFL = CFL_MAX;

    for (i = 0; i < nNode; i++) {
        DeltaT[i] = cVolume[i] * CFL / DeltaT[i];
        MinDeltaT = MIN(MinDeltaT, DeltaT[i]);
        MaxDeltaT = MAX(MaxDeltaT, DeltaT[i]);
    }

    // Check if Global Time Stepping is Required
    if (TimeAccuracy == SOLVER_TIME_ACCURACY_GLOBAL) {
        denom = DeltaT[0];
        for (i = 0; i < nNode; i++)
            denom = MIN(denom, DeltaT[i]);
        for (i = 0; i < nNode; i++)
            DeltaT[i] = denom;
    }
}

