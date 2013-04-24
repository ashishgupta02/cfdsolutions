/*******************************************************************************
 * File:        Jacobian.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _JACOBIAN_H
#define	_JACOBIAN_H

// Finite Difference Method
void FDJac_Init(void);
void FDJac_Finalize(void);
void FDJac_Reset(void);
void Compute_Jacobian_FiniteDifference(int AddTime, int Iteration);

// Roe Jacobian Methods
void RoeJac_Init(void);
void RoeJac_Finalize(void);
void RoeJac_Reset(void);
void Compute_Jacobian_Roe_Exact(int AddTime, int Iteration);
void Compute_Jacobian_Roe_Approximate(int AddTime, int Iteration);

// HLLC Jacobian Methods
void HLLCJac_Init(void);
void HLLCJac_Finalize(void);
void HLLCJac_Reset(void);
void Compute_Jacobian_HLLC_Exact(int AddTime, int Iteration);
void Compute_Jacobian_HLLC_Approximate(int AddTime, int Iteration);

// AUSM Jacobian Methods
void AUSMJac_Init(void);
void AUSMJac_Finalize(void);
void AUSMJac_Reset(void);
void Compute_Jacobian_AUSM_Exact(int AddTime, int Iteration);
void Compute_Jacobian_AUSM_Approximate(int AddTime, int Iteration);

// VanLeer Jacobian Methods
void VanLeerJac_Init(void);
void VanLeerJac_Finalize(void);
void VanLeerJac_Reset(void);
void Compute_Jacobian_VanLeer_Exact(int AddTime, int Iteration);
void Compute_Jacobian_VanLeer_Approximate(int AddTime, int Iteration);

// LDFSS Jacobian Methods
void LDFSSJac_Init(void);
void LDFSSJac_Finalize(void);
void LDFSSJac_Reset(void);
void Compute_Jacobian_LDFSS_Exact(int AddTime, int Iteration);
void Compute_Jacobian_LDFSS_Approximate(int AddTime, int Iteration);

// Osher Jacobian Methods
void OsherJac_Init(void);
void OsherJac_Finalize(void);
void OsherJac_Reset(void);
void Compute_Jacobian_Osher_Exact(int AddTime, int Iteration);
void Compute_Jacobian_Osher_Approximate(int AddTime, int Iteration);

// StegerWarming Jacobian Methods
void StegerWarmingJac_Init(void);
void StegerWarmingJac_Finalize(void);
void StegerWarmingJac_Reset(void);
void Compute_Jacobian_StegerWarming_Exact(int AddTime, int Iteration);
void Compute_Jacobian_StegerWarming_Approximate(int AddTime, int Iteration);

// JST Jacobian Methods
void JSTJac_Init(void);
void JSTJac_Finalize(void);
void JSTJac_Reset(void);
void Compute_Jacobian_JST_Exact(int AddTime, int Iteration);
void Compute_Jacobian_JST_Approximate(int AddTime, int Iteration);

// Main Jacobian Functions
void Jacobian_Init(void);
void Jacobian_Finalize(void);
void Compute_Jacobian(int AddTime, int Iteration);

#endif	/* _JACOBIAN_H */

