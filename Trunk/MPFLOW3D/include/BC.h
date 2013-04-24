/*******************************************************************************
 * File:        BC.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _BC_H
#define	_BC_H

// Boundary Conditions
void BC_Init(void);
void BC_Finalize(void);
void BC_Reset(void);

void Initialize_Boundary_Condition(void);
void Apply_Boundary_Condition(int Iteration);
void Apply_Boundary_Condition(int BEdgeID, int Iteration);

#endif	/* _BC_H */

