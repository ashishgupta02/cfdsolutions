/*******************************************************************************
 * File:        Euler2D_Solver_Finite_Difference.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EULER2D_SOLVER_FINITE_DIFFERENCE_H
#define	_EULER2D_SOLVER_FINITE_DIFFERENCE_H

typedef struct {
    double Resi_Old[4];
    double Resi_New[4];
} FD_NODE;

#endif	/* _EULER2D_SOLVER_FINITE_DIFFERENCE_H */

