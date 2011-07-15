/*******************************************************************************
 * File:        Euler2D_Solver_Roe.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _EULER2D_SOLVER_ROE_H
#define	_EULER2D_SOLVER_ROE_H

#include "MC.h"
#include "Euler2D_Mesh.h"
#include "Euler2D_Solver.h"

class Euler2D_Solver_Roe:
    virtual public Euler2D_Mesh,
    virtual public Euler2D_Solver {
protected:
    // Data
public:
    Euler2D_Solver_Roe();
    virtual ~Euler2D_Solver_Roe();
    void Get_Solver_Inputs_Roe(const char* FileName);
    void Solver_Prepare();
    void Solve();
    void Solver_Finalize();
protected:
    // Initialize the ROE Solver
    void Initialize_Solver_Roe();

    // Initialize the Solution
    void Initialize_Solution();

    // Compute Gauss Gradients
    void Compute_Gauss_Gradient();

    // Create and Initialize CRS Block Matrix
    void Create_CRS_SolverBlockMatrix();

    // Compute Residuals across Median Dual
    void Compute_Residual();

    // Compute Boundary Component of Residual
    void Compute_Boundary_Residual();

    // Compute Local Time Stepping deltaT
    void Compute_DeltaTime(int Iteration);

    // Compute and Fill CRS Matrix
    void Compute_CRS_SolverBlockMatrix(int AddTime);

    // Compute Boundary Component for CRS Matrix
    void Compute_Boundary_CRS_SolverBlockMatrix();

    // Update the Solution
    void Update_Solution();

    // Compute RMS
    double Compute_RMS();

    // Compute Finite Difference Jacobian
    void Compute_FD_Jacobian(int Iteration);
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_ROE_H */

