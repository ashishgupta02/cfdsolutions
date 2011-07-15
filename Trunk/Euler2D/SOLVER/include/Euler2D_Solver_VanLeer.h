/*******************************************************************************
 * File:        Euler2D_Solver_VanLeer.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _EULER2D_SOLVER_VANLEER_H
#define	_EULER2D_SOLVER_VANLEER_H

#include "MC.h"
#include "Euler2D_Mesh.h"
#include "Euler2D_Solver.h"

class Euler2D_Solver_VanLeer:
    virtual public Euler2D_Mesh,
    virtual public Euler2D_Solver {
protected:
    // Data
public:
    Euler2D_Solver_VanLeer();
    virtual ~Euler2D_Solver_VanLeer();
    void Get_Solver_Inputs_VanLeer(const char* FileName);
    void Solver_Prepare();
    void Solve();
    void Solver_Finalize();
protected:
    // Initialize the VanLeer Solver
    void Initialize_Solver_VanLeer();
    
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

    // Get the Fluxes: F, F+, F-
    void Compute_Flux(double *Q, double uNx, double uNy, double Mag, double *Flux);
    void Compute_Fplus(double *Q, double uNx, double uNy, double Mag, double *Fplus);
    void Compute_Fminus(double *Q, double uNx, double uNy, double Mag, double *Fminus);
    void Compute_Fplus_Manual(double *Q, double uNx, double uNy, double Mag, double *Fplus);
    void Compute_Fminus_Manual(double *Q, double uNx, double uNy, double Mag, double *Fminus);

    // Get the Flux Jacobians: A, A+, A-
    void Compute_FluxJacobian(double *Q, double uNx, double uNy, double Mag, double *Flux, double **A);
    void Compute_FluxJplus(double *Q, double uNx, double uNy, double Mag, double *Fplus, double **Ap);
    void Compute_FluxJminus(double *Q, double uNx, double uNy, double Mag, double *Fminus, double **Am);
    void Compute_FluxJplus_Manual(double *Q, double uNx, double uNy, double Mag, double *Fplus, double **Ap);
    void Compute_FluxJminus_Manual(double *Q, double uNx, double uNy, double Mag, double *Fminus, double **Am);
    
    // Get Van Leer Flux and Jacobians
    void Compute_Flux_Jacobians(double *Q, double uNx, double uNy, double Mag, double *Flux, double *Fplus, double *Fminus, double **A, double **Ap, double **Am);
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_VANLEER_H */

