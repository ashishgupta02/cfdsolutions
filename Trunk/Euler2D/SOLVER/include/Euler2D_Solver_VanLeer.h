/* 
 * File:   Euler2D_Solver_VanLeer.h
 * Author: Ashish Gupta
 *
 * Created on March 21, 2010, 4:31 PM
 */

#ifndef _EULER2D_SOLVER_VANLEER_H
#define	_EULER2D_SOLVER_VANLEER_H

#include "MC.h"
#include "Euler2D_Mesh.h"

class Euler2D_Solver_VanLeer: virtual public Euler2D_Mesh {
protected:
    // Iterations
    int NIteration;
    int InnerNIteration;
    double Relaxation;
    double RMS[4];
    // Constants
    double Gamma;
    // Reference Conditions
    double Ref_Mach;
    double Ref_Alpha;
    double Ref_Pressure;
    // CFL Number
    int CFL_Ramp;
    double CFL_Min;
    double CFL_Max;
    // Local Time
    double *DeltaT;
    // Matrix Computation Compressed Row Storage
    MC_CRS BlockMatrix;
public:
    Euler2D_Solver_VanLeer();
    virtual ~Euler2D_Solver_VanLeer();
    void Solve();
    void Get_Reference_Conditions();
protected:
    // Initialize the Solution
    void Initialize_Solution();

    // Compute Residuals across Median Dual
    void Compute_Residual();

    // Compute Boundary Component of Residual
    void Compute_Boundary_Residual();

    // Compute Local Time Stepping deltaT
    void Compute_DeltaTime(int Iteration);

    // Create and Initialize CRS Block Matrix
    void Create_CRS_BlockMatrix();

    // Compute and Fill CRS Matrix
    void Compute_CRS_BlockMatrix();

    // Compute Boundary Component for CRS Matrix
    void Compute_Boundary_CRS_BlockMatrix();

    // Update the Solution
    void Update_Solution();

    // Compute RMS
    double Compute_RMS();

    // Get the Fluxes: F, F+, F-
    void Compute_Flux(double *Q, double uNx, double uNy, double Mag, double *Flux);
    void Compute_Fplus(double *Q, double uNx, double uNy, double Mag, double *Fplus);
    void Compute_Fminus(double *Q, double uNx, double uNy, double Mag, double *Fminus);

    // Get the Flux Jacobians: A, A+, A-
    void Compute_FluxJacobian(double *Q, double uNx, double uNy, double Mag, double *Flux, double **A);
    void Compute_FluxJplus(double *Q, double uNx, double uNy, double Mag, double *Fplus, double **Ap);
    void Compute_FluxJminus(double *Q, double uNx, double uNy, double Mag, double *Fminus, double **Am);
    
    // Get Van Leer Flux and Jacobians
    void Compute_Flux_Jacobians(double *Q, double uNx, double uNy, double Mag, double *Flux, double *Fplus, double *Fminus, double **A, double **Ap, double **Am);
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_VANLEER_H */
