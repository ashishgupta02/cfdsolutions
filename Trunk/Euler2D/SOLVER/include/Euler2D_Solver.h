/*******************************************************************************
 * File:        Euler2D_Solver.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _EULER2D_SOLVER_H
#define	_EULER2D_SOLVER_H

#include <string>

#include "MC.h"
#include "Euler2D_Solver_Finite_Difference.h"

class Euler2D_Solver {
protected:
    int SolverBlockSize;
    int SolverVectorSize;   // No of Nodes
    int DoSolverValidation; // 0: No Validation, 1: Validation
    // Iterations
    int Order;
    int NIteration;
    int InnerNIteration;
    int DoRestart;
    double Relaxation;
    double RMS[4];
    double RMS_Res;
    // Constants
    double Gamma;
    // Reference Conditions
    double Ref_Mach;
    double Ref_Alpha;
    double Ref_Pressure;
    double Ref_Length;
    // CFL Number
    int CFL_Ramp;
    double CFL_Min;
    double CFL_Max;
    // Local Time
    double *DeltaT;
    // Gradients
    double **Qx;
    double **Qy;
    // Matrix Computation Compressed Row Storage
    MC_CRS SolverBlockMatrix;
    // Finite Difference Jacobian
    int     FDIteration;
    int     FDNodeID;
    int     FDPertQ;
    double  FDEpsilon;
    double  FDDR_DQ[4];
    FD_NODE *FDNode;

    // Input and Outputs
    std::string WKAMeshFileName;
    std::string InputRestartFileName;
    std::string OutputRestartFileName;
    std::string SLKMeshFileName;
    std::string VTKSolutionFileName;
public:
    Euler2D_Solver();
    virtual ~Euler2D_Solver();
    virtual void Solver_Prepare()=0;
    virtual void Solve()=0;
    virtual void Solver_Finalize()=0;
protected:
    // Get the Generic Solver Inputs
    void Get_Solver_Inputs(const char* FileName);
    
    // Initialize the Solver
    void Initialize_Solver(int InputBlockSize, int InputVectorSize);

    // Create and Initialize CRS Block Matrix
    virtual void Create_CRS_SolverBlockMatrix()=0;

    // Compute Residuals across Median Dual
    virtual void Compute_Residual()=0;

    // Compute Boundary Component of Residual
    virtual void Compute_Boundary_Residual()=0;

    // Compute Local Time Stepping deltaT
    virtual void Compute_DeltaTime(int Iteration)=0;

    // Compute and Fill CRS Matrix
    virtual void Compute_CRS_SolverBlockMatrix(int AddTime)=0;

    // Compute Boundary Component for CRS Matrix
    virtual void Compute_Boundary_CRS_SolverBlockMatrix()=0;

    // Update the Solution
    virtual void Update_Solution()=0;

    // Compute RMS
    virtual double Compute_RMS()=0;

    // Compute Finite Difference Jacobian
    virtual void Compute_FD_Jacobian(int Iteration)=0;
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_H */

