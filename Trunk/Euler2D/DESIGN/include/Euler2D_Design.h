/*
 * File:   Euler2D_Solver_Design.h
 * Author: Ashish Gupta
 *
 * Created on October 5, 2010, 4:53 PM
 */

#ifndef _EULER2D_SOLVER_DESIGN_H
#define	_EULER2D_SOLVER_DESIGN_H

#include "Euler2D_Solver_VanLeer.h"

class Euler2D_Solver_Design: virtual public Euler2D_Solver_VanLeer {
protected:
    double *dXdBeta;
    double *dYdBeta;
    double **dQdBeta;
    double **dRdX_dXdBeta;
    double **dIdQ;
    double dIdBeta;
    double dIdQ_dQdBeta;
    double dIdX_dXdBeta;
    // Matrix Computation CRS
    MC_CRS DesignBlockMatrix;
    // Cost Related Variables
    double Ref_CoeffLift;
    double Ref_CoeffDrag;
    double CoeffLift;
    double CoeffDrag;
    double Lift;
    double Drag;
    double I;
public:
    Euler2D_Solver_Design();
    virtual ~Euler2D_Solver_Design();
    void Design();
protected:
    // Initialize the Design
    void Initialize_Design();

    // Create and Initialize CRS Design Block Matrix
    void Create_CRS_DesignBlockMatrix();

    void Read_dXdBeta(const char* FileName);
    void Write_dXdBeta(const char* FileName);
    void Write_dQdBeta(const char* FileName);
    void Write_dRdX_dXdBeta(const char* FileName);

    // Compute (dR/dX)*(dX/dBeta)
    void Compute_dRdX_dXdBeta();
    void Compute_Boundary_dRdX_dXdBeta();
    void Compute_dFplusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FplusdBeta);
    void Compute_dFminusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FminusdBeta);
    void Verify_dRdX_dXdBeta();

    // Compute dQdBeta
    void Compute_dQdBeta();
    void Verify_dQdBeta();

    // Compute Cost: Coefficient of Lift and Drag Based
    void Compute_Cost();
    void Compute_dIdQ_dQdBeta();
    void Compute_dIdX_dXdBeta();
    void Compute_dIdBeta();
    
    // Back Up Code
    void Compute_dIdQ_bak();
    void Compute_Cost_bak();
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_DESIGN_H */
