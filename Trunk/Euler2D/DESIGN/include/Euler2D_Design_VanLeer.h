/*******************************************************************************
 * File:        Euler2D_Design_VanLeer.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EULER2D_DESIGN_VANLEER_H
#define	_EULER2D_DESIGN_VANLEER_H

#include "Euler2D_Design.h"
#include "Euler2D_Solver_VanLeer.h"

class Euler2D_Design_VanLeer:
        virtual public Euler2D_Solver_VanLeer,
        virtual public Euler2D_Design {
protected:
    // Data
public:
    Euler2D_Design_VanLeer();
    virtual ~Euler2D_Design_VanLeer();
    void Get_Design_Inputs_VanLeer(const char* FileName);
    void Design_Prepare();
    void Design();
    void Design_Finalize();
protected:
    // Initialize the Design
    void Initialize_Design_VanLeer();
    void Compute_dFplusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FplusdBeta);
    void Compute_dFminusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FminusdBeta);

    // Design Solver Cost
    void Design_Cost();
    // Design Solver: Direct/Forward Mode and Adjoint
    void Design_Direct();
    void Design_Adjoint();

    // Update the Mesh with Perturbation and Smooth
    void Update_Mesh_LinearElasticSmooth();

    // Update the Solution After Mesh Smooth
    void Compute_New_Solution();
    
    // Common Functions to Direct and Adjoint
    // --- START
    // Compute Cost: Coefficient of Lift and Drag Based
    void Compute_Cost();
    // Compute (dR/dX)*(dX/dBeta)
    void Compute_dRdX_dXdBeta();
    void Compute_Boundary_dRdX_dXdBeta();
    void Verify_dRdX_dXdBeta();
    void Compute_dIdX_dXdBeta(int iDesignVariable);
    void Compute_LinearElasticSmooth_dXdBeta(int cNode);
    // --- END
    
    // Compute Direct/Forward Mode Functions
    // Compute dQdBeta
    void Compute_Direct_dQdBeta();
    void Verify_Direct_dQdBeta();
    void Compute_Direct_dIdQ_dQdBeta(int iDesignVariable);
    void Compute_Direct_dIdBeta(int iDesignVariable);
    
    // Compute Adjoint: Functions
    void Compute_Adjoint_dIdQ();
    void Compute_Adjoint_Lambda();
    void Compute_Adjoint_dIdBeta(int iDesignVariable);
private:
    void Init();
};

#endif	/* _EULER2D_DESIGN_VANLEER_H */

