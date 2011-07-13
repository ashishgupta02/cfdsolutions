/*******************************************************************************
 * File:        Euler2D_Design.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _EULER2D_DESIGN_H
#define	_EULER2D_DESIGN_H

#include "MC.h"

class Euler2D_Design {
protected:
    int DesignAction; // 0: Cost 1: Gradient
    int BlockSize;
    int VectorSize;   // No of Nodes
    int DoValidation; // 0: No Validation, 1: Validation
    // 0: Coeff Lift, 1: Coeff Drag, 2: Coeff(Lift + Drag), 3: Pressure, 4: Quad Pressure
    int CostType; 
    int IsAdjoint;    // 0: Forward Mode, 1: Adjoint Mode
    int NDesignVariable;
    int *DesignVariable;
    int *DesignVariableSide; // 0: X-Coord, 1: Y-Coord
    double *DesignVariableValue;
    double *DesignVariableMin;
    double *DesignVariableMax;
    double *dXdBeta;
    double *dYdBeta;
    double **dQdBeta;
    double **dRdX_dXdBeta;
    double **dIdQ;
    double **Lambda;
    double *dIdBeta;
    double *dIdQ_dQdBeta;
    double *dIdX_dXdBeta;
    // Matrix Computation CRS
    MC_CRS DesignBlockMatrix;
    // Cost Related Variables
    double Ref_CoeffLift;
    double Ref_CoeffDrag;
    double Ref_Lift;
    double Ref_Drag;
    double Ref_Moment;
    double Weight_Lift;
    double Weight_Drag;
    double Weight_Moment;
    double Weight_Cp;
    double Lift;
    double Drag;
    double Moment;
    double CoeffLift;
    double CoeffDrag;
    double I;
public:
    Euler2D_Design();
    virtual ~Euler2D_Design();
    virtual void Design_Prepare()=0;
    virtual void Design()=0;
    virtual void Design_Finalize()=0;
protected:
    // Get the Generic Design Inputs
    void Get_Design_Inputs(const char* FileName);
    
    // Initialize the Design
    void Initialize_Design(int InputBlockSize, int InputVectorSize);

    // Create and Initialize CRS Design Block Matrix
    void Create_CRS_DesignBlockMatrix(MC_CRS *Object);

    void Read_DesignFile(const char* FileName);
    void Write_DesignFile(const char* FileName);
    void Read_MeshSensitivity_dXdBeta(const char* FileName);
    void Write_MeshSensitivity_dXdBeta(const char* FileName);
    void Write_dQdBeta(const char* FileName);
    void Write_dRdX_dXdBeta(const char* FileName);
    // Read - Write Boundary Field
    void Write_Boundary_Field(const char* FileName, int Size,
                        int *NodeMap, double *CoordX,
                        double *CoordY, double *FieldData);
    void Read_Boundary_Field(const char* FileName, int *Size,
                        int **NodeMap, double **CoordX,
                        double **CoordY, double **FieldData);
    
    // Design Solver Cost
    virtual void Design_Cost()=0;
    // Design Solver Gradient: Direct/Forward Mode and Adjoint
    virtual void Design_Direct()=0;
    virtual void Design_Adjoint()=0;

    // Common Functions to Direct and Adjoint
    // --- START
    // Compute Cost: Coefficient of Lift and Drag Based
    virtual void Compute_Cost()=0;
    // Compute (dR/dX)*(dX/dBeta)
    virtual void Compute_dRdX_dXdBeta()=0;
    virtual void Compute_Boundary_dRdX_dXdBeta()=0;
    virtual void Verify_dRdX_dXdBeta()=0;
    virtual void Compute_dIdX_dXdBeta(int iDesignVariable)=0;
    virtual void Compute_LinearElasticSmooth_dXdBeta(int cNode)=0;
    // --- END

    // Compute Direct/Forward Mode Functions
    // Compute dQdBeta
    virtual void Compute_Direct_dQdBeta()=0;
    virtual void Verify_Direct_dQdBeta()=0;
    virtual void Compute_Direct_dIdQ_dQdBeta(int iDesignVariable)=0;
    virtual void Compute_Direct_dIdBeta(int iDesignVariable)=0;

    // Compute Adjoint Functions
    virtual void Compute_Adjoint_dIdQ()=0;
    virtual void Compute_Adjoint_Lambda()=0;
    virtual void Compute_Adjoint_dIdBeta(int iDesignVariable)=0;
private:
    void Init();
};

#endif	/* _EULER2D_DESIGN_H */
