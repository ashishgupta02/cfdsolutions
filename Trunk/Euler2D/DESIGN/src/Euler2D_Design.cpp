/*******************************************************************************
 * File:        Euler2D_Design.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits.h>

#include "Utils.h"
#include "List.h"
#include "Euler2D_Design.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Design::Euler2D_Design() {
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Init() {
    DesignAction                 = 0;
    BlockSize                    = 0;
    VectorSize                   = 0;
    DoValidation                 = 0;
    CostType                     = 0;
    IsAdjoint                    = 0;
    NDesignVariable              = 0;
    DesignVariable               = NULL;
    DesignVariableSide           = NULL;
    DesignVariableValue          = NULL;
    DesignVariableMin            = NULL;
    DesignVariableMax            = NULL;
    dXdBeta                      = NULL;
    dYdBeta                      = NULL;
    dQdBeta                      = NULL;
    dRdX_dXdBeta                 = NULL;
    dIdQ                         = NULL;
    Lambda                       = NULL;
    dIdBeta                      = NULL;
    dIdQ_dQdBeta                 = NULL;
    dIdX_dXdBeta                 = NULL;
    DesignBlockMatrix.VectorSize = 0;
    DesignBlockMatrix.BlockSize  = 0;
    DesignBlockMatrix.CRSSize    = 0;
    DesignBlockMatrix.A          = NULL;
    DesignBlockMatrix.B          = NULL;
    DesignBlockMatrix.IA         = NULL;
    DesignBlockMatrix.IAU        = NULL;
    DesignBlockMatrix.JA         = NULL;
    DesignBlockMatrix.X          = NULL;
    Ref_CoeffLift                = 0.0;
    Ref_CoeffDrag                = 0.0;
    Ref_Lift                     = 0.0;
    Ref_Drag                     = 0.0;
    Ref_Moment                   = 0.0;
    Weight_Lift                  = 0.0;
    Weight_Drag                  = 0.0;
    Weight_Moment                = 0.0;
    Weight_Cp                    = 0.0;
    Lift                         = 0.0;
    Drag                         = 0.0;
    Moment                       = 0.0;
    CoeffLift                    = 0.0;
    CoeffDrag                    = 0.0;
    I                            = 0.0;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Design::~Euler2D_Design() {
    int i, j;

    if (DesignBlockMatrix.A != NULL) {
        for (i = 0; i < DesignBlockMatrix.CRSSize; i++) {
            if (DesignBlockMatrix.A[i] != NULL) {
                for (j = 0; j < DesignBlockMatrix.BlockSize; j++) {
                    if (DesignBlockMatrix.A[i][j] != NULL)
                        free(DesignBlockMatrix.A[i][j]);
                }
                free(DesignBlockMatrix.A[i]);
            }
        }
        free(DesignBlockMatrix.A);
    }
    if (DesignBlockMatrix.B != NULL) {
        for (i = 0; i < DesignBlockMatrix.VectorSize ; i++) {
            if (DesignBlockMatrix.B[i] != NULL)
                free(DesignBlockMatrix.B[i]);
        }
        free(DesignBlockMatrix.B);
    }
    if (DesignBlockMatrix.X != NULL) {
        for (i = 0; i < DesignBlockMatrix.VectorSize ; i++) {
            if (DesignBlockMatrix.X[i] != NULL)
                free(DesignBlockMatrix.X[i]);
        }
        free(DesignBlockMatrix.X);
    }
    
    // Reset the Shared Memory Block with Solver Object
    DesignBlockMatrix.IA  = NULL;
    DesignBlockMatrix.IAU = NULL;
    DesignBlockMatrix.JA  = NULL;
    
    // Free the memory for design variables
    if (dXdBeta != NULL)
        delete [] dXdBeta;
    if (dYdBeta != NULL)
        delete [] dYdBeta;

    if (dRdX_dXdBeta != NULL) {
        if (VectorSize > 0) {
            for (i = 0; i < VectorSize; i++) {
                if (dRdX_dXdBeta[i] != NULL)
                    delete[] dRdX_dXdBeta[i];
            }
            delete [] dRdX_dXdBeta;
        }
    }

    if (dQdBeta != NULL) {
        if (VectorSize > 0) {
            for (i = 0; i < VectorSize; i++) {
                if (dQdBeta[i] != NULL)
                    delete[] dQdBeta[i];
            }
            delete [] dQdBeta;
        }
    }

    if (dIdQ != NULL) {
        if (VectorSize > 0) {
            for (i = 0; i < VectorSize; i++) {
                if (dIdQ[i] != NULL)
                    delete[] dIdQ[i];
            }
            delete [] dIdQ;
        }
    }

    if (Lambda != NULL) {
        if (VectorSize > 0) {
            for (i = 0; i < VectorSize; i++) {
                if (Lambda[i] != NULL)
                    delete[] Lambda[i];
            }
            delete [] Lambda;
        }
    }

    if (dIdBeta != NULL)
        delete [] dIdBeta;
    if (dIdQ_dQdBeta != NULL)
        delete [] dIdQ_dQdBeta;
    if (dIdX_dXdBeta != NULL)
        delete [] dIdX_dXdBeta;
    if (DesignVariable != NULL)
        delete[] DesignVariable;
    if (DesignVariableValue != NULL)
        delete[] DesignVariableValue;
    if (DesignVariableSide != NULL)
        delete[] DesignVariableSide;
    if (DesignVariableMin != NULL)
        delete[] DesignVariableMin;
    if (DesignVariableMax != NULL)
        delete[] DesignVariableMax;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Get_Design_Inputs(const char* FileName) {
    FILE *fp;
    int iret, itmp;
    double dtmp;
    char dumstring[257];
    char dump[257];
    char *cdum;

    // Get the Generic Design Inputs
    // Open Input File
    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Get_Design_Inputs: Unable to Open Input File %s", FileName);

    // START -- Solver Input
    // Input Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);

    // Input Restart File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);

    // Output Restart File Name
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);

    // SLK Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);

    // VTK Based Solution File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);

    // Get the reference Conditions
    // Gamma
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);

    // Mach Number
    // Input Ref Mach Number
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);

    // Get Angle of Attack
    // Input Angle of Attack (deg)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);

    // Get the order of solver
    // Solver Order (1/2)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);

    // Get CFL Numbers
    // Input CFL Min
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);
    // Input CFL Max
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);
    // Input CFL Ramp
    iret = fscanf(fp, "\n");
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);

    // Get the Inner and Outer loop Conditions
    // Input Outer Iterations (0 = Convergence)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);
    // Input Linear Solver Iterations (0 = Convergence)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);
    // Input Linear Solver Relaxation Inner (0.5 < r < 1.0)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);

    // Restart File
    // Restart Capability (0=No, 1=Yes, 2=Create)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);

    // Finite Difference Jacobian
    // Finite Difference Jacobian (0=No, 1=Yes)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &DoValidation);
    // Finite Difference Node ID (-1 = All Nodes)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);
    // Finite Difference Perturbation Q
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);
    // Finite Difference Epsilon
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &dtmp);
    // Finite Difference Check Iteration
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &itmp);
    // END -- Solver Input

    // START -- Design Input
    // Design Compute Cost or Gradient (0: Cost, 1: Gradient)
    DesignAction = 0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &DesignAction);
    // Design Gradient Type (0: Direct, 1: Adjoint)
    IsAdjoint = 0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &IsAdjoint);
    // Get the Cost Type
    // 0: Coeff Lift, 1: Coeff Drag, 2: Coeff(Lift + Drag) 3: Pressure
    CostType = 0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &CostType);
    // Get Reference Coefficient of Lift
    Ref_CoeffLift = 0.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Ref_CoeffLift);
    // Get Reference Coefficient of Drag
    Ref_CoeffDrag = 0.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Ref_CoeffDrag);
    // END -- Design Input

#ifdef VERBOSE
    printf("Design Action       = %d\n",  DesignAction);
    printf("Gradient Type       = %d\n",  IsAdjoint);
    printf("Cost Type           = %d\n",  CostType);
    printf("Ref_CoeffLift       = %14.10e\n",  Ref_CoeffLift);
    printf("Ref_CoeffDrag       = %14.10e\n",  Ref_CoeffDrag);
    printf("-----------------------------------------------------------------------------\n");
#endif
    
    // Close Input File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Initialize_Design(int InputBlockSize, int InputVectorSize) {
    int i, j;

    BlockSize   = InputBlockSize;
    VectorSize  = InputVectorSize;
    // Allocate the memory for design variables
    dXdBeta      = new double[VectorSize];
    dYdBeta      = new double[VectorSize];
    dRdX_dXdBeta = new double*[VectorSize];
    dQdBeta      = new double*[VectorSize];
    dIdQ         = new double*[VectorSize];
    Lambda       = new double*[VectorSize];
    for (i = 0; i < VectorSize; i++) {
        dXdBeta[i]      = 0.0;
        dYdBeta[i]      = 0.0;
        dRdX_dXdBeta[i] = new double[BlockSize];
        dQdBeta[i]      = new double[BlockSize];
        dIdQ[i]         = new double[BlockSize];
        Lambda[i]       = new double[BlockSize];
        for (j = 0; j < BlockSize; j++) {
            dRdX_dXdBeta[i][j] = 0.0;
            dQdBeta[i][j]      = 0.0;
            dIdQ[i][j]         = 0.0;
            Lambda[i][j]       = 0.0;
        }
    }

//    Write_MeshSensitivity_dXdBeta("dXdBeta.q");

    // Read the dX/dBeta and Get No of Design Variables
//    Read_MeshSensitivity_dXdBeta("dXdBeta.q");

    // Read the Design File
    Read_DesignFile("Design.data");
    
    // Allocate Mememory to store design Sensitivities
    if (NDesignVariable <= 0)
        error("Initialize_Design: No of Design Variable %d", NDesignVariable);
    dIdQ_dQdBeta = new double[NDesignVariable];
    dIdX_dXdBeta = new double[NDesignVariable];
    for (i = 0; i < NDesignVariable; i++) {
        dIdQ_dQdBeta[i] = 0.0;
        dIdX_dXdBeta[i] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Create_CRS_DesignBlockMatrix(MC_CRS *Object) {
    int i, j, k;

    // Get the value from Solver Block Matrix
    DesignBlockMatrix.VectorSize = Object->VectorSize;
    DesignBlockMatrix.BlockSize  = Object->BlockSize;
    DesignBlockMatrix.CRSSize    = Object->CRSSize;
    // Share the memory with Solver Block Matrix
    DesignBlockMatrix.IA         = Object->IA;
    DesignBlockMatrix.IAU        = Object->IAU;
    DesignBlockMatrix.JA         = Object->JA;

    // Allocate Memory for Design Specific Computations
    // Allocate Memory for CRS Matrix
    DesignBlockMatrix.A = NULL;
    DesignBlockMatrix.A = (double ***) malloc (DesignBlockMatrix.CRSSize*sizeof(double**));
#ifdef DEBUG
    if (DesignBlockMatrix.A == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < DesignBlockMatrix.CRSSize; i++) {
        DesignBlockMatrix.A[i] = NULL;
        DesignBlockMatrix.A[i] = (double **) malloc (DesignBlockMatrix.BlockSize*sizeof(double*));
#ifdef DEBUG
        if (DesignBlockMatrix.A[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 2");
#endif
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++) {
            DesignBlockMatrix.A[i][j] = NULL;
            DesignBlockMatrix.A[i][j] = (double *) malloc (DesignBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
            if (DesignBlockMatrix.A[i] == NULL)
                error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 3");
#endif
            for (k = 0; k < DesignBlockMatrix.BlockSize; k++)
                DesignBlockMatrix.A[i][j][k] = 0.0;
        }
    }

    // Allocate Memory of RHS
    DesignBlockMatrix.B = NULL;
    DesignBlockMatrix.B = (double **) malloc (DesignBlockMatrix.VectorSize*sizeof(double*));
#ifdef DEBUG
    if (DesignBlockMatrix.B == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 4");
#endif
    for (i = 0; i < DesignBlockMatrix.VectorSize; i++) {
        DesignBlockMatrix.B[i] = NULL;
        DesignBlockMatrix.B[i] = (double *) malloc (DesignBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
        if (DesignBlockMatrix.B[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 5");
#endif
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++)
            DesignBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    DesignBlockMatrix.X = NULL;
    DesignBlockMatrix.X = (double **) malloc (DesignBlockMatrix.VectorSize*sizeof(double*));
#ifdef DEBUG
    if (DesignBlockMatrix.X == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 6");
#endif
    for (i = 0; i < DesignBlockMatrix.VectorSize; i++) {
        DesignBlockMatrix.X[i] = NULL;
        DesignBlockMatrix.X[i] = (double *) malloc (DesignBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
        if (DesignBlockMatrix.X[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 7");
#endif
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++)
            DesignBlockMatrix.X[i][j] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Read_DesignFile(const char* FileName) {
    char   dumstring[257];
    int    *side;
    int    *point;
    double *value;
    double *minv;
    double *maxv;
    double *grad;
    FILE   *fp;
    char   *cdum;
    int    iret, idum;
    List   unode;
    
    // Open Mesh File
    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Read_DesignFile: Unable to Open Design File %s", FileName);

    // Read Number of Design Variables
    iret = fscanf(fp, "%d", &NDesignVariable);

    // Allocate Memory
    side  = (int *)    calloc(NDesignVariable, sizeof(int));
    point = (int *)    calloc(NDesignVariable, sizeof(int));
    value = (double *) calloc(NDesignVariable, sizeof(double));
    minv  = (double *) calloc(NDesignVariable, sizeof(double));
    maxv  = (double *) calloc(NDesignVariable, sizeof(double));
    grad  = (double *) calloc(NDesignVariable, sizeof(double));

    // Read the Design Values and constraints
    for (int i = 0; i < NDesignVariable; i++)
        iret = fscanf(fp, "%d %s %d %d %lf %lf %lf", &idum, dumstring,
                &side[i], &point[i], &value[i], &minv[i], &maxv[i]);

    // Get the Cost
    iret = fscanf(fp, "%lf", &I);
    
    // Read the Gradients
    for (int i = 0; i < NDesignVariable; i++)
        iret = fscanf(fp, "%d %lf", &idum, &grad[i]);
    
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf %lf %lf %lf", &Weight_Lift, &Weight_Drag, &Weight_Moment, &Weight_Cp);

    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf %lf %lf", &Ref_Lift, &Ref_Drag, &Ref_Moment);

    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf %lf %lf", &Lift, &Drag, &Moment);

    fclose(fp);
    
    // Get the Nodes to be updated or Actual Design Variables
    for (int i = 0; i < NDesignVariable; i++) {
        if (minv[i] < maxv[i])
            unode.Check_List(i);
    }

    // Get the Real Design Variables
    NDesignVariable = unode.max;
    if (DesignVariable != NULL)
        delete[] DesignVariable;
    if (DesignVariableValue != NULL)
        delete[] DesignVariableValue;
    if (DesignVariableSide != NULL)
        delete[] DesignVariableSide;
    if (DesignVariableMin != NULL)
        delete[] DesignVariableMin;
    if (DesignVariableMax != NULL)
        delete[] DesignVariableMax;
    if (dIdBeta != NULL)
        delete[] dIdBeta;
    // Allocate the Memory
    DesignVariable      = new int[NDesignVariable];
    DesignVariableSide  = new int[NDesignVariable];
    DesignVariableValue = new double[NDesignVariable];
    DesignVariableMin   = new double[NDesignVariable];
    DesignVariableMax   = new double[NDesignVariable];
    dIdBeta             = new double[NDesignVariable];
    for (int i = 0; i < NDesignVariable; i++) {
        DesignVariable[i]      = point[unode.list[i]];
        DesignVariableValue[i] = value[unode.list[i]];
        DesignVariableSide[i]  = side[unode.list[i]];
        DesignVariableMin[i]   = minv[unode.list[i]];
        DesignVariableMax[i]   = maxv[unode.list[i]];
        dIdBeta[i]             = grad[unode.list[i]];
    }
    
    // Free Memory
    free(side);
    free(point);
    free(value);
    free(minv);
    free(maxv);
    free(grad);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Write_DesignFile(const char* FileName) {
    char   dumstring[257];
    FILE   *fp;
    
    // Open Mesh File
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_DesignFile: Unable to Write Design File %s", FileName);

    // Write Number of Design Variables
    fprintf(fp, " %6d\n", NDesignVariable);

    // Set the Dummy String
    strcpy(dumstring, "controlGridAirfoil.input");
    
    // Write the Design Values and constraints
    for (int i = 0; i < NDesignVariable; i++)
        fprintf(fp, " %6d %25s %6d %6d %19.10E %19.10E %19.10E\n", i, dumstring,
                DesignVariableSide[i], DesignVariable[i], DesignVariableValue[i],
                DesignVariableMin[i], DesignVariableMax[i]);

    // Write the Cost
    fprintf(fp, " %19.10E\n", I);

    // Write the Gradients
    for (int i = 0; i < NDesignVariable; i++)
        fprintf(fp, " %6d %19.10E\n", i, dIdBeta[i]);

    fprintf(fp, " Functions: Omega1(lift)  Omega2(drag)  Omega3(moment)  Omega4(Cp)\n");
    fprintf(fp, " %19.10E %19.10E %19.10E %19.10E\n", Weight_Lift, Weight_Drag, Weight_Moment, Weight_Cp);
    fprintf(fp, " Target values:  lift    drag    moment\n");
    fprintf(fp, " %19.10E %19.10E %19.10E\n", Ref_Lift, Ref_Drag, Ref_Moment);
    fprintf(fp, " Current values: lift    drag    moment\n");
    fprintf(fp, " %19.10E %19.10E %19.10E\n", Lift, Drag, Moment);
    
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Read_MeshSensitivity_dXdBeta(const char* FileName) {
    FILE *fp;
    int i, itmp, iret;

    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Read_dXdBeta: Unable to Open dXdBeta File %s", FileName);

    NDesignVariable = 0;
    for (i = 0; i < VectorSize; i++) {
        itmp = 0;
        iret = fscanf(fp, "%d %lf %lf", &itmp, &dXdBeta[i], &dYdBeta[i]);
        if (itmp !=0 )
            NDesignVariable++;
    }
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Write_MeshSensitivity_dXdBeta(const char* FileName) {
    FILE *fp;
    int i, itmp;
    
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dXdBeta: Unable to Open Mesh File %s", FileName);
    
    itmp = 0;
    for (i = 0; i < VectorSize; i++)
        fprintf(fp, "%22d %22.15e %22.15e\n", itmp, dXdBeta[i], dYdBeta[i]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Write_dQdBeta(const char* FileName) {
    FILE *fp;
    int i;
    
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dQdBeta: Unable to Open Mesh File %s", FileName);

    for (i = 0; i < VectorSize; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n",
                dQdBeta[i][0], dQdBeta[i][1], dQdBeta[i][2], dQdBeta[i][3]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Write_dRdX_dXdBeta(const char* FileName) {
    FILE *fp;
    int i;

    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dRdX_dXdBeta: Unable to Open Mesh File %s", FileName);

    for (i = 0; i < VectorSize; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n",
                dRdX_dXdBeta[i][0], dRdX_dXdBeta[i][1], dRdX_dXdBeta[i][2], dRdX_dXdBeta[i][3]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Write_Boundary_Field(const char* FileName, int Size,
        int* NodeMap, double *CoordX, double *CoordY, double* FieldData) {
    FILE *fp;
    int i;

    // Some Error Checks
    if ((Size <= 0)|| (NodeMap == NULL) || (CoordX == NULL)
            || (CoordY == NULL) || (FieldData == NULL))
        error("Write_Boundary_Field: Invalid Inputs");
    
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_Boundary_Field: Unable to Write File %s", FileName);

    // Number of Entries
    fprintf(fp,"%6d\n", Size);

    // NodeID and Field Data
    for (i = 0; i < Size; i++) {
        fprintf(fp, "%6d  %22.15e %22.15e %22.15e\n",
                NodeMap[i], CoordX[i], CoordY[i], FieldData[i]);
    }

    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design::Read_Boundary_Field(const char* FileName, int* Size,
        int** NodeMap, double **CoordX, double **CoordY, double** FieldData) {
    FILE   *fp;
    int    *Map;
    double *FData;
    double *XValue;
    double *YValue;
    int i, iret;

    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Read_Boundary_Field: Unable to Open File %s", FileName);

    *Size = 0;
    iret = fscanf(fp, "%d", Size);
    
    if (*Size <= 0) {
        warn("Read_Boundary_Field: No Data Field Found in %s", FileName);
        return;
    }

    // Allocate Some Memory
    Map    = (int *) malloc((*Size)*sizeof(int));
    FData  = (double *) malloc((*Size)*sizeof(double));
    XValue = (double *) malloc((*Size)*sizeof(double));
    YValue = (double *) malloc((*Size)*sizeof(double));

    // Now Get the NodeID and Field Data
    for (i = 0; i < *Size; i++) {
        iret = fscanf(fp, "%d %lf %lf %lf", &Map[i], &XValue[i], &YValue[i], &FData[i]);
    }

    *NodeMap   = Map;
    *FieldData = FData;
    *CoordX    = XValue;
    *CoordY    = YValue;

    Map    = NULL;
    FData  = NULL;
    XValue = NULL;
    YValue = NULL;
    fclose(fp);
}

