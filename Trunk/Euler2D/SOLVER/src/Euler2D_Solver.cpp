/*******************************************************************************
 * File:        Euler2D_Solver.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <limits.h>

#include "Utils.h"
#include "Euler2D_Solver.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver::Euler2D_Solver() {
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver::Init() {
    SolverBlockSize              = 0;
    SolverVectorSize             = 0;
    DoSolverValidation           = 0;
    Order                        = 2;
    NIteration                   = 0;
    InnerNIteration              = 0;
    DoRestart                    = 0;
    Relaxation                   = 1.0;
    RMS[0]                       = DBL_MIN;
    RMS[1]                       = DBL_MIN;
    RMS[2]                       = DBL_MIN;
    RMS[3]                       = DBL_MIN;
    RMS_Res                      = DBL_MIN;
    Gamma                        = 1.4;
    Ref_Mach                     = 0.0;
    Ref_Alpha                    = 0.0;
    Ref_Pressure                 = 0.0;
    Ref_Length                   = 0.0;
    CFL_Ramp                     = 0;
    CFL_Min                      = 0.0;
    CFL_Max                      = 0.0;
    SolverBlockMatrix.VectorSize = 0;
    SolverBlockMatrix.BlockSize  = 0;
    SolverBlockMatrix.CRSSize    = 0;
    SolverBlockMatrix.A          = NULL;
    SolverBlockMatrix.B          = NULL;
    SolverBlockMatrix.IA         = NULL;
    SolverBlockMatrix.IAU        = NULL;
    SolverBlockMatrix.JA         = NULL;
    SolverBlockMatrix.X          = NULL;
    DeltaT                       = NULL;
    Qx                           = NULL;
    Qy                           = NULL;
    FDIteration                  = 0;
    FDNodeID                     = -1;
    FDPertQ                      = -1;
    FDEpsilon                    = 1.0E-5;
    FDDR_DQ[0]                   = 0.0;
    FDDR_DQ[1]                   = 0.0;
    FDDR_DQ[2]                   = 0.0;
    FDDR_DQ[3]                   = 0.0;
    FDNode                       = NULL;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver::~Euler2D_Solver() {
    int i, j;

    if (SolverBlockMatrix.A != NULL) {
        for (i = 0; i < SolverBlockMatrix.CRSSize; i++) {
            if (SolverBlockMatrix.A[i] != NULL) {
                for (j = 0; j < SolverBlockMatrix.BlockSize; j++) {
                    if (SolverBlockMatrix.A[i][j] != NULL)
                        free(SolverBlockMatrix.A[i][j]);
                }
                free(SolverBlockMatrix.A[i]);
            }
        }
        free(SolverBlockMatrix.A);
    }
    if (SolverBlockMatrix.B != NULL) {
        for (i = 0; i < SolverBlockMatrix.VectorSize ; i++) {
            if (SolverBlockMatrix.B[i] != NULL)
                free(SolverBlockMatrix.B[i]);
        }
        free(SolverBlockMatrix.B);
    }
    if (SolverBlockMatrix.X != NULL) {
        for (i = 0; i < SolverBlockMatrix.VectorSize ; i++) {
            if (SolverBlockMatrix.X[i] != NULL)
                free(SolverBlockMatrix.X[i]);
        }
        free(SolverBlockMatrix.X);
    }
    if (SolverBlockMatrix.IA != NULL)
        delete[] SolverBlockMatrix.IA;
    if (SolverBlockMatrix.IAU != NULL)
        delete[] SolverBlockMatrix.IAU;
    if (SolverBlockMatrix.JA != NULL)
        delete[] SolverBlockMatrix.JA;

    if (DeltaT != NULL)
        delete[] DeltaT;

    if (Qx != NULL) {
        for (i = 0; i < SolverVectorSize; i++) {
            if (Qx[i] != NULL)
                delete[] Qx[i];
        }
        delete[] Qx;
    }

    if (Qy != NULL) {
        for (i = 0; i < SolverVectorSize; i++) {
            if (Qy[i] != NULL)
                delete[] Qy[i];
        }
        delete[] Qy;
    }

    if (FDNode != NULL)
        delete[] FDNode;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver::Get_Solver_Inputs(const char* FileName) {
    FILE *fp;
    int iret;
    char dumstring[257];
    char dump[257];
    char *cdum;

    // Open Input File
    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Get_Solver_Inputs: Unable to Open Input File %s", FileName);

    // Input Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);
    WKAMeshFileName.append(dump);

    // Input Restart File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);
    InputRestartFileName.append(dump);

    // Output Restart File Name
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);
    OutputRestartFileName.append(dump);

    // SLK Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);
    SLKMeshFileName.append(dump);

    // VTK Based Solution File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%s\n", dump);
    VTKSolutionFileName.append(dump);

    // Get the reference Conditions
    // Gamma
    Gamma  = 1.4;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Gamma);
    Ref_Pressure    = 1.0/Gamma;
    Ref_Length      = 1.0;

    // Mach Number
    // Input Ref Mach Number
    Ref_Mach        = 0.5;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Ref_Mach);

    // Get Angle of Attack
    // Input Angle of Attack (deg)
    Ref_Alpha       = 0.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Ref_Alpha);
    Ref_Alpha       = Ref_Alpha * M_PI / 180.0;

    // Get the order of solver
    // Solver Order (1/2)
    Order = 2;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &Order);

    if ((Order < 1) || (Order > 2))
        error("Get_Reference_Conditions: %s\n", "Invalid Solver Order");

    // Get CFL Numbers
    // Input CFL Min
    CFL_Min         = 1.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &CFL_Min);
    // Input CFL Max
    CFL_Max         = 20.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &CFL_Max);
    // Input CFL Ramp
    CFL_Ramp        = 100;
    iret = fscanf(fp, "\n");
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &CFL_Ramp);

    // Get the Inner and Outer loop Conditions
    // Input Outer Iterations (0 = Convergence)
    NIteration      = 1000;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &NIteration);
    if (NIteration == 0)
        NIteration = INT_MAX/10;
    // Input Linear Solver Iterations (0 = Convergence)
    InnerNIteration = 20;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &InnerNIteration);
    // Input Linear Solver Relaxation Inner (0.5 < r < 1.0)
    Relaxation      = 1.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%lf\n", &Relaxation);

    // Restart File
    // Restart Capability (0=No, 1=Yes, 2=Create)
    DoRestart = 0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &DoRestart);

    // Finite Difference Jacobian
    FDNodeID = 0;
    DoSolverValidation = 0;
    // Finite Difference Jacobian (0=No, 1=Yes)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 256, fp);
    iret = fscanf(fp, "%d\n", &DoSolverValidation);
    if (DoSolverValidation == 1) {
        // Finite Difference Node ID (-1 = All Nodes)
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 256, fp);
        iret = fscanf(fp, "%d\n", &FDNodeID);
        // Finite Difference Perturbation Q
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 256, fp);
        iret = fscanf(fp, "%d\n", &FDPertQ);
        // Finite Difference Epsilon
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 256, fp);
        iret = fscanf(fp, "%lf\n", &FDEpsilon);
        // Finite Difference Check Iteration
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 256, fp);
        iret = fscanf(fp, "%d\n", &FDIteration);
    } else {
        FDNodeID = -2;
    }

#ifdef VERBOSE
    printf("Input Mesh File     = %s\n",  WKAMeshFileName.c_str());
    printf("Input Restart File  = %s\n",  InputRestartFileName.c_str());
    printf("Output Restart File = %s\n",  OutputRestartFileName.c_str());
    printf("Output Mesh File    = %s\n",  SLKMeshFileName.c_str());
    printf("Solution File       = %s\n",  VTKSolutionFileName.c_str());
    printf("Gamma               = %lf\n", Gamma);
    printf("Reference Mach      = %lf\n", Ref_Mach);
    printf("Reference Alpha     = %lf\n", (180.0/M_PI)*Ref_Alpha);
    printf("Solver Order        = %d\n",  Order);
    printf("CFL Minimum         = %lf\n", CFL_Min);
    printf("CFL Maximum         = %lf\n", CFL_Max);
    printf("CFL Ramp            = %d\n",  CFL_Ramp);
    printf("Solver Iteration    = %d\n",  NIteration);
    printf("LS Solver Iteration = %d\n",  InnerNIteration);
    printf("LS Relaxation       = %lf\n", Relaxation);
    printf("Restart             = %d\n",  DoRestart);
    printf("FD Node ID          = %d\n",  FDNodeID);
    printf("FD Pertuabation Q   = %d\n",  FDPertQ);
    printf("FD Epsilon          = %lf\n", FDEpsilon);
    printf("FD Check Iteration  = %d\n",  FDIteration);
    printf("-----------------------------------------------------------------------------\n");
#endif
    
    // Close Input File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver::Initialize_Solver(int InputBlockSize, int InputVectorSize) {
    int i, iNode;

    // Set the Dimensions
    SolverBlockSize  = InputBlockSize;
    SolverVectorSize = InputVectorSize;
    
    // Allocate the Memory to store Local Time Stepping
    DeltaT = new double[SolverVectorSize];
#ifdef DEBUG
    if (DeltaT == NULL)
        error("Initialize_Solver: %s\n", "Error Allocating Memory 1");
#endif

    if (Order == 2) {
        // Allocate Memory to store gradient of conserved variables
        Qx = new double*[SolverVectorSize];
        Qy = new double*[SolverVectorSize];
#ifdef DEBUG
        if ((Qx == NULL) || (Qy == NULL))
            error("Initialize_Solver: %s\n", "Error Allocating Memory 2");
#endif
        for (iNode = 0; iNode < SolverVectorSize; iNode++) {
            Qx[iNode] = NULL;
            Qy[iNode] = NULL;
            Qx[iNode] = new double[SolverBlockSize];
            Qy[iNode] = new double[SolverBlockSize];
#ifdef DEBUG
            if ((Qx[iNode] == NULL) || (Qy[iNode] == NULL))
                error("Initialize_Solver: %s\n", "Error Allocating Memory 3");
#endif
            // Initialize
            for (i = 0; i < SolverBlockSize; i++) {
                Qx[iNode][i] = 0.0;
                Qy[iNode][i] = 0.0;
            }
        }
    }

    // Allocate the Finite Difference Jacobian
    FDNode = (FD_NODE *) malloc(sizeof(FD_NODE)*SolverVectorSize);
#ifdef DEBUG
    if (FDNode == NULL)
        error("Initialize_Solver: %s\n", "Error Allocating Memory 4");
#endif
    // Initialize the Finite Difference Jacobian Variables
    for (iNode = 0; iNode < SolverVectorSize; iNode++) {
        for (i = 0; i < SolverBlockSize; i++) {
            FDNode[iNode].Resi_Old[i] = 0.0;
            FDNode[iNode].Resi_New[i] = 0.0;
        }
    }
}

