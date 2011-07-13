/*******************************************************************************
 * File:        Euler2D_Solver_StegerWarming.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <limits.h>

#include "Utils.h"
#include "Euler2D_Solver_StegerWarming.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_StegerWarming::Euler2D_Solver_StegerWarming() {
#ifdef VERBOSE
    printf("============================================================\n");
    printf("      Euler2D : Steger Warming Flux Vector Splitting        \n");
    printf("============================================================\n");
#endif
    
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Init() {
    // Data Initialization
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_StegerWarming::~Euler2D_Solver_StegerWarming() {
    // Free the Resource Used
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Solver_Prepare() {
    // Read Mesh File
    WKA_MeshReader(WKAMeshFileName.c_str());
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Solver_Finalize() {
    // Write SLK Mesh file
    SLK_MeshWriter(SLKMeshFileName.c_str());
    // Write VTK Solution File
    Write_VTK_Unstructured_File(VTKSolutionFileName.c_str());
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Get_Solver_Inputs_StegerWarming(const char* FileName) {
    // Get the Generic Solver Inputs
    Get_Solver_Inputs(FileName);
    // Now Get the Steger Warming Solver Inputs
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Initialize_Solver_StegerWarming() {
    // Initialize the Solver Data Field
    Initialize_Solver(4, mesh.nnodes);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Solve() {
    // Algorithm
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Initialize_Solution() {
    int i, iNode;
    double Q[4];

    // Calculte Q
    Q[0] = 1.0;
    Q[1] = Q[0]*Ref_Mach*cos(Ref_Alpha);
    Q[2] = Q[0]*Ref_Mach*sin(Ref_Alpha);
    Q[3] = (Ref_Pressure/(Gamma - 1.0)) + 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0]));

    // Initialize the Solution Field
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        node[iNode].area = 0.0;
        for (i = 0; i < 4; i++) {
            node[iNode].Q[i] = Q[i];
            node[iNode].Resi[i] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_Gauss_Gradient() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_Residual() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_Boundary_Residual() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_FD_Jacobian(int Iteration) {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_DeltaTime(int Iteration) {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Create_CRS_SolverBlockMatrix() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_CRS_SolverBlockMatrix(int AddTime) {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Compute_Boundary_CRS_SolverBlockMatrix() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Update_Solution() {
    // TODO
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Solver_StegerWarming::Compute_RMS() {
    // TODO
    return 0.0;
}

