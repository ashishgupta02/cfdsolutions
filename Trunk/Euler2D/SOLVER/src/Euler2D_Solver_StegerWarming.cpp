/* 
 * File:   Euler2D_Solver_StegerWarming.cpp
 * Author: Ashish Gupta
 * 
 * Created on March 21, 2010, 11:10 PM
 */

#include <stdlib.h>

#include "Utils.h"
#include "Euler2D_Solver_StegerWarming.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_StegerWarming::Euler2D_Solver_StegerWarming() {
    printf("============================================================\n");
    printf("      Euler2D : Steger Warming Flux Vector Splitting        \n");
    printf("============================================================\n");

    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Init() {
    Gamma        = 1.4;
    Ref_Mach     = 0.0;
    Ref_Alpha    = 0.0;
    Ref_Pressure = 0.0;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_StegerWarming::~Euler2D_Solver_StegerWarming() {

}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_StegerWarming::Solve() {
    // Get the reference Conditions
    Ref_Mach     = 0.5;
    Ref_Alpha    = 0.0 * M_PI / 180.0;
    Ref_Pressure = 1.0/Gamma;
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

