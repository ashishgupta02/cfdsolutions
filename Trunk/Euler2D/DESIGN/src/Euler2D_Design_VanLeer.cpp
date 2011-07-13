/*******************************************************************************
 * File:        Euler2D_Design_VanLeer.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits.h>

#include "List.h"
#include "Utils.h"
#include "Mesh.h"
#include "Euler2D_Design.h"
#include "Euler2D_Design_VanLeer.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Design_VanLeer::Euler2D_Design_VanLeer() {
#ifdef VERBOSE
    printf("=============================================================================\n");
    printf("      Euler2D : Van Leer Computational Design                                \n");
    printf("=============================================================================\n");
#endif

    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Init() {
    // Data Initialization
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Design_VanLeer::~Euler2D_Design_VanLeer() {
    // Free the Resource Used
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Get_Design_Inputs_VanLeer(const char* FileName) {
    // Get the Solver Inputs First
    Get_Solver_Inputs_VanLeer(FileName);
    // Get the Generic Design Inputs
    Get_Design_Inputs(FileName);
    // Now Get the VanLeer Design Inputs
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design_Prepare() {
    // Prepare the Solver
    Solver_Prepare();
    // Now Prepare the Design System
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design() {

    // Create and Initialize Solver and Design Data
    Initialize_Design_VanLeer();
    
    switch (DesignAction) {
        // Compute Cost
        case 0:
            Design_Cost();
            break;
        // Compute Gradient
        case 1:
            // Forward Mode
            if (IsAdjoint == 0) {
                Design_Direct();
#ifdef VERBOSE
                printf("DesignV       dIdQ_dQdBeta       dIdX_dXdBeta            dIdBeta\n");
                for (int iDV = 0; iDV < NDesignVariable; iDV++)
                    printf("%7d %18.10E %18.10E %18.10E\n", iDV,
                            dIdQ_dQdBeta[iDV], dIdX_dXdBeta[iDV], dIdBeta[iDV]);
#endif
            // Adjoint Mode
            } else {
                Design_Adjoint();
#ifdef VERBOSE
                printf("DesignV  -(dRdX_dXdBeta)Lamda    dIdX_dXdBeta            dIdBeta\n");
                for (int iDV = 0; iDV < NDesignVariable; iDV++)
                    printf("%7d %18.10E %18.10E %18.10E\n", iDV,
                            (dIdBeta[iDV] - dIdX_dXdBeta[iDV]), dIdX_dXdBeta[iDV], dIdBeta[iDV]);
#endif
            }
            break;
        // Dummy Do Nothing
        default:
            break;
    }
    
#ifdef VERBOSE
    printf("Cl = %lf, Cd = %lf, L = %lf, D = %lf, I = %lf\n",
            CoeffLift, CoeffDrag, Lift, Drag, I);
#endif
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design_Finalize() {
    // Write the Design File
    Write_DesignFile("Design.data");
    // Write New Mesh File
    WKA_MeshWriter("New.mesh");
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Initialize_Design_VanLeer() {
    // START -- Solver Initialization
    // Initialize the Solver Data Field
    Initialize_Solver_VanLeer();
    // END   -- Solver Initialization

    // START -- Design Initialization
    // Initialize the Design Data Field
    Initialize_Design(4, mesh.nnodes);

    // Create CRS Design Matrix
    Create_CRS_DesignBlockMatrix(&SolverBlockMatrix);
    // END   -- Design Initialization
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Update_Mesh_LinearElasticSmooth() {
    int iNode, icount, iter;
    double relax;
    double *U;
    double *V;

#ifdef VERBOSE
    info("Linear Elastic Mesh Smooth Starts");
#endif

    // Initialize
    U = V = NULL;
    // Create the Purturbation Array
    U = new double[mesh.nnodes];
    V = new double[mesh.nnodes];

    icount = 0;
    // Right Now Only Purturbation in Y is done
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        if ((icount < NDesignVariable) && (DesignVariable[icount] == iNode)) {
            U[iNode] = 0.0;
            V[iNode] = DesignVariableValue[icount] - node[iNode].y;
            icount++;
        } else {
            U[iNode] = 0.0;
            V[iNode] = 0.0;
        }
    }

    relax = 0.5;
    iter  = 0;
    LESmooth.Mesh_Smoother_Prepare();
    LESmooth.Initialize_Mesh_Smoother("domain_baseline.mesh", &DesignBlockMatrix);
//    LESmooth.Initialize_Mesh_Smoother(mesh, cell, edge, node, boundaryEdge,
//            boundaryNode, BNTag, &DesignBlockMatrix);
    LESmooth.Mesh_Smoother(iter, relax, U, V);
    LESmooth.Mesh_Smoother_Finalize();

    // Apply the Purturbations
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        node[iNode].x += U[iNode];
        node[iNode].y += V[iNode];
    }

    // Update the Geometric Properties
    Compute_Geometric_Properties();
    
    // Delete Memory
    delete []U;
    delete []V;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_New_Solution() {
    int iter, AddTime;
    double max_rms, lrms;

#ifdef VERBOSE
    info("Solver Computation Starts");
    printf("-----------------------------------------------------------------------------\n");
    printf(" Iter        LRMS     RMS_RHO    RMS_RHOU    RMS_RHOV       RMS_E     RMS_RES\n");
    printf("-----------------------------------------------------------------------------\n");
#endif
    
    for (iter = 0; iter < NIteration; iter++) {
        // Compute Local Time Stepping 
        Compute_DeltaTime(iter);

        // Compute Residials
        Compute_Gauss_Gradient();
        Compute_Residual();
        Compute_Boundary_Residual();
  
        // Compute and Fill CRS Matrix
        AddTime = 1;
        Compute_CRS_SolverBlockMatrix(AddTime);
        Compute_Boundary_CRS_SolverBlockMatrix();

        // Solve for Solution
        lrms = MC_Iterative_Block_Jacobi_CRS(InnerNIteration, Relaxation, SolverBlockMatrix);
//        lrms = MC_Iterative_Block_LU_Jacobi_CRS(InnerNIteration, 0, SolverBlockMatrix);

        // Update the Solution
        Update_Solution();

        // Compute RMS
        max_rms = Compute_RMS();

#ifdef VERBOSE
        printf("%5d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n", iter+1, lrms, RMS[0], RMS[1], RMS[2], RMS[3], RMS_Res);
#endif
        
        if ((max_rms < (DBL_EPSILON*10.0))|| ((iter+1) == NIteration))
            iter = NIteration + 1;

        if (isnan(max_rms)) {
            info("Solve: Outer Iteration: NAN Encountered ! - Abort");
            iter = NIteration + 1;
        }

        if (isnan(lrms)) {
            info("Solve: Liner Solver Iteration: NAN Encountered ! - Abort");
            iter = NIteration + 1;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design_Cost() {
    // Snap the mesh to new location and Smooth
    Update_Mesh_LinearElasticSmooth();

    // Compute New Solution
    Compute_New_Solution();
    
    // Now Compute Cost
    Compute_Cost();
    
#ifndef DEBUG
    Write_VTK_Unstructured_File("Debug.vtu");
#endif
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design_Direct() {
    int iDV, AddTime;
    
    // Compute Cost
    Compute_Cost();

    for (iDV = 0; iDV < NDesignVariable; iDV++) {
        // Perturb the Node
        // --- TODO

        // Compute the Expensive Geometric Properties
        Compute_Geometric_Properties();

        // Compute Local Time Stepping
        CFL_Ramp += 1;
        Compute_DeltaTime(CFL_Ramp);

        // Compute Residuals
        Compute_Gauss_Gradient();
        Compute_Residual();
        Compute_Boundary_Residual();

        // Compute and Fill CRS Solver Matrix: dR/dQ
        AddTime = 0;
        Compute_CRS_SolverBlockMatrix(AddTime);
        Compute_Boundary_CRS_SolverBlockMatrix();

        // Compute (dR/dX)*(dX/dBeta)
        Compute_dRdX_dXdBeta();
        Compute_Boundary_dRdX_dXdBeta();
        // Perform Validation
        if (DoValidation) {
            Verify_dRdX_dXdBeta();
            Write_dRdX_dXdBeta("dRdX_dXdBeta.q");
        }
        // Compute dQ/dBeta
        Compute_Direct_dQdBeta();
        // Perform Validation
        if (DoValidation) {
            Verify_Direct_dQdBeta();
            Write_dQdBeta("dQdBeta.q");
        }
        // Compute Cost Solution Sensitivity
        Compute_Direct_dIdQ_dQdBeta(iDV);
        // Compute Cost Mesh Sensitivity
        Compute_dIdX_dXdBeta(iDV);
        // Compute Forward Mode Gradient
        Compute_Direct_dIdBeta(iDV);

        // UnPerturb the Node
        // --- TODO
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Design_Adjoint() {
    int iDV, AddTime;

    // Snap the mesh to new location and Smooth
    Update_Mesh_LinearElasticSmooth();

    // Compute New Solution
    Compute_New_Solution();
    
    // Compute Cost
    Compute_Cost();

#ifndef DEBUG
    Write_VTK_Unstructured_File("Debug.vtu");
#endif
    
    // Compute Local Time Stepping
    CFL_Ramp += 1;
    Compute_DeltaTime(CFL_Ramp);

    // Compute Residuals
    Compute_Gauss_Gradient();
    Compute_Residual();
    Compute_Boundary_Residual();

    // Compute and Fill CRS Solver Matrix: dR/dQ
    AddTime = 0;
    Compute_CRS_SolverBlockMatrix(AddTime);
    Compute_Boundary_CRS_SolverBlockMatrix();

    Compute_Adjoint_dIdQ();
    Compute_Adjoint_Lambda();

    // Loop Over all Design Variables
    for (iDV = 0; iDV < NDesignVariable; iDV++) {
        // Compute the Linear Elastic Smooth dXdBeta
        Compute_LinearElasticSmooth_dXdBeta(DesignVariable[iDV]);
        
        // Compute (dR/dX)*(dX/dBeta)
        Compute_dRdX_dXdBeta();
        Compute_Boundary_dRdX_dXdBeta();
        
        // Perform Validation
        if (DoValidation) {
            FDNodeID = DesignVariable[iDV];
            Verify_dRdX_dXdBeta();
            Write_dRdX_dXdBeta("dRdX_dXdBeta.q");
        }
        // Compute Cost Mesh Sensitivity
        Compute_dIdX_dXdBeta(iDV);
        // Compute Adjoint Gradient
        Compute_Adjoint_dIdBeta(iDV);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_dRdX_dXdBeta() {
    int i, iNode, iCell;
    int n1, n2, n3;
    double FplusdBeta[4], FminusdBeta[4];
    double x[3],   xc,   x12,   x23,   x31,   nxc12,   nxc23,   nxc31,   unxc12,   unxc23,   unxc31;
    double x_b[3], xc_b, x12_b, x23_b, x31_b, nxc12_b, nxc23_b, nxc31_b, unxc12_b, unxc23_b, unxc31_b;
    double y[3],   yc,   y12,   y23,   y31,   nyc12,   nyc23,   nyc31,   unyc12,   unyc23,   unyc31;
    double y_b[3], yc_b, y12_b, y23_b, y31_b, nyc12_b, nyc23_b, nyc31_b, unyc12_b, unyc23_b, unyc31_b;
    double magc12,   magc23,   magc31;
    double magc12_b, magc23_b, magc31_b;

    // Initialize
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++)
            dRdX_dXdBeta[iNode][i] = 0.0;
    }

    // Loop over cells
    for (iCell = 0; iCell < mesh.inside; iCell++) {
        n1 = cell[iCell].node1;
        n2 = cell[iCell].node2;
        n3 = cell[iCell].node3;

        // Get the coordinates
        x[0] = node[n1].x;
        x[1] = node[n2].x;
        x[2] = node[n3].x;
        y[0] = node[n1].y;
        y[1] = node[n2].y;
        y[2] = node[n3].y;

        // Get the Coordinate Derivates
        x_b[0] = dXdBeta[n1];
        x_b[1] = dXdBeta[n2];
        x_b[2] = dXdBeta[n3];
        y_b[0] = dYdBeta[n1];
        y_b[1] = dYdBeta[n2];
        y_b[2] = dYdBeta[n3];

        // Compute the Centriod of Cell
        xc = (x[0] + x[1] + x[2]) / 3.0;
        yc = (y[0] + y[1] + y[2]) / 3.0;

        // Compute the Centriod of Cell Derivates
        xc_b = (x_b[0] + x_b[1] + x_b[2]) / 3.0;
        yc_b = (y_b[0] + y_b[1] + y_b[2]) / 3.0;

        // Compute the Edge Centers
        x12 = (x[0] + x[1]) / 2.0;
        y12 = (y[0] + y[1]) / 2.0;
        x23 = (x[1] + x[2]) / 2.0;
        y23 = (y[1] + y[2]) / 2.0;
        x31 = (x[2] + x[0]) / 2.0;
        y31 = (y[2] + y[0]) / 2.0;

        // Compute the Edge Centers Derivates
        x12_b = (x_b[0] + x_b[1]) / 2.0;
        y12_b = (y_b[0] + y_b[1]) / 2.0;
        x23_b = (x_b[1] + x_b[2]) / 2.0;
        y23_b = (y_b[1] + y_b[2]) / 2.0;
        x31_b = (x_b[2] + x_b[0]) / 2.0;
        y31_b = (y_b[2] + y_b[0]) / 2.0;

        // Compute the Normal, Unit Normal to Median Dual Edge and its Magnitude
        nxc12  = +(yc - y12);
        nyc12  = -(xc - x12);
        magc12 = sqrt(nxc12*nxc12 + nyc12*nyc12);
        unxc12 = nxc12/magc12;
        unyc12 = nyc12/magc12;

        nxc23  = +(yc - y23);
        nyc23  = -(xc - x23);
        magc23 = sqrt(nxc23*nxc23 + nyc23*nyc23);
        unxc23 = nxc23/magc23;
        unyc23 = nyc23/magc23;

        nxc31  = +(yc - y31);
        nyc31  = -(xc - x31);
        magc31 = sqrt(nxc31*nxc31 + nyc31*nyc31);
        unxc31 = nxc31/magc31;
        unyc31 = nyc31/magc31;

        // Compute the Normal, Unit Normal to Median Dual Edge and its Magnitude Derivates
        nxc12_b  = +(yc_b - y12_b);
        nyc12_b  = -(xc_b - x12_b);
        magc12_b = (nxc12*nxc12_b + nyc12*nyc12_b)/magc12;
        unxc12_b = (nxc12_b*magc12 - magc12_b*nxc12)/(magc12*magc12);
        unyc12_b = (nyc12_b*magc12 - magc12_b*nyc12)/(magc12*magc12);

        nxc23_b  = +(yc_b - y23_b);
        nyc23_b  = -(xc_b - x23_b);
        magc23_b = (nxc23*nxc23_b + nyc23*nyc23_b)/magc23;
        unxc23_b = (nxc23_b*magc23 - magc23_b*nxc23)/(magc23*magc23);
        unyc23_b = (nyc23_b*magc23 - magc23_b*nyc23)/(magc23*magc23);

        nxc31_b  = +(yc_b - y31_b);
        nyc31_b  = -(xc_b - x31_b);
        magc31_b = (nxc31*nxc31_b + nyc31*nyc31_b)/magc31;
        unxc31_b = (nxc31_b*magc31 - magc31_b*nxc31)/(magc31*magc31);
        unyc31_b = (nyc31_b*magc31 - magc31_b*nyc31)/(magc31*magc31);

        // ********************  Edge 1-2  *******************************
        Compute_dFplusdX_dXdBeta(node[n1].Q, unxc12, unyc12, magc12, unxc12_b, unyc12_b, magc12_b, FplusdBeta);
        Compute_dFminusdX_dXdBeta(node[n2].Q, unxc12, unyc12, magc12, unxc12_b, unyc12_b, magc12_b, FminusdBeta);
        for (i = 0; i < 4; i++) {
            dRdX_dXdBeta[n1][i] += (FplusdBeta[i] + FminusdBeta[i]);
            dRdX_dXdBeta[n2][i] -= (FplusdBeta[i] + FminusdBeta[i]);
        }

        // ********************  Edge 2-3  *******************************
        Compute_dFplusdX_dXdBeta(node[n2].Q, unxc23, unyc23, magc23, unxc23_b, unyc23_b, magc23_b, FplusdBeta);
        Compute_dFminusdX_dXdBeta(node[n3].Q, unxc23, unyc23, magc23, unxc23_b, unyc23_b, magc23_b, FminusdBeta);
        for (i = 0; i < 4; i++) {
            dRdX_dXdBeta[n2][i] += (FplusdBeta[i] + FminusdBeta[i]);
            dRdX_dXdBeta[n3][i] -= (FplusdBeta[i] + FminusdBeta[i]);
        }

        // ********************  Edge 3-1  *******************************
        Compute_dFplusdX_dXdBeta(node[n3].Q, unxc31, unyc31, magc31, unxc31_b, unyc31_b, magc31_b, FplusdBeta);
        Compute_dFminusdX_dXdBeta(node[n1].Q, unxc31, unyc31, magc31, unxc31_b, unyc31_b, magc31_b, FminusdBeta);
        for (i = 0; i < 4; i++) {
            dRdX_dXdBeta[n3][i] += (FplusdBeta[i] + FminusdBeta[i]);
            dRdX_dXdBeta[n1][i] -= (FplusdBeta[i] + FminusdBeta[i]);
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Boundary_dRdX_dXdBeta() {
    int i, n1, n2, iBEdge, iedge;
    double Pressure;
    double x_b[2], Nx_b;
    double y_b[2], Ny_b;

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[iBEdge].bcType < 1000) || (boundaryEdge[iBEdge].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Get the Coordinate Derivates
        x_b[0] = dXdBeta[n1];
        x_b[1] = dXdBeta[n2];
        y_b[0] = dYdBeta[n1];
        y_b[1] = dYdBeta[n2];

        // Compute the Normal Derivatives
        Nx_b  = +0.5*(y_b[1] - y_b[0]);
        Ny_b  = -0.5*(x_b[1] - x_b[0]);

        // dResi/dBeta at Node1
        Pressure = (Gamma - 1.0)*(node[n1].Q[3] - 0.5*node[n1].Q[0]*((node[n1].Q[1]/node[n1].Q[0])*(node[n1].Q[1]/node[n1].Q[0])
                + (node[n1].Q[2]/node[n1].Q[0])*(node[n1].Q[2]/node[n1].Q[0])));
        dRdX_dXdBeta[n1][1] += Pressure*Nx_b;
        dRdX_dXdBeta[n1][2] += Pressure*Ny_b;

        // dResi/dBeta at Node2
        Pressure = (Gamma - 1.0)*(node[n2].Q[3] - 0.5*node[n2].Q[0]*((node[n2].Q[1]/node[n2].Q[0])*(node[n2].Q[1]/node[n2].Q[0])
                + (node[n2].Q[2]/node[n2].Q[0])*(node[n2].Q[2]/node[n2].Q[0])));
        dRdX_dXdBeta[n2][1] += Pressure*Nx_b;
        dRdX_dXdBeta[n2][2] += Pressure*Ny_b;
    }

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 2
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 2: 2000-2999
        if ((boundaryEdge[iBEdge].bcType < 2000) || (boundaryEdge[iBEdge].bcType > 2999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Get the Coordinate Derivates
        x_b[0] = dXdBeta[n1];
        x_b[1] = dXdBeta[n2];
        y_b[0] = dYdBeta[n1];
        y_b[1] = dYdBeta[n2];

        // Compute the Normal Derivatives
        Nx_b  = +0.5*(y_b[1] - y_b[0]);
        Ny_b  = -0.5*(x_b[1] - x_b[0]);

        // dResi/dBeta at Node1
        Pressure = (Gamma - 1.0)*(node[n1].Q[3] - 0.5*node[n1].Q[0]*((node[n1].Q[1]/node[n1].Q[0])*(node[n1].Q[1]/node[n1].Q[0])
                + (node[n1].Q[2]/node[n1].Q[0])*(node[n1].Q[2]/node[n1].Q[0])));
        dRdX_dXdBeta[n1][1] += Pressure*Nx_b;
        dRdX_dXdBeta[n1][2] += Pressure*Ny_b;

        // dResi/dBeta at Node2
        Pressure = (Gamma - 1.0)*(node[n2].Q[3] - 0.5*node[n2].Q[0]*((node[n2].Q[1]/node[n2].Q[0])*(node[n2].Q[1]/node[n2].Q[0])
                + (node[n2].Q[2]/node[n2].Q[0])*(node[n2].Q[2]/node[n2].Q[0])));
        dRdX_dXdBeta[n2][1] += Pressure*Nx_b;
        dRdX_dXdBeta[n2][2] += Pressure*Ny_b;
    }

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Freestream Boundary: BCTag = 3
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 3: 3000-3999 : Free Stream
        if ((boundaryEdge[iBEdge].bcType < 3000) || (boundaryEdge[iBEdge].bcType > 3999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // dResi/dBeta at Node1 and Node2
        for (i = 0; i < 4; i++) {
            dRdX_dXdBeta[n1][i] = 0.0;
            dRdX_dXdBeta[n2][i] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_dFplusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FplusdBeta) {
    double C, Ubar, Mbar, Pressure;
    double Ubar_b, Fplus1;
    double Var1, Var2, Var3;

    Pressure = (Gamma - 1.0)*(Q[3] - 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0])));
    C        = sqrt(Gamma*Pressure/Q[0]);
    Ubar     = ((Q[1]/Q[0])*uNx) + ((Q[2]/Q[0])*uNy);
    Ubar_b   = ((Q[1]/Q[0])*uNx_b) + ((Q[2]/Q[0])*uNy_b);
    Mbar     = Ubar/C;

    if (fabs(Mbar) <= 1.0) {
        Fplus1 = Mag*0.25*Q[0]*C*(Mbar + 1.0)*(Mbar + 1.0);

        Var1 = 0.5*Mag*Q[0]*(Mbar + 1.0)*Ubar_b;
        Var2 = 0.25*C*Q[0]*(Mbar + 1.0)*(Mbar + 1.0)*Mag_b;
        FplusdBeta[0] = Var1 + Var2;

        Var1 = FplusdBeta[0]*((uNx/Gamma)*(-Ubar + 2.0*C) + Q[1]/Q[0]);
        Var2 = Fplus1*(uNx_b/ Gamma)*(-Ubar + 2.0*C);
        Var3 = Fplus1*(uNx/Gamma)*(-Ubar_b);
        FplusdBeta[1] = Var1 + Var2 + Var3;

        Var1 = FplusdBeta[0]*((uNy/Gamma)*(-Ubar + 2.0*C) + Q[2]/Q[0]);
        Var2 = Fplus1*(uNy_b/Gamma)*(-Ubar + 2.0*C);
        Var3 = Fplus1*(uNy/Gamma)*(-Ubar_b);
        FplusdBeta[2] = Var1 + Var2 + Var3;

        Var1 = FplusdBeta[0]*((((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0))
                + (0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(Q[0]*Q[0])));
        Var2 = Fplus1*(((-2.0*(Gamma - 1.0)*Ubar*Ubar_b) + (2.0*(Gamma - 1.0)*Ubar_b*C))/(Gamma*Gamma - 1.0));
        FplusdBeta[3] = Var1 + Var2;
    } else if (Mbar > 1.0) {
        FplusdBeta[0] = Q[0]*(Mag_b*Ubar + Mag*Ubar_b);
        FplusdBeta[1] = Q[1]*(Mag_b*Ubar + Mag*Ubar_b) + Pressure*(Mag_b*uNx + Mag*uNx_b);
        FplusdBeta[2] = Q[2]*(Mag_b*Ubar + Mag*Ubar_b) + Pressure*(Mag_b*uNy + Mag*uNy_b);
        FplusdBeta[3] = (Q[3] + Pressure)*(Mag_b*Ubar + Mag*Ubar_b);
    } else {
        FplusdBeta[0] = 0.0;
        FplusdBeta[1] = 0.0;
        FplusdBeta[2] = 0.0;
        FplusdBeta[3] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_dFminusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
            double uNx_b, double uNy_b, double Mag_b, double *FminusdBeta) {
    double C, Ubar, Mbar, Pressure;
    double Ubar_b, Fminus1;
    double Var1, Var2, Var3;

    Pressure = (Gamma - 1.0)*(Q[3] - 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0])));
    C        = sqrt(Gamma*Pressure/Q[0]);
    Ubar     = ((Q[1]/Q[0])*uNx) + ((Q[2]/Q[0])*uNy);
    Ubar_b   = ((Q[1]/Q[0])*uNx_b) + ((Q[2]/Q[0])*uNy_b);
    Mbar     = Ubar/C;

    if (fabs(Mbar) <= 1.0) {
        Fminus1 = -Mag*0.25*Q[0]*C*(Mbar - 1.0)*(Mbar - 1.0);

        Var1 = -0.5*Mag*Q[0]*(Mbar - 1.0)*Ubar_b;
        Var2 = -0.25*C*Q[0]*(Mbar - 1.0)*(Mbar - 1.0)*Mag_b;
        FminusdBeta[0] = Var1 + Var2;

        Var1 = FminusdBeta[0]*((uNx/Gamma)*(-Ubar - 2.0*C) + Q[1]/Q[0]);
        Var2 = Fminus1*(uNx_b/ Gamma)*(-Ubar - 2.0*C);
        Var3 = Fminus1*(uNx/Gamma)*(-Ubar_b);
        FminusdBeta[1] = Var1 + Var2 + Var3;

        Var1 = FminusdBeta[0]*((uNy/Gamma)*(-Ubar - 2.0*C) + Q[2]/Q[0]);
        Var2 = Fminus1*(uNy_b/Gamma)*(-Ubar - 2.0*C);
        Var3 = Fminus1*(uNy/Gamma)*(-Ubar_b);
        FminusdBeta[2] = Var1 + Var2 + Var3;

        Var1 = FminusdBeta[0]*((((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0))
                + (0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(Q[0]*Q[0])));
        Var2 = Fminus1*(((-2.0*(Gamma - 1.0)*Ubar*Ubar_b) - (2.0*(Gamma - 1.0)*Ubar_b*C))/(Gamma*Gamma - 1.0));
        FminusdBeta[3] = Var1 + Var2;

    } else if (Mbar < -1.0) {
        FminusdBeta[0] = Q[0]*(Mag_b*Ubar + Mag*Ubar_b);
        FminusdBeta[1] = Q[1]*(Mag_b*Ubar + Mag*Ubar_b) + Pressure*(Mag_b*uNx + Mag*uNx_b);
        FminusdBeta[2] = Q[2]*(Mag_b*Ubar + Mag*Ubar_b) + Pressure*(Mag_b*uNy + Mag*uNy_b);
        FminusdBeta[3] = (Q[3] + Pressure)*(Mag_b*Ubar + Mag*Ubar_b);
    } else {
        FminusdBeta[0] = 0.0;
        FminusdBeta[1] = 0.0;
        FminusdBeta[2] = 0.0;
        FminusdBeta[3] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Verify_dRdX_dXdBeta() {
    int i, iNode;
    double R_Up[4], R_Dn[4];
    double FD_dRdX_dXdBeta[4];

    if (Order != 1)
        return;

    // No FD Validation computation
    if (FDNodeID < 0 || FDNodeID >= mesh.nnodes)
        return;

    // Store the Old Residual
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++)
            FDNode[iNode].Resi_Old[i] = node[iNode].Resi[i];
    }

    // Now Purturb the Node Up with Epsilon
    node[FDNodeID].y += FDEpsilon;
    // Now Compute the New Residual
    Compute_Geometric_Properties();
    Compute_Residual();
    Compute_Boundary_Residual();

    //  R_Up
    for (i = 0; i < 4; i++)
        R_Up[i] = node[FDNodeID].Resi[i];

    // Now Purturb the Node down with Epsilon
    node[FDNodeID].y -= 2.0*FDEpsilon;
    // Now Compute the New Residual
    Compute_Geometric_Properties();
    Compute_Residual();
    Compute_Boundary_Residual();

    //  R_Dn
    for (i = 0; i < 4; i++)
        R_Dn[i] = node[FDNodeID].Resi[i];

    // Now compute the Finite Diffence dR/dBeta
    for (i = 0; i < 4; i++)
        FD_dRdX_dXdBeta[i] = (R_Up[i] - R_Dn[i])/(2.0*FDEpsilon);

    printf("Finite Difference dRdX_dXdBeta :\n %15.9e %15.9e %15.9e %15.9e \n",
            FD_dRdX_dXdBeta[0], FD_dRdX_dXdBeta[1], FD_dRdX_dXdBeta[2], FD_dRdX_dXdBeta[3]);
    printf("Analytic Difference dRdX_dXdBeta:\n %15.9e %15.9e %15.9e %15.9e \n",
            dRdX_dXdBeta[FDNodeID][0], dRdX_dXdBeta[FDNodeID][1], dRdX_dXdBeta[FDNodeID][2], dRdX_dXdBeta[FDNodeID][3]);

    // Restore back the Changes
    node[FDNodeID].y += FDEpsilon;
    Compute_Geometric_Properties();
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++)
            node[iNode].Resi[i] = FDNode[iNode].Resi_Old[i];
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Direct_dQdBeta() {
    int i, j, k;
    int iNode;
    double lrms;

    // Update the CRS Design Matrix with solver dR/dQ and
    for (i = 0; i < DesignBlockMatrix.CRSSize; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++)
                DesignBlockMatrix.A[i][j][k] = SolverBlockMatrix.A[i][j][k];
        }
    }

    // Update the CRS Design Matrix Vector B with dR/dBeta
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++) {
            DesignBlockMatrix.B[iNode][j] = -dRdX_dXdBeta[iNode][j];
        }
    }

    // Solve for dQdBeta
    InnerNIteration = 0;
    Relaxation = 0.9;
    lrms = MC_Iterative_Block_Jacobi_CRS(InnerNIteration, Relaxation, DesignBlockMatrix);
//    lrms = MC_Iterative_Block_LU_Jacobi_CRS(InnerNIteration, 2, DesignBlockMatrix);

#ifdef VERBOSE
    info("dIdQ LRMS = %10.5e", lrms);
#endif

    // Get the Solution
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++)
            dQdBeta[iNode][j] = DesignBlockMatrix.X[iNode][j];
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Verify_Direct_dQdBeta() {
    double FD_dQdBeta[4];
    double Q1[4], Q2[4];
    int i, iNode, var;
    size_t sdum;
    FILE *fp;

    if (Order != 1)
        return;

    // No FD Validation computation
    if (FDNodeID >= mesh.nnodes)
        return;

    // Create Solver1 and Converge
    {
        Euler2D_Solver_VanLeer Solver1;
        Solver1.Get_Solver_Inputs_VanLeer("Euler_Solver1.inp");
        Solver1.Solver_Prepare();
        Solver1.Solve();
        Solver1.Solver_Finalize();
    };
    // Create Solver2 and Converge
    {
        Euler2D_Solver_VanLeer Solver2;
        Solver2.Get_Solver_Inputs_VanLeer("Euler_Solver2.inp");
        Solver2.Solver_Prepare();
        Solver2.Solve();
        Solver2.Solver_Finalize();
    };

    // Read the Solver1 Q for the Perticular Node
    // Open Mesh File
    if ((fp = fopen("OutputRestart1.q", "rb")) == (FILE *) NULL)
        error("Verify_dQdBeta: Unable to Open Mesh File OutputRestart1.q");
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != mesh.nnodes)
        error("Verify_dQdBeta: Mismatched in Restart and Mesh File OutputRestart1.q");
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    for (iNode = 0; iNode < FDNodeID+1; iNode++) {
        sdum = fread(&Q1[0], sizeof (double), 1, fp);
        sdum = fread(&Q1[1], sizeof (double), 1, fp);
        sdum = fread(&Q1[2], sizeof (double), 1, fp);
        sdum = fread(&Q1[3], sizeof (double), 1, fp);
    }
    fclose(fp);

    // Read the Solver2 Q for the Perticular Node
    // Open Mesh File
    if ((fp = fopen("OutputRestart2.q", "rb")) == (FILE *) NULL)
        error("Verify_dQdBeta: Unable to Open Mesh File OutputRestart2.q");
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != mesh.nnodes)
        error("Verify_dQdBeta: Mismatched in Restart and Mesh File OutputRestart2.q");
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    for (iNode = 0; iNode < FDNodeID+1; iNode++) {
        sdum = fread(&Q2[0], sizeof (double), 1, fp);
        sdum = fread(&Q2[1], sizeof (double), 1, fp);
        sdum = fread(&Q2[2], sizeof (double), 1, fp);
        sdum = fread(&Q2[3], sizeof (double), 1, fp);
    }
    fclose(fp);

    // Compute Finite Difference dQ/dBeta
    for (i = 0; i < 4; i++)
        FD_dQdBeta[i] = (Q1[i] - Q2[i])/(2.0*FDEpsilon);

    printf("Verification of dQ/dBeta at Node %d \n", FDNodeID);
    printf("Finite Difference:\n %15.9e %15.9e %15.9e %15.9e \n",
            FD_dQdBeta[0], FD_dQdBeta[1], FD_dQdBeta[2], FD_dQdBeta[3]);
    printf("Analytic Difference:\n %15.9e %15.9e %15.9e %15.9e \n",
            dQdBeta[FDNodeID][0], dQdBeta[FDNodeID][1], dQdBeta[FDNodeID][2], dQdBeta[FDNodeID][3]);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Direct_dIdQ_dQdBeta(int iDesignVariable) {
    int  i, iBEdge, iedge;
    int bn[2];
    double x[2], Nx;
    double y[2], Ny;
    double Var1, Var2;
    double Rho,   U,   V,   E,   P[2],   Pavg;
    double Rho_b, U_b, V_b, E_b, P_b[2], Pavg_b;
    double Lift_b, Drag_b, CoeffLift_b, CoeffDrag_b;

    dIdQ_dQdBeta[iDesignVariable] = 0.0;
    Lift_b = 0.0;
    Drag_b = 0.0;

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[iBEdge].bcType < 1000) || (boundaryEdge[iBEdge].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        bn[0] = edge[iedge].node1;
        bn[1] = edge[iedge].node2;

        // Get the Coordinates
        x[0]  = node[bn[0]].x;
        y[0]  = node[bn[0]].y;
        x[1]  = node[bn[1]].x;
        y[1]  = node[bn[1]].y;

        // Compute the Normal
        Nx    = +(y[1] - y[0]);
        Ny    = -(x[1] - x[0]);

        // Compute the partial derivative of P => dPdBeta for Node1 & Node2
        for (i = 0; i < 2; i++) {
            // Calculate Rho and Partial derivatives
            Rho    = node[bn[i]].Q[0];
            Rho_b = dQdBeta[bn[i]][0];

            // Calculate U and Partial derivatives
            U   = node[bn[i]].Q[1]/node[bn[i]].Q[0];
            U_b = (Rho*dQdBeta[bn[i]][1] - Rho_b*node[bn[i]].Q[1])/(Rho * Rho);

            // Calculate V and Partial derivatives
            V   = node[bn[i]].Q[2]/node[bn[i]].Q[0];
            V_b = (Rho*dQdBeta[bn[i]][2] - Rho_b*node[bn[i]].Q[2])/(Rho * Rho);

            // Calculate E and Partial derivatives
            E   = node[bn[i]].Q[3];
            E_b = dQdBeta[bn[i]][3];

            // Calculate P and Partial derivatives
            P[i]   = (Gamma - 1.0) * (E - 0.5 * Rho * (U * U + V * V));
            Var1   = 0.5*Rho_b*(U*U + V*V);
            Var2   = 0.5*Rho*(2.0*U*U_b + 2.0*V*V_b);
            P_b[i] = (Gamma - 1.0)*(E_b - Var1 - Var2);
        }

        // Get the Average Pressure and its derivates
        Pavg   = 0.5*(P[0] + P[1]);
        Pavg_b = 0.5*(P_b[0] + P_b[1]);

        // Compute Lift and its derivatives
        Lift_b += Pavg_b*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny);

        // Compute Drag and its derivatives
        Drag_b += Pavg_b*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny);
    }

    // Compute Coefficient of Lift and Drag and its derivates
    CoeffLift_b = Lift_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag_b = Drag_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute (dI/dQ)*(dQ/dBeta)
    dIdQ_dQdBeta[iDesignVariable] = (CoeffLift - Ref_CoeffLift)*CoeffLift_b
            + (CoeffDrag - Ref_CoeffDrag)*CoeffDrag_b;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_dIdX_dXdBeta(int iDesignVariable) {
    int  i, n1, n2, iBEdge, iedge, Size;
    int bn[2];
    double x1, x2, Nx, x_b[2], Nx_b, uNx, uNx_b, Lx, Lx_b;
    double y1, y2, Ny, y_b[2], Ny_b, uNy, uNy_b, Ly, Ly_b;
    double MagN, MagN_b, MagL, MagL_b;
    double Rho,   U,   V,   E,   P[2],   Pavg;
    double Lift_b, Drag_b, CoeffLift_b, CoeffDrag_b;
    List bnodeRef;
    int    *NodeMapRef = NULL;
    double *PressureMapRef = NULL;
    double *CoordXRef = NULL;
    double *CoordYRef = NULL;
    double PressureCorrection_b;

    n1 = n2 = -1;
    Lift_b = 0.0;
    Drag_b = 0.0;
    PressureCorrection_b = 0.0;

    if ((CostType == 3) || (CostType == 4)) {
        // Read the Reference Boundary Pressure Distribution
        Read_Boundary_Field("Pressure_Baseline.dat", &Size, &NodeMapRef,
                &CoordXRef, &CoordYRef, &PressureMapRef);
        for (i = 0; i < Size; i++)
            bnodeRef.Check_List(NodeMapRef[i]);
    }

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[iBEdge].bcType < 1000) || (boundaryEdge[iBEdge].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        bn[0] = edge[iedge].node1;
        bn[1] = edge[iedge].node2;

        // Get the Coordinate
        x1 = node[bn[0]].x;
        y1 = node[bn[0]].y;
        x2 = node[bn[1]].x;
        y2 = node[bn[1]].y;

        // Get the Coordinate Derivates
        x_b[0] = dXdBeta[bn[0]];
        y_b[0] = dYdBeta[bn[0]];
        x_b[1] = dXdBeta[bn[1]];
        y_b[1] = dYdBeta[bn[1]];

        // Compute the Length
        Lx   = (x2 - x1);
        Ly   = (y2 - y1);
        MagL = sqrt(Lx*Lx + Ly*Ly);

        // Compute the Length Derivative
        Lx_b   = x_b[1] - x_b[0];
        Ly_b   = y_b[1] - y_b[0];
        MagL_b = (Lx*Lx_b + Ly*Ly_b)/MagL;
        
        // Compute the Normals
        Nx   = +(y2 - y1);
        Ny   = -(x2 - x1);
        MagN = sqrt(Nx*Nx + Ny*Ny);
        uNx  = Nx/MagN;
        uNy  = Ny/MagN;
        
        // Compute the Unit Normal Derivatives
        Nx_b   = +(y_b[1] - y_b[0]);
        Ny_b   = -(x_b[1] - x_b[0]);
        MagN_b = (Nx*Nx_b + Ny*Ny_b)/MagN;
        uNx_b  = (Nx_b*MagN - MagN_b*Nx)/(MagN*MagN);
        uNy_b  = (Ny_b*MagN - MagN_b*Ny)/(MagN*MagN);

        // Compute the Pressure for Node1 & Node2
        for (i = 0; i < 2; i++) {
            // Calculate Rho
            Rho    = node[bn[i]].Q[0];
            // Calculate U
            U   = node[bn[i]].Q[1]/node[bn[i]].Q[0];
            // Calculate V
            V   = node[bn[i]].Q[2]/node[bn[i]].Q[0];
            // Calculate E
            E   = node[bn[i]].Q[3];
            // Calculate P
            P[i]   = (Gamma - 1.0) * (E - 0.5 * Rho * (U * U + V * V));
        }

        // Get the Average Pressure
        Pavg   = 0.5*(P[0] + P[1]);

        // Compute Lift derivatives => (dLift/dX)*(dX/dBeta)
        Lift_b += Pavg*(-sin(Ref_Alpha)*uNx_b + cos(Ref_Alpha)*uNy_b)*MagN
                + Pavg*(-sin(Ref_Alpha)*uNx   + cos(Ref_Alpha)*uNy)*MagN_b;

        // Compute Drag derivatives => (dDrag/dX)*(dX/dBeta)
        Drag_b += Pavg*( cos(Ref_Alpha)*uNx_b + sin(Ref_Alpha)*uNy_b)*MagN
                + Pavg*( cos(Ref_Alpha)*uNx   + sin(Ref_Alpha)*uNy)*MagN_b;
        
        // Compute Pressure Correction Derivative => (dPcorr/dX)*(dX/dBeta)
        if ((CostType == 3) || (CostType == 4)) {
            n1 = bnodeRef.Index(bn[0]);
            n2 = bnodeRef.Index(bn[1]);
        }
        // Linear Implementation
        if (CostType == 3) {
            PressureCorrection_b += (P[0] - PressureMapRef[n1])*MagL_b;
            PressureCorrection_b += (P[1] - PressureMapRef[n2])*MagL_b;
        }
        // Quadratic Implementation
        if (CostType == 4) {
            PressureCorrection_b += 0.5*(P[0] - PressureMapRef[n1])*(P[0] - PressureMapRef[n1])*MagL_b;
            PressureCorrection_b += 0.5*(P[1] - PressureMapRef[n2])*(P[1] - PressureMapRef[n2])*MagL_b;
        }
    }

    // Compute Coefficient of Lift and Drag derivates
    CoeffLift_b = Lift_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag_b = Drag_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute (dI/dX)*(dX/dBeta)
    dIdX_dXdBeta[iDesignVariable] = 0.0;
    // Lift Only
    if (CostType == 0)
        dIdX_dXdBeta[iDesignVariable] = (CoeffLift - Ref_CoeffLift)*CoeffLift_b;

    // Drag Only
    if (CostType == 1)
        dIdX_dXdBeta[iDesignVariable] = (CoeffDrag - Ref_CoeffDrag)*CoeffDrag_b;

    // Lift and Drag Only
    if (CostType == 2)
        dIdX_dXdBeta[iDesignVariable] = (CoeffLift - Ref_CoeffLift)*CoeffLift_b
            + (CoeffDrag - Ref_CoeffDrag)*CoeffDrag_b;

    // Pressure
    if ((CostType == 3) || (CostType == 4)) {
        dIdX_dXdBeta[iDesignVariable] = PressureCorrection_b;
    
        // Free Memory
        if (NodeMapRef != NULL)
            free(NodeMapRef);
        if (PressureMapRef != NULL)
            free(PressureMapRef);
        if (CoordXRef != NULL)
            free(CoordXRef);
        if (CoordYRef != NULL)
            free(CoordYRef);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Direct_dIdBeta(int iDesignVariable) {
    dIdBeta[iDesignVariable] = dIdQ_dQdBeta[iDesignVariable]
                                + dIdX_dXdBeta[iDesignVariable];
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Cost() {
    int  i, n1, n2, iBEdge, iedge, Size;
    int bn[2];
    double x[2], Nx;
    double y[2], Ny;
    double Rho,   U,   V,   E,   P[2],   Pavg, MagN;
    List bnode;
    List bnodeRef;
    int    *NodeMap = NULL;
    int    *NodeMapRef = NULL;
    double *PressureMap = NULL;
    double *PressureMapRef = NULL;
    double *CoordXRef = NULL;
    double *CoordYRef = NULL;
    double *CoordX = NULL;
    double *CoordY = NULL;
    double PressureCorrection;

    n1 = n2 = -1;
    Lift = 0.0;
    Drag = 0.0;
    PressureCorrection = 0.0;

    if ((CostType == 3) || (CostType == 4)) {
        // Read the Reference Boundary Pressure Distribution
        Read_Boundary_Field("Pressure_Baseline.dat", &Size, &NodeMapRef,
                &CoordXRef, &CoordYRef, &PressureMapRef);
        for (i = 0; i < Size; i++)
            bnodeRef.Check_List(NodeMapRef[i]);
    }
    
    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[iBEdge].bcType < 1000) || (boundaryEdge[iBEdge].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        bn[0] = edge[iedge].node1;
        bn[1] = edge[iedge].node2;

        // Add to Boundary Node List
        bnode.Check_List(bn[0]);
        bnode.Check_List(bn[1]);
        
        // Get the Coordinates
        x[0]  = node[bn[0]].x;
        y[0]  = node[bn[0]].y;
        x[1]  = node[bn[1]].x;
        y[1]  = node[bn[1]].y;

        // Compute the Normal
        Nx    = +(y[1] - y[0]);
        Ny    = -(x[1] - x[0]);
        MagN  = sqrt(Nx*Nx + Ny*Ny);

        // Compute the Pressure for Node1 & Node2
        for (i = 0; i < 2; i++) {
            // Calculate Rho
            Rho    = node[bn[i]].Q[0];
            // Calculate U
            U   = node[bn[i]].Q[1]/node[bn[i]].Q[0];
            // Calculate V
            V   = node[bn[i]].Q[2]/node[bn[i]].Q[0];
            // Calculate E
            E   = node[bn[i]].Q[3];
            // Calculate P
            P[i]   = (Gamma - 1.0) * (E - 0.5 * Rho * (U * U + V * V));
        }

        // Get the Average Pressure
        Pavg   = 0.5*(P[0] + P[1]);

        // Compute Lift and Drag
        Lift += (-Pavg*sin(Ref_Alpha)*Nx + Pavg*cos(Ref_Alpha)*Ny);
        Drag += ( Pavg*cos(Ref_Alpha)*Nx + Pavg*sin(Ref_Alpha)*Ny);

        // Pressure
        if ((CostType == 3) || (CostType == 4)) {
            // Compute Pressure Correction
            n1 = bnodeRef.Index(bn[0]);
            n2 = bnodeRef.Index(bn[1]);
        }
        
        // Linear Implementation
        if (CostType == 3) {
            PressureCorrection += (P[0] - PressureMapRef[n1])*MagN;
            PressureCorrection += (P[1] - PressureMapRef[n2])*MagN;
        }

        // Quadratic Implementation
        if (CostType == 4) {
            PressureCorrection += 0.5*(P[0] - PressureMapRef[n1])*(P[0] - PressureMapRef[n1])*MagN;
            PressureCorrection += 0.5*(P[1] - PressureMapRef[n2])*(P[1] - PressureMapRef[n2])*MagN;
        }
    }
    
    // Compute Coefficient of Lift and Drag and its derivates
    CoeffLift   = Lift/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag   = Drag/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    
    // Compute Cost Function
    I = 0.0;

    // Lift Only
    if (CostType == 0)
        I = 0.5*(CoeffLift - Ref_CoeffLift)*(CoeffLift - Ref_CoeffLift);

    // Drag Only
    if (CostType == 1)
        I = 0.5*(CoeffDrag - Ref_CoeffDrag)*(CoeffDrag - Ref_CoeffDrag);

    // Lift and Drag
    if (CostType == 2)
        I = 0.5*(CoeffLift - Ref_CoeffLift)*(CoeffLift - Ref_CoeffLift)
            + 0.5*(CoeffDrag - Ref_CoeffDrag)*(CoeffDrag - Ref_CoeffDrag);

    // Pressure
    if ((CostType == 3) || (CostType == 4)) {
        I = PressureCorrection;
        // Write Boundary Pressure Map
        NodeMap     = (int *) malloc(bnode.max*sizeof(int));
        PressureMap = (double *) malloc(bnode.max*sizeof(double));
        CoordX      = (double *) malloc(bnode.max*sizeof(double));
        CoordY      = (double *) malloc(bnode.max*sizeof(double));
        for (i = 0; i < bnode.max; i++) {
            NodeMap[i] = bnode.list[i];
            // Get the coordinates
            CoordX[i] = node[bnode.list[i]].x;
            CoordY[i] = node[bnode.list[i]].y;
            // Calculate Rho
            Rho = node[bnode.list[i]].Q[0];
            // Calculate U
            U   = node[bnode.list[i]].Q[1]/node[bnode.list[i]].Q[0];
            // Calculate V
            V   = node[bnode.list[i]].Q[2]/node[bnode.list[i]].Q[0];
            // Calculate E
            E   = node[bnode.list[i]].Q[3];
            // Calculate P
            PressureMap[i] = (Gamma - 1.0) * (E - 0.5 * Rho * (U * U + V * V));
        }
        Write_Boundary_Field("Pressure.dat", bnode.max, NodeMap, CoordX, CoordY, PressureMap);

        // Free Memory
        if (NodeMapRef != NULL)
            free(NodeMapRef);
        if (PressureMapRef != NULL)
            free(PressureMapRef);
        if (NodeMap != NULL)
            free(NodeMap);
        if (PressureMap != NULL)
            free(PressureMap);
        if (CoordXRef != NULL)
            free(CoordXRef);
        if (CoordYRef != NULL)
            free(CoordYRef);
        if (CoordX != NULL)
            free(CoordX);
        if (CoordY != NULL)
            free(CoordY);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Adjoint_dIdQ() {
    int  i, j, iNode, iBEdge, iedge, Size;
    int bn[2], pbn[2];
    double Nx, Ny, Lx, Ly, MagL;
    double x[2], y[2], Var1, Var2;
    double U, V;
    double Rho, P[2], E, Pavg;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    // I_Q    => Cost Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], E_Q[4], P_Q[2][4], Pavg_Q[2][4], I_Q[4];
    List bnodeRef;
    int    *NodeMapRef = NULL;
    double *PressureMapRef = NULL;
    double *CoordXRef = NULL;
    double *CoordYRef = NULL;

    // Initialize
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < 4; j++)
            dIdQ[iNode][j] = 0.0;
    }

    // Read the Reference Boundary Pressure Distribution
    if ((CostType == 3) || (CostType == 4)) {
        Read_Boundary_Field("Pressure_Baseline.dat", &Size, &NodeMapRef,
                &CoordXRef, &CoordYRef, &PressureMapRef);
        for (i = 0; i < Size; i++)
            bnodeRef.Check_List(NodeMapRef[i]);
    }
    
    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (iBEdge = 0; iBEdge < mesh.nbedges; iBEdge++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[iBEdge].bcType < 1000) || (boundaryEdge[iBEdge].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[iBEdge].edgeNumber;
        bn[0] = edge[iedge].node1;
        bn[1] = edge[iedge].node2;

        // Get the coordiantes
        x[0]  = node[bn[0]].x;
        y[0]  = node[bn[0]].y;
        x[1]  = node[bn[1]].x;
        y[1]  = node[bn[1]].y;

        // Compute Length
        Lx    = x[1] - x[0];
        Ly    = y[1] - y[0];
        MagL  = sqrt(Lx*Lx + Ly*Ly);

        // Compute the Normal
        Nx    = +(y[1] - y[0]);
        Ny    = -(x[1] - x[0]);
        
        // Compute the partial derivative of P => dPdQ for Node1 & Node2
        for (i = 0; i < 2; i++) {
            // Calculate Rho and Partial derivatives
            Rho = node[bn[i]].Q[0];
            Rho_Q[0] = 1.0;
            Rho_Q[1] = 0.0;
            Rho_Q[2] = 0.0;
            Rho_Q[3] = 0.0;

            // Calculate U and Partial derivatives
            U = node[bn[i]].Q[1]/node[bn[i]].Q[0];
            U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[bn[i]].Q[1])) / (Rho * Rho);
            U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * node[bn[i]].Q[1])) / (Rho * Rho);
            U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * node[bn[i]].Q[1])) / (Rho * Rho);
            U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[bn[i]].Q[1])) / (Rho * Rho);

            // Calculate V and Partial derivatives
            V = node[bn[i]].Q[2]/node[bn[i]].Q[0];
            V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[bn[i]].Q[2])) / (Rho * Rho);
            V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * node[bn[i]].Q[2])) / (Rho * Rho);
            V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * node[bn[i]].Q[2])) / (Rho * Rho);
            V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[bn[i]].Q[2])) / (Rho * Rho);

            // Calculate E and Partial derivatives
            E = node[bn[i]].Q[3];
            E_Q[0] = 0.0;
            E_Q[1] = 0.0;
            E_Q[2] = 0.0;
            E_Q[3] = 1.0;

            // Calculate P and Partial derivatives
            P[i] = (Gamma - 1.0) * (node[bn[i]].Q[3] - 0.5 * Rho * (U * U + V * V));
            for (j = 0; j < 4; j++) {
                Var1 = 0.5*Rho_Q[j]*(U*U + V*V);
                Var2 = 0.5*Rho*(2.0*U*U_Q[j] + 2.0*V*V_Q[i]);
                P_Q[i][j] = (Gamma - 1.0) * (E_Q[j] - Var1 - Var2);
            }
        }

        // Get the Average Pressure and its derivates
        Pavg = 0.5*(P[0] + P[1]);
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 4; j++) {
                Pavg_Q[i][j] = 0.5*P_Q[i][j];
            }
        }

        // Get Pressure Correction Nodes
        if ((CostType == 3) || (CostType == 4)) {
            pbn[0] = bnodeRef.Index(bn[0]);
            pbn[1] = bnodeRef.Index(bn[1]);
        }

        Var1 = (2.0/(Gamma*Ref_Mach*Ref_Mach));
        // Computing dI/dQ = dCl/dQ + dCd/dQ + dPcor/dQ
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 4; j++) {
                // Lift - Contribution
                if (CostType == 0)
                    I_Q[j] = (CoeffLift - Ref_CoeffLift)*(Pavg_Q[i][j]*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny))/Var1;

                // Drag - Contribution
                if (CostType == 1)
                    I_Q[j] = (CoeffDrag - Ref_CoeffDrag)*(Pavg_Q[i][j]*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny))/Var1;

                // Lift and Drag - Contribution
                if (CostType == 2)
                    I_Q[j] = (CoeffLift - Ref_CoeffLift)*(Pavg_Q[i][j]*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny))/Var1
                            + (CoeffDrag - Ref_CoeffDrag)*(Pavg_Q[i][j]*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny))/Var1;

                // Linear Pressure - Contribution
                if (CostType == 3)
                    I_Q[j] = P_Q[i][j]*MagL;

                // Quadratic Pressure - Contribution
                if (CostType == 4)
                    I_Q[j] = (P[i] - PressureMapRef[pbn[i]])*P_Q[i][j]*MagL;
                
                // Total
                dIdQ[bn[i]][j] += I_Q[j];
            }
        }
    }

    if ((CostType == 3) || (CostType == 4)) {
        // Free Memory
        if (NodeMapRef != NULL)
            free(NodeMapRef);
        if (PressureMapRef != NULL)
            free(PressureMapRef);
        if (CoordXRef != NULL)
            free(CoordXRef);
        if (CoordYRef != NULL)
            free(CoordYRef);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Adjoint_Lambda() {
    int i, j, k, count;
    int iNode;
    double lrms;
    int *connect = NULL;
    
    // Update the CRS Design Matrix with solver trans[dR/dQ]
    // Allocate memory to store the transpose locations
    connect = (int *) malloc(DesignBlockMatrix.CRSSize*sizeof(int));
    // Step - 1 Get the transpose location
    count = 0;
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < DesignBlockMatrix.CRSSize; i++) {
            if (DesignBlockMatrix.JA[i] == iNode) {
                connect[count] = i;
                count++;
            }
        }
    }
    
    // Step - 2: Copy the transpose to new location
    for (i = 0; i < DesignBlockMatrix.CRSSize; i++) {
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++) {
            for (k = 0; k < DesignBlockMatrix.BlockSize; k++)
                DesignBlockMatrix.A[i][k][j] = SolverBlockMatrix.A[connect[i]][j][k];
        }
    }

    // Update the CRS Design Matrix Vector B with dR/dBeta
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++) {
            DesignBlockMatrix.B[iNode][j] = dIdQ[iNode][j];
        }
    }
    
    // Solve for dQdBeta
    InnerNIteration = 0;
    Relaxation = 0.9;
    lrms = MC_Iterative_Block_Jacobi_CRS(InnerNIteration, Relaxation, DesignBlockMatrix);
//    lrms = MC_Iterative_Block_LU_Jacobi_CRS(InnerNIteration, 2, DesignBlockMatrix);

#ifdef VERBOSE
    info("Lambda LRMS = %10.5e", lrms);
#endif
    // Get the Solution
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.BlockSize; j++)
            Lambda[iNode][j] = DesignBlockMatrix.X[iNode][j];
    }
    
    // Free the Memory
    free(connect);
    connect = NULL;   
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_Adjoint_dIdBeta(int iDesignVariable) {
    int iNode, j;
    double Var1, Var2;

    dIdBeta[iDesignVariable] = 0.0;
    // Compute the Vector, Vector Dot Product
    Var1 = 0.0;
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        Var2 = 0.0;
        for (j = 0; j < BlockSize; j++)
            Var2 += dRdX_dXdBeta[iNode][j] * Lambda[iNode][j];
        Var1 += Var2;
    }

    dIdBeta[iDesignVariable] = -Var1 + dIdX_dXdBeta[iDesignVariable];
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Compute_LinearElasticSmooth_dXdBeta(int cNode) {
    int iNode, iter;
    double relax;

    // Initialize
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        dXdBeta[iNode] = 0.0;
        dYdBeta[iNode] = 0.0;
    }

    // Right Now Only Purturbation in Y is done
    dXdBeta[cNode] = 0.0;
    dYdBeta[cNode] = 1.0;
    
    relax = 0.5;
    iter  = 0;
    LESmooth.Mesh_Smoother_Prepare();
    LESmooth.Initialize_Mesh_Smoother("domain_baseline.mesh", &DesignBlockMatrix);
//    LESmooth.Initialize_Mesh_Smoother(mesh, cell, edge, node, boundaryEdge,
//            boundaryNode, BNTag, &DesignBlockMatrix);
    LESmooth.Mesh_Smoother(iter, relax, dXdBeta, dYdBeta);
    LESmooth.Mesh_Smoother_Finalize();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Design_VanLeer::Write_Boundary_DesignFile(const char* FileName) {
    char   dumstring[130];
    FILE   *fp;

    // Open Mesh File
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_DesignFile: Unable to Write Design File %s", FileName);

    // Write Number of Design Variables
    fprintf(fp, " %6d\n", mesh.nbnodes);

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

