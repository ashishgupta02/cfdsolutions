/*
 * File:   Euler2D_Solver_Design.cpp
 * Author: Ashish Gupta
 *
 * Created on October 5, 2010, 4:53 PM
 */

#include <stdlib.h>
#include <iostream>
#include <limits.h>

#include "Utils.h"
#include "Euler2D_Design.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_Design::Euler2D_Solver_Design() {
    printf("=============================================================================\n");
    printf("      Euler2D : Computational Design                                         \n");
    printf("=============================================================================\n");

    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Init() {
    dXdBeta                      = NULL;
    dYdBeta                      = NULL;
    dQdBeta                      = NULL;
    dRdX_dXdBeta                 = NULL;
    dIdQ                         = NULL;
    dIdBeta                      = 0.0;
    dIdQ_dQdBeta                 = 0.0;
    dIdX_dXdBeta                 = 0.0;
    DesignBlockMatrix.nCOL       = 0;
    DesignBlockMatrix.nROW       = 0;
    DesignBlockMatrix.Block_nCol = 0;
    DesignBlockMatrix.Block_nRow = 0;
    DesignBlockMatrix.DIM        = 0;
    DesignBlockMatrix.A          = NULL;
    DesignBlockMatrix.B          = NULL;
    DesignBlockMatrix.IA         = NULL;
    DesignBlockMatrix.IAU        = NULL;
    DesignBlockMatrix.JA         = NULL;
    DesignBlockMatrix.X          = NULL;
    Ref_CoeffLift                = 0.0;
    Ref_CoeffDrag                = 0.0;
    CoeffLift                    = 0.0;
    CoeffDrag                    = 0.0;
    Lift                         = 0.0;
    Drag                         = 0.0;
    I                            = 0.0;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_Design::~Euler2D_Solver_Design() {

}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Design() {
    int AddTime;

    // Extract Connectivity
    WKA_ExtractCells();
    Tag_Boundary_Nodes();
    
    // Compute the Expensive Geometric Properties
    Compute_Geometric_Properties();

    // Initialize the solution
    Initialize_Solution();

    // Read Restart File
    Read_RestartFile(InputRestartFileName.c_str());

    // Create CRS Solver Matrix
    Create_CRS_BlockMatrix();

    // Initialize the Design
    Initialize_Design();

//    Write_dXdBeta("dXdBeta.q");
//    exit(1);
    
    // Create CRS Design Matrix
    Create_CRS_DesignBlockMatrix();

    // Compute Local Time Stepping
    CFL_Ramp += 1;
    Compute_DeltaTime(CFL_Ramp);

    // Compute Residuals
    Compute_Gauss_Gradient();
    Compute_Residual();
    Compute_Boundary_Residual();

    // Compute and Fill CRS Solver Matrix
    AddTime = 0;
    Compute_CRS_BlockMatrix(AddTime);
    Compute_Boundary_CRS_BlockMatrix();

    // Read the dX/dBeta
    Read_dXdBeta("dXdBeta.q");

    // Compute (dR/dX)*(dX/dBeta)
    Compute_dRdX_dXdBeta();
    Verify_dRdX_dXdBeta();
    Write_dRdX_dXdBeta("dRdX_dXdBeta.q");
    
    // Compute dQ/dBeta
    Compute_dQdBeta();
    Verify_dQdBeta();
    Write_dQdBeta("dQdBeta.q");

    // Compute Cost and its derivates
    Compute_Cost();
    Compute_dIdQ_dQdBeta();
    Compute_dIdX_dXdBeta();
    Compute_dIdBeta();
    printf("Cl = %lf, Cd = %lf, L = %lf, D = %lf, I = %lf\n",
            CoeffLift, CoeffDrag, Lift, Drag, I);
    printf("dIdQ_dQdBeta = %lf, dIdX_dXdBeta = %lf, dIdBeta = %lf\n", dIdQ_dQdBeta, dIdX_dXdBeta, dIdBeta);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Initialize_Design() {
    int i, j;

    // Allocate the memory for design variables
    dXdBeta      = new double[mesh.nnodes];
    dYdBeta      = new double[mesh.nnodes];
    dRdX_dXdBeta = new double*[mesh.nnodes];
    dQdBeta      = new double*[mesh.nnodes];
    dIdQ         = new double*[mesh.nnodes];
    for (i = 0; i < mesh.nnodes; i++) {
        dXdBeta[i]      = 0.0;
        dYdBeta[i]      = 0.0;
        dRdX_dXdBeta[i] = new double[4];
        dQdBeta[i]      = new double[4];
        dIdQ[i]         = new double[4];
        for (j = 0; j < 4; j++) {
            dRdX_dXdBeta[i][j] = 0.0;
            dQdBeta[i][j]      = 0.0;
            dIdQ[i][j]         = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Create_CRS_DesignBlockMatrix() {
    int i, j, k;

    // Get the value from Solver Block Matrix
    DesignBlockMatrix.nROW       = BlockMatrix.nROW;
    DesignBlockMatrix.nCOL       = BlockMatrix.nCOL;
    DesignBlockMatrix.Block_nRow = BlockMatrix.Block_nRow;
    DesignBlockMatrix.Block_nCol = BlockMatrix.Block_nCol;
    DesignBlockMatrix.DIM        = BlockMatrix.DIM;
    // Share the memory with Solver Block Matrix
    DesignBlockMatrix.IA         = BlockMatrix.IA;
    DesignBlockMatrix.IAU        = BlockMatrix.IAU;
    DesignBlockMatrix.JA         = BlockMatrix.JA;

    // Allocate Memory for Design Specific Computations
    // Allocate Memory for CRS Matrix
    DesignBlockMatrix.A = NULL;
    DesignBlockMatrix.A = (double ***) malloc (DesignBlockMatrix.DIM*sizeof(double**));
#ifdef DEBUG
    if (DesignBlockMatrix.A == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < DesignBlockMatrix.DIM; i++) {
        DesignBlockMatrix.A[i] = NULL;
        DesignBlockMatrix.A[i] = (double **) malloc (DesignBlockMatrix.Block_nRow*sizeof(double*));
#ifdef DEBUG
        if (DesignBlockMatrix.A[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 2");
#endif
        for (j = 0; j < DesignBlockMatrix.Block_nRow; j++) {
            DesignBlockMatrix.A[i][j] = NULL;
            DesignBlockMatrix.A[i][j] = (double *) malloc (DesignBlockMatrix.Block_nCol*sizeof(double));
#ifdef DEBUG
            if (DesignBlockMatrix.A[i] == NULL)
                error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 3");
#endif
            for (k = 0; k < DesignBlockMatrix.Block_nCol; k++)
                DesignBlockMatrix.A[i][j][k] = 0.0;
        }
    }

    // Allocate Memory of RHS
    DesignBlockMatrix.B = NULL;
    DesignBlockMatrix.B = (double **) malloc (DesignBlockMatrix.nCOL*sizeof(double*));
#ifdef DEBUG
    if (DesignBlockMatrix.B == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 4");
#endif
    for (i = 0; i < DesignBlockMatrix.nCOL; i++) {
        DesignBlockMatrix.B[i] = NULL;
        DesignBlockMatrix.B[i] = (double *) malloc (DesignBlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (DesignBlockMatrix.B[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 5");
#endif
        for (j = 0; j < DesignBlockMatrix.Block_nRow; j++)
            DesignBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    DesignBlockMatrix.X = NULL;
    DesignBlockMatrix.X = (double **) malloc (DesignBlockMatrix.nCOL*sizeof(double*));
#ifdef DEBUG
    if (DesignBlockMatrix.X == NULL)
        error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 6");
#endif
    for (i = 0; i < DesignBlockMatrix.nCOL; i++) {
        DesignBlockMatrix.X[i] = NULL;
        DesignBlockMatrix.X[i] = (double *) malloc (DesignBlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (DesignBlockMatrix.X[i] == NULL)
            error("Create_CRS_DesignBlockMatrix: %s\n", "Error Allocating Memory 7");
#endif
        for (j = 0; j < DesignBlockMatrix.Block_nRow; j++)
            DesignBlockMatrix.X[i][j] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Read_dXdBeta(const char* FileName) {
    FILE *fp;
    int i;
    double tempvar1, tempvar2;

    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Read_dXdBeta: Unable to Open dXdBeta File %s", FileName);

    for (i = 0; i < mesh.nnodes; i++)
        fscanf(fp, "%lf %lf %lf %lf", &tempvar1, &tempvar2, &dXdBeta[i], &dYdBeta[i]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Write_dXdBeta(const char* FileName) {
    FILE *fp;
    int i;
    double tempvar1, tempvar2;

    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dXdBeta: Unable to Open Mesh File %s", FileName);

    tempvar1 = 0.0;
    tempvar2 = 0.0;
    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n", tempvar1, tempvar2, dXdBeta[i], dYdBeta[i]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Write_dQdBeta(const char* FileName) {
    FILE *fp;
    int i;
    
    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dQdBeta: Unable to Open Mesh File %s", FileName);

    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n",
                dQdBeta[i][0], dQdBeta[i][1], dQdBeta[i][2], dQdBeta[i][3]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Write_dRdX_dXdBeta(const char* FileName) {
    FILE *fp;
    int i;

    if ((fp = fopen(FileName, "w")) == (FILE *) NULL)
        error("Write_dRdX_dXdBeta: Unable to Open Mesh File %s", FileName);

    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n",
                dRdX_dXdBeta[i][0], dRdX_dXdBeta[i][1], dRdX_dXdBeta[i][2], dRdX_dXdBeta[i][3]);
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_dRdX_dXdBeta() {
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
void Euler2D_Solver_Design::Compute_Boundary_dRdX_dXdBeta() {
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
void Euler2D_Solver_Design::Compute_dFplusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
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
void Euler2D_Solver_Design::Compute_dFminusdX_dXdBeta(double *Q, double uNx, double uNy, double Mag,
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
void Euler2D_Solver_Design::Verify_dRdX_dXdBeta() {
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
void Euler2D_Solver_Design::Compute_dQdBeta() {
    int i, j, k;
    int iNode;
    double lrms;

    // Update the CRS Design Matrix with solver dR/dQ and
    for (i = 0; i < DesignBlockMatrix.DIM; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++)
                DesignBlockMatrix.A[i][j][k] = BlockMatrix.A[i][j][k];
        }
    }

    // Update the CRS Design Matrix Vector B with dR/dBeta
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.Block_nRow; j++) {
            DesignBlockMatrix.B[iNode][j] = -dRdX_dXdBeta[iNode][j];
        }
    }

    // Solve for dQdBeta
    InnerNIteration = 0;
    Relaxation = 0.9;
    lrms = MC_Iterative_Block_Jacobi_CRS(InnerNIteration, Relaxation, DesignBlockMatrix);
//    lrms = MC_Iterative_Block_LU_Jacobi_CRS(InnerNIteration, 2, DesignBlockMatrix);

    printf("LRMS = %10.5e\n", lrms);

    // Get the Solution
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < DesignBlockMatrix.Block_nRow; j++)
            dQdBeta[iNode][j] = DesignBlockMatrix.X[iNode][j];
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Verify_dQdBeta() {
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
        Solver1.Get_Solver_Inputs("Euler_Solver1.inp");
        Solver1.Solver_Prepare();
        Solver1.Solve();
        Solver1.Solver_Finalize();
    };
    // Create Solver2 and Converge
    {
        Euler2D_Solver_VanLeer Solver2;
        Solver2.Get_Solver_Inputs("Euler_Solver2.inp");
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
void Euler2D_Solver_Design::Compute_dIdQ_dQdBeta() {
    int  i, iBEdge, iedge;
    int bn[2];
    double x[2], Nx;
    double y[2], Ny;
    double Var1, Var2;
    double Rho,   U,   V,   E,   P[2],   Pavg,   Cp;
    double Rho_b, U_b, V_b, E_b, P_b[2], Pavg_b, Cp_b;
    double Lift_b, Drag_b, CoeffLift_b, CoeffDrag_b;

    dIdQ_dQdBeta = 0.0;
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
        Nx    = +0.5*(y[1] - y[0]);
        Ny    = -0.5*(x[1] - x[0]);

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

        // Get Coefficent of Pressure and its derivates
        Var1 = (2.0/(Gamma*Ref_Mach*Ref_Mach));
        Cp   = Var1 * (Pavg/Ref_Pressure - 1.0);
        Cp_b = Var1 * Pavg_b/Ref_Pressure;
        
        // Compute Lift and its derivatives
        Lift_b += Cp_b*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny);

        // Compute Drag and its derivatives
        Drag_b += Cp_b*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny);
    }

    // Compute Coefficient of Lift and Drag and its derivates
    CoeffLift_b = Lift_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag_b = Drag_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute (dI/dQ)*(dQ/dBeta)
    dIdQ_dQdBeta = (CoeffLift - Ref_CoeffLift)*CoeffLift_b;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_dIdX_dXdBeta() {
    int  i, iBEdge, iedge;
    int bn[2];
    double x_b[2], Nx_b;
    double y_b[2], Ny_b;
    double Var1;
    double Rho,   U,   V,   E,   P[2],   Pavg,   Cp;
    double Lift_b, Drag_b, CoeffLift_b, CoeffDrag_b;

    dIdX_dXdBeta = 0.0;
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

        // Get the Coordinate Derivates
        x_b[0] = dXdBeta[bn[0]];
        y_b[0] = dYdBeta[bn[0]];
        x_b[1] = dXdBeta[bn[1]];
        y_b[1] = dYdBeta[bn[1]];

        // Compute the Normal Derivatives
        Nx_b  = +0.5*(y_b[1] - y_b[0]);
        Ny_b  = -0.5*(x_b[1] - x_b[0]);

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

        // Get Coefficent of Pressure and its derivates
        Var1 = (2.0/(Gamma*Ref_Mach*Ref_Mach));
        Cp   = Var1 * (Pavg/Ref_Pressure - 1.0);

        // Compute Lift derivatives => (dLift/dX)*(dX/dBeta)
        Lift_b += Cp*(-sin(Ref_Alpha)*Nx_b + cos(Ref_Alpha)*Ny_b);

        // Compute Drag derivatives => (dDrag/dX)*(dX/dBeta)
        Drag_b += Cp*( cos(Ref_Alpha)*Nx_b + sin(Ref_Alpha)*Ny_b);
    }

    // Compute Coefficient of Lift and Drag derivates
    CoeffLift_b = Lift_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag_b = Drag_b/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute (dI/dX)*(dX/dBeta)
    dIdX_dXdBeta = (CoeffLift - Ref_CoeffLift)*CoeffLift_b;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_dIdBeta() {
    dIdBeta = dIdQ_dQdBeta + dIdX_dXdBeta;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_Cost() {
    int  i, iBEdge, iedge;
    int bn[2];
    double x[2], Nx;
    double y[2], Ny;
    double Var1;
    double Rho,   U,   V,   E,   P[2],   Pavg,   Cp;
    
    Lift = 0.0;
    Drag = 0.0;

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
        Nx    = +0.5*(y[1] - y[0]);
        Ny    = -0.5*(x[1] - x[0]);

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

        // Get Coefficent of Pressure and its derivates
        Var1 = (2.0/(Gamma*Ref_Mach*Ref_Mach));
        Cp   = Var1 * (Pavg/Ref_Pressure - 1.0);

        // Compute Lift
        Lift   += Cp*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny);

        // Compute Drag
        Drag   += Cp*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny);
    }

    // Compute Coefficient of Lift and Drag and its derivates
    CoeffLift   = Lift/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag   = Drag/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute Cost Function
    I = 0.5*(CoeffLift - Ref_CoeffLift)*(CoeffLift - Ref_CoeffLift);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_dIdQ_bak() {
    int  i, j, iNode, iBEdge, iedge;
    int bn[2];
    double Nx, Ny;
    double Var1, Var2;
    double U, V;
    double Rho, P[2], E, Pavg, Cp;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], E_Q[4], P_Q[2][4], Pavg_Q[2][4], Cp_Q[2][4];

    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < 4; j++)
            dIdQ[iNode][j] = 0.0;
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

        // Get the Normal - Components
        Nx = edge[iedge].unx*edge[iedge].mag;
        Ny = edge[iedge].uny*edge[iedge].mag;

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

        // Get Coefficent of Pressure and its derivates
        Var1 = (2.0/(Gamma*Ref_Mach*Ref_Mach));
        Cp = Var1 * (Pavg/Ref_Pressure - 1.0);
        Var1 = Var1/Ref_Pressure;
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 4; j++) {
                Cp_Q[i][j] = Var1*Pavg_Q[i][j];
            }
        }

        // Compute Lift and Drag
        Lift += Cp*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny);
        Drag += Cp*( cos(Ref_Alpha)*Nx + sin(Ref_Alpha)*Ny);

        for (i = 0; i < 2; i++) {
            for (j = 0; j < 4; j++) {
                dIdQ[bn[i]][j] += Cp_Q[i][j]*(-sin(Ref_Alpha)*Nx + cos(Ref_Alpha)*Ny);
            }
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_Design::Compute_Cost_bak() {
    int i, n1, n2, iedge;
    double Pressure, Nx, Ny;
    double P1, P2;

    Lift = 0.0;
    Drag = 0.0;

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[i].bcType < 1000) || (boundaryEdge[i].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[i].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Get the Normal - Components
        Nx = edge[iedge].unx*edge[iedge].mag;
        Ny = edge[iedge].uny*edge[iedge].mag;

        // Pressure at Node1
        P1 = (Gamma - 1.0)*(node[n1].Q[3] - 0.5*node[n1].Q[0]*((node[n1].Q[1]/node[n1].Q[0])*(node[n1].Q[1]/node[n1].Q[0])
                + (node[n1].Q[2]/node[n1].Q[0])*(node[n1].Q[2]/node[n1].Q[0])));

        // Pressure at Node2
        P2 = (Gamma - 1.0)*(node[n2].Q[3] - 0.5*node[n2].Q[0]*((node[n2].Q[1]/node[n2].Q[0])*(node[n2].Q[1]/node[n2].Q[0])
                + (node[n2].Q[2]/node[n2].Q[0])*(node[n2].Q[2]/node[n2].Q[0])));

        // Get the Average Pressure
        Pressure = 0.5*(P1 + P2);

        // Compute Lift and Drag
        Lift += (-Pressure*sin(Ref_Alpha)*Nx + Pressure*cos(Ref_Alpha)*Ny);
        Drag += ( Pressure*cos(Ref_Alpha)*Nx + Pressure*sin(Ref_Alpha)*Ny);
    }
    CoeffLift = Lift/(0.5*Ref_Mach*Ref_Mach*Ref_Length);
    CoeffDrag = Drag/(0.5*Ref_Mach*Ref_Mach*Ref_Length);

    // Compute the value of Cost
    I = 0.5*(CoeffLift - Ref_CoeffLift)*(CoeffLift - Ref_CoeffLift);
}

