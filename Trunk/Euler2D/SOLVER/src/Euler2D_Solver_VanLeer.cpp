/* 
 * File:   Euler2D_Solver_VanLeer.cpp
 * Author: Ashish Gupta
 * 
 * Created on March 21, 2010, 4:31 PM
 */

#include <stdlib.h>
#include <iostream>
#include <limits.h>

#include "Utils.h"
#include "Euler2D_Solver_VanLeer.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_VanLeer::Euler2D_Solver_VanLeer() {
    printf("=============================================================================\n");
    printf("      Euler2D : Van Leer Flux Vector Splitting                               \n");
    printf("=============================================================================\n");
    
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Init() {
    Order                  = 2;
    NIteration             = 0;
    InnerNIteration        = 0;
    Relaxation             = 1.0;
    RMS[0]                 = DBL_MIN;
    RMS[1]                 = DBL_MIN;
    RMS[2]                 = DBL_MIN;
    RMS[3]                 = DBL_MIN;
    RMS_Res                = DBL_MIN;
    Gamma                  = 1.4;
    Ref_Mach               = 0.0;
    Ref_Alpha              = 0.0;
    Ref_Pressure           = 0.0;
    Ref_Length             = 0.0;
    CFL_Ramp               = 0;
    CFL_Min                = 0.0;
    CFL_Max                = 0.0;
    BlockMatrix.nCOL       = 0;
    BlockMatrix.nROW       = 0;
    BlockMatrix.Block_nCol = 0;
    BlockMatrix.Block_nRow = 0;
    BlockMatrix.DIM        = 0;
    BlockMatrix.A          = NULL;
    BlockMatrix.B          = NULL;
    BlockMatrix.IA         = NULL;
    BlockMatrix.IAU        = NULL;
    BlockMatrix.JA         = NULL;
    BlockMatrix.X          = NULL;
    DeltaT                 = NULL;
    Qx                     = NULL;
    Qy                     = NULL;
    FDIteration            = 0;
    FDNodeID               = -1;
    FDPertQ                = -1;
    FDEpsilon              = 1.0E-5;
    FDDR_DQ[0]             = 0.0;
    FDDR_DQ[1]             = 0.0;
    FDDR_DQ[2]             = 0.0;
    FDDR_DQ[3]             = 0.0;
    FDNode                 = NULL;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Solver_VanLeer::~Euler2D_Solver_VanLeer() {
    int i, j;
    
    if (BlockMatrix.A != NULL) {
        for (i = 0; i < BlockMatrix.DIM; i++) {
            if (BlockMatrix.A[i] != NULL) {
                for (j = 0; j < BlockMatrix.Block_nRow; j++) {
                    if (BlockMatrix.A[i][j] != NULL)
                        free(BlockMatrix.A[i][j]);
                }
                free(BlockMatrix.A[i]);
            }
        }
        free(BlockMatrix.A);
    }
    if (BlockMatrix.B != NULL) {
        for (i = 0; i < BlockMatrix.nCOL ; i++) {
            if (BlockMatrix.B[i] != NULL)
                free(BlockMatrix.B[i]);
        }
        free(BlockMatrix.B);
    }
    if (BlockMatrix.X != NULL) {
        for (i = 0; i < BlockMatrix.nCOL ; i++) {
            if (BlockMatrix.X[i] != NULL)
                free(BlockMatrix.X[i]);
        }
        free(BlockMatrix.X);
    }
    if (BlockMatrix.IA != NULL)
        delete[] BlockMatrix.IA;
    if (BlockMatrix.IAU != NULL)
        delete[] BlockMatrix.IAU;
    if (BlockMatrix.JA != NULL)
        delete[] BlockMatrix.JA;

    if (DeltaT != NULL)
        delete[] DeltaT;

    if (Qx != NULL) {
        for (i = 0; i < mesh.nnodes; i++) {
            if (Qx[i] != NULL)
                delete[] Qx[i];
        }
        delete[] Qx;
    }

    if (Qy != NULL) {
        for (i = 0; i < mesh.nnodes; i++) {
            if (Qy[i] != NULL)
                delete[] Qy[i];
        }
        delete[] Qy;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Solver_Prepare() {
    // Read Mesh File
    WKA_MeshReader(WKAMeshFileName.c_str());
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Solver_Finalize() {
    // Write SLK Mesh file
    SLK_MeshWriter(SLKMeshFileName.c_str());
    // Write VTK Solution File
    Write_VTK_Unstructured_File(VTKSolutionFileName.c_str());
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Get_Solver_Inputs(const char* FileName) {
    FILE *fp;
    int iret;
    char dumstring[257];
    char dump[257];
    char *cdum;

    // Open Input File
    if ((fp = fopen(FileName, "r")) == (FILE *) NULL)
        error("Unable to Open Input File %s", FileName);

    // Input Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%s\n", dump);
    WKAMeshFileName.append(dump);

    // Input Restart File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%s\n", dump);
    InputRestartFileName.append(dump);

    // Output Restart File Name
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%s\n", dump);
    OutputRestartFileName.append(dump);

    // SLK Mesh File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%s\n", dump);
    SLKMeshFileName.append(dump);

    // VTK Based Solution File
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%s\n", dump);
    VTKSolutionFileName.append(dump);
    
    // Get the reference Conditions
    // Gamma
    Gamma  = 1.4;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &Gamma);
    Ref_Pressure    = 1.0/Gamma;
    Ref_Length      = 1.0;
    
    // Mach Number
    // Input Ref Mach Number
    Ref_Mach        = 0.5;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &Ref_Mach);

    // Get Angle of Attack
    // Input Angle of Attack (deg)
    Ref_Alpha       = 0.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &Ref_Alpha);
    Ref_Alpha       = Ref_Alpha * M_PI / 180.0;

    // Get the order of solver
    // Solver Order (1/2)
    Order = 2;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &Order);
    
    if ((Order < 1) || (Order > 2))
        error("Get_Reference_Conditions: %s\n", "Invalid Solver Order");
    
    // Get CFL Numbers
    // Input CFL Min
    CFL_Min         = 1.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &CFL_Min);
    // Input CFL Max
    CFL_Max         = 20.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &CFL_Max);
    // Input CFL Ramp
    CFL_Ramp        = 100;
    iret = fscanf(fp, "\n");
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &CFL_Ramp);

    // Get the Inner and Outer loop Conditions
    // Input Outer Iterations (0 = Convergence)
    NIteration      = 1000;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &NIteration);
    if (NIteration == 0)
        NIteration = INT_MAX/10;
    // Input Linear Solver Iterations (0 = Convergence)
    InnerNIteration = 20;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &InnerNIteration);
    // Input Linear Solver Relaxation Inner (0.5 < r < 1.0)
    Relaxation      = 1.0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%lf\n", &Relaxation);

    // Restart File
    // Restart Capability (0=No, 1=Yes, 2=Create)
    restart = 0;
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &restart);

    // Finite Difference Jacobian
    FDNodeID = 0;
    // Finite Difference Jacobian (0=No, 1=Yes)
    iret = fscanf(fp, "\n");
    cdum = fgets(dumstring, 100, fp);
    iret = fscanf(fp, "%d\n", &FDNodeID);
    if (FDNodeID == 1) {
        // Finite Difference Node ID (-1 = All Nodes)
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 100, fp);
        iret = fscanf(fp, "%d\n", &FDNodeID);
        // Finite Difference Perturbation Q
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 100, fp);
        iret = fscanf(fp, "%d\n", &FDPertQ);
        // Finite Difference Epsilon
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 100, fp);
        iret = fscanf(fp, "%lf\n", &FDEpsilon);
        // Finite Difference Check Iteration
        iret = fscanf(fp, "\n");
        cdum = fgets(dumstring, 100, fp);
        iret = fscanf(fp, "%d\n", &FDIteration);
    } else {
        FDNodeID = -2;
    }

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
    printf("Restart             = %d\n",  restart);
    printf("FD Node ID          = %d\n",  FDNodeID);
    printf("FD Pertuabation Q   = %d\n",  FDPertQ);
    printf("FD Epsilon          = %lf\n", FDEpsilon);
    printf("FD Check Iteration  = %d\n",  FDIteration);
    printf("-----------------------------------------------------------------------------\n");
    // Close Input File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Solve() {
    int iter, AddTime;
    double max_rms, lrms;
    
    // Extract Connectivity
    WKA_ExtractCells();
    Tag_Boundary_Nodes();
    
    // Compute the Expensive Geometric Properties
    Compute_Geometric_Properties();
    
    // Initialize the solution
    Initialize_Solution();

    // Read Restart File
    if (restart == 1) {
        Read_RestartFile(InputRestartFileName.c_str());
//        // Purturb Solution Q
//        for (int iNode = 0; iNode < mesh.nnodes; iNode++) {
//            for (int j = 0; j < 4; j++)
//                node[iNode].Q[j] = node[iNode].Q[j] + 0.000001;
//        }
    }

    // Create CRS Matrix
    Create_CRS_BlockMatrix();

    info("Euler2D Computation Starts");
    printf("-----------------------------------------------------------------------------\n");
    printf(" Iter        LRMS     RMS_RHO    RMS_RHOU    RMS_RHOV       RMS_E     RMS_RES\n");
    printf("-----------------------------------------------------------------------------\n");

    for (iter = 0; iter < NIteration; iter++) {
        // Compute Local Time Stepping
        Compute_DeltaTime(iter);

        // Compute Residials
        Compute_Gauss_Gradient();
        Compute_Residual();
        Compute_Boundary_Residual();
        
        // Compute and Fill CRS Matrix
        AddTime = 1;
        Compute_CRS_BlockMatrix(AddTime);
        Compute_Boundary_CRS_BlockMatrix();

        // Compute Finite Difference Jacobian
        Compute_FD_Jacobian(iter);

#ifdef DEBUG
        int inode, irow, ibrow, ibcol;
        for (inode = 0; inode < mesh.nnodes; inode++) {
            printf("------------------------------------------------------------\n");
            printf("Diagonal = %d \t", BlockMatrix.JA[BlockMatrix.IAU[inode]]);
            printf("Row Nodes = { ");
            for (irow = BlockMatrix.IA[inode]; irow < BlockMatrix.IA[inode+1]; irow++)
                printf("%5d ",BlockMatrix.JA[irow]);
            printf("}\n");

            for (irow = BlockMatrix.IA[inode]; irow < BlockMatrix.IA[inode+1]; irow++) {
                printf("Node = %d\n", BlockMatrix.JA[irow]);
                for (ibrow = 0; ibrow < 4; ibrow++) {
                    for (ibcol = 0; ibcol < 4; ibcol++)
                        printf("%22.15e ", BlockMatrix.A[irow][ibrow][ibcol]);
                    printf("\n");
                }
                printf("\n");
            }
        }
        printf("------------------------------------------------------------\n");
#endif
        // Solve for Solution
        lrms = MC_Iterative_Block_Jacobi_CRS(InnerNIteration, Relaxation, BlockMatrix);

        // Update the Solution
        Update_Solution();

        // Compute RMS
        max_rms = Compute_RMS();

        printf("%5d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n", iter+1, lrms, RMS[0], RMS[1], RMS[2], RMS[3], RMS_Res);
        
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
    printf("-----------------------------------------------------------------------------\n");
    if (restart)
        Write_RestartFile(OutputRestartFileName.c_str());
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Initialize_Solution() {
    int i, iNode;
    double Q[4];

    // Allocate the Memory to store Local Time Stepping
    DeltaT = new double[mesh.nnodes];
#ifdef DEBUG
    if (DeltaT == NULL)
        error("Initialize_Solution: %s\n", "Error Allocating Memory 1");
#endif

    if (Order == 2) {
        // Allocate Memory to store gradient of conserved variables
        Qx = new double*[mesh.nnodes];
        Qy = new double*[mesh.nnodes];
#ifdef DEBUG
        if ((Qx == NULL) || (Qy == NULL))
            error("Initialize_Solution: %s\n", "Error Allocating Memory 2");
#endif
        for (iNode = 0; iNode < mesh.nnodes; iNode++) {
            Qx[iNode] = NULL;
            Qy[iNode] = NULL;
            Qx[iNode] = new double[4];
            Qy[iNode] = new double[4];
#ifdef DEBUG
            if ((Qx[iNode] == NULL) || (Qy[iNode] == NULL))
                error("Initialize_Solution: %s\n", "Error Allocating Memory 3");
#endif
            // Initialize
            for (i = 0; i < 4; i++) {
                Qx[iNode][i] = 0.0;
                Qy[iNode][i] = 0.0;
            }
        }
    }
    
    // Calculte Q
    Q[0] = 1.0;
    Q[1] = Q[0]*Ref_Mach*cos(Ref_Alpha);
    Q[2] = Q[0]*Ref_Mach*sin(Ref_Alpha);
    Q[3] = (Ref_Pressure/(Gamma - 1.0)) + 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0]));
    
    // Initialize the Solution Field
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        DeltaT[iNode] = 0.0;
        for (i = 0; i < 4; i++) {
            node[iNode].Q[i] = Q[i];
            node[iNode].Resi[i] = 0.0;
        }
    }

    // Allocate the Finite Difference Jacobian
    FDNode = (FD_NODE *) malloc(sizeof(FD_NODE)*mesh.nnodes);
#ifdef DEBUG
    if (FDNode == NULL)
        error("Initialize_Solution: %s\n", "Error Allocating Memory 4");
#endif
    // Initialize the Finite Difference Jacobian Variables
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++) {
            FDNode[iNode].Resi_Old[i] = 0.0;
            FDNode[iNode].Resi_New[i] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Gauss_Gradient() {
    int i, n1, n2, n3, iNode, iCell;
    double Q1, Q2, Q3, Q4;
    double tx, ty;
    double Partx[4], Party[4];
    
    if (Order != 2)
        return;

    // Initialize
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++) {
            Qx[iNode][i] = 0.0;
            Qy[iNode][i] = 0.0;
        }
    }

    for (iCell = 0; iCell  < mesh.inside; iCell++) {
        n1 = cell[iCell].node1;
        n2 = cell[iCell].node2;
        n3 = cell[iCell].node3;

        // Node 1
        Q1 = 0.5*(node[n2].Q[0] + node[n3].Q[0]);
        Q2 = 0.5*(node[n2].Q[1] + node[n3].Q[1]);
        Q3 = 0.5*(node[n2].Q[2] + node[n3].Q[2]);
        Q4 = 0.5*(node[n2].Q[3] + node[n3].Q[3]);

        tx = cell[iCell].unx23*cell[iCell].mag23;
        ty = cell[iCell].uny23*cell[iCell].mag23;

        Partx[0] = Q1*tx;
        Partx[1] = Q2*tx;
        Partx[2] = Q3*tx;
        Partx[3] = Q4*tx;
        
        Party[0] = Q1*ty;
        Party[1] = Q2*ty;
        Party[2] = Q3*ty;
        Party[3] = Q4*ty;
        
        // Node 2
        Q1 = 0.5*(node[n3].Q[0] + node[n1].Q[0]);
        Q2 = 0.5*(node[n3].Q[1] + node[n1].Q[1]);
        Q3 = 0.5*(node[n3].Q[2] + node[n1].Q[2]);
        Q4 = 0.5*(node[n3].Q[3] + node[n1].Q[3]);

        tx = cell[iCell].unx31*cell[iCell].mag31;
        ty = cell[iCell].uny31*cell[iCell].mag31;

        Partx[0] += Q1*tx;
        Partx[1] += Q2*tx;
        Partx[2] += Q3*tx;
        Partx[3] += Q4*tx;

        Party[0] += Q1*ty;
        Party[1] += Q2*ty;
        Party[2] += Q3*ty;
        Party[3] += Q4*ty;

        // Node 3
        Q1 = 0.5*(node[n1].Q[0] + node[n2].Q[0]);
        Q2 = 0.5*(node[n1].Q[1] + node[n2].Q[1]);
        Q3 = 0.5*(node[n1].Q[2] + node[n2].Q[2]);
        Q4 = 0.5*(node[n1].Q[3] + node[n2].Q[3]);

        tx = cell[iCell].unx12*cell[iCell].mag12;
        ty = cell[iCell].uny12*cell[iCell].mag12;

        Partx[0] += Q1*tx;
        Partx[1] += Q2*tx;
        Partx[2] += Q3*tx;
        Partx[3] += Q4*tx;

        Party[0] += Q1*ty;
        Party[1] += Q2*ty;
        Party[2] += Q3*ty;
        Party[3] += Q4*ty;

        // Finally the Gradients
        // Node 1
        Qx[n1][0] += Partx[0];
        Qx[n1][1] += Partx[1];
        Qx[n1][2] += Partx[2];
        Qx[n1][3] += Partx[3];

        Qy[n1][0] += Party[0];
        Qy[n1][1] += Party[1];
        Qy[n1][2] += Party[2];
        Qy[n1][3] += Party[3];

        // Node 2
        Qx[n2][0] += Partx[0];
        Qx[n2][1] += Partx[1];
        Qx[n2][2] += Partx[2];
        Qx[n2][3] += Partx[3];

        Qy[n2][0] += Party[0];
        Qy[n2][1] += Party[1];
        Qy[n2][2] += Party[2];
        Qy[n2][3] += Party[3];

        // Node 3
        Qx[n3][0] += Partx[0];
        Qx[n3][1] += Partx[1];
        Qx[n3][2] += Partx[2];
        Qx[n3][3] += Partx[3];

        Qy[n3][0] += Party[0];
        Qy[n3][1] += Party[1];
        Qy[n3][2] += Party[2];
        Qy[n3][3] += Party[3];
    }
    
    // Finally Division by Area
    // Median Dual Area is 1/3 times the area enclosed by connected triangles
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++) {
            Qx[iNode][i] /= 3.0*node[iNode].area;
            Qy[iNode][i] /= 3.0*node[iNode].area;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Residual() {
    int i, n1, n2, n3, iNode, iCell;
    double Flux[4], Fplus[4], Fminus[4];
    
    // Initialize the Residual
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++)
            node[iNode].Resi[i] = 0.0;
    }

    // First Order
    if (Order == 1) {
        // Compute Residual of All Nodes
        for (iCell = 0; iCell < mesh.inside; iCell++) {
            n1 = cell[iCell].node1;
            n2 = cell[iCell].node2;
            n3 = cell[iCell].node3;

            // Get the Flux Across the Median Dual Line
            // Edge 1-2
            Compute_Fplus(node[n1].Q, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fplus);
            Compute_Fminus(node[n2].Q, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n1].Resi[i] += Flux[i];
                node[n2].Resi[i] -= Flux[i];
            }

            // Edge 2-3
            Compute_Fplus(node[n2].Q, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fplus);
            Compute_Fminus(node[n3].Q, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n2].Resi[i] += Flux[i];
                node[n3].Resi[i] -= Flux[i];
            }

            // Edge 3-1
            Compute_Fplus(node[n3].Q, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fplus);
            Compute_Fminus(node[n1].Q, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n3].Resi[i] += Flux[i];
                node[n1].Resi[i] -= Flux[i];
            }
        }
    }

    // Second Order
    if (Order == 2) {
        double QL[4], QR[4];
        double xf, yf;
        // Compute Residual of All Nodes
        for (iCell = 0; iCell < mesh.inside; iCell++) {
            n1 = cell[iCell].node1;
            n2 = cell[iCell].node2;
            n3 = cell[iCell].node3;

            // Get the Flux Across the Median Dual Line
            // Edge 1-2
            // Get QLeft and QRight
            xf = 0.5*cell[iCell].xc + 0.25*(node[n1].x + node[n2].x);
            yf = 0.5*cell[iCell].yc + 0.25*(node[n1].y + node[n2].y);
            for (i = 0; i < 4; i++) {
                QL[i] = node[n1].Q[i] + Qx[n1][i]*(xf - node[n1].x) + Qy[n1][i]*(yf - node[n1].y);
                QR[i] = node[n2].Q[i] + Qx[n2][i]*(xf - node[n2].x) + Qy[n2][i]*(yf - node[n2].y);
            }
            
            Compute_Fplus(QL, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fplus);
            Compute_Fminus(QR, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n1].Resi[i] += Flux[i];
                node[n2].Resi[i] -= Flux[i];
            }

            // Edge 2-3
            // Get QLeft and QRight
            xf = 0.5*cell[iCell].xc + 0.25*(node[n2].x + node[n3].x);
            yf = 0.5*cell[iCell].yc + 0.25*(node[n2].y + node[n3].y);
            for (i = 0; i < 4; i++) {
                QL[i] = node[n2].Q[i] + Qx[n2][i]*(xf - node[n2].x) + Qy[n2][i]*(yf - node[n2].y);
                QR[i] = node[n3].Q[i] + Qx[n3][i]*(xf - node[n3].x) + Qy[n3][i]*(yf - node[n3].y);
            }

            Compute_Fplus(QL, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fplus);
            Compute_Fminus(QR, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n2].Resi[i] += Flux[i];
                node[n3].Resi[i] -= Flux[i];
            }

            // Edge 3-1
            // Get QLeft and QRight
            xf = 0.5*cell[iCell].xc + 0.25*(node[n3].x + node[n3].x);
            yf = 0.5*cell[iCell].yc + 0.25*(node[n1].y + node[n1].y);
            for (i = 0; i < 4; i++) {
                QL[i] = node[n3].Q[i] + Qx[n3][i]*(xf - node[n3].x) + Qy[n3][i]*(yf - node[n3].y);
                QR[i] = node[n1].Q[i] + Qx[n1][i]*(xf - node[n1].x) + Qy[n1][i]*(yf - node[n1].y);
            }

            Compute_Fplus(QL, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fplus);
            Compute_Fminus(QR, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fminus);
            for (i = 0; i < 4; i++) {
                Flux[i] = Fplus[i] + Fminus[i];
                node[n3].Resi[i] += Flux[i];
                node[n1].Resi[i] -= Flux[i];
            }
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Boundary_Residual() {
    int i, j, n1, n2, iedge;
    double Pressure, Nx, Ny;

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Freestream Boundary: BCTag = 3
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 3: 3000-3999 : Free Stream
        if ((boundaryEdge[i].bcType < 3000) || (boundaryEdge[i].bcType > 3999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[i].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Residual
        for (j = 0; j < 4; j++) {
            node[n1].Resi[j] = 0.0;
            node[n2].Resi[j] = 0.0;
        }
    }
    
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
        Nx = 0.5*edge[iedge].unx*edge[iedge].mag;
        Ny = 0.5*edge[iedge].uny*edge[iedge].mag;

        // Residual at Node1
        Pressure = (Gamma - 1.0)*(node[n1].Q[3] - 0.5*node[n1].Q[0]*((node[n1].Q[1]/node[n1].Q[0])*(node[n1].Q[1]/node[n1].Q[0])
                + (node[n1].Q[2]/node[n1].Q[0])*(node[n1].Q[2]/node[n1].Q[0])));
        node[n1].Resi[1] += Pressure*Nx;
        node[n1].Resi[2] += Pressure*Ny;

        // Residual at Node2
        Pressure = (Gamma - 1.0)*(node[n2].Q[3] - 0.5*node[n2].Q[0]*((node[n2].Q[1]/node[n2].Q[0])*(node[n2].Q[1]/node[n2].Q[0])
                + (node[n2].Q[2]/node[n2].Q[0])*(node[n2].Q[2]/node[n2].Q[0])));
        node[n2].Resi[1] += Pressure*Nx;
        node[n2].Resi[2] += Pressure*Ny;
    }

    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 2
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 2: 2000-2999
        if ((boundaryEdge[i].bcType < 2000) || (boundaryEdge[i].bcType > 2999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[i].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Get the Normal - Components
        Nx = 0.5*edge[iedge].unx*edge[iedge].mag;
        Ny = 0.5*edge[iedge].uny*edge[iedge].mag;

        // Residual at Node1
        Pressure = (Gamma - 1.0)*(node[n1].Q[3] - 0.5*node[n1].Q[0]*((node[n1].Q[1]/node[n1].Q[0])*(node[n1].Q[1]/node[n1].Q[0])
                + (node[n1].Q[2]/node[n1].Q[0])*(node[n1].Q[2]/node[n1].Q[0])));
        node[n1].Resi[1] += Pressure*Nx;
        node[n1].Resi[2] += Pressure*Ny;

        // Residual at Node2
        Pressure = (Gamma - 1.0)*(node[n2].Q[3] - 0.5*node[n2].Q[0]*((node[n2].Q[1]/node[n2].Q[0])*(node[n2].Q[1]/node[n2].Q[0])
                + (node[n2].Q[2]/node[n2].Q[0])*(node[n2].Q[2]/node[n2].Q[0])));
        node[n2].Resi[1] += Pressure*Nx;
        node[n2].Resi[2] += Pressure*Ny;
    }

//    // Calculate the Resi by Closing the Control Volume for Boundary Nodes
//    // Freestream Boundary: BCTag = 3
//    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
//    // Hence leads to CCW oriented boundary triangles
//    for (i = 0; i < mesh.nbedges; i++) {
//        // BCType 3: 3000-3999 : Free Stream
//        if ((boundaryEdge[i].bcType < 3000) || (boundaryEdge[i].bcType > 3999))
//            continue;
//
//        // Get The Nodes of Edge
//        iedge = boundaryEdge[i].edgeNumber;
//        n1 = edge[iedge].node1;
//        n2 = edge[iedge].node2;
//
//        // Residual
//        for (j = 0; j < 4; j++) {
//            node[n1].Resi[j] = 0.0;
//            node[n2].Resi[j] = 0.0;
//        }
//    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_FD_Jacobian(int Iteration) {
    int i, j, iNode, iNode_start, iNode_end;
    double DR_DQ[4];
    double **Data = NULL;

    if (Order != 1)
        return;
    
    if (Iteration != FDIteration)
        return;

    // No FD Validation computation
    if (FDNodeID <= -2 || FDNodeID >= mesh.nnodes)
        return;
    
    // Store the Old Residual
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (i = 0; i < 4; i++)
            FDNode[iNode].Resi_Old[i] = node[iNode].Resi[i];
    }

    if (FDNodeID < 0) {
        iNode_start = 0;
        iNode_end   = mesh.nnodes;
        Data = (double **) malloc(4*sizeof(double*));
        for (i = 0; i < 4; i++) {
            Data[i] = (double *) malloc(mesh.nnodes*sizeof(double));
            for (j = 0; j < mesh.nnodes; j++)
                Data[i][j] = 0.0;
        }
    } else {
        iNode_start = FDNodeID;
        iNode_end   = FDNodeID + 1;
    }
    
    for (iNode = iNode_start; iNode < iNode_end; iNode++) {
        // Now Purturb the Desired Node Particular Q
        node[iNode].Q[FDPertQ] += FDEpsilon;

        // Now Compute the New Residual
        Compute_Residual();
        Compute_Boundary_Residual();

        // Now Store the New Residual
        for (i = 0; i < mesh.nnodes; i++) {
            for (j = 0; j < 4; j++)
                FDNode[i].Resi_New[j] = node[i].Resi[j];
        }

        // Compute the Finite Difference DR_DQ
        for (i = 0; i < 4; i++)
            FDDR_DQ[i] = (FDNode[iNode].Resi_New[i] - FDNode[iNode].Resi_Old[i])/FDEpsilon;

        // Get the Analytical DR_DQ
        // Get the diagonal location
        int idgn = BlockMatrix.IAU[iNode];
        for (i = 0; i < 4; i++)
            DR_DQ[i] = BlockMatrix.A[idgn][i][FDPertQ];
        // Subtract I/DT
        DR_DQ[FDPertQ] -= node[iNode].area/DeltaT[iNode];

        if (FDNodeID < 0) {
            for (i = 0; i < 4; i++)
                Data[i][iNode] = fabs(FDDR_DQ[i] - DR_DQ[i]);
        } else {
            printf("Finite Difference:\n %15.9e %15.9e %15.9e %15.9e \n",
                    FDDR_DQ[0], FDDR_DQ[1], FDDR_DQ[2], FDDR_DQ[3]);
            printf("Analytic Difference:\n %15.9e %15.9e %15.9e %15.9e \n",
                    DR_DQ[0], DR_DQ[1], DR_DQ[2], DR_DQ[3]);
        }
        
        // Restore back the Old Residual and Desired Node Particular Q
        for (i = 0; i < mesh.nnodes; i++) {
            for (j = 0; j < 4; j++)
                node[i].Resi[j] = FDNode[i].Resi_Old[j];
        }
        node[iNode].Q[FDPertQ] -= FDEpsilon;
    }

    if (FDNodeID < 0) {
        Write_VTK_Unstructured_DebugFile("DRdQ_Valid1.vtu", Data[0]);
        Write_VTK_Unstructured_DebugFile("DRdQ_Valid2.vtu", Data[1]);
        Write_VTK_Unstructured_DebugFile("DRdQ_Valid3.vtu", Data[2]);
        Write_VTK_Unstructured_DebugFile("DRdQ_Valid4.vtu", Data[3]);
        for (i = 0; i < 4; i++)
            free(Data[i]);
        free(Data);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_DeltaTime(int Iteration) {
    int i, n1, n2, n3, iNode, iCell, iEdge;
    double U1, V1, Rho1, C1, P1, E1;
    double U2, V2, Rho2, C2, P2, E2;
    double U3, V3, Rho3, C3, P3, E3;
    double U, V, C, Vn, CFL;

    // Initialize the DeltaT
    for (iNode = 0; iNode < mesh.nnodes; iNode++)
        DeltaT[iNode] = 0.0;

    // Compute Local Time Stepping of All Nodes
    for (iCell = 0; iCell < mesh.inside; iCell++) {
        n1 = cell[iCell].node1;
        n2 = cell[iCell].node2;
        n3 = cell[iCell].node3;

        Rho1 = node[n1].Q[0];
        U1   = node[n1].Q[1]/Rho1;
        V1   = node[n1].Q[2]/Rho1;
        E1   = node[n1].Q[3];
        P1   = (Gamma - 1.0)*(E1 - 0.5*Rho1*(U1*U1 + V1*V1));
        C1   = sqrt(Gamma*P1/Rho1);

        Rho2 = node[n2].Q[0];
        U2   = node[n2].Q[1]/Rho2;
        V2   = node[n2].Q[2]/Rho2;
        E2   = node[n2].Q[3];
        P2   = (Gamma - 1.0)*(E2 - 0.5*Rho2*(U2*U2 + V2*V2));
        C2   = sqrt(Gamma*P2/Rho2);

        Rho3 = node[n3].Q[0];
        U3   = node[n3].Q[1]/Rho3;
        V3   = node[n3].Q[2]/Rho3;
        E3   = node[n3].Q[3];
        P3   = (Gamma - 1.0)*(E3 - 0.5*Rho3*(U3*U3 + V3*V3));
        C3   = sqrt(Gamma*P3/Rho3);

        // Get the Flux Across the Median Dual Line
        // Edge 1-2
        U = 0.5*(U1+U2);
        V = 0.5*(V1+V2);
        C = 0.5*(C1+C2);
        Vn = (fabs(U*cell[iCell].unxc12 + V*cell[iCell].unyc12) + C)*cell[iCell].magc12;

        DeltaT[n1] += Vn;
        DeltaT[n2] += Vn;

        // Edge 2-3
        U = 0.5*(U2+U3);
        V = 0.5*(V2+V3);
        C = 0.5*(C2+C3);
        Vn = (fabs(U*cell[iCell].unxc23 + V*cell[iCell].unyc23) + C)*cell[iCell].magc23;

        DeltaT[n2] += Vn;
        DeltaT[n3] += Vn;

        // Edge 3-1
        U = 0.5*(U3+U1);
        V = 0.5*(V3+V1);
        C = 0.5*(C3+C1);
        Vn = (fabs(U*cell[iCell].unxc31 + V*cell[iCell].unyc31) + C)*cell[iCell].magc31;
        
        DeltaT[n3] += Vn;
        DeltaT[n1] += Vn;
    }

    // Calculate the Local Time Stepping by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 1
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[i].bcType < 1000) || (boundaryEdge[i].bcType > 1999))
            continue;

        // Get The Nodes of Edge
        iEdge = boundaryEdge[i].edgeNumber;
        n1 = edge[iEdge].node1;
        n2 = edge[iEdge].node2;

        Rho1 = node[n1].Q[0];
        U1   = node[n1].Q[1]/Rho1;
        V1   = node[n1].Q[2]/Rho1;
        E1   = node[n1].Q[3];
        P1   = (Gamma - 1.0)*(E1 - 0.5*Rho1*(U1*U1 + V1*V1));
        C1   = sqrt(Gamma*P1/Rho1);

        Rho2 = node[n2].Q[0];
        U2   = node[n2].Q[1]/Rho2;
        V2   = node[n2].Q[2]/Rho2;
        E2   = node[n2].Q[3];
        P2   = (Gamma - 1.0)*(E2 - 0.5*Rho2*(U2*U2 + V2*V2));
        C2   = sqrt(Gamma*P2/Rho2);

        // Node1
        U = (5.0/6.0)*U1+ (1.0/6.0)*U2;
        V = (5.0/6.0)*V1+ (1.0/6.0)*V2;
        C = (5.0/6.0)*C1+ (1.0/6.0)*C2;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;
        
        DeltaT[n1] += Vn;
        
        // Node2
        U = (5.0/6.0)*U2+ (1.0/6.0)*U1;
        V = (5.0/6.0)*V2+ (1.0/6.0)*V1;
        C = (5.0/6.0)*C2+ (1.0/6.0)*C1;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;
        
        DeltaT[n2] += Vn;
    }

    // Calculate the Local Time Stepping by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 2
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 2: 2000-2999
        if ((boundaryEdge[i].bcType < 2000) || (boundaryEdge[i].bcType > 2999))
            continue;

        // Get The Nodes of Edge
        iEdge = boundaryEdge[i].edgeNumber;
        n1 = edge[iEdge].node1;
        n2 = edge[iEdge].node2;
        
        Rho1 = node[n1].Q[0];
        U1   = node[n1].Q[1]/Rho1;
        V1   = node[n1].Q[2]/Rho1;
        E1   = node[n1].Q[3];
        P1   = (Gamma - 1.0)*(E1 - 0.5*Rho1*(U1*U1 + V1*V1));
        C1   = sqrt(Gamma*P1/Rho1);

        Rho2 = node[n2].Q[0];
        U2   = node[n2].Q[1]/Rho2;
        V2   = node[n2].Q[2]/Rho2;
        E2   = node[n2].Q[3];
        P2   = (Gamma - 1.0)*(E2 - 0.5*Rho2*(U2*U2 + V2*V2));
        C2   = sqrt(Gamma*P2/Rho2);

        // Node1
        U = (5.0/6.0)*U1+ (1.0/6.0)*U2;
        V = (5.0/6.0)*V1+ (1.0/6.0)*V2;
        C = (5.0/6.0)*C1+ (1.0/6.0)*C2;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;

        DeltaT[n1] += Vn;

        // Node2
        U = (5.0/6.0)*U2+ (1.0/6.0)*U1;
        V = (5.0/6.0)*V2+ (1.0/6.0)*V1;
        C = (5.0/6.0)*C2+ (1.0/6.0)*C1;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;

        DeltaT[n2] += Vn;
    }

    // Calculate the Local Time Stepping by Closing the Control Volume for Boundary Nodes
    // Freestream Boundary: BCTag = 3
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 3: 3000-3999 : Free Stream
        if ((boundaryEdge[i].bcType < 3000) || (boundaryEdge[i].bcType > 3999))
            continue;

        // Get The Nodes of Edge
        iEdge = boundaryEdge[i].edgeNumber;
        n1 = edge[iEdge].node1;
        n2 = edge[iEdge].node2;

        Rho1 = node[n1].Q[0];
        U1   = node[n1].Q[1]/Rho1;
        V1   = node[n1].Q[2]/Rho1;
        E1   = node[n1].Q[3];
        P1   = (Gamma - 1.0)*(E1 - 0.5*Rho1*(U1*U1 + V1*V1));
        C1   = sqrt(Gamma*P1/Rho1);

        Rho2 = node[n2].Q[0];
        U2   = node[n2].Q[1]/Rho2;
        V2   = node[n2].Q[2]/Rho2;
        E2   = node[n2].Q[3];
        P2   = (Gamma - 1.0)*(E2 - 0.5*Rho2*(U2*U2 + V2*V2));
        C2   = sqrt(Gamma*P2/Rho2);

        // Node1
        U = (5.0/6.0)*U1+ (1.0/6.0)*U2;
        V = (5.0/6.0)*V1+ (1.0/6.0)*V2;
        C = (5.0/6.0)*C1+ (1.0/6.0)*C2;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;

        DeltaT[n1] += Vn;
        
        // Node2
        U = (5.0/6.0)*U2+ (1.0/6.0)*U1;
        V = (5.0/6.0)*V2+ (1.0/6.0)*V1;
        C = (5.0/6.0)*C2+ (1.0/6.0)*C1;
        Vn = 0.5*(fabs(U*edge[iEdge].unx + V*edge[iEdge].uny) + C)*edge[iEdge].mag;

        DeltaT[n2] += Vn;
    }

    // Finally Compute the Local Time Stepping
    if (CFL_Ramp > 1) {
        if (Iteration < CFL_Ramp)
            CFL = CFL_Min + (CFL_Max - CFL_Min)*(((double)Iteration)/((double)(CFL_Ramp-1)));
        else
            CFL = CFL_Max;
    } else
        CFL = CFL_Max;
    
    for (iNode = 0; iNode < mesh.nnodes; iNode++)
        DeltaT[iNode] = CFL*node[iNode].area/DeltaT[iNode];
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Create_CRS_BlockMatrix() {
    int i, j, jstart, jend, k, ksave;
    int node1, node2, index, index1, index2;
    int min, minsave;
    int *degree = NULL;

    BlockMatrix.nROW = mesh.nnodes;
    BlockMatrix.nCOL = mesh.nnodes;
    BlockMatrix.Block_nRow = 4;
    BlockMatrix.Block_nCol = 4;
    
    // Allocate Memory to Find Node2Node Connectivity Degree
    degree = new int[mesh.nnodes];
#ifdef DEBUG
    if (degree == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < mesh.nnodes; i++)
        degree[i] = 0;

    for (i = 0; i < mesh.nedges; i++) {
        degree[edge[i].node1] += 1;
        degree[edge[i].node2] += 1;
    }
    
    // Allocate Memory of IA Array to Store Start and End Location of Row
    BlockMatrix.IA = NULL;
    BlockMatrix.IA = new int[mesh.nnodes+1];
#ifdef DEBUG
    if (BlockMatrix.IA == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 2");
#endif

    // Start Filling the Row Location and
    // Get the No of Non Zero Entries for BlockMatrix
    BlockMatrix.DIM = 0;
    BlockMatrix.IA[0] = 0;
    for (i = 0; i < mesh.nnodes; i++) {
        BlockMatrix.DIM += degree[i] + 1;
        BlockMatrix.IA[i+1] = BlockMatrix.IA[i] + degree[i] + 1;
    }
    
    // Allocate Memory to IAU Array to store Diagonal Location
    BlockMatrix.IAU = NULL;
    BlockMatrix.IAU = new int[mesh.nnodes];
#ifdef DEBUG
    if (BlockMatrix.IAU == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 3");
#endif

    // Allocate Memory to JA Array to store location of Non Zero Entries
    BlockMatrix.JA = NULL;
    BlockMatrix.JA = new int[BlockMatrix.DIM];
#ifdef DEBUG
    if (BlockMatrix.JA == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 4");
#endif

    for (i = 0; i < mesh.nnodes; i++) {
        index = BlockMatrix.IA[i];
        BlockMatrix.JA[index] = i;
        degree[i] = 1;
    }
    
    for (i = 0; i < mesh.nedges; i++) {
        node1 = edge[i].node1;
        node2 = edge[i].node2;

        index1 = BlockMatrix.IA[node1] + degree[node1];
        index2 = BlockMatrix.IA[node2] + degree[node2];
        
        degree[node1] += 1;
        degree[node2] += 1;

        BlockMatrix.JA[index1] = node2;
        BlockMatrix.JA[index2] = node1;
    }

    /* Now Sort JA and Find IAU */
    for (i = 0; i < mesh.nnodes; i++) {
        jstart = BlockMatrix.IA[i];
        jend   = BlockMatrix.IA[i + 1];
        for (j = jstart; j < jend; j++) {
            min = BlockMatrix.JA[j];
            minsave = BlockMatrix.JA[j];
            ksave = j;
            for (k = j + 1; k < jend; k++) {
                if (BlockMatrix.JA[k] < min) {
                    min = BlockMatrix.JA[k];
                    ksave = k;
                }
            }
            BlockMatrix.JA[j] = min;
            BlockMatrix.JA[ksave] = minsave;
            if (BlockMatrix.JA[j] == i)
                BlockMatrix.IAU[i] = j;
        }
    }

    // Allocate Memory for CRS Matrix
    BlockMatrix.A = NULL;
    BlockMatrix.A = (double ***) malloc (BlockMatrix.DIM*sizeof(double**));
#ifdef DEBUG
    if (BlockMatrix.A == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 5");
#endif
    for (i = 0; i < BlockMatrix.DIM; i++) {
        BlockMatrix.A[i] = NULL;
        BlockMatrix.A[i] = (double **) malloc (BlockMatrix.Block_nRow*sizeof(double*));
#ifdef DEBUG
        if (BlockMatrix.A[i] == NULL)
            error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 6");
#endif
        for (j = 0; j < BlockMatrix.Block_nRow; j++) {
            BlockMatrix.A[i][j] = NULL;
            BlockMatrix.A[i][j] = (double *) malloc (BlockMatrix.Block_nCol*sizeof(double));
#ifdef DEBUG
            if (BlockMatrix.A[i] == NULL)
                error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 7");
#endif
            for (k = 0; k < BlockMatrix.Block_nCol; k++)
                BlockMatrix.A[i][j][k] = 0.0;
        }
    }

    // Allocate Memory of RHS
    BlockMatrix.B = NULL;
    BlockMatrix.B = (double **) malloc (BlockMatrix.nCOL*sizeof(double*));
#ifdef DEBUG
    if (BlockMatrix.B == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 8");
#endif
    for (i = 0; i < BlockMatrix.nCOL; i++) {
        BlockMatrix.B[i] = NULL;
        BlockMatrix.B[i] = (double *) malloc (BlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (BlockMatrix.B[i] == NULL)
            error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 9");
#endif
        for (j = 0; j < BlockMatrix.Block_nRow; j++)
            BlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    BlockMatrix.X = NULL;
    BlockMatrix.X = (double **) malloc (BlockMatrix.nCOL*sizeof(double*));
#ifdef DEBUG
    if (BlockMatrix.X == NULL)
        error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 10");
#endif
    for (i = 0; i < BlockMatrix.nCOL; i++) {
        BlockMatrix.X[i] = NULL;
        BlockMatrix.X[i] = (double *) malloc (BlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (BlockMatrix.X[i] == NULL)
            error("Create_CRS_BlockMatrix: %s\n", "Error Allocating Memory 11");
#endif
        for (j = 0; j < BlockMatrix.Block_nRow; j++)
            BlockMatrix.X[i][j] = 0.0;
    }
    
    // Free Memory
    if (degree != NULL) {
        delete[] degree;
        degree = NULL;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_CRS_BlockMatrix(int AddTime) {
    int iNode, i, j, k, idgn, idgn1, idgn2, idgn3;
    int ofdgn1[2], ofdgn2[2], ofdgn3[2];
    int n1, n2, n3, iCell;
    double Fplus[4], Fminus[4];
    double **Ap,**Am;

    // Initialize the CRS Matrix
    for (i = 0; i < BlockMatrix.DIM; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++)
                BlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Allocate Memory to Store Jacobians
    // Aplus
    Ap = NULL;
    Ap = (double **) malloc(4*sizeof(double*));
#ifdef DEBUG
    if (Ap == NULL)
        error("Compute_CRS_BlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < 4; i++) {
        Ap[i] = NULL;
        Ap[i] = (double *) malloc(4*sizeof(double));
#ifdef DEBUG
        if (Ap[i] == NULL)
            error("Compute_CRS_BlockMatrix: %s\n", "Error Allocating Memory 2");
#endif
        for (j = 0; j < 4; j++)
            Ap[i][j] = 0.0;
    }

    // Aminus
    Am = NULL;
    Am = (double **) malloc(4*sizeof(double*));
#ifdef DEBUG
    if (Am == NULL)
        error("Compute_CRS_BlockMatrix: %s\n", "Error Allocating Memory 3");
#endif
    for (i = 0; i < 4; i++) {
        Am[i] = NULL;
        Am[i] = (double *) malloc(4*sizeof(double));
#ifdef DEBUG
        if (Am[i] == NULL)
            error("Compute_CRS_BlockMatrix: %s\n", "Error Allocating Memory 4");
#endif
        for (j = 0; j < 4; j++)
            Am[i][j] = 0.0;
    }

    // Copy the Residuals to Block Matrix which is B
    // And Copy I/DeltaT
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        for (j = 0; j < BlockMatrix.Block_nRow; j++) {
            BlockMatrix.B[iNode][j] = -node[iNode].Resi[j];
            if (AddTime == 1) {
                // Get the diagonal location
                idgn = BlockMatrix.IAU[iNode];
                for (k = 0; k < BlockMatrix.Block_nCol; k++)
                    if (k == j)
                        BlockMatrix.A[idgn][j][k] = node[iNode].area/DeltaT[iNode];
            }
        }
    }

    // Get the Diagonal Component and Off-Diagonal Terms
    for (iCell = 0; iCell < mesh.inside; iCell++) {
        n1 = cell[iCell].node1;
        n2 = cell[iCell].node2;
        n3 = cell[iCell].node3;
        
        // Get the diagonal locations
        idgn1 = BlockMatrix.IAU[n1];
        idgn2 = BlockMatrix.IAU[n2];
        idgn3 = BlockMatrix.IAU[n3];

        // Get the Off-Diagonal Locations
        // Node1: ofdgn1[0]-> n2; ofdgn1[1]->n3
        // Node2: ofdgn2[0]-> n3; ofdgn2[1]->n1
        // Node3: ofdgn3[0]-> n1; ofdgn3[1]->n2
        ofdgn1[0] = ofdgn1[1] = -1;
        for (i = BlockMatrix.IA[n1]; i < BlockMatrix.IA[n1+1]; i++) {
            if (BlockMatrix.JA[i] == n2)
                ofdgn1[0] = i;
            if (BlockMatrix.JA[i] == n3)
                ofdgn1[1] = i;
        }
        ofdgn2[0] = ofdgn2[1] = -1;
        for (i = BlockMatrix.IA[n2]; i < BlockMatrix.IA[n2+1]; i++) {
            if (BlockMatrix.JA[i] == n3)
                ofdgn2[0] = i;
            if (BlockMatrix.JA[i] == n1)
                ofdgn2[1] = i;
        }
        ofdgn3[0] = ofdgn3[1] = -1;
        for (i = BlockMatrix.IA[n3]; i < BlockMatrix.IA[n3+1]; i++) {
            if (BlockMatrix.JA[i] == n1)
                ofdgn3[0] = i;
            if (BlockMatrix.JA[i] == n2)
                ofdgn3[1] = i;
        }

        // Get the Flux Across the Median Dual Line
        // ********************  Edge 1-2  *******************************
        // Initialize
        for (i = 0; i < 4; i++) {
            Fplus[i]  = 0.0;
            Fminus[i] = 0.0;
            for (j = 0; j < 4; j++) {
                Ap[i][j] = 0.0;
                Am[i][j] = 0.0;
            }
        }
        // Compute Split Flux Jacobians
        Compute_FluxJplus(node[n1].Q, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fplus, Ap);
        Compute_FluxJminus(node[n2].Q, cell[iCell].unxc12, cell[iCell].unyc12, cell[iCell].magc12, Fminus, Am);
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                // Diagonal
                BlockMatrix.A[idgn1][i][j] += Ap[i][j];
                BlockMatrix.A[idgn2][i][j] -= Am[i][j];
                // Off-Diagonal
                BlockMatrix.A[ofdgn1[0]][i][j] += Am[i][j];
                BlockMatrix.A[ofdgn2[1]][i][j] -= Ap[i][j];
            }
        }
        
        // ********************  Edge 2-3  *******************************
        // Initialize
        for (i = 0; i < 4; i++) {
            Fplus[i]  = 0.0;
            Fminus[i] = 0.0;
            for (j = 0; j < 4; j++) {
                Ap[i][j] = 0.0;
                Am[i][j] = 0.0;
            }
        }
        // Compute Split Flux Jacobians
        Compute_FluxJplus(node[n2].Q, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fplus, Ap);
        Compute_FluxJminus(node[n3].Q, cell[iCell].unxc23, cell[iCell].unyc23, cell[iCell].magc23, Fminus, Am);
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                // Diagonal
                BlockMatrix.A[idgn2][i][j] += Ap[i][j];
                BlockMatrix.A[idgn3][i][j] -= Am[i][j];
                // Off-Diagonal
                BlockMatrix.A[ofdgn2[0]][i][j] += Am[i][j];
                BlockMatrix.A[ofdgn3[1]][i][j] -= Ap[i][j];
            }
        }

        // ********************  Edge 3-1  *******************************
        // Initialize
        for (i = 0; i < 4; i++) {
            Fplus[i]  = 0.0;
            Fminus[i] = 0.0;
            for (j = 0; j < 4; j++) {
                Ap[i][j] = 0.0;
                Am[i][j] = 0.0;
            }
        }
        // Compute Split Flux Jacobians
        Compute_FluxJplus(node[n3].Q, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fplus, Ap);
        Compute_FluxJminus(node[n1].Q, cell[iCell].unxc31, cell[iCell].unyc31, cell[iCell].magc31, Fminus, Am);
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                // Diagonal
                BlockMatrix.A[idgn3][i][j] += Ap[i][j];
                BlockMatrix.A[idgn1][i][j] -= Am[i][j];
                // Off-Diagonal
                BlockMatrix.A[ofdgn3[0]][i][j] += Am[i][j];
                BlockMatrix.A[ofdgn1[1]][i][j] -= Ap[i][j];
            }
        }
    }

    // Free the Memory
    if (Ap != NULL) {
        for (i = 0; i < 4; i++) {
            if (Ap[i] != NULL)
                free(Ap[i]);
        }
        free(Ap);
    }

    if (Am != NULL) {
        for (i = 0; i < 4; i++) {
            if (Am[i] != NULL)
                free(Am[i]);
        }
        free(Am);
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Boundary_CRS_BlockMatrix() {
    int i, j, n1, n2, idgn, iedge;
    double Nx, Ny;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], E_Q[4], P_Q[4];
    double Rho, U, V, E, P;
    double Var1, Var2;

    // Calculate the Jacobian by Closing the Control Volume for Boundary Nodes
    // Freestream Boundary: BCTag = 3
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 3: 3000-3999 : Free Stream
        if ((boundaryEdge[i].bcType < 3000) || (boundaryEdge[i].bcType > 3999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[i].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // ************************ Node1 **********************************
        int k, l, jstart, jend;
        jstart = BlockMatrix.IA[n1];
        jend   = BlockMatrix.IA[n1+1];
        for (j = jstart; j < jend; j++) {
            for (k = 0; k < 4; k++) {
                for (l = 0; l < 4; l++)
                    BlockMatrix.A[j][k][l] = 0.0;
            }
        }
        idgn = BlockMatrix.IAU[n1];
        for (k = 0; k < 4; k++) {
            for (l = 0; l < 4; l++)
                if (l == k)
                    BlockMatrix.A[idgn][k][l] = 1.0;
        }

        // ************************ Node2 **********************************
        jstart = BlockMatrix.IA[n2];
        jend   = BlockMatrix.IA[n2+1];
        for (j = jstart; j < jend; j++) {
            for (k = 0; k < 4; k++) {
                for (l = 0; l < 4; l++)
                    BlockMatrix.A[j][k][l] = 0.0;
            }
        }
        idgn = BlockMatrix.IAU[n2];
        for (k = 0; k < 4; k++) {
            for (l = 0; l < 4; l++)
                if (l == k)
                    BlockMatrix.A[idgn][k][l] = 1.0;
        }
    }
    
    // Calculate the Jacobian by Closing the Control Volume for Boundary Nodes
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
        Nx = 0.5*edge[iedge].unx*edge[iedge].mag;
        Ny = 0.5*edge[iedge].uny*edge[iedge].mag;

        // ************************ Node1 **********************************
        // Calculate Rho and Partial derivatives
        Rho = node[n1].Q[0];
        Rho_Q[0] = 1.0;
        Rho_Q[1] = 0.0;
        Rho_Q[2] = 0.0;
        Rho_Q[3] = 0.0;

        // Calculate U and Partial derivatives
        U = node[n1].Q[1]/Rho;
        U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n1].Q[1])) / (Rho * Rho);

        // Calculate V and Partial derivatives
        V = node[n1].Q[2]/Rho;
        V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n1].Q[2])) / (Rho * Rho);

        // Calculate E and Partial derivatives
        E = node[n1].Q[3];
        E_Q[0] = 0.0;
        E_Q[1] = 0.0;
        E_Q[2] = 0.0;
        E_Q[3] = 1.0;

        // Calculate P and Partial derivatives
        P = (Gamma - 1.0)*(E - 0.5*Rho*(U*U + V*V));
        for (j = 0; j < 4; j++) {
            Var1 = 0.5*Rho_Q[j]*(U*U + V*V);
            Var2 = 0.5*Rho*(2.0*U*U_Q[j] + 2.0*V*V_Q[j]);
            P_Q[j] = (Gamma - 1.0) * (E_Q[j] - Var1 - Var2);
        }
        
        // Get the diagonal term
        idgn = BlockMatrix.IAU[n1];
        for (j = 0; j < 4; j++) {
            BlockMatrix.A[idgn][1][j] += Nx*P_Q[j];
            BlockMatrix.A[idgn][2][j] += Ny*P_Q[j];
        }
        
        // ************************ Node2 **********************************
        // Calculate Rho and Partial derivatives
        Rho = node[n2].Q[0];
        Rho_Q[0] = 1.0;
        Rho_Q[1] = 0.0;
        Rho_Q[2] = 0.0;
        Rho_Q[3] = 0.0;

        // Calculate U and Partial derivatives
        U = node[n2].Q[1]/Rho;
        U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n2].Q[1])) / (Rho * Rho);

        // Calculate V and Partial derivatives
        V = node[n2].Q[2]/Rho;
        V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n2].Q[2])) / (Rho * Rho);

        // Calculate E and Partial derivatives
        E = node[n2].Q[3];
        E_Q[0] = 0.0;
        E_Q[1] = 0.0;
        E_Q[2] = 0.0;
        E_Q[3] = 1.0;

        // Calculate P and Partial derivatives
        P = (Gamma - 1.0)*(E - 0.5*Rho*(U*U + V*V));
        for (j = 0; j < 4; j++) {
            Var1 = 0.5*Rho_Q[j]*(U*U + V*V);
            Var2 = 0.5*Rho*(2.0*U*U_Q[j] + 2.0*V*V_Q[j]);
            P_Q[j] = (Gamma - 1.0) * (E_Q[j] - Var1 - Var2);
        }

        // Get the diagonal term
        idgn = BlockMatrix.IAU[n2];
        for (j = 0; j < 4; j++) {
            BlockMatrix.A[idgn][1][j] += Nx*P_Q[j];
            BlockMatrix.A[idgn][2][j] += Ny*P_Q[j];
        }
    }
    
    // Calculate the Jacobian by Closing the Control Volume for Boundary Nodes
    // Solid Boundary: BCTag = 2
    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
    // Hence leads to CCW oriented boundary triangles
    for (i = 0; i < mesh.nbedges; i++) {
        // BCType 2: 2000-2999
        if ((boundaryEdge[i].bcType < 2000) || (boundaryEdge[i].bcType > 2999))
            continue;

        // Get The Nodes of Edge
        iedge = boundaryEdge[i].edgeNumber;
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        // Get the Normal - Components
        Nx = 0.5*edge[iedge].unx*edge[iedge].mag;
        Ny = 0.5*edge[iedge].uny*edge[iedge].mag;

        // ************************ Node1 **********************************
        // Calculate Rho and Partial derivatives
        Rho = node[n1].Q[0];
        Rho_Q[0] = 1.0;
        Rho_Q[1] = 0.0;
        Rho_Q[2] = 0.0;
        Rho_Q[3] = 0.0;

        // Calculate U and Partial derivatives
        U = node[n1].Q[1]/Rho;
        U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * node[n1].Q[1])) / (Rho * Rho);
        U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n1].Q[1])) / (Rho * Rho);

        // Calculate V and Partial derivatives
        V = node[n1].Q[2]/Rho;
        V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * node[n1].Q[2])) / (Rho * Rho);
        V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n1].Q[2])) / (Rho * Rho);

        // Calculate E and Partial derivatives
        E = node[n1].Q[3];
        E_Q[0] = 0.0;
        E_Q[1] = 0.0;
        E_Q[2] = 0.0;
        E_Q[3] = 1.0;

        // Calculate P and Partial derivatives
        P = (Gamma - 1.0)*(E - 0.5*Rho*(U*U + V*V));
        for (j = 0; j < 4; j++) {
            Var1 = 0.5*Rho_Q[j]*(U*U + V*V);
            Var2 = 0.5*Rho*(2.0*U*U_Q[j] + 2.0*V*V_Q[j]);
            P_Q[j] = (Gamma - 1.0) * (E_Q[j] - Var1 - Var2);
        }

        // Get the diagonal term
        idgn = BlockMatrix.IAU[n1];
        for (j = 0; j < 4; j++) {
            BlockMatrix.A[idgn][1][j] += Nx*P_Q[j];
            BlockMatrix.A[idgn][2][j] += Ny*P_Q[j];
        }

        // ************************ Node2 **********************************
        // Calculate Rho and Partial derivatives
        Rho = node[n2].Q[0];
        Rho_Q[0] = 1.0;
        Rho_Q[1] = 0.0;
        Rho_Q[2] = 0.0;
        Rho_Q[3] = 0.0;

        // Calculate U and Partial derivatives
        U = node[n2].Q[1]/Rho;
        U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * node[n2].Q[1])) / (Rho * Rho);
        U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n2].Q[1])) / (Rho * Rho);

        // Calculate V and Partial derivatives
        V = node[n2].Q[2]/Rho;
        V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * node[n2].Q[2])) / (Rho * Rho);
        V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * node[n2].Q[2])) / (Rho * Rho);

        // Calculate E and Partial derivatives
        E = node[n2].Q[3];
        E_Q[0] = 0.0;
        E_Q[1] = 0.0;
        E_Q[2] = 0.0;
        E_Q[3] = 1.0;

        // Calculate P and Partial derivatives
        P = (Gamma - 1.0)*(E - 0.5*Rho*(U*U + V*V));
        for (j = 0; j < 4; j++) {
            Var1 = 0.5*Rho_Q[j]*(U*U + V*V);
            Var2 = 0.5*Rho*(2.0*U*U_Q[j] + 2.0*V*V_Q[j]);
            P_Q[j] = (Gamma - 1.0) * (E_Q[j] - Var1 - Var2);
        }

        // Get the diagonal term
        idgn = BlockMatrix.IAU[n2];
        for (j = 0; j < 4; j++) {
            BlockMatrix.A[idgn][1][j] += Nx*P_Q[j];
            BlockMatrix.A[idgn][2][j] += Ny*P_Q[j];
        }
    }

//    // Calculate the Jacobian by Closing the Control Volume for Boundary Nodes
//    // Freestream Boundary: BCTag = 3
//    // Boundary Edge is such that Node1->Node2 direction mesh inside is on left
//    // Hence leads to CCW oriented boundary triangles
//    for (i = 0; i < mesh.nbedges; i++) {
//        // BCType 3: 3000-3999 : Free Stream
//        if ((boundaryEdge[i].bcType < 3000) || (boundaryEdge[i].bcType > 3999))
//            continue;
//
//        // Get The Nodes of Edge
//        iedge = boundaryEdge[i].edgeNumber;
//        n1 = edge[iedge].node1;
//        n2 = edge[iedge].node2;
//
//        // ************************ Node1 **********************************
//        int k, l, jstart, jend;
//        jstart = BlockMatrix.IA[n1];
//        jend   = BlockMatrix.IA[n1+1];
//        for (j = jstart; j < jend; j++) {
//            for (k = 0; k < 4; k++) {
//                for (l = 0; l < 4; l++)
//                    BlockMatrix.A[j][k][l] = 0.0;
//            }
//        }
//        idgn = BlockMatrix.IAU[n1];
//        for (k = 0; k < 4; k++) {
//            for (l = 0; l < 4; l++)
//                if (l == k)
//                    BlockMatrix.A[idgn][k][l] = 1.0;
//        }
//
//        // ************************ Node2 **********************************
//        jstart = BlockMatrix.IA[n2];
//        jend   = BlockMatrix.IA[n2+1];
//        for (j = jstart; j < jend; j++) {
//            for (k = 0; k < 4; k++) {
//                for (l = 0; l < 4; l++)
//                    BlockMatrix.A[j][k][l] = 0.0;
//            }
//        }
//        idgn = BlockMatrix.IAU[n2];
//        for (k = 0; k < 4; k++) {
//            for (l = 0; l < 4; l++)
//                if (l == k)
//                    BlockMatrix.A[idgn][k][l] = 1.0;
//        }
//    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Update_Solution() {
    int i, j;

    // DelataQ = Q(n+1) - Q(n)
    for (i = 0; i < mesh.nnodes; i++) {
        for (j = 0; j < BlockMatrix.Block_nRow; j++)
            node[i].Q[j] += BlockMatrix.X[i][j];
    }
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Solver_VanLeer::Compute_RMS() {
    int i, j;
    double max_rms = DBL_MIN;
    
    RMS[0] = RMS[1] = RMS[2] = RMS[3] = RMS_Res = 0.0;
    // DelataQ = Q(n+1) - Q(n)
    for (i = 0; i < mesh.nnodes; i++) {
        for (j = 0; j < 4; j++) {
            RMS[j] += BlockMatrix.X[i][j]*BlockMatrix.X[i][j];
            RMS_Res += node[i].Resi[j]*node[i].Resi[j];
        }
    }

    for (i = 0; i < 4; i++) {
        RMS[i] = sqrt(RMS[i])/mesh.nnodes;
        max_rms = MAX(max_rms, RMS[i]);
    }
    RMS_Res = sqrt(RMS_Res)/mesh.nnodes;
    max_rms = MAX(max_rms, RMS_Res);
    return max_rms;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Flux(double *Q, double uNx, double uNy, double Mag, double *Flux) {
    double Ubar, Pressure;
    
    Pressure = (Gamma - 1.0)*(Q[3] - 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0])));
    Ubar = ((Q[1]/Q[0])*uNx) + ((Q[2]/Q[0])*uNy);

    Flux[0] = Mag*Q[0]*Ubar;
    Flux[1] = Mag*(Q[1]*Ubar + uNx*Pressure);
    Flux[2] = Mag*(Q[2]*Ubar + uNy*Pressure);
    Flux[3] = Mag*(Q[3] + Pressure)*Ubar;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Fplus(double *Q, double uNx, double uNy, double Mag, double *Fplus) {
    double C, Ubar, Mbar, Pressure;
    
    Pressure = (Gamma - 1.0)*(Q[3] - 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0])));
    C = sqrt(Gamma*Pressure/Q[0]);
    Ubar = ((Q[1]/Q[0])*uNx) + ((Q[2]/Q[0])*uNy);
    Mbar = Ubar/C;

    if (fabs(Mbar) <= 1.0) {
        Fplus[0] = Mag*0.25*Q[0]*C*(Mbar + 1.0)*(Mbar + 1.0);
        Fplus[1] = Fplus[0]*((uNx/Gamma)*(-Ubar + 2.0*C) + Q[1]/Q[0]);
        Fplus[2] = Fplus[0]*((uNy/Gamma)*(-Ubar + 2.0*C) + Q[2]/Q[0]);
        Fplus[3] = Fplus[0]*((((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0))
                + (0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(Q[0]*Q[0])));
    } else if (Mbar >= 1.0) {
        Fplus[0] = Mag*Q[0]*Ubar;
        Fplus[1] = Mag*(Q[1]*Ubar + uNx*Pressure);
        Fplus[2] = Mag*(Q[2]*Ubar + uNy*Pressure);
        Fplus[3] = Mag*(Q[3] + Pressure)*Ubar;
    } else {
        Fplus[0] = 0.0;
        Fplus[1] = 0.0;
        Fplus[2] = 0.0;
        Fplus[3] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Fminus(double *Q, double uNx, double uNy, double Mag, double *Fminus) {
    double C, Ubar, Mbar, Pressure;
    
    Pressure = (Gamma - 1.0)*(Q[3] - 0.5*Q[0]*((Q[1]/Q[0])*(Q[1]/Q[0]) + (Q[2]/Q[0])*(Q[2]/Q[0])));
    C = sqrt(Gamma*Pressure/Q[0]);
    Ubar = ((Q[1]/Q[0])*uNx) + ((Q[2]/Q[0])*uNy);
    Mbar = Ubar/C;

    if (fabs(Mbar) <= 1.0) {
        Fminus[0] = -Mag*0.25*Q[0]*C*(Mbar - 1.0)*(Mbar - 1.0);
        Fminus[1] = Fminus[0]*((uNx/Gamma)*(-Ubar - 2.0*C) + Q[1]/Q[0]);
        Fminus[2] = Fminus[0]*((uNy/Gamma)*(-Ubar - 2.0*C) + Q[2]/Q[0]);
        Fminus[3] = Fminus[0]*((((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0))
                + (0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(Q[0]*Q[0])));
    } else if (Mbar <= -1.0) {
        Fminus[0] = Mag*Q[0]*Ubar;
        Fminus[1] = Mag*(Q[1]*Ubar + uNx*Pressure);
        Fminus[2] = Mag*(Q[2]*Ubar + uNy*Pressure);
        Fminus[3] = Mag*(Q[3] + Pressure)*Ubar;
    } else {
        Fminus[0] = 0.0;
        Fminus[1] = 0.0;
        Fminus[2] = 0.0;
        Fminus[3] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_FluxJacobian(double *Q, double uNx, double uNy, double Mag, double *Flux, double **A) {
    int i, j;
    double Var1, Var2;
    double U, V, Ubar;
    double Rho, P, E;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // Ubar_Q => Velocity Flux Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], Ubar_Q[4], E_Q[4], P_Q[4];
    
    // Calculate Rho and Partial derivatives
    Rho = Q[0];
    Rho_Q[0] = 1.0;
    Rho_Q[1] = 0.0;
    Rho_Q[2] = 0.0;
    Rho_Q[3] = 0.0;

    // Calculate U and Partial derivatives
    U = Q[1]/Q[0];
    U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[1])) / (Rho * Rho);
    U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * Q[1])) / (Rho * Rho);
    U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * Q[1])) / (Rho * Rho);
    U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[1])) / (Rho * Rho);

    // Calculate V and Partial derivatives
    V = Q[2]/Q[0];
    V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[2])) / (Rho * Rho);
    V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * Q[2])) / (Rho * Rho);
    V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * Q[2])) / (Rho * Rho);
    V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[2])) / (Rho * Rho);

    // Calculate Ubar and Partial derivatives
    Ubar = U * uNx + V * uNy;
    for (i = 0; i < 4; i++)
        Ubar_Q[i] = U_Q[i] * uNx + V_Q[i] * uNy;

    // Calculate E and Partial derivatives
    E = Q[3];
    E_Q[0] = 0.0;
    E_Q[1] = 0.0;
    E_Q[2] = 0.0;
    E_Q[3] = 1.0;

    // Calculate P and Partial derivatives
    P = (Gamma - 1.0) * (Q[3] - 0.5 * Rho * (U * U + V * V));
    for (i = 0; i < 4; i++) {
        Var1 = 0.5*Rho_Q[i]*(U*U + V*V);
        Var2 = 0.5*Rho*(2.0*U*U_Q[i] + 2.0*V*V_Q[i]);
        P_Q[i] = (Gamma - 1.0) * (E_Q[i] - Var1 - Var2);
    }

    // Calculate the Flux Vector
    Flux[0] = (Rho * Ubar);
    Flux[1] = (Rho * U * Ubar + uNx * P);
    Flux[2] = (Rho * V * Ubar + uNy * P);
    Flux[3] = ((E + P) * Ubar);

    for (i = 0; i < 4; i++)
        Flux[i] = Flux[i] * Mag;

    // Calculate the Flux Jacobian
    for (i = 0; i < 4; i++) {
        // Row 1
        A[0][i] = (Rho*Ubar_Q[i]) + (Rho_Q[i]*Ubar);
        // Row 2
        A[1][i] = (Rho*U*Ubar_Q[i]) + (Rho*U_Q[i]*Ubar) + (Rho_Q[i]*U*Ubar) + uNx*P_Q[i];
        // Row 3
        A[2][i] = (Rho*V*Ubar_Q[i]) + (Rho*V_Q[i]*Ubar) + (Rho_Q[i]*V*Ubar) + uNy*P_Q[i];
        // Row 4
        A[3][i] = ((E_Q[i] + P_Q[i])*Ubar) + ((E + P)*Ubar_Q[i]);
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            A[i][j] = A[i][j] * Mag;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_FluxJplus(double *Q, double uNx, double uNy, double Mag, double *Fplus, double **Ap) {
    int i, j;
    double Var1, Var2, Var3;
    double U, V, Ubar, Mbar;
    double Rho, C, P, E;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // Ubar_Q => Velocity Flux Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    // C_Q    => Speed of Sound Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], Ubar_Q[4], E_Q[4], P_Q[4], C_Q[4];
    
    // Calculate Rho and Partial derivatives
    Rho = Q[0];
    Rho_Q[0] = 1.0;
    Rho_Q[1] = 0.0;
    Rho_Q[2] = 0.0;
    Rho_Q[3] = 0.0;

    // Calculate U and Partial derivatives
    U = Q[1]/Q[0];
    U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[1])) / (Rho * Rho);
    U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * Q[1])) / (Rho * Rho);
    U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * Q[1])) / (Rho * Rho);
    U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[1])) / (Rho * Rho);

    // Calculate V and Partial derivatives
    V = Q[2]/Q[0];
    V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[2])) / (Rho * Rho);
    V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * Q[2])) / (Rho * Rho);
    V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * Q[2])) / (Rho * Rho);
    V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[2])) / (Rho * Rho);

    // Calculate Ubar and Partial derivatives
    Ubar = U * uNx + V * uNy;
    for (i = 0; i < 4; i++)
        Ubar_Q[i] = U_Q[i] * uNx + V_Q[i] * uNy;

    // Calculate E and Partial derivatives
    E = Q[3];
    E_Q[0] = 0.0;
    E_Q[1] = 0.0;
    E_Q[2] = 0.0;
    E_Q[3] = 1.0;

    // Calculate P and Partial derivatives
    P = (Gamma - 1.0) * (Q[3] - 0.5 * Rho * (U * U + V * V));
    for (i = 0; i < 4; i++) {
        Var1 = 0.5*Rho_Q[i]*(U*U + V*V);
        Var2 = 0.5*Rho*(2.0*U*U_Q[i] + 2.0*V*V_Q[i]);
        P_Q[i] = (Gamma - 1.0) * (E_Q[i] - Var1 - Var2);
    }

    // Calculate C and Partial derivatives
    C = sqrt(Gamma * (P / Rho));
    Var1 = sqrt(Gamma * (P / Rho));
    for (i = 0; i < 4; i++) {
        Var2 = (Gamma * Rho * P_Q[i] - Gamma * P * Rho_Q[i]) / (Rho * Rho);
        C_Q[i] = 0.5 * (Var2 / Var1);
    }

    Mbar = Ubar/C;

    // Calculate the Van Leer Split Flux and Splited Jacobians
    if (fabs(Mbar) <= 1.0) {
        // Compute Van-Leer Fplus
        Fplus[0] = 0.25*Rho*C*(Mbar + 1.0)*(Mbar + 1.0);
        Fplus[1] = Fplus[0]*((uNx/Gamma)*(-Ubar + 2.0*C) + U);
        Fplus[2] = Fplus[0]*((uNy/Gamma)*(-Ubar + 2.0*C) + V);
        Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0);
        Fplus[3] = Fplus[0]*(Var1 + 0.5*(U*U + V*V));

        // Compute Van-Leer Aplus
        for (i = 0; i < 4; i++) {
            // Row 1
            Var1 = Mbar + 1.0;
            Var2 = (C*Ubar_Q[i] - Ubar*C_Q[i])/(C*C);
            Ap[0][i] = 0.25*((Rho_Q[i]*C*Var1*Var1) + (Rho*C_Q[i]*Var1*Var1) + (2.0*Rho*C*Var1*Var2));
            // Row 2
            Var1 = ((uNx*(-Ubar + 2.0*C))/Gamma) + U;
            Var2 = ((uNx*(-Ubar_Q[i] + 2.0*C_Q[i]))/Gamma) + U_Q[i];
            Ap[1][i] = Ap[0][i]*Var1 + Fplus[0]*Var2;
            // Row 3
            Var1 = ((uNy*(-Ubar + 2.0*C))/Gamma) + V;
            Var2 = ((uNy*(-Ubar_Q[i] + 2.0*C_Q[i]))/Gamma) + V_Q[i];
            Ap[2][i] = Ap[0][i]*Var1 + Fplus[0]*Var2;
            // Row 4
            Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0) + 0.5*(U*U + V*V);
            Var2 = (-2.0*(Gamma - 1.0)*Ubar*Ubar_Q[i]) + (2.0*(Gamma - 1.0)*Ubar_Q[i]*C) + (2.0*(Gamma - 1.0)*Ubar*C_Q[i]) + 4.0*C*C_Q[i];
            Var3 = (Var2/(Gamma*Gamma - 1.0)) + U*U_Q[i] + V*V_Q[i];
            Ap[3][i] = Ap[0][i]*Var1 + Fplus[0]*Var3;
        }

        for (i = 0; i < 4; i++) {
            Fplus[i] = Fplus[i] * Mag;
            for (j = 0; j < 4; j++)
                Ap[i][j] = Ap[i][j] * Mag;
        }
    } else if (Mbar >= 1.0) {
        // Compute Van-Leer Fplus
        Fplus[0] = (Rho * Ubar);
        Fplus[1] = (Rho * U * Ubar + uNx * P);
        Fplus[2] = (Rho * V * Ubar + uNy * P);
        Fplus[3] = ((E + P) * Ubar);

        // Compute Van-Leer Aplus
        for (i = 0; i < 4; i++) {
            // Row 1
            Ap[0][i] = (Rho*Ubar_Q[i]) + (Rho_Q[i]*Ubar);
            // Row 2
            Ap[1][i] = (Rho*U*Ubar_Q[i]) + (Rho*U_Q[i]*Ubar) + (Rho_Q[i]*U*Ubar) + uNx*P_Q[i];
            // Row 3
            Ap[2][i] = (Rho*V*Ubar_Q[i]) + (Rho*V_Q[i]*Ubar) + (Rho_Q[i]*V*Ubar) + uNy*P_Q[i];
            // Row 4
            Ap[3][i] = ((E_Q[i] + P_Q[i])*Ubar) + ((E + P)*Ubar_Q[i]);
        }

        for (i = 0; i < 4; i++) {
            Fplus[i] = Fplus[i] * Mag;
            for (j = 0; j < 4; j++)
                Ap[i][j] = Ap[i][j] * Mag;
        }
    } else {
        for (i = 0; i < 4; i++) {
            Fplus[i] = 0.0;
            for (j = 0; j < 4; j++)
                Ap[i][j] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_FluxJminus(double *Q, double uNx, double uNy, double Mag, double *Fminus, double **Am) {
    int i, j;
    double Var1, Var2, Var3;
    double U, V, Ubar, Mbar;
    double Rho, C, P, E;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // Ubar_Q => Velocity Flux Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    // C_Q    => Speed of Sound Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], Ubar_Q[4], E_Q[4], P_Q[4], C_Q[4];

    // Calculate Rho and Partial derivatives
    Rho = Q[0];
    Rho_Q[0] = 1.0;
    Rho_Q[1] = 0.0;
    Rho_Q[2] = 0.0;
    Rho_Q[3] = 0.0;

    // Calculate U and Partial derivatives
    U = Q[1]/Q[0];
    U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[1])) / (Rho * Rho);
    U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * Q[1])) / (Rho * Rho);
    U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * Q[1])) / (Rho * Rho);
    U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[1])) / (Rho * Rho);

    // Calculate V and Partial derivatives
    V = Q[2]/Q[0];
    V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[2])) / (Rho * Rho);
    V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * Q[2])) / (Rho * Rho);
    V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * Q[2])) / (Rho * Rho);
    V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[2])) / (Rho * Rho);

    // Calculate Ubar and Partial derivatives
    Ubar = U * uNx + V * uNy;
    for (i = 0; i < 4; i++)
        Ubar_Q[i] = U_Q[i] * uNx + V_Q[i] * uNy;

    // Calculate E and Partial derivatives
    E = Q[3];
    E_Q[0] = 0.0;
    E_Q[1] = 0.0;
    E_Q[2] = 0.0;
    E_Q[3] = 1.0;

    // Calculate P and Partial derivatives
    P = (Gamma - 1.0) * (Q[3] - 0.5 * Rho * (U * U + V * V));
    for (i = 0; i < 4; i++) {
        Var1 = 0.5*Rho_Q[i]*(U*U + V*V);
        Var2 = 0.5*Rho*(2.0*U*U_Q[i] + 2.0*V*V_Q[i]);
        P_Q[i] = (Gamma - 1.0) * (E_Q[i] - Var1 - Var2);
    }

    // Calculate C and Partial derivatives
    C = sqrt(Gamma * (P / Rho));
    Var1 = sqrt(Gamma * (P / Rho));
    for (i = 0; i < 4; i++) {
        Var2 = (Gamma * Rho * P_Q[i] - Gamma * P * Rho_Q[i]) / (Rho * Rho);
        C_Q[i] = 0.5 * (Var2 / Var1);
    }

    Mbar = Ubar/C;

    // Calculate the Van Leer Split Flux and Splited Jacobians
    if (fabs(Mbar) <= 1.0) {
        // Compute Van-Leer Fminus
        Fminus[0] = -0.25*Rho*C*(Mbar - 1.0)*(Mbar - 1.0);
        Fminus[1] = Fminus[0]*((uNx/Gamma)*(-Ubar - 2.0*C) + U);
        Fminus[2] = Fminus[0]*((uNy/Gamma)*(-Ubar - 2.0*C) + V);
        Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0);
        Fminus[3] = Fminus[0]*(Var1 + 0.5*(U*U + V*V));

        // Compute Van-Leer Aminus
        for (i = 0; i < 4; i++) {
            // Row 1
            Var1 = Mbar - 1.0;
            Var2 = (C *Ubar_Q[i] - Ubar*C_Q[i])/(C*C);
            Am[0][i] = -0.25*((Rho_Q[i]*C*Var1*Var1) + (Rho*C_Q[i]*Var1*Var1) + (2.0*Rho*C*Var1*Var2));
            // Row 2
            Var1 = ((uNx*(-Ubar - 2.0*C))/Gamma) + U;
            Var2 = ((uNx*(-Ubar_Q[i] - 2.0*C_Q[i]))/Gamma) + U_Q[i];
            Am[1][i] = Am[0][i]*Var1 + Fminus[0]*Var2;
            // Row 3
            Var1 = ((uNy*(-Ubar - 2.0*C))/Gamma) + V;
            Var2 = ((uNy*(-Ubar_Q[i] - 2.0*C_Q[i]))/Gamma) + V_Q[i];
            Am[2][i] = Am[0][i]*Var1 + Fminus[0]*Var2;
            // Row 4
            Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0) + 0.5*(U*U + V*V);
            Var2 = (-2.0*(Gamma - 1.0)*Ubar*Ubar_Q[i]) - (2.0*(Gamma - 1.0)*Ubar_Q[i]*C) - (2.0*(Gamma - 1.0)*Ubar*C_Q[i]) + 4.0*C*C_Q[i];
            Var3 = (Var2/(Gamma*Gamma - 1.0)) + U*U_Q[i] + V*V_Q[i];
            Am[3][i] = Am[0][i]*Var1 + Fminus[0]*Var3;
        }

        for (i = 0; i < 4; i++) {
            Fminus[i] = Fminus[i] * Mag;
            for (j = 0; j < 4; j++)
                Am[i][j] = Am[i][j] * Mag;
        }
    } else if (Mbar <= -1.0) {
        // Compute Van-Leer Fminus
        Fminus[0] = (Rho * Ubar);
        Fminus[1] = (Rho * U * Ubar + uNx * P);
        Fminus[2] = (Rho * V * Ubar + uNy * P);
        Fminus[3] = ((E + P) * Ubar);

        // Compute Van-Leer Aminus
        for (i = 0; i < 4; i++) {
            // Row 1
            Am[0][i] = (Rho*Ubar_Q[i]) + (Rho_Q[i]*Ubar);
            // Row 2
            Am[1][i] = (Rho*U*Ubar_Q[i]) + (Rho*U_Q[i]*Ubar) + (Rho_Q[i]*U*Ubar) + uNx*P_Q[i];
            // Row 3
            Am[2][i] = (Rho*V*Ubar_Q[i]) + (Rho*V_Q[i]*Ubar) + (Rho_Q[i]*V*Ubar) + uNy*P_Q[i];
            // Row 4
            Am[3][i] = ((E_Q[i] + P_Q[i])*Ubar) + ((E + P)*Ubar_Q[i]);
        }

        for (i = 0; i < 4; i++) {
            Fminus[i] = Fminus[i] * Mag;
            for (j = 0; j < 4; j++)
                Am[i][j] = Am[i][j] * Mag;
        }
    } else {
        for (i = 0; i < 4; i++) {
            Fminus[i] = 0.0;
            for (j = 0; j < 4; j++)
                Am[i][j] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Solver_VanLeer::Compute_Flux_Jacobians(double *Q, double uNx, double uNy, double Mag,
        double *Flux, double *Fplus, double *Fminus, double **A, double **Ap, double **Am) {
    int i, j;
    double Var1, Var2, Var3;
    double U, V, Ubar, Mbar;
    double Rho, C, P, E;
    // Variables to Store Partial Derivatives with Conserved Variables
    // Rho_Q  => Density Derivatives
    // U_Q    => U Velocity Derivatives
    // V_Q    => V Velocity Derivatives
    // Ubar_Q => Velocity Flux Derivatives
    // E_Q    => Energy Derivatives
    // P_Q    => Pressure Derivatives
    // C_Q    => Speed of Sound Derivatives
    double Rho_Q[4], U_Q[4], V_Q[4], Ubar_Q[4], E_Q[4], P_Q[4], C_Q[4];

    // Calculate Rho and Partial derivatives
    Rho = Q[0];
    Rho_Q[0] = 1.0;
    Rho_Q[1] = 0.0;
    Rho_Q[2] = 0.0;
    Rho_Q[3] = 0.0;

    // Calculate U and Partial derivatives
    U = Q[1]/Q[0];
    U_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[1])) / (Rho * Rho);
    U_Q[1] = ((Rho * 1.0)-(Rho_Q[1] * Q[1])) / (Rho * Rho);
    U_Q[2] = ((Rho * 0.0)-(Rho_Q[2] * Q[1])) / (Rho * Rho);
    U_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[1])) / (Rho * Rho);

    // Calculate V and Partial derivatives
    V = Q[2]/Q[0];
    V_Q[0] = ((Rho * 0.0)-(Rho_Q[0] * Q[2])) / (Rho * Rho);
    V_Q[1] = ((Rho * 0.0)-(Rho_Q[1] * Q[2])) / (Rho * Rho);
    V_Q[2] = ((Rho * 1.0)-(Rho_Q[2] * Q[2])) / (Rho * Rho);
    V_Q[3] = ((Rho * 0.0)-(Rho_Q[3] * Q[2])) / (Rho * Rho);

    // Calculate Ubar and Partial derivatives
    Ubar = U * uNx + V * uNy;
    for (i = 0; i < 4; i++)
        Ubar_Q[i] = U_Q[i] * uNx + V_Q[i] * uNy;

    // Calculate E and Partial derivatives
    E = Q[3];
    E_Q[0] = 0.0;
    E_Q[1] = 0.0;
    E_Q[2] = 0.0;
    E_Q[3] = 1.0;

    // Calculate P and Partial derivatives
    P = (Gamma - 1.0) * (Q[3] - 0.5 * Rho * (U * U + V * V));
    for (i = 0; i < 4; i++) {
        Var1 = 0.5*Rho_Q[i]*(U*U + V*V);
        Var2 = 0.5*Rho*(2.0*U*U_Q[i] + 2.0*V*V_Q[i]);
        P_Q[i] = (Gamma - 1.0) * (E_Q[i] - Var1 - Var2);
    }

    // Calculate C and Partial derivatives
    C = sqrt(Gamma * (P / Rho));
    Var1 = sqrt(Gamma * (P / Rho));
    for (i = 0; i < 4; i++) {
        Var2 = (Gamma * Rho * P_Q[i] - Gamma * P * Rho_Q[i]) / (Rho * Rho);
        C_Q[i] = 0.5 * (Var2 / Var1);
    }

    Mbar = Ubar/C;

    // Calculate the Flux Vector
    Flux[0] = (Rho * Ubar);
    Flux[1] = (Rho * U * Ubar + uNx * P);
    Flux[2] = (Rho * V * Ubar + uNy * P);
    Flux[3] = ((E + P) * Ubar);

    for (i = 0; i < 4; i++)
        Flux[i] = Flux[i] * Mag;

    // Calculate the Flux Jacobian
    for (i = 0; i < 4; i++) {
        // Row 1
        A[0][i] = (Rho*Ubar_Q[i]) + (Rho_Q[i]*Ubar);
        // Row 2
        A[1][i] = (Rho*U*Ubar_Q[i]) + (Rho*U_Q[i]*Ubar) + (Rho_Q[i]*U*Ubar) + uNx*P_Q[i];
        // Row 3
        A[2][i] = (Rho*V*Ubar_Q[i]) + (Rho*V_Q[i]*Ubar) + (Rho_Q[i]*V*Ubar) + uNy*P_Q[i];
        // Row 4
        A[3][i] = ((E_Q[i] + P_Q[i])*Ubar) + ((E + P)*Ubar_Q[i]);
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            A[i][j] = A[i][j] * Mag;
    }

    // Calculate the Van Leer Split Flux and Splited Jacobians
    if (Mbar > 1.0) {
        for (i = 0; i < 4; i++) {
            Fplus[i] = Flux[i];
            Fminus[i] = 0.0;
            for (j = 0; j < 4; j++) {
                Ap[i][j] = A[i][j];
                Am[i][j] = 0.0;
            }
        }
    } else if (Mbar < -1.0) {
        for (i = 0; i < 4; i++) {
            Fplus[i] = 0.0;
            Fminus[i] = Flux[i];
            for (j = 0; j < 4; j++) {
                Ap[i][j] = 0.0;
                Am[i][j] = A[i][j];
            }
        }
    } else {
        // Compute Van-Leer Fplus
        Fplus[0] = 0.25*Rho*C*(Mbar + 1.0)*(Mbar + 1.0);
        Fplus[1] = Fplus[0]*((uNx/Gamma)*(-Ubar + 2.0*C) + U);
        Fplus[2] = Fplus[0]*((uNy/Gamma)*(-Ubar + 2.0*C) + V);
        Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0);
        Fplus[3] = Fplus[0]*(Var1 + 0.5*(U*U + V*V));

        // Compute Van-Leer Aplus
        for (i = 0; i < 4; i++) {
            // Row 1
            Var1 = Mbar + 1.0;
            Var2 = (C*Ubar_Q[i] - Ubar*C_Q[i])/(C*C);
            Ap[0][i] = 0.25*((Rho_Q[i]*C*Var1*Var1) + (Rho*C_Q[i]*Var1*Var1) + (2.0*Rho*C*Var1*Var2));
            // Row 2
            Var1 = ((uNx*(-Ubar + 2.0*C))/Gamma) + U;
            Var2 = ((uNx*(-Ubar_Q[i] + 2.0*C_Q[i]))/Gamma) + U_Q[i];
            Ap[1][i] = Ap[0][i]*Var1 + Fplus[0]*Var2;
            // Row 3
            Var1 = ((uNy*(-Ubar + 2.0*C))/Gamma) + V;
            Var2 = ((uNy*(-Ubar_Q[i] + 2.0*C_Q[i]))/Gamma) + V_Q[i];
            Ap[2][i] = Ap[0][i]*Var1 + Fplus[0]*Var2;
            // Row 4
            Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) + (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0) + 0.5*(U*U + V*V);
            Var2 = (-2.0*(Gamma - 1.0)*Ubar*Ubar_Q[i]) + (2.0*(Gamma - 1.0)*Ubar_Q[i]*C) + (2.0*(Gamma - 1.0)*Ubar*C_Q[i]) + 4.0*C*C_Q[i];
            Var3 = (Var2/(Gamma*Gamma - 1.0)) + U*U_Q[i] + V*V_Q[i];
            Ap[3][i] = Ap[0][i]*Var1 + Fplus[0]*Var3;
        }

        for (i = 0; i < 4; i++) {
            Fplus[i] = Fplus[i] * Mag;
            for (j = 0; j < 4; j++)
                Ap[i][j] = Ap[i][j] * Mag;
        }

        // Compute Van-Leer Fminus
        Fminus[0] = -0.25*Rho*C*(Mbar - 1.0)*(Mbar - 1.0);
        Fminus[1] = Fminus[0]*((uNx/Gamma)*(-Ubar - 2.0*C) + U);
        Fminus[2] = Fminus[0]*((uNy/Gamma)*(-Ubar - 2.0*C) + V);
        Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0);
        Fminus[3] = Fminus[0]*(Var1 + 0.5*(U*U + V*V));


        // Compute Van-Leer Aminus
        for (i = 0; i < 4; i++) {
            // Row 1
            Var1 = Mbar - 1.0;
            Var2 = (C *Ubar_Q[i] - Ubar*C_Q[i])/(C*C);
            Am[0][i] = -0.25*((Rho_Q[i]*C*Var1*Var1) + (Rho*C_Q[i]*Var1*Var1) + (2.0*Rho*C*Var1*Var2));
            // Row 2
            Var1 = ((uNx*(-Ubar - 2.0*C))/Gamma) + U;
            Var2 = ((uNx*(-Ubar_Q[i] - 2.0*C_Q[i]))/Gamma) + U_Q[i];
            Am[1][i] = Am[0][i]*Var1 + Fminus[0]*Var2;
            // Row 3
            Var1 = ((uNy*(-Ubar - 2.0*C))/Gamma) + V;
            Var2 = ((uNy*(-Ubar_Q[i] - 2.0*C_Q[i]))/Gamma) + V_Q[i];
            Am[2][i] = Am[0][i]*Var1 + Fminus[0]*Var2;
            // Row 4
            Var1 = ((-(Gamma - 1.0)*Ubar*Ubar) - (2.0*(Gamma - 1.0)*Ubar*C) + 2.0*C*C)/(Gamma*Gamma - 1.0) + 0.5*(U*U + V*V);
            Var2 = (-2.0*(Gamma - 1.0)*Ubar*Ubar_Q[i]) - (2.0*(Gamma - 1.0)*Ubar_Q[i]*C) - (2.0*(Gamma - 1.0)*Ubar*C_Q[i]) + 4.0*C*C_Q[i];
            Var3 = (Var2/(Gamma*Gamma - 1.0)) + U*U_Q[i] + V*V_Q[i];
            Am[3][i] = Am[0][i]*Var1 + Fminus[0]*Var3;
        }

        for (i = 0; i < 4; i++) {
            Fminus[i] = Fminus[i] * Mag;
            for (j = 0; j < 4; j++)
                Am[i][j] = Am[i][j] * Mag;
        }
    }
}

