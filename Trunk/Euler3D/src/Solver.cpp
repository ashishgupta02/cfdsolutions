/*******************************************************************************
 * File:        Solver.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "RestartIO.h"
#include "MeshIO.h"
#include "Commons.h"
#include "Solver.h"
#include "MC.h"
#include "Gradient.h"
#include "CompressibleUtils.h"
#include "Residual_Smoothing.h"
#include "Material.h"
#include "DebugSolver.h"

int SolverIteration;

// Local Time
double *DeltaT;
double MinDeltaT;
double MaxDeltaT;

// Conservative or Primitive Variables
double *Q1;
double *Q2;
double *Q3;
double *Q4;
double *Q5;

// Conservative or Primitive Variables Gradients
double *Q1x;
double *Q1y;
double *Q1z;

double *Q2x;
double *Q2y;
double *Q2z;

double *Q3x;
double *Q3y;
double *Q3z;

double *Q4x;
double *Q4y;
double *Q4z;

double *Q5x;
double *Q5y;
double *Q5z;

// Residuals (Convective and Dissipative)
double *Res1;
double *Res2;
double *Res3;
double *Res4;
double *Res5;

double *Res1_Diss;
double *Res2_Diss;
double *Res3_Diss;
double *Res4_Diss;
double *Res5_Diss;

// Limiters
double *Limiter_Phi1;
double *Limiter_Phi2;
double *Limiter_Phi3;
double *Limiter_Phi4;
double *Limiter_Phi5;

//RMS
double RMS[5];
double RMS_Res;

// Min and Max EigenValue Value
double MinEigenLamda1;
double MaxEigenLamda1;
double MinEigenLamda4;
double MaxEigenLamda4;
double MinEigenLamda5;
double MaxEigenLamda5;

// Min and Max Precondition Variable
double *PrecondSigma;
double MinPrecondSigma;
double MaxPrecondSigma;

// MC CRS Matrix
MC_CRS SolverBlockMatrix;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Init(void) {
    // Local Time
    DeltaT          = NULL;
    MinDeltaT       = DBL_MAX;
    MaxDeltaT       = DBL_MIN;

    // Conservative or Primitive Variables
    Q1              = NULL;
    Q2              = NULL;
    Q3              = NULL;
    Q4              = NULL;
    Q5              = NULL;

    // Conservative or Primitive Variables Gradients
    Q1x             = NULL;
    Q1y             = NULL;
    Q1z             = NULL;

    Q2x             = NULL;
    Q2y             = NULL;
    Q2z             = NULL;

    Q3x             = NULL;
    Q3y             = NULL;
    Q3z             = NULL;

    Q4x             = NULL;
    Q4y             = NULL;
    Q4z             = NULL;

    Q5x             = NULL;
    Q5y             = NULL;
    Q5z             = NULL;
    
    // Residuals (Convective and Dissipative)
    Res1            = NULL;
    Res2            = NULL;
    Res3            = NULL;
    Res4            = NULL;
    Res5            = NULL;

    Res1_Diss       = NULL;
    Res2_Diss       = NULL;
    Res3_Diss       = NULL;
    Res4_Diss       = NULL;
    Res5_Diss       = NULL;
    
    // Limiters
    Limiter_Phi1    = NULL;
    Limiter_Phi2    = NULL;
    Limiter_Phi3    = NULL;
    Limiter_Phi4    = NULL;
    Limiter_Phi5    = NULL;
    
    // RMS
    RMS[0]          = 0.0;
    RMS[1]          = 0.0;
    RMS[2]          = 0.0;
    RMS[3]          = 0.0;
    RMS[4]          = 0.0;
    RMS_Res         = 0.0;
    
    // Min and Max EigenValue Value
    MinEigenLamda1   = DBL_MAX;
    MaxEigenLamda1   = DBL_MIN;
    MinEigenLamda4   = DBL_MAX;
    MaxEigenLamda4   = DBL_MIN;
    MinEigenLamda5   = DBL_MAX;
    MaxEigenLamda5   = DBL_MIN;
    
    // Min and Max Precondition Variable
    PrecondSigma    = NULL;
    MinPrecondSigma = DBL_MAX;
    MaxPrecondSigma = DBL_MIN;
    
    // MC CRS Matrix
    SolverBlockMatrix.Block_nRow = 0;
    SolverBlockMatrix.Block_nCol = 0;
    SolverBlockMatrix.nROW       = 0;
    SolverBlockMatrix.nCOL       = 0;
    SolverBlockMatrix.DIM        = 0;
    SolverBlockMatrix.A          = NULL;
    SolverBlockMatrix.B          = NULL;
    SolverBlockMatrix.X          = NULL;
    SolverBlockMatrix.IA         = NULL;
    SolverBlockMatrix.IAU        = NULL;
    SolverBlockMatrix.JA         = NULL;
    
    // RMS Writer Initialization
    RMS_Writer_Init();
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Finalize(void) {
    // Conservative or Primitive Variables
    if (Q1 != NULL)
        delete[] Q1;
    if (Q2 != NULL)
        delete[] Q2;
    if (Q3 != NULL)
        delete[] Q3;
    if (Q4 != NULL)
        delete[] Q4;
    if (Q5 != NULL)
        delete[] Q5;
    Q1 = NULL;
    Q2 = NULL;
    Q3 = NULL;
    Q4 = NULL;
    Q5 = NULL;

     // Conservative or Primitive Variables Gradient
    if (Q1x != NULL)
        delete[] Q1x;
    if (Q1y != NULL)
        delete[] Q1y;
    if (Q1z != NULL)
        delete[] Q1z;
    Q1x = NULL;
    Q1y = NULL;
    Q1z = NULL;

    if (Q2x != NULL)
        delete[] Q2x;
    if (Q2y != NULL)
        delete[] Q2y;
    if (Q2z != NULL)
        delete[] Q2z;
    Q2x = NULL;
    Q2y = NULL;
    Q2z = NULL;

    if (Q3x != NULL)
        delete[] Q3x;
    if (Q3y != NULL)
        delete[] Q3y;
    if (Q3z != NULL)
        delete[] Q3z;
    Q3x = NULL;
    Q3y = NULL;
    Q3z = NULL;
    
    if (Q4x != NULL)
        delete[] Q4x;
    if (Q4y != NULL)
        delete[] Q4y;
    if (Q4z != NULL)
        delete[] Q4z;
    Q4x = NULL;
    Q4y = NULL;
    Q4z = NULL;

    if (Q5x != NULL)
        delete[] Q5x;
    if (Q5y != NULL)
        delete[] Q5y;
    if (Q5z != NULL)
        delete[] Q5z;
    Q5x = NULL;
    Q5y = NULL;
    Q5z = NULL;

    // Limiters
    if (Limiter_Phi1 != NULL)
        delete[] Limiter_Phi1;
    if (Limiter_Phi2 != NULL)
        delete[] Limiter_Phi2;
    if (Limiter_Phi3 != NULL)
        delete[] Limiter_Phi3;
    if (Limiter_Phi4 != NULL)
        delete[] Limiter_Phi4;
    if (Limiter_Phi5 != NULL)
        delete[] Limiter_Phi5;
    Limiter_Phi1 = NULL;
    Limiter_Phi2 = NULL;
    Limiter_Phi3 = NULL;
    Limiter_Phi4 = NULL;
    Limiter_Phi5 = NULL;
    
    // Residuals (Convective and Dissipative)
    if (Res1 != NULL)
        delete[] Res1;
    if (Res2 != NULL)
        delete[] Res2;
    if (Res3 != NULL)
        delete[] Res3;
    if (Res4 != NULL)
        delete[] Res4;
    if (Res5 != NULL)
        delete[] Res5;
    Res1 = NULL;
    Res2 = NULL;
    Res3 = NULL;
    Res4 = NULL;
    Res5 = NULL;

    if (Res1_Diss != NULL)
        delete[] Res1_Diss;
    if (Res2_Diss != NULL)
        delete[] Res2_Diss;
    if (Res3_Diss != NULL)
        delete[] Res3_Diss;
    if (Res4_Diss != NULL)
        delete[] Res4_Diss;
    if (Res5_Diss != NULL)
        delete[] Res5_Diss;
    Res1_Diss = NULL;
    Res2_Diss = NULL;
    Res3_Diss = NULL;
    Res4_Diss = NULL;
    Res5_Diss = NULL;
    
    // Local Time
    if (DeltaT != NULL)
        delete[] DeltaT;
    DeltaT = NULL;
    
    // Precondition Variable
    if (PrecondSigma != NULL)
        delete[] PrecondSigma;
    PrecondSigma = NULL;
    
    // Finalize Gradient Infrastructure
    if (SolverOrder == SOLVER_ORDER_SECOND)
        Gradient_Finalize();
    
    // Check if Implicit Scheme
    if (SolverScheme == SOLVER_SCHEME_IMPLICIT)
        Delete_CRS_SolverBlockMatrix();
    
    // RMS Writer Finalize
    RMS_Writer_Finalize();
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Set Initial Conditions
//------------------------------------------------------------------------------
void Solver_Set_Initial_Conditions(void) {

    // Allocate Memory to Store Conservative or Primitive Variables
    Q1 = new double[nNode + nBNode];
    Q2 = new double[nNode + nBNode];
    Q3 = new double[nNode + nBNode];
    Q4 = new double[nNode + nBNode];
    Q5 = new double[nNode + nBNode];

    // Compute Far Field Conditions
    ComputeFarFieldCondition(0);
    
    // Initialize the variable with reference conditions
    switch (VariableType) {
        case VARIABLE_CONSERVATIVE:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_Rho*Inf_U;
                Q3[i] = Inf_Rho*Inf_V;
                Q4[i] = Inf_Rho*Inf_W;
                Q5[i] = Inf_Rho*Inf_Et;
            }
            break;
        case VARIABLE_PRIMITIVE_PUT:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Pressure;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Temperature;
            }
            break;
        case VARIABLE_PRIMITIVE_RUP:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Pressure;
            }
            break;
         case VARIABLE_PRIMITIVE_RUT:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Temperature;
            }
            break;
        default:
            error("Solver_Set_Initial_Conditions: Undefined Variable Type - %d", VariableType);
            break;
    }
    

    // Allocate Memory to Store Conservative or Primitive Variables Gradients
    if (SolverOrder == SOLVER_ORDER_SECOND) {
        // Initialize Gradient Infrastructure
        Gradient_Init();

        Q1x = new double[nNode];
        Q1y = new double[nNode];
        Q1z = new double[nNode];
        Gradient_Add_Function(Q1, Q1x, Q1y, Q1z, nNode);
        
        Q2x = new double[nNode];
        Q2y = new double[nNode];
        Q2z = new double[nNode];
        Gradient_Add_Function(Q2, Q2x, Q2y, Q2z, nNode);

        Q3x = new double[nNode];
        Q3y = new double[nNode];
        Q3z = new double[nNode];
        Gradient_Add_Function(Q3, Q3x, Q3y, Q3z, nNode);

        Q4x = new double[nNode];
        Q4y = new double[nNode];
        Q4z = new double[nNode];
        Gradient_Add_Function(Q4, Q4x, Q4y, Q4z, nNode);

        Q5x = new double[nNode];
        Q5y = new double[nNode];
        Q5z = new double[nNode];
        Gradient_Add_Function(Q5, Q5x, Q5y, Q5z, nNode);
        
        // Check if Limiter Memory is Required
        if (LimiterMethod != LIMITER_METHOD_NONE) {
            Limiter_Phi1 = new double[nNode];
            Limiter_Phi2 = new double[nNode];
            Limiter_Phi3 = new double[nNode];
            Limiter_Phi4 = new double[nNode];
            Limiter_Phi5 = new double[nNode];
        }
    }
    
    // Allocate  Memory to Store Residuals (Convective and Dissipative)
    Res1 = new double[nNode];
    Res2 = new double[nNode];
    Res3 = new double[nNode];
    Res4 = new double[nNode];
    Res5 = new double[nNode];
    for (int i = 0; i < nNode; i++) {
        Res1[i] = 0.0;
        Res2[i] = 0.0;
        Res3[i] = 0.0;
        Res4[i] = 0.0;
        Res5[i] = 0.0;
    }
    
    Res1_Diss = new double[nNode];
    Res2_Diss = new double[nNode];
    Res3_Diss = new double[nNode];
    Res4_Diss = new double[nNode];
    Res5_Diss = new double[nNode];
    for (int i = 0; i < nNode; i++) {
        Res1_Diss[i] = 0.0;
        Res2_Diss[i] = 0.0;
        Res3_Diss[i] = 0.0;
        Res4_Diss[i] = 0.0;
        Res5_Diss[i] = 0.0;
    }
    
    // Min and Max EigenValue Value
    MinEigenLamda1   = DBL_MAX;
    MaxEigenLamda1   = DBL_MIN;
    MinEigenLamda4   = DBL_MAX;
    MaxEigenLamda4   = DBL_MIN;
    MinEigenLamda5   = DBL_MAX;
    MaxEigenLamda5   = DBL_MIN;
    
    // Min and Max Precondition Variable
    if (PrecondMethod != PRECOND_METHOD_NONE) {
        PrecondSigma    = new double[nNode];
        for (int i = 0; i < nNode; i++)
            PrecondSigma[i] = 0.0;
        MinPrecondSigma = DBL_MAX;
        MaxPrecondSigma = DBL_MIN;
    } else {
        MinPrecondSigma = 1.0;
        MaxPrecondSigma = 1.0;
    }
    
    // Allocate Memory to Store Time Step
    DeltaT = new double[nNode];
    // Initialize
    for (int i = 0; i < nNode; i++)
        DeltaT[i] = 0.0;
    MinDeltaT = DBL_MAX;
    MaxDeltaT = DBL_MIN;
    
    // Set Boundary Conditions
    Initialize_Boundary_Condition();
    
    // Check if Implicit Scheme
    if (SolverScheme == SOLVER_SCHEME_IMPLICIT)
        Create_CRS_SolverBlockMatrix();
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
int Solver(void) {
    int rvalue = EXIT_FAILURE;
    
    // Initialize the Solver Scheme Data Structure
    switch (FluxScheme) {
        case FLUX_SCHEME_ROE: // Roe
            Roe_Init_New();
            break;
        case FLUX_SCHEME_HLLC: // HLLC
            HLLC_Init();
            break;
        case FLUX_SCHEME_AUSM: // AUSM
            AUSM_Init();
            break;
        case FLUX_SCHEME_VANLEER: // Van Leer
            VanLeer_Init();
            break;
        case FLUX_SCHEME_LDFSS: // LDFSS
            LDFSS_Init();
            break;
        case FLUX_SCHEME_OSHER: // Osher
            Osher_Init();
            break;
        case FLUX_SCHEME_STEGERWARMING: // Steger Warming
            StegerWarming_Init();
            break;
        default:
            error("Solver: Invalid Flux Scheme - %d", FluxScheme);
            break;
    }
    
    // Check if Solution Restart is Requested
    if (RestartInput)
        Restart_Reader(RestartInputFilename);
    
    // Check if Residual Smoothing is Requested
    if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
        Residual_Smoothing_Init();
    
    // Select the Solver Method Type
    switch (SolverMethod) {
        case SOLVER_METHOD_STEADY:
            switch (SolverScheme) {
                case SOLVER_SCHEME_EXPLICIT:
                    rvalue = Solver_Steady_Explicit();
                    break;
                case SOLVER_SCHEME_IMPLICIT:
                    rvalue = Solver_Steady_Implicit();
                    break;
                default:
                    error("Solver: Invalid Solver Scheme - %d -1", SolverScheme);
                    break;
            }
            break;
        case SOLVER_METHOD_UNSTEADY:
            switch (SolverScheme) {
                case SOLVER_SCHEME_EXPLICIT:
                    rvalue = Solver_Unsteady_Explicit();
                    break;
                case SOLVER_SCHEME_IMPLICIT:
                    rvalue = Solver_Unsteady_Implicit();
                    break;
                default:
                    error("Solver: Invalid Solver Scheme - %d -2", SolverScheme);
                    break;
            }
            break;
        default:
            error("Solver: Invalid Solver Method - %d", SolverMethod);
            break;
    }
    
    // Finalize Residual Smoothing Data Structure
    if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
        Residual_Smoothing_Finalize();
    
    // Finalize the Solver Scheme Data Structure
    switch (FluxScheme) {
        case FLUX_SCHEME_ROE: // Roe
            Roe_Finalize_New();
            break;
        case FLUX_SCHEME_HLLC: // HLLC
            HLLC_Finalize();
            break;
        case FLUX_SCHEME_AUSM: // AUSM
            AUSM_Finalize();
            break;
        case FLUX_SCHEME_VANLEER: // Van Leer
            VanLeer_Finalize();
            break;
        case FLUX_SCHEME_LDFSS: // LDFSS
            LDFSS_Finalize();
            break;
        case FLUX_SCHEME_OSHER: // Osher
            Osher_Finalize();
            break;
        case FLUX_SCHEME_STEGERWARMING: // Steger Warming
            StegerWarming_Finalize();
            break;
        default:
            error("Solver: Invalid Flux Scheme - %d", FluxScheme);
            break;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Create_CRS_SolverBlockMatrix(void) {
    int i, j, jstart, jend, k, ksave;
    int degree, index, min, minsave;
    
    SolverBlockMatrix.nROW       = nNode;
    SolverBlockMatrix.nCOL       = nNode;
    SolverBlockMatrix.Block_nRow = NEQUATIONS;
    SolverBlockMatrix.Block_nCol = NEQUATIONS;
    
    // Allocate Memory of IA Array to Store Start and End Location of Row
    SolverBlockMatrix.IA = new int[nNode+1];
    
    // Start Filling the Row Location and
    // Get the No of Non Zero Entries for SolverBlockMatrix
    SolverBlockMatrix.DIM = 0;
    SolverBlockMatrix.IA[0] = 0;
    for (i = 0; i < nNode; i++) {
        degree = crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i];
        SolverBlockMatrix.DIM += degree + 1;
        SolverBlockMatrix.IA[i+1] = SolverBlockMatrix.IA[i] + degree + 1;
    }
    
    // Allocate Memory to IAU Array to store Diagonal Location
    SolverBlockMatrix.IAU = new int[nNode];
    
    // Allocate Memory to JA Array to store location of Non Zero Entries
    SolverBlockMatrix.JA = new int[SolverBlockMatrix.DIM];
    // Get the values for JA
    for (i = 0; i < nNode; i++) {
        // Make the first start entry as main node id
        index = SolverBlockMatrix.IA[i];
        SolverBlockMatrix.JA[index] = i;
        
        // Now add the nodes connected to main node
        for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
            index++;
            SolverBlockMatrix.JA[index] = crs_JA_Node2Node[j];
        }
    }
    
    /* Now Sort JA and Find IAU */
    // This step is necessary for computation of Transpose in Design Code - Adjoint
    for (i = 0; i < nNode; i++) {
        jstart = SolverBlockMatrix.IA[i];
        jend   = SolverBlockMatrix.IA[i + 1];
        for (j = jstart; j < jend; j++) {
            min = SolverBlockMatrix.JA[j];
            minsave = SolverBlockMatrix.JA[j];
            ksave = j;
            for (k = j + 1; k < jend; k++) {
                if (SolverBlockMatrix.JA[k] < min) {
                    min = SolverBlockMatrix.JA[k];
                    ksave = k;
                }
            }
            SolverBlockMatrix.JA[j] = min;
            SolverBlockMatrix.JA[ksave] = minsave;
            if (SolverBlockMatrix.JA[j] == i)
                SolverBlockMatrix.IAU[i] = j;
        }
    }
    
    // Allocate Memory for CRS Matrix
    SolverBlockMatrix.A = (double ***) malloc (SolverBlockMatrix.DIM*sizeof(double**));
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        SolverBlockMatrix.A[i] = NULL;
        SolverBlockMatrix.A[i] = (double **) malloc (SolverBlockMatrix.Block_nRow*sizeof(double*));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            SolverBlockMatrix.A[i][j] = NULL;
            SolverBlockMatrix.A[i][j] = (double *) malloc (SolverBlockMatrix.Block_nCol*sizeof(double));
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Allocate Memory of RHS
    SolverBlockMatrix.B = (double **) malloc (SolverBlockMatrix.nROW*sizeof(double*));
    for (i = 0; i < SolverBlockMatrix.nROW; i++) {
        SolverBlockMatrix.B[i] = NULL;
        SolverBlockMatrix.B[i] = (double *) malloc (SolverBlockMatrix.Block_nRow*sizeof(double));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++)
            SolverBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    SolverBlockMatrix.X = (double **) malloc (SolverBlockMatrix.nROW*sizeof(double*));
    for (i = 0; i < SolverBlockMatrix.nROW; i++) {
        SolverBlockMatrix.X[i] = NULL;
        SolverBlockMatrix.X[i] = (double *) malloc (SolverBlockMatrix.Block_nRow*sizeof(double));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++)
            SolverBlockMatrix.X[i][j] = 0.0;
    }
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Delete_CRS_SolverBlockMatrix(void) {
    int i, j;
    
    if (SolverBlockMatrix.A != NULL) {
        for (i = 0; i < SolverBlockMatrix.DIM; i++) {
            if (SolverBlockMatrix.A[i] != NULL) {
                for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                    if (SolverBlockMatrix.A[i][j] != NULL)
                        free(SolverBlockMatrix.A[i][j]);
                }
                free(SolverBlockMatrix.A[i]);
            }
        }
        free(SolverBlockMatrix.A);
    }
    
    if (SolverBlockMatrix.B != NULL) {
        for (i = 0; i < SolverBlockMatrix.nROW; i++) {
            if (SolverBlockMatrix.B[i] != NULL)
                free(SolverBlockMatrix.B[i]);
        }
        free(SolverBlockMatrix.B);
    }
    
    if (SolverBlockMatrix.X != NULL) {
        for (i = 0; i < SolverBlockMatrix.nROW; i++) {
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
}

