/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "Solver.h"

// Linear Solver Params
int    Order;
int    NIteration;
int    InnerNIteration;
double Relaxation;

// Reference Conditions
double Ref_Mach;
double Ref_Alpha;
double Ref_Pressure;

// Free Stream Conditions
double Inf_Rho;
double Inf_U;
double Inf_V;
double Inf_W;
double Inf_Et;
double Inf_Pressure;

// Constants
double Gamma;

// CFL Conditons
int    CFL_Ramp;
double CFL_MAX;
double CFL_MIN;
double CFL;

// Local Time
double *DeltaT;

// Conservative Variables
double *Q1;
double *Q2;
double *Q3;
double *Q4;
double *Q5;

// Residuals
double *Res1;
double *Res2;
double *Res3;
double *Res4;
double *Res5;

//RMS
double RMS[5];
double RMS_Res;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Init(void) {
    // Linear Solver Params
    Order           = 0;
    NIteration      = 0;
    InnerNIteration = 0;
    Relaxation      = 0.0;

    // Reference Conditions
    Ref_Mach        = 0.0;
    Ref_Alpha       = 0.0;
    Ref_Pressure    = 0.0;

    // Free Stream Conditions
    Inf_Rho         = 0.0;
    Inf_U           = 0.0;
    Inf_V           = 0.0;
    Inf_W           = 0.0;
    Inf_Et          = 0.0;
    Inf_Pressure    = 0.0;
    
    // Constants
    Gamma           = 0.0;

    // CFL Conditons
    CFL_Ramp        = 0;
    CFL_MAX         = 0.0;
    CFL_MIN         = 0.0;
    CFL             = 0.0;

    // Local Time
    DeltaT          = NULL;

    // Conservative Variables
    Q1              = NULL;
    Q2              = NULL;
    Q3              = NULL;
    Q4              = NULL;
    Q5              = NULL;

    // Residuals
    Res1            = NULL;
    Res2            = NULL;
    Res3            = NULL;
    Res4            = NULL;
    Res5            = NULL;

    // RMS
    RMS[0]          = 0.0;
    RMS[1]          = 0.0;
    RMS[2]          = 0.0;
    RMS[3]          = 0.0;
    RMS[4]          = 0.0;
    RMS_Res         = 0.0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Finalize(void) {
    // Conservative Variables
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

    // Residuals
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

    // Local Time
    if (DeltaT != NULL)
        delete[] DeltaT;
    DeltaT = NULL;
    
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Read_Params(const char *filename) {
    FILE *fp;
    int bdim = 256;
    char buff[256];
    char *dummy;

    if ((fp = fopen(filename, "r")) == NULL)
        error("Solver_Read_Params: Unable to Read Parameter File - %s", filename);
    
    // Number of Boundaries
    int nb;
    // Read Number of Boundaries
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &nb);

    // Allocate memory to store boundary type and read
    bndType = new int[nBC];
    for (int i = 0; i < nb; i++) {
        dummy = fgets(buff, bdim, fp);
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &bndType[i]);
    }

    // Get the Order
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Order);

    // Get Number of Outer Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &NIteration);

    // Get Number of Inner Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &InnerNIteration);
    
    // Get the Relaxation Factor
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Relaxation);

    // Read Gamma
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Gamma);

    // Read Mach
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Mach);

    // Read Pressure
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Pressure);

    // Read Alpha
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Alpha);

    // Read CFL Ramp
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &CFL_Ramp);

    // Read CFL MIN
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MIN);

    // Read CFL MAX
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MAX);

    CFL             = CFL_MIN;
    Ref_Alpha       = Ref_Alpha * M_PI / 180.0;
    
    // Close file
    fclose(fp);
}

//------------------------------------------------------------------------------
//! Set Initial Conditions
//------------------------------------------------------------------------------
void Solver_Set_Initial_Conditions() {

    // Allocate Memory to Store Conservative Variables
    Q1 = new double[nNode + nBNode];
    Q2 = new double[nNode + nBNode];
    Q3 = new double[nNode + nBNode];
    Q4 = new double[nNode + nBNode];
    Q5 = new double[nNode + nBNode];

    // Intialize the variable with reference conditions
    for (int i = 0; i < (nNode + nBNode); i++) {
        Q1[i] = 1.0;
        Q2[i] = Q1[i]*Ref_Mach*cos(Ref_Alpha);
        Q3[i] = Q1[i]*Ref_Mach*sin(Ref_Alpha);
        Q4[i] = 0.0;
        Q5[i] = 1.0 / (Gamma * (Gamma - 1.0)) + 0.5 *(Q2[i]*Q2[i] + Q3[i]*Q3[i] + Q4[i]*Q4[i])/Q1[i];
    }
    
    // Allocate  Memory to Store Residuals
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
    
    // Allocate Memory to Store Time Step
    DeltaT = new double[nNode];
    // Initialize
    for (int i = 0; i < nNode; i++)
        DeltaT[i] = 0.0;

    // Set Boundary Conditions
    Initialize_Boundary_Condition();

    // Set the Free Stream Conditions
    Inf_Rho      = 1.0;
    Inf_U        = Inf_Rho*Ref_Mach*cos(Ref_Alpha);
    Inf_V        = Inf_Rho*Ref_Mach*sin(Ref_Alpha);
    Inf_W        = 0.0;
    Inf_Et       = 1.0 / (Gamma * (Gamma - 1.0)) + 0.5 *(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W)/Inf_Rho;
    Inf_Pressure = (Gamma - 1.0) * Inf_Rho * (Inf_Et - 0.5 * (Inf_U * Inf_U + Inf_V * Inf_V + Inf_W * Inf_W));
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Solve() {
    printf("=============================================================================\n");
    printf("-----------------------------------------------------------------------------\n");
    printf(" Iter        RMS_RHO    RMS_RHOU    RMS_RHOV    RMS_RHOW    RMS_E     RMS_RES\n");
    printf("-----------------------------------------------------------------------------\n");
    for (int iter = 0; iter < NIteration; iter++) {
        // Reset Residuals and DeltaT
        for (int i = 0; i < nNode; i++) {
            Res1[i]   = 0.0;
            Res2[i]   = 0.0;
            Res3[i]   = 0.0;
            Res4[i]   = 0.0;
            Res5[i]   = 0.0;
            DeltaT[i] = 0.0;
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition();

        // Compute Residuals
        Compute_Residual();

        // Compute Local Time Stepping
        Compute_DeltaT(iter);

        // Update Conservative Variables
        for (int i = 0; i < nNode; i++) {
            Q1[i] -= (DeltaT[i] / cVolume[i]) * Res1[i];
            Q2[i] -= (DeltaT[i] / cVolume[i]) * Res2[i];
            Q3[i] -= (DeltaT[i] / cVolume[i]) * Res3[i];
            Q4[i] -= (DeltaT[i] / cVolume[i]) * Res4[i];
            Q5[i] -= (DeltaT[i] / cVolume[i]) * Res5[i];
        }
        
        // Compute RMS
        RMS[0] = RMS[1] = RMS[2] = RMS[3] = RMS[4] = 0.0;
        for (int i = 0; i < nNode; i++) {
            RMS[0] += Res1[i]*Res1[i];
            RMS[1] += Res2[i]*Res2[i];
            RMS[2] += Res3[i]*Res3[i];
            RMS[3] += Res4[i]*Res4[i];
            RMS[4] += Res5[i]*Res5[i];
        }
        RMS_Res = RMS[0] + RMS[1] + RMS[2] + RMS[3] + RMS[4];
        RMS_Res = sqrt(RMS_Res/(5.0 * (double)nNode));
        RMS[0] = sqrt(RMS[0]/(double)nNode);
        RMS[1] = sqrt(RMS[1]/(double)nNode);
        RMS[2] = sqrt(RMS[2]/(double)nNode);
        RMS[3] = sqrt(RMS[3]/(double)nNode);
        RMS[4] = sqrt(RMS[4]/(double)nNode);

        printf("%5d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",
                iter+1, RMS[0], RMS[1], RMS[2], RMS[3], RMS[4], RMS_Res);

        if (isnan(RMS_Res)) {
            info("Solve: NAN Encountered ! - Abort");
            iter = NIteration + 1;
        }

        if ((RMS_Res < (DBL_EPSILON*10.0))|| ((iter+1) == NIteration))
            iter = NIteration + 1;
    }
}

