/*******************************************************************************
 * File:        SolverDB.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#include "SolverDB.h"
#include "Commons.h"
#include "Solver.h"
#include "Material.h"
#include "Utils.h"
#include "EOS.h"
#include "Material.h"

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
CSolverDBCore::CSolverDBCore() {
    // Initialize
    dvRho           = 0.0;
    dvRhoL          = 0.0;
    dvRhoV          = 0.0;
    dvPressure      = 0.0;
    dvTemperature   = 0.0;
    dvVelocity_U    = 0.0;
    dvVelocity_V    = 0.0;
    dvVelocity_W    = 0.0;
    dvQ2            = 0.0;
    dvSpeedSound    = 0.0;
    dvMach          = 0.0;
    dvTotalEnergy   = 0.0;
    dvTotalEnthalpy = 0.0;
    dvDPDRho        = 0.0;
    dvDPDT          = 0.0;
    dvDRhoDT        = 0.0;
    dvDEDRho_P      = 0.0;
    dvDEDRho_T      = 0.0;
    dvDEDP_Rho      = 0.0;
    dvDEDT_Rho      = 0.0;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
CSolverDBCore::~CSolverDBCore() {
    // Nothing to do
}

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
CSolverDB::CSolverDB() {
    for (int i = 0; i < NEQUATIONS; i++) {
        daQ[i]           = 0.0;
        daQx[i]          = 0.0;
        daQy[i]          = 0.0;
        daQz[i]          = 0.0;
        Residual_Conv[i] = 0.0;
        Residual_Diss[i] = 0.0;
        Limiter_Phi[i]   = 0.0;
    }
    ivNodeID  = -1;
    ivBndFlag = -1;
    DeltaT    = 0.0;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
CSolverDB::~CSolverDB() {
    // Nothing to do
}

//------------------------------------------------------------------------------
//! Set the Boundary Flag
//------------------------------------------------------------------------------
void CSolverDB::Set_BoundaryFlag(int BndFlag) {
    ivBndFlag = BndFlag;
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties
//------------------------------------------------------------------------------
void CSolverDB::Compute_Properties() {
    double tmp;
    double daProperty[EOS_THERM_DIM];
    // Get the Q's
    daQ[0] = Q1[ivNodeID];
    daQ[1] = Q2[ivNodeID];
    daQ[2] = Q3[ivNodeID];
    daQ[3] = Q4[ivNodeID];
    daQ[4] = Q5[ivNodeID];
    Material_Get_Properties(daQ, daProperty);

    // Update the values
    Properties.dvRho           = daProperty[ 0];
    Properties.dvRhoL          = daProperty[ 1];
    Properties.dvRhoV          = daProperty[ 2];
    Properties.dvPressure      = daProperty[ 3];
    Properties.dvTemperature   = daProperty[ 4];
    Properties.dvVelocity_U    = daProperty[ 5];
    Properties.dvVelocity_V    = daProperty[ 6];
    Properties.dvVelocity_W    = daProperty[ 7];
    Properties.dvQ2            = daProperty[ 8];
    Properties.dvSpeedSound    = daProperty[ 9];
    Properties.dvMach          = daProperty[10];
    Properties.dvTotalEnthalpy = daProperty[14];
    Properties.dvTotalEnergy   = daProperty[15];

//    // Fix the Pressure Bounce
//    tmp = Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W;
//    tmp = Inf_Rho*tmp;
//    Properties.dvPressure = MIN(Properties.dvPressure, Inf_Pressure+tmp);
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties
//------------------------------------------------------------------------------
void CSolverDB::Compute_Extended_Properties() {
    double tmp;
    double daProperty[EOS_EX_THERM_DIM];
    // Get the Q's
    daQ[0] = Q1[ivNodeID];
    daQ[1] = Q2[ivNodeID];
    daQ[2] = Q3[ivNodeID];
    daQ[3] = Q4[ivNodeID];
    daQ[4] = Q5[ivNodeID];
    Material_Get_Extended_Properties(daQ, daProperty);

    // Update the values
    Properties.dvRho           = daProperty[ 0];
    Properties.dvRhoL          = daProperty[ 1];
    Properties.dvRhoV          = daProperty[ 2];
    Properties.dvPressure      = daProperty[ 3];
    Properties.dvTemperature   = daProperty[ 4];
    Properties.dvVelocity_U    = daProperty[ 5];
    Properties.dvVelocity_V    = daProperty[ 6];
    Properties.dvVelocity_W    = daProperty[ 7];
    Properties.dvQ2            = daProperty[ 8];
    Properties.dvSpeedSound    = daProperty[ 9];
    Properties.dvMach          = daProperty[10];
    Properties.dvTotalEnthalpy = daProperty[14];
    Properties.dvTotalEnergy   = daProperty[15];
    Properties.dvDPDRho        = daProperty[18];
    Properties.dvDPDT          = daProperty[19];
    Properties.dvDRhoDT        = daProperty[20];
    Properties.dvDEDT_Rho      = daProperty[31];
    Properties.dvDEDRho_T      = daProperty[33];
    Properties.dvDEDRho_P      = daProperty[34];
    Properties.dvDEDP_Rho      = daProperty[36];

//    // Fix the Pressure Bounce
//    tmp = Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W;
//    tmp = Inf_Rho*tmp;
//    Properties.dvPressure = MIN(Properties.dvPressure, Inf_Pressure+tmp);
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties: First Order
//------------------------------------------------------------------------------
void CSolverDB::Get_Properties(double &Rho, double &RhoL, double &RhoV, double &Pressure,
        double &Temperature, double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy) {

    Rho           = Properties.dvRho;
    RhoL          = Properties.dvRhoL;
    RhoV          = Properties.dvRhoV;
    Pressure      = Properties.dvPressure;
    Temperature   = Properties.dvTemperature;
    Velocity_U    = Properties.dvVelocity_U;
    Velocity_V    = Properties.dvVelocity_V;
    Velocity_W    = Properties.dvVelocity_W;
    Q2            = Properties.dvQ2;
    SpeedSound    = Properties.dvSpeedSound;
    Mach          = Properties.dvMach;
    TotalEnergy   = Properties.dvTotalEnergy;
    TotalEnthalpy = Properties.dvTotalEnthalpy;
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties: First Order
//------------------------------------------------------------------------------
void CSolverDB::Get_Extended_Properties(double &Rho, double &RhoL, double &RhoV, double &Pressure,
        double &Temperature, double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy,
        double &DPDRho, double &DPDT, double &DRhoDT, double &DEDRho_P,
        double &DEDRho_T, double &DEDP_Rho, double &DEDT_Rho) {

    Rho           = Properties.dvRho;
    RhoL          = Properties.dvRhoL;
    RhoV          = Properties.dvRhoV;
    Pressure      = Properties.dvPressure;
    Temperature   = Properties.dvTemperature;
    Velocity_U    = Properties.dvVelocity_U;
    Velocity_V    = Properties.dvVelocity_V;
    Velocity_W    = Properties.dvVelocity_W;
    Q2            = Properties.dvQ2;
    SpeedSound    = Properties.dvSpeedSound;
    Mach          = Properties.dvMach;
    TotalEnergy   = Properties.dvTotalEnergy;
    TotalEnthalpy = Properties.dvTotalEnthalpy;
    DPDRho        = Properties.dvDPDRho;
    DPDT          = Properties.dvDPDT;
    DRhoDT        = Properties.dvDRhoDT;
    DEDRho_P      = Properties.dvDEDRho_P;
    DEDRho_T      = Properties.dvDEDRho_T;
    DEDP_Rho      = Properties.dvDEDP_Rho;
    DEDT_Rho      = Properties.dvDEDT_Rho;
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties: Q's are Input
//------------------------------------------------------------------------------
void CSolverDB::Get_Recomputed_Properties(double *dpVariableIn, double &Rho, double &RhoL, double &RhoV,
    double &Pressure, double &Temperature,
    double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
    double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy) {

    double tmp;
    double daProperty[EOS_THERM_DIM];

    // Recompute the Material Properties
    Material_Get_Properties(dpVariableIn, daProperty);

    // Update the values
    Rho           = daProperty[ 0];
    RhoL          = daProperty[ 1];
    RhoV          = daProperty[ 2];
    Pressure      = daProperty[ 3];
    Temperature   = daProperty[ 4];
    Velocity_U    = daProperty[ 5];
    Velocity_V    = daProperty[ 6];
    Velocity_W    = daProperty[ 7];
    Q2            = daProperty[ 8];
    SpeedSound    = daProperty[ 9];
    Mach          = daProperty[10];
    TotalEnthalpy = daProperty[14];
    TotalEnergy   = daProperty[15];

//    // Fix the Pressure Bounce
//    tmp = Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W;
//    tmp = Inf_Rho*tmp;
//    Pressure = MIN(Pressure, Inf_Pressure+tmp);
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties: Q's are Input
//------------------------------------------------------------------------------
void CSolverDB::Get_Recomputed_Extended_Properties(double *dpVariableIn, double &Rho, double &RhoL, double &RhoV,
    double &Pressure, double &Temperature,
    double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
    double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy,
    double &DPDRho, double &DPDT, double &DRhoDT, double &DEDRho_P,
    double &DEDRho_T, double &DEDP_Rho, double &DEDT_Rho) {

    double tmp;
    double daProperty[EOS_EX_THERM_DIM];

    // Recompute the Material Extended Properties
    Material_Get_Extended_Properties(dpVariableIn, daProperty);

    // Update the values
    Rho           = daProperty[ 0];
    RhoL          = daProperty[ 1];
    RhoV          = daProperty[ 2];
    Pressure      = daProperty[ 3];
    Temperature   = daProperty[ 4];
    Velocity_U    = daProperty[ 5];
    Velocity_V    = daProperty[ 6];
    Velocity_W    = daProperty[ 7];
    Q2            = daProperty[ 8];
    SpeedSound    = daProperty[ 9];
    Mach          = daProperty[10];
    TotalEnthalpy = daProperty[14];
    TotalEnergy   = daProperty[15];
    DPDRho        = daProperty[18];
    DPDT          = daProperty[19];
    DRhoDT        = daProperty[20];
    DEDT_Rho      = daProperty[31];
    DEDRho_T      = daProperty[33];
    DEDRho_P      = daProperty[34];
    DEDP_Rho      = daProperty[36];

//    // Fix the Pressure Bounce
//    tmp = Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W;
//    tmp = Inf_Rho*tmp;
//    Pressure = MIN(Pressure, Inf_Pressure+tmp);
}

//------------------------------------------------------------------------------
//! Get All Densities
//------------------------------------------------------------------------------
void CSolverDB::Get_Density_All(double* dpDensity) {
    dpDensity[0] = Properties.dvRho;
    dpDensity[1] = Properties.dvRhoL;
    dpDensity[2] = Properties.dvRhoV;
}

//------------------------------------------------------------------------------
//! Get Fluid Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Density() {
    return Properties.dvRho;
}

//------------------------------------------------------------------------------
//! Get Fluid Liquid Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Density_Liquid() {
    return Properties.dvRhoL;
}

//------------------------------------------------------------------------------
//! Get Fluid Vapor Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Density_Vapor() {
    return Properties.dvRhoV;
}

//------------------------------------------------------------------------------
//! Get Fluid Velocity_U
//------------------------------------------------------------------------------
double CSolverDB::Get_Velocity_U() {
    return Properties.dvVelocity_U;
}

//------------------------------------------------------------------------------
//! Get Fluid Velocity_V
//------------------------------------------------------------------------------
double CSolverDB::Get_Velocity_V() {
    return Properties.dvVelocity_V;
}

//------------------------------------------------------------------------------
//! Get Fluid Velocity_W
//------------------------------------------------------------------------------
double CSolverDB::Get_Velocity_W() {
    return Properties.dvVelocity_W;
}

//------------------------------------------------------------------------------
//! Get Fluid Pressure
//------------------------------------------------------------------------------
double CSolverDB::Get_Pressure() {
    return Properties.dvPressure;
}

//------------------------------------------------------------------------------
//! Get Fluid Temperature
//------------------------------------------------------------------------------
double CSolverDB::Get_Temperature() {
    return Properties.dvTemperature;
}

//------------------------------------------------------------------------------
//! Get Fluid Mach
//------------------------------------------------------------------------------
double CSolverDB::Get_Mach() {
    return Properties.dvMach;
}

//------------------------------------------------------------------------------
//! Get Fluid Total Energy
//------------------------------------------------------------------------------
double CSolverDB::Get_TotalEnergy() {
    return Properties.dvTotalEnergy;
}

//------------------------------------------------------------------------------
//! Get Fluid Speed of Sound
//------------------------------------------------------------------------------
double CSolverDB::Get_SpeedSound() {
    return Properties.dvSpeedSound;
}

//------------------------------------------------------------------------------
//! Get All Densities
//------------------------------------------------------------------------------
void CSolverDB::Get_Recomputed_Density_All(double *dpVariableIn, double* dpDensity) {
    // Recompute the Material All Densities
    Material_Get_Density_All(dpVariableIn, dpDensity);
}

//------------------------------------------------------------------------------
//! Get Fluid Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Density(double *dpVariableIn) {
    // Recompute the Material Density
    return Material_Get_Density(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Liquid Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Density_Liquid(double *dpVariableIn) {
    // Recompute the Material Liquid Density
    return Material_Get_Density_Liquid(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Vapor Density
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Density_Vapor(double *dpVariableIn) {
    // Recompute the Material Vapor Density
    return Material_Get_Density_Vapor(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Pressure
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Pressure(double *dpVariableIn) {
    // Recompute the Material Pressure
    return Material_Get_Pressure(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Temperature
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Temperature(double *dpVariableIn) {
    // Recompute the Material Temperature
    return Material_Get_Temperature(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Mach
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_Mach(double *dpVariableIn) {
    // Recompute the Material Mach
    return Material_Get_Mach(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Total Energy
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_TotalEnergy(double *dpVariableIn) {
    // Recompute the Material Total Energy
    return Material_Get_TotalEnergy(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Get Fluid Speed of Sound
//------------------------------------------------------------------------------
double CSolverDB::Get_Recomputed_SpeedSound(double *dpVariableIn) {
    // Recompute the Material Speed of Sound
    return Material_Get_SpeedSound(dpVariableIn);
}

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
CSolver::CSolver() {
    ivNodeDBDim       = 0;
    CpNodeDB          = NULL;
    FluxRecomputeFlag = FALSE;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
CSolver::~CSolver() {
    // Nothing to do
    if (CpNodeDB != NULL)
        delete[] CpNodeDB;
}

//------------------------------------------------------------------------------
//! Create
//------------------------------------------------------------------------------
void CSolver::Create(int ivNodes) {
    // Create and Initialize the Node Data Base
    ivNodeDBDim = ivNodes;
    CpNodeDB    = new CSolverDB[ivNodes];

    for (int i = 0; i < ivNodes; i++) {
        CpNodeDB[i].ivNodeID  = i;
        CpNodeDB[i].ivBndFlag = -2; // Default All to Interior Node
    }
}

//------------------------------------------------------------------------------
//! Set the Boundary Flags 
//  Flag = -2: Interior Physical Node
//  Flag = -1: Boundary Physical Node
//  Flag = BC_TYPE Ghost Nodes
//------------------------------------------------------------------------------
void CSolver::Set_BoundaryFlag () {
    int node_L, node_R;

    // Initialize all Nodes to Interior
    for (int i = 0; i < ivNodeDBDim; i++)
        CpNodeDB[i].ivBndFlag = -2;

    // Loop over all Boundary Edges and Set the Appropriate Flags
    for (int iBEdge = 0; iBEdge < nBEdge; iBEdge++) {
        node_L = bndEdge[iBEdge].node[0]; // Physical Node
        node_R = bndEdge[iBEdge].node[1]; // Ghost Node
        CpNodeDB[node_L].ivBndFlag = -1;
        CpNodeDB[node_R].ivBndFlag = bndEdge[iBEdge].type;
    }
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties of Specific Node
//------------------------------------------------------------------------------
void CSolver::Compute_Properties (int ivNode) {
    // Compute the Property of Particular Node
    CpNodeDB[ivNode].Compute_Properties();
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties of Specific Node
//------------------------------------------------------------------------------
void CSolver::Compute_Extended_Properties (int ivNode) {
    // Compute the Property of Particular Node
    CpNodeDB[ivNode].Compute_Extended_Properties();
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties of Range of Nodes
//------------------------------------------------------------------------------
void CSolver::Compute_Bulk_Properties (int ivNodeFrom, int ivNodeTo) {
    // Compute the Property of Range of Nodes
    for (int ivNode = ivNodeFrom; ivNode < ivNodeTo; ivNode++)
        CpNodeDB[ivNode].Compute_Properties();
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties of Range of Nodes
//------------------------------------------------------------------------------
void CSolver::Compute_Bulk_Extended_Properties (int ivNodeFrom, int ivNodeTo) {
    // Compute the Extended Property of Range of Nodes
    for (int ivNode = ivNodeFrom; ivNode < ivNodeTo; ivNode++)
        CpNodeDB[ivNode].Compute_Extended_Properties();
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties of Physical Boundary Nodes
//------------------------------------------------------------------------------
void CSolver::Compute_Boundary_Nodes_Properties () {
    // Compute the Property
    for (int ivNode = 0; ivNode < ivNodeDBDim; ivNode++) {
        // Only Physical Boundary Nodes
        if (CpNodeDB[ivNode].ivBndFlag == -1)
            CpNodeDB[ivNode].Compute_Properties();
    }
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Extended Properties of Physical Boundary Nodes
//------------------------------------------------------------------------------
void CSolver::Compute_Boundary_Nodes_Extended_Properties () {
    // Compute the Extended Property
    for (int ivNode = 0; ivNode < ivNodeDBDim; ivNode++) {
        // Only Physical Boundary Nodes
        if (CpNodeDB[ivNode].ivBndFlag == -1)
            CpNodeDB[ivNode].Compute_Extended_Properties();
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CSolver::Smooth_Solution() {
    double Coeff1, Coeff2;
    int nid;
    double *Q_Old = NULL;
    double *Q_New = NULL;
    double *Qswap = NULL;

    Q_Old = new double [NEQUATIONS*nNode];
    Q_New = new double [NEQUATIONS*nNode];

    // Get the Old Q's
    for (int inode = 0; inode < nNode; inode++) {
        Q_Old[NEQUATIONS*inode + 0] = Q1[inode];
        Q_Old[NEQUATIONS*inode + 1] = Q2[inode];
        Q_Old[NEQUATIONS*inode + 2] = Q3[inode];
        Q_Old[NEQUATIONS*inode + 3] = Q4[inode];
        Q_Old[NEQUATIONS*inode + 4] = Q5[inode];
    }

    Coeff1 = SolutionSmoothRelaxation;
    for (int iSmooth = 0; iSmooth < SolutionSmoothNIteration; iSmooth++) {
        for (int inode = 0; inode < nNode; inode++) {
            Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
            // First Term
            for (int j = 0; j < NEQUATIONS; j++)
                Q_New[NEQUATIONS*inode + j] = Coeff1*Q_Old[NEQUATIONS*inode + j];
            // Second Term
            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                nid = crs_JA_Node2Node[i];
                for (int j = 0; j < NEQUATIONS; j++)
                    Q_New[NEQUATIONS*inode + j] += Coeff2*Q_Old[NEQUATIONS*nid + j];
            }
        }
        // Swap the New to Old
        Qswap = Q_New;
        Q_New = Q_Old;
        Q_Old = Qswap; // Old is Now New solution
        Qswap = NULL;
    }

    // Smooth Solution
    for (int inode = 0; inode < nNode; inode++) {
        // Dissipative Term
        Q1[inode] = Q_Old[NEQUATIONS*inode + 0];
        Q2[inode] = Q_Old[NEQUATIONS*inode + 1];
        Q3[inode] = Q_Old[NEQUATIONS*inode + 2];
        Q4[inode] = Q_Old[NEQUATIONS*inode + 3];
        Q5[inode] = Q_Old[NEQUATIONS*inode + 4];
    }

    // Free Memory
    delete[] Q_Old;
    delete[] Q_New;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CSolver::Smooth_Solution_Gradient() {
    double Coeff1, Coeff2;
    int nid;
    double *Qx_Old = NULL;
    double *Qx_New = NULL;
    double *Qy_Old = NULL;
    double *Qy_New = NULL;
    double *Qz_Old = NULL;
    double *Qz_New = NULL;
    double *Qswap = NULL;

    Qx_Old = new double [NEQUATIONS*nNode];
    Qx_New = new double [NEQUATIONS*nNode];
    Qy_Old = new double [NEQUATIONS*nNode];
    Qy_New = new double [NEQUATIONS*nNode];
    Qz_Old = new double [NEQUATIONS*nNode];
    Qz_New = new double [NEQUATIONS*nNode];
    
    // Get the Old Q's
    for (int inode = 0; inode < nNode; inode++) {
        // X
        Qx_Old[NEQUATIONS*inode + 0] = Q1x[inode];
        Qx_Old[NEQUATIONS*inode + 1] = Q2x[inode];
        Qx_Old[NEQUATIONS*inode + 2] = Q3x[inode];
        Qx_Old[NEQUATIONS*inode + 3] = Q4x[inode];
        Qx_Old[NEQUATIONS*inode + 4] = Q5x[inode];

        // Y
        Qy_Old[NEQUATIONS*inode + 0] = Q1y[inode];
        Qy_Old[NEQUATIONS*inode + 1] = Q2y[inode];
        Qy_Old[NEQUATIONS*inode + 2] = Q3y[inode];
        Qy_Old[NEQUATIONS*inode + 3] = Q4y[inode];
        Qy_Old[NEQUATIONS*inode + 4] = Q5y[inode];

        // Z
        Qz_Old[NEQUATIONS*inode + 0] = Q1z[inode];
        Qz_Old[NEQUATIONS*inode + 1] = Q2z[inode];
        Qz_Old[NEQUATIONS*inode + 2] = Q3z[inode];
        Qz_Old[NEQUATIONS*inode + 3] = Q4z[inode];
        Qz_Old[NEQUATIONS*inode + 4] = Q5z[inode];
    }
    
    Coeff1 = SolutionSmoothRelaxation;
    for (int iSmooth = 0; iSmooth < SolutionSmoothNIteration; iSmooth++) {
        for (int inode = 0; inode < nNode; inode++) {
            Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
            // First Term
            for (int j = 0; j < NEQUATIONS; j++) {
                Qx_New[NEQUATIONS*inode + j] = Coeff1*Qx_Old[NEQUATIONS*inode + j];
                Qy_New[NEQUATIONS*inode + j] = Coeff1*Qy_Old[NEQUATIONS*inode + j];
                Qz_New[NEQUATIONS*inode + j] = Coeff1*Qz_Old[NEQUATIONS*inode + j];
            }
            // Second Term
            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                nid = crs_JA_Node2Node[i];
                for (int j = 0; j < NEQUATIONS; j++) {
                    Qx_New[NEQUATIONS*inode + j] += Coeff2*Qx_Old[NEQUATIONS*nid + j];
                    Qy_New[NEQUATIONS*inode + j] += Coeff2*Qy_Old[NEQUATIONS*nid + j];
                    Qz_New[NEQUATIONS*inode + j] += Coeff2*Qz_Old[NEQUATIONS*nid + j];
                }
            }
        }
        // Swap the New to Old
        // X
        Qswap  = Qx_New;
        Qx_New = Qx_Old;
        Qx_Old = Qswap; // Old is Now New solution
        Qswap  = NULL;

        // Y
        Qswap  = Qy_New;
        Qy_New = Qy_Old;
        Qy_Old = Qswap; // Old is Now New solution
        Qswap  = NULL;

        // Z
        Qswap  = Qz_New;
        Qz_New = Qz_Old;
        Qz_Old = Qswap; // Old is Now New solution
        Qswap  = NULL;
    }

    // Smooth Solution
    for (int inode = 0; inode < nNode; inode++) {
        // X
        Q1x[inode] = Qx_Old[NEQUATIONS*inode + 0];
        Q2x[inode] = Qx_Old[NEQUATIONS*inode + 1];
        Q3x[inode] = Qx_Old[NEQUATIONS*inode + 2];
        Q4x[inode] = Qx_Old[NEQUATIONS*inode + 3];
        Q5x[inode] = Qx_Old[NEQUATIONS*inode + 4];

        // Y
        Q1y[inode] = Qy_Old[NEQUATIONS*inode + 0];
        Q2y[inode] = Qy_Old[NEQUATIONS*inode + 1];
        Q3y[inode] = Qy_Old[NEQUATIONS*inode + 2];
        Q4y[inode] = Qy_Old[NEQUATIONS*inode + 3];
        Q5y[inode] = Qy_Old[NEQUATIONS*inode + 4];

        // Z
        Q1z[inode] = Qz_Old[NEQUATIONS*inode + 0];
        Q2z[inode] = Qz_Old[NEQUATIONS*inode + 1];
        Q3z[inode] = Qz_Old[NEQUATIONS*inode + 2];
        Q4z[inode] = Qz_Old[NEQUATIONS*inode + 3];
        Q5z[inode] = Qz_Old[NEQUATIONS*inode + 4];
    }
    
    // Free Memory
    delete[] Qx_Old;
    delete[] Qx_New;
    delete[] Qy_Old;
    delete[] Qy_New;
    delete[] Qz_Old;
    delete[] Qz_New;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CSolver::Smooth_Stagnation_Solution() {
    double Coeff1, Coeff2;
    int nid;
    int    *SmoothNode = NULL;
    double *Q_Old = NULL;
    double *Q_New = NULL;
    double *Qswap = NULL;

    SmoothNode = new int [nNode];
    Q_Old = new double [NEQUATIONS*nNode];
    Q_New = new double [NEQUATIONS*nNode];

    // Get the Old Q's and Set the Flags
    Coeff2 = PrecondGlobalMach/10.0;
    for (int inode = 0; inode < nNode; inode++) {
        // Compute the Properties
        CpNodeDB[inode].Compute_Properties();
        // Get the Mach number
        Coeff1 = CpNodeDB[inode].Get_Mach();
        if (Coeff1 <= Coeff2)
            SmoothNode[inode] = 1;
        else
            SmoothNode[inode] = 0;
        Q_Old[NEQUATIONS*inode + 0] = Q1[inode];
        Q_Old[NEQUATIONS*inode + 1] = Q2[inode];
        Q_Old[NEQUATIONS*inode + 2] = Q3[inode];
        Q_Old[NEQUATIONS*inode + 3] = Q4[inode];
        Q_Old[NEQUATIONS*inode + 4] = Q5[inode];
    }

    Coeff1 = SolutionSmoothRelaxation;
    for (int iSmooth = 0; iSmooth < SolutionSmoothNIteration; iSmooth++) {
        for (int inode = 0; inode < nNode; inode++) {
            if (SmoothNode[inode] == 1) {
                Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                // First Term
                for (int j = 0; j < NEQUATIONS; j++)
                    Q_New[NEQUATIONS*inode + j] = Coeff1*Q_Old[NEQUATIONS*inode + j];
                // Second Term
                for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                    nid = crs_JA_Node2Node[i];
                    for (int j = 0; j < NEQUATIONS; j++)
                        Q_New[NEQUATIONS*inode + j] += Coeff2*Q_Old[NEQUATIONS*nid + j];
                }
            } else {
                for (int j = 0; j < NEQUATIONS; j++)
                    Q_New[NEQUATIONS*inode + j] = Q_Old[NEQUATIONS*inode + j];
            }
        }
        // Swap the New to Old
        Qswap = Q_New;
        Q_New = Q_Old;
        Q_Old = Qswap; // Old is Now New solution
        Qswap = NULL;
    }

    // Smooth Solution
    for (int inode = 0; inode < nNode; inode++) {
        // Dissipative Term
        Q1[inode] = Q_Old[NEQUATIONS*inode + 0];
//        Q2[inode] = Q_Old[NEQUATIONS*inode + 1];
//        Q3[inode] = Q_Old[NEQUATIONS*inode + 2];
//        Q4[inode] = Q_Old[NEQUATIONS*inode + 3];
        Q5[inode] = Q_Old[NEQUATIONS*inode + 4];
    }

    // Free Memory
    delete[] Q_Old;
    delete[] Q_New;
    delete[] SmoothNode;
}
