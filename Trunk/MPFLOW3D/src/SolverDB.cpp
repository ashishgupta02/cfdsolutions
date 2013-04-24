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
    // Get the Q's 
    daQ[0] = Q1[ivNodeID];
    daQ[1] = Q2[ivNodeID];
    daQ[2] = Q3[ivNodeID];
    daQ[3] = Q4[ivNodeID];
    daQ[4] = Q5[ivNodeID];
    Material_Get_ControlVolume_Properties(daQ, PropertiesFirst.dvRho, PropertiesFirst.dvPressure, PropertiesFirst.dvTemperature, 
            PropertiesFirst.dvVelocity_U, PropertiesFirst.dvVelocity_V, PropertiesFirst.dvVelocity_W, PropertiesFirst.dvQ2, 
            PropertiesFirst.dvSpeedSound, PropertiesFirst.dvMach, PropertiesFirst.dvTotalEnergy, PropertiesFirst.dvTotalEnthalpy);
    // Check if Second Order compute is required
    if (ivNodeID < nNode) {
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            
        }
    }
}

//------------------------------------------------------------------------------
//! Compute the Thermodynamic Properties
//------------------------------------------------------------------------------
void CSolverDB::Get_Properties(int ivOrder, double &Rho, double &Pressure, double &Temperature,
        double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy) {
    
    switch (ivOrder) {
        case SOLVER_ORDER_FIRST:
            Rho           = PropertiesFirst.dvRho;
            Pressure      = PropertiesFirst.dvPressure;
            Temperature   = PropertiesFirst.dvTemperature;
            Velocity_U    = PropertiesFirst.dvVelocity_U;
            Velocity_V    = PropertiesFirst.dvVelocity_V;
            Velocity_W    = PropertiesFirst.dvVelocity_W;
            Q2            = PropertiesFirst.dvQ2;
            SpeedSound    = PropertiesFirst.dvSpeedSound;
            Mach          = PropertiesFirst.dvMach;
            TotalEnergy   = PropertiesFirst.dvTotalEnergy;
            TotalEnthalpy = PropertiesFirst.dvTotalEnthalpy;
            break;
        case SOLVER_ORDER_SECOND:
            Rho           = PropertiesSecond.dvRho;
            Pressure      = PropertiesSecond.dvPressure;
            Temperature   = PropertiesSecond.dvTemperature;
            Velocity_U    = PropertiesSecond.dvVelocity_U;
            Velocity_V    = PropertiesSecond.dvVelocity_V;
            Velocity_W    = PropertiesSecond.dvVelocity_W;
            Q2            = PropertiesSecond.dvQ2;
            SpeedSound    = PropertiesSecond.dvSpeedSound;
            Mach          = PropertiesSecond.dvMach;
            TotalEnergy   = PropertiesSecond.dvTotalEnergy;
            TotalEnthalpy = PropertiesSecond.dvTotalEnthalpy;
            break;
        default:
            error("CSolverDB::Get_Properties: Invalid Solver Order - %d", ivOrder);
            break;
    }
}

