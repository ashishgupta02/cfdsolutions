/*******************************************************************************
 * File:        SolverDB.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _SOLVER_DB_H
#define	_SOLVER_DB_H

#include "SolverParameters.h"

using namespace std;

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
class CSolverDBCore {
// Public Data - Attributes
public:
    double dvRho;
    double dvRhoL;
    double dvRhoV;
    double dvPressure;
    double dvTemperature;
    double dvVelocity_U;
    double dvVelocity_V;
    double dvVelocity_W;
    double dvQ2;
    double dvSpeedSound;
    double dvMach;
    double dvTotalEnergy; 
    double dvTotalEnthalpy;
// Constructors, Destructor and Operators
public:
    CSolverDBCore();
    virtual ~CSolverDBCore();
};

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
class CSolverDB {
// Protected Data - Attributes
protected:
    int    ivNodeID;
    int    ivBndFlag;
    double  daQ[NEQUATIONS];
    double daQx[NEQUATIONS];
    double daQy[NEQUATIONS];
    double daQz[NEQUATIONS];
    double Residual_Conv[NEQUATIONS];
    double Residual_Diss[NEQUATIONS];
    double Limiter_Phi[NEQUATIONS];
    double DeltaT;
    
    CSolverDBCore PropertiesFirst;
    CSolverDBCore PropertiesSecond;
// Constructors, Destructor and Operators
public:
    CSolverDB();
    virtual ~CSolverDB();
// Public - Functions
public:
    void Set_BoundaryFlag(int BndFlag);
    void Compute_Properties();
    void Get_Properties(int ivOrder, double &Rho, double &Pressure, double &Temperature,
            double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
            double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy);
};


#endif	/* _SOLVER_DB_H */

