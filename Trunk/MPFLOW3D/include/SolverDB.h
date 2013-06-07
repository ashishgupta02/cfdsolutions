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
    double dvDPDRho;
    double dvDPDT;
    double dvDRhoDT;
    double dvDEDRho_P;
    double dvDEDRho_T;
    double dvDEDP_Rho;
    double dvDEDT_Rho;
// Constructors, Destructor and Operators
public:
    CSolverDBCore();
    virtual ~CSolverDBCore();
};

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
class CSolverDB {
// Public Data - Attributes
public:
    int    ivNodeID;
    int    ivBndFlag;
// Protected Data - Attributes
protected:
    double  daQ[NEQUATIONS];
    double daQx[NEQUATIONS];
    double daQy[NEQUATIONS];
    double daQz[NEQUATIONS];
    double Residual_Conv[NEQUATIONS];
    double Residual_Diss[NEQUATIONS];
    double Limiter_Phi[NEQUATIONS];
    double DeltaT;
    CSolverDBCore Properties;
    
// Constructors, Destructor and Operators
public:
    CSolverDB();
    virtual ~CSolverDB();
// Public - Functions
public:
    void Set_BoundaryFlag(int BndFlag);
    void Compute_Properties();
    void Compute_Extended_Properties();
    void Get_Properties(double &Rho, double &RhoL, double &RhoV, double &Pressure,
            double &Temperature, double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
            double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy);
    void Get_Extended_Properties(double &Rho, double &RhoL, double &RhoV, double &Pressure,
            double &Temperature, double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
            double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy,
            double &DPDRho, double &DPDT, double &DRhoDT, double &DEDRho_P, 
            double &DEDRho_T, double &DEDP_Rho, double &DEDT_Rho);
    void Get_Recomputed_Properties(double *dpVariableIn, double &Rho, double &RhoL, double &RhoV,
            double &Pressure, double &Temperature,
            double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
            double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy);
    void Get_Recomputed_Extended_Properties(double *dpVariableIn, double &Rho, double &RhoL, double &RhoV,
            double &Pressure, double &Temperature,
            double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
            double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy,
            double &DPDRho, double &DPDT, double &DRhoDT, double &DEDRho_P, 
            double &DEDRho_T, double &DEDP_Rho, double &DEDT_Rho);
    
    void   Get_Density_All(double *dpDensity);
    double Get_Density();
    double Get_Density_Liquid();
    double Get_Density_Vapor();
    double Get_Velocity_U();
    double Get_Velocity_V();
    double Get_Velocity_W();
    double Get_Pressure();
    double Get_Temperature();
    double Get_Mach();
    double Get_TotalEnergy();
    double Get_SpeedSound();
    
    void   Get_Recomputed_Density_All(double *dpVariableIn, double *dpDensity);
    double Get_Recomputed_Density(double *dpVariableIn);
    double Get_Recomputed_Density_Liquid(double *dpVariableIn);
    double Get_Recomputed_Density_Vapor(double *dpVariableIn);
    double Get_Recomputed_Pressure(double *dpVariableIn);
    double Get_Recomputed_Temperature(double *dpVariableIn);
    double Get_Recomputed_Mach(double *dpVariableIn);
    double Get_Recomputed_TotalEnergy(double *dpVariableIn);
    double Get_Recomputed_SpeedSound(double *dpVariableIn);
};

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
class CSolver {
protected:
    int ivNodeDBDim;
public:
    //! Constructors and Destructors
    CSolver();
    ~CSolver();
    //! Attributes
    CSolverDB *CpNodeDB;
    bool FluxRecomputeFlag;
    // Operations
    void Create(int ivNodes);
    void Set_BoundaryFlag(void);
    void Compute_Properties(int ivNode);
    void Compute_Extended_Properties(int ivNode);
    void Compute_Bulk_Properties(int ivNodeFrom, int ivNodeTo);
    void Compute_Bulk_Extended_Properties(int ivNodeFrom, int ivNodeTo);
    void Compute_Boundary_Nodes_Properties(void);
    void Compute_Boundary_Nodes_Extended_Properties(void);
    void Smooth_Solution(void);
    void Smooth_Stagnation_Solution(void);
};

#endif	/* _SOLVER_DB_H */

