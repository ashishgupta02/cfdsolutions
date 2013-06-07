/*******************************************************************************
 * File:        Material.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _MATERIAL_H
#define	_MATERIAL_H

// Material and Non Dimensionalization Type
extern int    MaterialType;
extern int    MaterialCompType;
extern int    MaterialEOS_IO_Type;
extern char   MaterialName[256];
extern double MaterialPropertyLimits[9];

// Dimensional Properties (SI Units)
extern double Ref_Rho;
extern double Ref_Pressure;
extern double Ref_Temperature;
extern double Ref_Velocity;
extern double Ref_TotalEnergy;
extern double Ref_Length;
extern double Ref_Time;
extern double Ref_Mach;

// Other Reference Properties
extern double Gamma;

// Far Field Conditions
extern double Inf_Rho;
extern double Inf_U;
extern double Inf_V;
extern double Inf_W;
extern double Inf_Et;
extern double Inf_Pressure;
extern double Inf_Temperature;
extern double Inf_SpeedSound;
extern double Inf_Mach;
extern double Inf_Alpha;
extern double Inf_Beta;
// Infinity Mach Ramping
extern int    Inf_Mach_Ramp;
extern double Inf_Mach_MAX;
extern double Inf_Mach_MIN;

// Gauge Variables
extern double Gauge_Pressure;
extern double Outflow_Pressure;

// Material User Property Limits
extern double Limit_Min_Pressure;
extern double Limit_Max_Pressure;
extern double Limit_Min_Rho;
extern double Limit_Max_Rho;
extern double Limit_Min_Temperature;
extern double Limit_Max_Temperature;

// Material Property Functions
void Material_Init(void);
void Material_Set_Properties(void);

// Compute Far Field Condition with Mach Ramping
void Material_Set_InfinityCondition(int Iteration);

// Bound the Solution in Limits
void Material_Limit_Solution(void);

// Equation of State
void Material_Get_Properties(double *dpVariableIn, double *dpPropertyOut);

void Material_Get_Extended_Properties(double *dpVariableIn, double *dpPropertyOut);

void Material_Get_ControlVolume_Properties(double *dpVariableIn,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &total_energy, double &total_enthalpy);

void Material_Get_Face_Properties(double *dpVariableIn, double nx, double ny, double nz,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &ubar, double &total_energy, double &total_enthalpy);

void Material_Get_RUH_To_Q(double Rho, double Velocity_U, double Velocity_V, 
        double Velocity_W, double Enthalpy, double *Q);

void   Material_Get_Density_All(double *dpVariableIn, double *dpDensity);
double Material_Get_Density(double *dpVariableIn);
double Material_Get_Density_Liquid(double *dpVariableIn);
double Material_Get_Density_Vapor(double *dpVariableIn);
double Material_Get_Pressure(double *dpVariableIn);
double Material_Get_Temperature(double *dpVariableIn);
double Material_Get_Mach(double *dpVariableIn);
double Material_Get_TotalEnergy(double *dpVariableIn);
double Material_Get_SpeedSound(double *dpVariableIn);
double Material_Get_Quality(double *dpVariableIn);
double Material_Get_DH_SpeedSound(double dvDensity, double dvEnthalpy);
double Material_Get_DH_Pressure(double dvDensity, double dvEnthalpy);
double Material_Get_DH_Temperature(double dvDensity, double dvEnthalpy);

// Variable Transformations
void Material_Get_Transformation_Matrix(int ivNodeID, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);
void Material_Get_Transformation_Matrix(double *dpPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);
void Material_Get_Transformation_Matrix_Properties(double *dpExPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);

#endif	/* _MATERIAL_H */

