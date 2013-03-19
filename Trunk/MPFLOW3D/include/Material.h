/*******************************************************************************
 * File:        Material.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _MATERIAL_H
#define	_MATERIAL_H

#define UNIV_GAS_CONST          8314.47215

// Material and Non Dimensionalization Type
extern int    NonDimensionalMethod;
extern int    MaterialType;

// Dimensional Properties
extern double Ref_R;
extern double Ref_MolecularWeight;
extern double Ref_Length;
extern double Ref_Area;
extern double Ref_Time;
extern double Ref_Rho;
extern double Ref_Pressure;
extern double Ref_Temperature;
extern double Ref_Velocity;
extern double Ref_Mach;

// Other Reference Properties
extern double Gamma;
extern double Ref_Alpha;
extern double Ref_Beta;

// Non Dimensional Properties
extern double NonDim_R;
extern double NonDim_Time;
extern double NonDim_Rho;
extern double NonDim_Pressure;
extern double NonDim_Temperature;

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

// Gauge Variables
extern double Gauge_Pressure;
extern double Outflow_Pressure;

// Material Property Functions
void Material_Init(void);
void Material_Set_Properties(void);

// Compute Far Field Condition with Mach Ramping
void ComputeFarFieldCondition(int Iteration);

// Equation of State
void Compute_EOS_Variables_Face(double *Q, double nx, double ny, double nz,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &ubar, double &total_energy, double &total_enthalpy);
void Compute_EOS_Variables_ControlVolume(double *Q,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &total_energy, double &total_enthalpy);
void Compute_Q_From_EOS_Variables(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double Total_Enthalpy, double *Q);
double Get_Rho(double Pressure, double Temperature);
double Get_Pressure(double Rho, double Temperature);
double Get_Pressure(double *Q);
double Get_Temperature(double Rho, double Pressure);
double Get_Temperature(double *Q);
double Get_SpeedSound(double Rho, double Pressure);
double Get_SpeedSound(double Temperature);
double Get_SpeedSound(double Velocity_U, double Velocity_V, double Velocity_W, double Total_Enthalpy);
double Get_SpeedSoundSquare(double Velocity_U, double Velocity_V, double Velocity_W, double Total_Enthalpy);
double Get_TotalEnergy(double Rho, double Pressure, double Velocity_U, double Velocity_V, double Velocity_W);
double Get_Mach(double *Q);

// Variable Transformations
void Compute_Transformation_Conservative_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUP_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_Conservative_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUT_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_Conservative_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUS_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_Conservative_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUT_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUP_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUT_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUP_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUS_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUP_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUT_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUT_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUS_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUT_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUT_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_PUS_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_RUT_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);
void Compute_Transformation_Matrix(int VarType1, int VarType2, double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix);

#endif	/* _MATERIAL_H */

