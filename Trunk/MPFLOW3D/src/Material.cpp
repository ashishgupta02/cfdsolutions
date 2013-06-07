/*******************************************************************************
 * File:        Material.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Material.h"
#include "SolverParameters.h"
#include "Trim_Utils.h"
#include "EOS.h"
#include "Solver.h"
#include "Commons.h"

// Material Variables
int    MaterialType;
int    MaterialCompType;
int    MaterialEOS_IO_Type;
char   MaterialName[256];
double MaterialPropertyLimits[9];

// Dimensional Properties (SI Units)
double Ref_Rho;
double Ref_Pressure;
double Ref_Temperature;
double Ref_Velocity;
double Ref_TotalEnergy;
double Ref_Length;
double Ref_Time;
double Ref_Mach;

// Other Reference Properties
double Gamma;

// Far Field Conditions
double Inf_Rho;
double Inf_U;
double Inf_V;
double Inf_W;
double Inf_Et;
double Inf_Pressure;
double Inf_Temperature;
double Inf_SpeedSound;
double Inf_Mach;
double Inf_Alpha;
double Inf_Beta;
// Infinity Mach Ramping
int    Inf_Mach_Ramp;
double Inf_Mach_MAX;
double Inf_Mach_MIN;

// Gauge Variables
double Gauge_Pressure;
double Outflow_Pressure;

// Material User Limits
double Limit_Min_Pressure;
double Limit_Max_Pressure;
double Limit_Min_Rho;
double Limit_Max_Rho;
double Limit_Min_Temperature;
double Limit_Max_Temperature;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Init(void) {
    // Initialize
    MaterialType        = MATERIAL_TYPE_NONE;
    MaterialCompType    = MATERIAL_COMP_TYPE_NONE;
    MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_NONE;
    str_blank(MaterialName);
    for (int i = 0; i < 9; i++)
        MaterialPropertyLimits[i] = 0.0;
    
    // Dimensional Reference Properties
    Ref_Rho             = 0.0;
    Ref_Pressure        = 0.0;
    Ref_Temperature     = 0.0;
    Ref_Velocity        = 0.0;
    Ref_TotalEnergy     = 0.0;
    Ref_Length          = 0.0;
    Ref_Time            = 0.0;
    Ref_Mach            = 0.0;
    
    // Other Reference Properties
    Gamma               = 0.0;
    
    // Far Field Conditions
    Inf_Rho             = 0.0;
    Inf_U               = 0.0;
    Inf_V               = 0.0;
    Inf_W               = 0.0;
    Inf_Et              = 0.0;
    Inf_Pressure        = 0.0;
    Inf_Temperature     = 0.0;
    Inf_SpeedSound      = 0.0;
    Inf_Mach            = 0.0;
    Inf_Alpha           = 0.0;
    Inf_Beta            = 0.0;
    // Mach Ramping
    Inf_Mach_Ramp       = 0;
    Inf_Mach_MAX        = 0.0;
    Inf_Mach_MIN        = 0.0;
    
    // Gauge Variables
    Gauge_Pressure      = 0.0;
    Outflow_Pressure    = 0.0;
    
    // Material User Limits
    Limit_Min_Pressure    = 0.0;
    Limit_Max_Pressure    = 0.0;
    Limit_Min_Rho         = 0.0;
    Limit_Max_Rho         = 0.0;
    Limit_Min_Temperature = 0.0;
    Limit_Max_Temperature = 0.0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Set_Properties(void) {
    double Q[NEQUATIONS];
    double Density[3];
    double Prop[EOS_THERM_DIM];
    
    Gamma = 1.4;
    // Initialize and Set the EOS based on Material Type
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEALGAS:
            EOS_Init(EOS_MODEL_IDEALGAS);
            break;
        case MATERIAL_TYPE_NIST:
            EOS_Init(EOS_MODEL_NIST);
            break;
        default:
            error("Material_Set_Properties:1: Undefined Material Type - %d", MaterialType);
            break;
    }
    EOS_Set(MaterialName);
    
    // Set the EOS Reference Properties
    EOS_Set_Reference_Properties(Ref_Pressure, Ref_Temperature, Ref_Length);
    EOS_Print_Reference_Properties();
    // Get the Reference Properties
    EOS_Get_Reference_Properties(Prop);
    
    // Update the Final Reference Properties (SI Units)
    Ref_Rho             = Prop[ 0];
    Ref_Pressure        = Prop[ 1];
    Ref_Temperature     = Prop[ 2];
    Ref_Velocity        = Prop[ 3];
    Ref_Length          = Prop[ 4];
    Ref_Mach            = Prop[ 6];
    Ref_Time            = Prop[ 7];
    Ref_TotalEnergy     = Prop[12];
    
    // Compute the Material Property Limits for NIST Material
    if (MaterialType == MATERIAL_TYPE_NIST) {
        EOS_Get_Fluid_Information(Prop);
        MaterialPropertyLimits[0] = Prop[ 2]; // Critical Temperature
        MaterialPropertyLimits[1] = Prop[ 3]; // Critical Pressure
        MaterialPropertyLimits[2] = Prop[ 5]; // Critical Density
        MaterialPropertyLimits[3] = Prop[ 6]; // Triple Point Temperature
        MaterialPropertyLimits[4] = Prop[ 7]; // Normal Boiling Point Temperature
        MaterialPropertyLimits[5] = Prop[10]; // Minimum Temperature 
        MaterialPropertyLimits[6] = Prop[11]; // Maximum Temperature
        MaterialPropertyLimits[7] = Prop[12]; // Maximum Pressure
        MaterialPropertyLimits[8] = Prop[13]; // Maximum Density
        
        // Check and set the valid User Property Limits
        // Pressure Limits
        if ((Limit_Min_Pressure > 0.0) || (Limit_Max_Pressure > 0.0)) {
            if (Limit_Max_Pressure < MaterialPropertyLimits[7])
                MaterialPropertyLimits[7] = Limit_Max_Pressure;
        }
        // Density Limits
        if((Limit_Min_Rho > 0.0) || (Limit_Max_Rho > 0.0)) {
            if (Limit_Max_Rho < MaterialPropertyLimits[8])
                MaterialPropertyLimits[8] = Limit_Max_Rho;
        }
        // Temperature Limits
        if ((Limit_Min_Temperature > 0.0) || (Limit_Max_Temperature > 0.0)) {
            // Maximum Temperature
            if ((Limit_Max_Temperature < MaterialPropertyLimits[6]) &&
                    (Limit_Max_Temperature > MaterialPropertyLimits[5]))
                MaterialPropertyLimits[6] = Limit_Max_Temperature;
            // Minimum Temperature
            if ((Limit_Min_Temperature > MaterialPropertyLimits[5]) &&
                    (Limit_Min_Temperature < MaterialPropertyLimits[6]))
                MaterialPropertyLimits[5] = Limit_Min_Temperature;
        }
    }
    
    // Compute the Infinity Conditions (SI Units)
    EOS_Get_PT_Density(EOS_DIMENSIONAL_IO_D_D, Inf_Pressure, Inf_Temperature, Density);
    Inf_Rho = Density[0];
    Q[0]    = Inf_Rho;
    Q[1]    = Q[2] = Q[3] = 0.0;
    Q[4]    = Inf_Temperature;
    EOS_Get_Properties(EOS_DIMENSIONAL_IO_D_D, EOS_VARIABLE_RUT, Q, Prop);
    Inf_Alpha       = Inf_Alpha * M_PI / 180.0;
    Inf_Beta        = Inf_Beta * M_PI / 180.0;
    Inf_SpeedSound  = Prop[9];
    Inf_U           = Inf_Mach*Inf_SpeedSound*cos(Inf_Alpha)*cos(Inf_Beta);
    Inf_V           = Inf_Mach*Inf_SpeedSound*sin(Inf_Alpha);
    Inf_W           = Inf_Mach*Inf_SpeedSound*sin(Inf_Beta)*cos(Inf_Alpha);
    Inf_Et          = Prop[13] + 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
    
    // Print the Infinity Conditions
    info("Infinity Conditions:");
    printf("-----------------------------------------------------------------------------\n");
    info("Density_Inf -------------------------------: %15.6f kg/m^3", Inf_Rho);
    info("Pressure_Inf ------------------------------: %15.6f kPa",    Inf_Pressure/1000.0);
    info("Temperature_Inf ---------------------------: %15.6f K",      Inf_Temperature);
    info("Velocity_U_Inf ----------------------------: %15.6f m/s",    Inf_U);
    info("Velocity_V_Inf ----------------------------: %15.6f m/s",    Inf_V);
    info("Velocity_W_Inf ----------------------------: %15.6f m/s",    Inf_W);
    info("SpeedSound_Inf ----------------------------: %15.6f m/s",    Inf_SpeedSound);
    info("Mach_Inf ----------------------------------: %15.6f ",       Inf_Mach);
    info("TotalEnergy_Inf ---------------------------: %15.6f J/kg",   Inf_Et);
    printf("=============================================================================\n");
    
    // Check the Material Computation Type
    switch (MaterialCompType) {
        case MATERIAL_COMP_TYPE_NDIM:
            //--START Non-Dimensionalization 
            // Infinity Conditions
            Inf_Rho         /= Ref_Rho;
            Inf_Pressure    /= Ref_Pressure;
            Inf_Temperature /= Ref_Temperature;
            Inf_U           /= Ref_Velocity;
            Inf_V           /= Ref_Velocity;
            Inf_W           /= Ref_Velocity;
            Inf_SpeedSound  /= Ref_Velocity;
            Inf_Et          /= Ref_TotalEnergy;
            // Pressure Conditions
            Gauge_Pressure   /= Ref_Pressure;
            Outflow_Pressure /= Ref_Pressure;
            // Set Computations I/O Mode
            MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_ND_ND;
            // Print the Infinity Conditions
            info("Non-Dimensional - Infinity and Other Conditions:");
            printf("-----------------------------------------------------------------------------\n");
            info("Density_Inf -------------------------------: %15.6f", Inf_Rho);
            info("Pressure_Inf ------------------------------: %15.6f", Inf_Pressure);
            info("Temperature_Inf ---------------------------: %15.6f", Inf_Temperature);
            info("Velocity_U_Inf ----------------------------: %15.6f", Inf_U);
            info("Velocity_V_Inf ----------------------------: %15.6f", Inf_V);
            info("Velocity_W_Inf ----------------------------: %15.6f", Inf_W);
            info("SpeedSound_Inf ----------------------------: %15.6f", Inf_SpeedSound);
            info("Mach_Inf ----------------------------------: %15.6f", Inf_Mach);
            info("TotalEnergy_Inf ---------------------------: %15.6f", Inf_Et);
            info("Gauge_Pressure ----------------------------: %15.6f", Gauge_Pressure);
            info("Outflow_Pressure --------------------------: %15.6f", Outflow_Pressure);
            printf("=============================================================================\n");
            // Non-Dimensionalize the Material Property Limits for NIST Material
            if (MaterialType == MATERIAL_TYPE_NIST) {
                MaterialPropertyLimits[0] /= Ref_Temperature; // Critical Temperature
                MaterialPropertyLimits[1] /= Ref_Pressure; // Critical Pressure
                MaterialPropertyLimits[2] /= Ref_Rho; // Critical Density
                MaterialPropertyLimits[3] /= Ref_Temperature; // Triple Point Temperature
                MaterialPropertyLimits[4] /= Ref_Temperature; // Normal Boiling Point Temperature
                MaterialPropertyLimits[5] /= Ref_Temperature; // Minimum Temperature 
                MaterialPropertyLimits[6] /= Ref_Temperature; // Maximum Temperature
                MaterialPropertyLimits[7] /= Ref_Pressure; // Maximum Pressure
                MaterialPropertyLimits[8] /= Ref_Rho; // Maximum Density;
            }
            
            // Non-Dimensionalize the User Property Limits
            Limit_Min_Pressure    /= Ref_Pressure;
            Limit_Max_Pressure    /= Ref_Pressure;
            Limit_Min_Rho         /= Ref_Rho;
            Limit_Max_Rho         /= Ref_Rho;
            Limit_Min_Temperature /= Ref_Temperature;
            Limit_Max_Temperature /= Ref_Temperature;
            
            //--END Non-Dimensionalization
            break;
        case MATERIAL_COMP_TYPE_DIM:
            // Set Computations I/O Mode
            MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D;
            break;
        default:
            error("Material_Set_Properties:2: Undefined Material Computation Type - %d", MaterialCompType);
            break;
    }
    
    // Set All Pressure be Perturbations
    Inf_Pressure     -= Gauge_Pressure;
    Outflow_Pressure -= Gauge_Pressure;
    
    PrecondGlobalMach = Inf_Mach_MAX;
}

//------------------------------------------------------------------------------
//! Set Far Field Condition with Mach Ramping
//------------------------------------------------------------------------------
void Material_Set_InfinityCondition(int Iteration) {
    double tmpMach = 0.0;

    // Compute the Far Field Condition Ramping
    if ((Inf_Mach_Ramp > 1) && (Inf_Mach_MAX > Inf_Mach_MIN)) {
        if (Iteration < Inf_Mach_Ramp)
            tmpMach = Inf_Mach_MIN + (Inf_Mach_MAX - Inf_Mach_MIN)*(((double)Iteration)/((double)(Inf_Mach_Ramp-1)));
        else
            tmpMach = Inf_Mach_MAX;
    } else
        tmpMach = Inf_Mach_MAX;

    // Set the Far Field Conditions
    Inf_Mach = tmpMach;
    Inf_Et   = Inf_Et - 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W); 
    Inf_U    = Inf_Mach*Inf_SpeedSound*cos(Inf_Alpha)*cos(Inf_Beta);
    Inf_V    = Inf_Mach*Inf_SpeedSound*sin(Inf_Alpha);
    Inf_W    = Inf_Mach*Inf_SpeedSound*sin(Inf_Beta)*cos(Inf_Alpha);
    Inf_Et   = Inf_Et + 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
}

//------------------------------------------------------------------------------
//! Bound the Solution in Limits
//------------------------------------------------------------------------------
void Material_Limit_Solution(void) {
    double Coeff1, Coeff2;
    int nid;
    double DensityEps = 0.0, PressureEps = 0.0, TemperatureEps = 0.0;
    double SmoothVar1 = 0.0;
    double SmoothVar2 = 0.0;
    
    // Set the Material Eps 0.1 - 1% of the properties
    // Check the Material Computation Type
    switch (MaterialCompType) {
        case MATERIAL_COMP_TYPE_NDIM:
            DensityEps     = 0.0001*Inf_Rho/Ref_Rho;
            PressureEps    = 0.0001*Inf_Pressure/Ref_Pressure;
            TemperatureEps = 0.001*MaterialPropertyLimits[5];
            break;
        case MATERIAL_COMP_TYPE_DIM:
            DensityEps     = 0.0001*Inf_Rho;
            PressureEps    = 0.0001*Inf_Pressure;
            TemperatureEps = 0.001*MaterialPropertyLimits[5];
            break;
        default:
            error("Material_Limit_Solution:1: Undefined Material Computation Type - %d", MaterialCompType);
            break;
    }
    
    // Bound the Solution in Limits for NIST Materials
    if (MaterialType == MATERIAL_TYPE_NIST) {
        switch (VariableType) {
            // Conservative Variable Formulation
            case VARIABLE_CON:
                error("Material_Limit_Solution:2: No Support Yet - %d", VariableType);
                break;
            // Primitive Variable Formulation Density Velocity Pressure
            case VARIABLE_RUP:
                for (int inode = 0; inode < nNode; inode++) {
                    // Maximum Density
                    if (Q1[inode] >= MaterialPropertyLimits[8]) {
                        // Get the Average from the neighbor
                        SmoothVar1 = 0.0; // Density
                        SmoothVar2 = 0.0; // Pressure
                        for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                            nid = crs_JA_Node2Node[i];
                            SmoothVar1 += Q1[nid];
                            SmoothVar2 += Q5[nid];
                        }
                        SmoothVar1 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        SmoothVar2 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        // Check if density is back in limits
                        if (SmoothVar1 >= MaterialPropertyLimits[8])
                            Q1[inode] = 0.5*(SmoothVar1 + MaterialPropertyLimits[8] - DensityEps);
                        else
                            Q1[inode] = SmoothVar1;
                        
                        // Update Pressure with average
                        Q5[inode] = SmoothVar2;
                        
                        // Smooth the Density and Pressure
                        Coeff1 = SolutionSmoothRelaxation;
                        Coeff2 = ((1.0 - Coeff1)/((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode])));
                        for (int iSmooth = 0; iSmooth < 5; iSmooth++) {
                            // First Term
                            SmoothVar1 = Coeff1*Q1[inode];
                            SmoothVar2 = Coeff1*Q5[inode];
                            // Second Term
                            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                                nid = crs_JA_Node2Node[i];
                                SmoothVar1 += Coeff2*Q1[nid];
                                SmoothVar2 += Coeff2*Q5[nid];
                            }
                            Q1[inode] = SmoothVar1;
                            Q5[inode] = SmoothVar2;
                        }
                    }
                    // Maximum Pressure
                    if ((Q5[inode] + Gauge_Pressure) >= MaterialPropertyLimits[7]) {
                        // Get the Average from the neighbor
                        SmoothVar1 = 0.0; // Pressure
                        SmoothVar2 = 0.0; // Density
                        for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                            nid = crs_JA_Node2Node[i];
                            SmoothVar1 += Q5[nid];
                            SmoothVar2 += Q1[nid];
                        }
                        SmoothVar1 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        SmoothVar2 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        if ((SmoothVar1 + Gauge_Pressure) >= MaterialPropertyLimits[7])
                            Q5[inode] = 0.5*(SmoothVar1 + MaterialPropertyLimits[7] - PressureEps - Gauge_Pressure);
                        else
                            Q5[inode] = SmoothVar1;
                        
                        // Update Density with average
                        Q1[inode] = SmoothVar2;
                        
                        // Smooth the Pressure and Density
                        Coeff1 = SolutionSmoothRelaxation;
                        Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        for (int iSmooth = 0; iSmooth < 5; iSmooth++) {
                            // First Term
                            SmoothVar1 = Coeff1*Q5[inode];
                            SmoothVar2 = Coeff1*Q1[inode];
                            // Second Term
                            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                                nid = crs_JA_Node2Node[i];
                                SmoothVar1 += Coeff2*Q5[nid];
                                SmoothVar2 += Coeff2*Q1[nid];
                            }
                            Q5[inode] = SmoothVar1;
                            Q1[inode] = SmoothVar2;
                        }
                    }
                }
                break;
            // Primitive Variable Formulation Pressure Velocity Temperature
            case VARIABLE_PUT:
                error("Material_Limit_Solution:3: No Support Yet - %d", VariableType);
                break;
            // Primitive Variable Formulation Density Velocity Temperature
            case VARIABLE_RUT:
                for (int inode = 0; inode < nNode; inode++) {
                    // Maximum Density
                    if (Q1[inode] >= MaterialPropertyLimits[8]) {
                        // Get the Average from the neighbor
                        SmoothVar1 = 0.0; // Density
                        SmoothVar2 = 0.0; // Temperature
                        for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                            nid = crs_JA_Node2Node[i];
                            SmoothVar1 += Q1[nid];
                            SmoothVar2 += Q5[nid];
                        }
                        SmoothVar1 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        SmoothVar2 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        if (SmoothVar1 >= MaterialPropertyLimits[8])
                            Q1[inode] = 0.5*(SmoothVar1 + MaterialPropertyLimits[8] - DensityEps);
                        else
                            Q1[inode] = SmoothVar1;
                        
                        // Update Temperature with average
                        Q5[inode] = SmoothVar2;
                        
                        // Smooth the Density and Temperature
                        Coeff1 = SolutionSmoothRelaxation;
                        Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        for (int iSmooth = 0; iSmooth < 5; iSmooth++) {
                            // First Term
                            SmoothVar1 = Coeff1*Q1[inode];
                            SmoothVar2 = Coeff1*Q5[inode];
                            // Second Term
                            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                                nid = crs_JA_Node2Node[i];
                                SmoothVar1 += Coeff2*Q1[nid];
                                SmoothVar2 += Coeff2*Q5[nid];
                            }
                            Q1[inode] = SmoothVar1;
                            Q5[inode] = SmoothVar2;
                        }
                    }
                    // Temperature Minimum
                    if ((Q5[inode]) <= MaterialPropertyLimits[5]) {
                        // Get the Average from the neighbor
                        SmoothVar1 = 0.0; // Temperature
                        SmoothVar2 = 0.0; // Density
                        for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                            nid = crs_JA_Node2Node[i];
                            SmoothVar1 += Q5[nid];
                            SmoothVar2 += Q1[nid];
                        }
                        SmoothVar1 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        SmoothVar2 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        if (SmoothVar1 <= MaterialPropertyLimits[5])
                            Q5[inode] = 0.5*(SmoothVar1 + MaterialPropertyLimits[5] + TemperatureEps);
                        else
                            Q5[inode] = SmoothVar1;
                        
                        // Update Density with average
                        Q1[inode] = SmoothVar2;
                        
                        // Smooth the Temperature and Density
                        Coeff1 = SolutionSmoothRelaxation;
                        Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        for (int iSmooth = 0; iSmooth < 5; iSmooth++) {
                            // First Term
                            SmoothVar1 = Coeff1*Q5[inode];
                            SmoothVar2 = Coeff1*Q1[inode];
                            // Second Term
                            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                                nid = crs_JA_Node2Node[i];
                                SmoothVar1 += Coeff2*Q5[nid];
                                SmoothVar2 += Coeff2*Q1[nid];
                            }
                            Q5[inode] = SmoothVar1;
                            Q1[inode] = SmoothVar2;
                        }
                    }
                    // Temperature Maximum
                    if ((Q5[inode]) >= MaterialPropertyLimits[6]) {
                        // Get the Average from the neighbor
                        SmoothVar1 = 0.0; // Temperature
                        SmoothVar2 = 0.0; // Density
                        for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                            nid = crs_JA_Node2Node[i];
                            SmoothVar1 += Q5[nid];
                            SmoothVar2 += Q1[nid];
                        }
                        SmoothVar1 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        SmoothVar2 /= ((double)(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        if (SmoothVar1 >= MaterialPropertyLimits[6])
                            Q5[inode] = 0.5*(SmoothVar1 + MaterialPropertyLimits[6] - TemperatureEps);
                        else
                            Q5[inode] = SmoothVar1;
                        
                        // Update Density with average
                        Q1[inode] = SmoothVar2;
                        
                        // Smooth the Temperature
                        Coeff1 = SolutionSmoothRelaxation;
                        Coeff2 = ((1.0 - Coeff1)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                        for (int iSmooth = 0; iSmooth < 5; iSmooth++) {
                            // First Term
                            SmoothVar1 = Coeff1*Q5[inode];
                            SmoothVar2 = Coeff1*Q1[inode];
                            // Second Term
                            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                                nid = crs_JA_Node2Node[i];
                                SmoothVar1 += Coeff2*Q5[nid];
                                SmoothVar2 += Coeff2*Q1[nid];
                            }
                            Q5[inode] = SmoothVar1;
                            Q1[inode] = SmoothVar2;
                        }
                    }
                }
                break;
            default:
                error("Material_Limit_Solution:4: Undefined Variable Type - %d", VariableType);
                break;
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Equation of State Properties at Control Volume
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_Properties(double *dpVariableIn, double *dpPropertyOut) {
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable, dpPropertyOut);
            // Update the variables
            dpPropertyOut[ 0] = dpVariableIn[0]; // Density
            dpPropertyOut[ 5] = dpVariableIn[1]/dpVariableIn[0];
            dpPropertyOut[ 6] = dpVariableIn[2]/dpVariableIn[0];
            dpPropertyOut[ 7] = dpVariableIn[3]/dpVariableIn[0];
            dpPropertyOut[15] = dpVariableIn[4]/dpVariableIn[0];
            dpPropertyOut[ 3] = dpPropertyOut[3] - Gauge_Pressure;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 0] = dpVariableIn[0];
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 3] = dpVariableIn[4]; // Perturbation
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 3] = dpVariableIn[0]; // Perturbation
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 4] = dpVariableIn[4];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 0] = dpVariableIn[0];
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 4] = dpVariableIn[4];
            dpPropertyOut[ 3] = dpPropertyOut[ 3] - Gauge_Pressure;
            break;
        default:
            error("Material_Get_Properties:1: Undefined Variable Type - %d", VariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Equation of State Extended Properties at Control Volume
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_Extended_Properties(double *dpVariableIn, double *dpPropertyOut) {
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            EOS_Get_Extended_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable, dpPropertyOut);
            // Update the variables
            dpPropertyOut[ 0] = dpVariableIn[0]; // Density
            dpPropertyOut[ 5] = dpVariableIn[1]/dpVariableIn[0];
            dpPropertyOut[ 6] = dpVariableIn[2]/dpVariableIn[0];
            dpPropertyOut[ 7] = dpVariableIn[3]/dpVariableIn[0];
            dpPropertyOut[15] = dpVariableIn[4]/dpVariableIn[0];
            dpPropertyOut[ 3] = dpPropertyOut[3] - Gauge_Pressure;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            EOS_Get_Extended_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 0] = dpVariableIn[0];
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 3] = dpVariableIn[4]; // Perturbation
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            EOS_Get_Extended_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 3] = dpVariableIn[0]; // Perturbation
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 4] = dpVariableIn[4];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            EOS_Get_Extended_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable, dpPropertyOut);
            dpPropertyOut[ 0] = dpVariableIn[0];
            dpPropertyOut[ 5] = dpVariableIn[1];
            dpPropertyOut[ 6] = dpVariableIn[2];
            dpPropertyOut[ 7] = dpVariableIn[3];
            dpPropertyOut[ 4] = dpVariableIn[4];
            dpPropertyOut[ 3] = dpPropertyOut[ 3] - Gauge_Pressure;
            break;
        default:
            error("Material_Get_Extended_Properties:1: Undefined Variable Type - %d", VariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Equation of State Proterties at Control Volume
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_ControlVolume_Properties(double *dpVariableIn,
        double &Rho, double &Pressure, double &Temperature,
        double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy) {
    double Prop[EOS_THERM_DIM];
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable, Prop);
            // Update the variables
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1]/Rho;
            Velocity_V  = dpVariableIn[2]/Rho;
            Velocity_W  = dpVariableIn[3]/Rho;
            TotalEnergy = dpVariableIn[4]/Rho;
            Pressure    = Prop[ 3] - Gauge_Pressure;
            Temperature = Prop[ 4];
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable, Prop);
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Pressure    = dpVariableIn[4]; // Perturbation
            Temperature = Prop[ 4];
            TotalEnergy = Prop[15];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable, Prop);
            Pressure    = dpVariableIn[0]; // Perturbation
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Temperature = dpVariableIn[4];
            Rho         = Prop[ 0];
            TotalEnergy = Prop[15];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable, Prop);
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Temperature = dpVariableIn[4];
            Pressure    = Prop[ 3] - Gauge_Pressure;
            TotalEnergy = Prop[15];
            break;
        default:
            error("Material_Get_ControlVolume_Properties:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    // Update the Variables
    Q2            = Prop[ 8];
    SpeedSound    = Prop[ 9];
    Mach          = Prop[10];
    TotalEnthalpy = Prop[14];
}

//------------------------------------------------------------------------------
//! Compute Equation of State Variables at Face/Edge
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_Face_Properties(double *dpVariableIn, double nx, double ny, double nz,
        double &Rho, double &Pressure, double &Temperature,
        double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &Ubar, double &TotalEnergy, 
        double &TotalEnthalpy) {
    
    Material_Get_ControlVolume_Properties(dpVariableIn, Rho, Pressure, Temperature,
            Velocity_U, Velocity_V, Velocity_W, Q2, SpeedSound, Mach, 
            TotalEnergy, TotalEnthalpy);
    
    // Compute the Ubar
    Ubar = Velocity_U*nx + Velocity_V*ny + Velocity_W*nz;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Get_Density_All(double *dpVariableIn, double *dpDensity) {
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            EOS_Get_Density_All(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable, dpDensity);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            EOS_Get_Density_All(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable, dpDensity);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            EOS_Get_Density_All(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable, dpDensity);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            EOS_Get_Density_All(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable, dpDensity);
            break;
        default:
            error("Material_Get_Density_All:1: Undefined Variable Type - %d", VariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Density(double *dpVariableIn) {
    double Rho = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Density:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Rho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Density_Liquid(double *dpVariableIn) {
    double RhoL = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            RhoL = EOS_Get_Density_Liquid(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            RhoL = EOS_Get_Density_Liquid(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            RhoL = EOS_Get_Density_Liquid(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            RhoL = EOS_Get_Density_Liquid(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Density_Liquid:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return RhoL;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Density_Vapor(double *dpVariableIn) {
    double RhoV = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            RhoV = EOS_Get_Density_Vapor(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            RhoV = EOS_Get_Density_Vapor(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            RhoV = EOS_Get_Density_Vapor(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            RhoV = EOS_Get_Density_Vapor(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Density_Vapor:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return RhoV;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Pressure(double *dpVariableIn) {
    double Pressure     = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Pressure:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return (Pressure - Gauge_Pressure);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Temperature(double *dpVariableIn) {
    double Temperature = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Temperature:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Temperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Mach(double *dpVariableIn) {
    double Mach = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Mach:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Mach;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_TotalEnergy(double *dpVariableIn) {
    double TotalEnergy = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_TotalEnergy:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return TotalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_SpeedSound(double *dpVariableIn) {
    double SpeedSound = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_SpeedSound:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Quality(double *dpVariableIn) {
    double Quality = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Quality = EOS_Get_Quality(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Quality = EOS_Get_Quality(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Quality = EOS_Get_Quality(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Quality = EOS_Get_Quality(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Quality:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Quality;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_SpeedSound(double dvDensity, double dvEnthalpy) {
    double SpeedSound = 0.0;
    
    SpeedSound = EOS_Get_DH_SpeedSound(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_Pressure(double dvDensity, double dvEnthalpy) {
    double Pressure = 0.0;
    
    Pressure = EOS_Get_DH_Pressure(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return Pressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_Temperature(double dvDensity, double dvEnthalpy) {
    double Temperature = 0.0;
    
    Temperature = EOS_Get_DH_Temperature(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return Temperature;
}

//------------------------------------------------------------------------------
//! Compute Q's from EOS Variables
//------------------------------------------------------------------------------
void Material_Get_RUH_To_Q(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double Enthalpy, double *Q) {
    double Pressure, Q2;
    
    // Compute the EOS Based on Non Dimensionalization and Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Q2       = Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W;
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Rho;
            Q[1]     = Rho*Velocity_U;
            Q[2]     = Rho*Velocity_V;
            Q[3]     = Rho*Velocity_W;
            Q[4]     = Rho*(Enthalpy + 0.5*Q2) - Pressure;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Pressure - Gauge_Pressure;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Material_Get_DH_Temperature(Rho, Enthalpy);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Rho;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Pressure - Gauge_Pressure;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Q[0]     = Rho;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Material_Get_DH_Temperature(Rho, Enthalpy);
            break;
        default:
            error("Compute_Q_From_EOS_Variables:1: Undefined Variable Type - %d", VariableType);
            break;
    }
}

