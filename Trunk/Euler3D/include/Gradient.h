/*******************************************************************************
 * File:        Gradient.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef GRADIENT_H
#define	GRADIENT_H

typedef struct _GradientData {
    double r11;
    double r12;
    double r13;
    double r22;
    double ro22;
    double r23;
    double r33;
    double ro33;
} GradientData;

// Data
extern int NFunction;
// -1: Recompute, 0: Compute, 1: Computed
extern int FlagLSCoeff;
// 0: Unweighted, 1: Inverse of Distance
extern int LSWeight;
extern double **Function;
extern double **GradientX;
extern double **GradientY;
extern double **GradientZ;
extern GradientData *LSCoeff;

// Functions Definations
void Gradient_Init(void);
void Gradient_Finalize(void);
// Add the Function to the Gradient Computation List
int Gradient_Add_Function(double *NewFunction, double *NewGradientX,
        double *NewGradientY, double *NewGradientZ, int Size);
// WeightType: 0 => Unweighted 1 => Inverse Distance
void Compute_Least_Square_Gradient(int WeightType);

#endif	/* GRADIENT_H */

