/*******************************************************************************
 * File:        Roe_EntropyFix.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Entropy Fix
//------------------------------------------------------------------------------
void Roe_EntropyFix(double ubar_L, double c_L, double ubar_R, double c_R, double ubar, double c, double **Eigen) {
    double lambda_L, lambda_R, lambda;
    double epsilon;

    switch (EntropyFix) {
        case 1: // M J Kermani AIAA 2001-0083
            // Eigenvalue 1
            Eigen[0][0] = fabs(ubar);
            Eigen[1][1] = Eigen[0][0];
            Eigen[2][2] = Eigen[0][0];
            
            // Eigenvalue 4
            lambda_L = ubar_L + c_L;
            lambda_R = ubar_R + c_R;
            lambda   = ubar + c;
            epsilon = 4.0*MAX(0.0, MAX((lambda - lambda_L), (lambda_R - lambda)));
            if (fabs(lambda) < epsilon)
                Eigen[3][3] = (lambda*lambda + epsilon*epsilon)/(2.0*epsilon);
            else
                Eigen[3][3] = fabs(lambda);

            // Eigenvalue 5
            lambda_L = ubar_L - c_L;
            lambda_R = ubar_R - c_R;
            lambda   = ubar - c;
            epsilon = 4.0*MAX(0.0, MAX((lambda - lambda_L), (lambda_R - lambda)));
            if (fabs(lambda) < epsilon)
                Eigen[4][4] = (lambda*lambda + epsilon*epsilon)/(2.0*epsilon);
            else
                Eigen[4][4] = fabs(lambda);

            break;
        default: // Normal Eigenvalues
            Eigen[0][0] = fabs(ubar);
            Eigen[1][1] = Eigen[0][0];
            Eigen[2][2] = Eigen[0][0];
            Eigen[3][3] = fabs(ubar + c);
            Eigen[4][4] = fabs(ubar - c);
            break;
    }
}

