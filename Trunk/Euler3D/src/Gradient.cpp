/*******************************************************************************
 * File:        Gradient.cpp
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "Gradient.h"

int NFunction;
int FlagLSCoeff;
int LSWeight;
double **Function;
double **GradientX;
double **GradientY;
double **GradientZ;
GradientData *LSCoeff;

//------------------------------------------------------------------------------
//! Initialize the Gradient Compution Resources
//------------------------------------------------------------------------------
void Gradient_Init(void) {
    NFunction   = 0;
    FlagLSCoeff = 0;
    LSWeight    = 0;
    Function  = NULL;
    GradientX = NULL;
    GradientY = NULL;
    GradientZ = NULL;
    LSCoeff   = NULL;
}

//------------------------------------------------------------------------------
//! Free and Reset the Gradient Compution Resources
//------------------------------------------------------------------------------
void Gradient_Finalize(void) {
    // Free Memory
    if (Function != NULL)
        free(Function);
    if (GradientX != NULL)
        free(GradientX);
    if (GradientY != NULL)
        free(GradientY);
    if (GradientZ != NULL)
        free(GradientZ);
    if (LSCoeff != NULL)
        free(LSCoeff);
    
    // Reset the values
    NFunction   = 0;
    FlagLSCoeff = 0;
    LSWeight    = 0;
    Function  = NULL;
    GradientX = NULL;
    GradientY = NULL;
    GradientZ = NULL;
    LSCoeff   = NULL;
}

//------------------------------------------------------------------------------
//! Add the Function to the List for Computation of Gradient
//------------------------------------------------------------------------------
int Gradient_Add_Function(double *NewFunction, double *NewGradientX,
        double *NewGradientY, double *NewGradientZ, int Size) {
    int FunctionID;
    double **tmpFunction  = NULL;
    double **tmpGradientX = NULL;
    double **tmpGradientY = NULL;
    double **tmpGradientZ = NULL;

    // Perform Some Basic Checking
    if ((Size != nNode) || (NewFunction == NULL) || (NewGradientX == NULL)
            || (NewGradientY == NULL) || (NewGradientZ == NULL))
        error("Gradient_Add_Function: Invalid Array Size or Storage Location");

    // Check the size of Function Array
    if (NFunction > 0) {
        tmpFunction  = Function;
        tmpGradientX = GradientX;
        tmpGradientY = GradientY;
        tmpGradientZ = GradientZ;
        
        Function  = NULL;
        GradientX = NULL;
        GradientY = NULL;
        GradientZ = NULL;
        Function  = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientX = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientY = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientZ = (double **) malloc ((NFunction+1)*sizeof(double *));
        for (int i = 0; i < NFunction; i++) {
            Function[i]  = tmpFunction[i];
            GradientX[i] = tmpGradientX[i];
            GradientY[i] = tmpGradientY[i];
            GradientZ[i] = tmpGradientZ[i];
        }
        Function[NFunction]  = NewFunction;
        GradientX[NFunction] = NewGradientX;
        GradientY[NFunction] = NewGradientY;
        GradientZ[NFunction] = NewGradientZ;
        
        // Free the memory
        free(tmpFunction);
        free(tmpGradientX);
        free(tmpGradientY);
        free(tmpGradientZ);
        tmpFunction  = NULL;
        tmpGradientX = NULL;
        tmpGradientY = NULL;
        tmpGradientZ = NULL;
    } else {
        NFunction = 0;
        Function  = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientX = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientY = (double **) malloc ((NFunction+1)*sizeof(double *));
        GradientZ = (double **) malloc ((NFunction+1)*sizeof(double *));
        Function[NFunction]  = NewFunction;
        GradientX[NFunction] = NewGradientX;
        GradientY[NFunction] = NewGradientY;
        GradientZ[NFunction] = NewGradientZ;
    }

    FunctionID = NFunction;
    // Increment the Size
    NFunction++;
    return FunctionID;
}

//------------------------------------------------------------------------------
//! Compute Least Squares Gradient Coefficients
//------------------------------------------------------------------------------
static void Compute_Least_Square_Gradient_Coefficients() {
    int i, iNode, nid;
    double dx, dy, dz, alpha;
    double dtmp;
    
    // Check if Memory is already allocated
    if (FlagLSCoeff == 0)
        LSCoeff = (GradientData *) malloc(nNode*sizeof(GradientData));
    
    // Compute the coefficients
    for (iNode = 0; iNode < nNode; iNode++) {
        LSCoeff[iNode].r11  = 0.0;
        LSCoeff[iNode].r12  = 0.0;
        LSCoeff[iNode].r13  = 0.0;
        LSCoeff[iNode].r22  = 0.0;
        LSCoeff[iNode].r23  = 0.0;
        LSCoeff[iNode].r33  = 0.0;
        LSCoeff[iNode].ro22 = 0.0;
        LSCoeff[iNode].ro33 = 0.0;
        
        // Compute r11, r12, r13
        for (i = crs_IA_Node2Node[iNode]; i < crs_IA_Node2Node[iNode+1]; i++) {
            nid = crs_JA_Node2Node[i];
            dx = coordXYZ[3*nid]     - coordXYZ[3*iNode];
            dy = coordXYZ[3*nid + 1] - coordXYZ[3*iNode + 1];
            dz = coordXYZ[3*nid + 2] - coordXYZ[3*iNode + 2];
            
            // Unweighted
            if (LSWeight == 0)
                alpha = 1.0;
            else // Inverse Distance
                alpha = 1.0/sqrt(dx*dx + dy*dy + dz*dz);

            dx = alpha*dx;
            dy = alpha*dy;
            dz = alpha*dz;

            LSCoeff[iNode].r11 += dx*dx;
            LSCoeff[iNode].r12 += dx*dy;
            LSCoeff[iNode].r13 += dx*dz;
        }

        // Compute r22, r23, ro22
        for (i = crs_IA_Node2Node[iNode]; i < crs_IA_Node2Node[iNode+1]; i++) {
            nid = crs_JA_Node2Node[i];
            dx = coordXYZ[3*nid]     - coordXYZ[3*iNode];
            dy = coordXYZ[3*nid + 1] - coordXYZ[3*iNode + 1];
            dz = coordXYZ[3*nid + 2] - coordXYZ[3*iNode + 2];

            // Unweighted
            if (LSWeight == 0)
                alpha = 1.0;
            else // Inverse Distance
                alpha = 1.0/sqrt(dx*dx + dy*dy + dz*dz);

            dx = alpha*dx;
            dy = alpha*dy;
            dz = alpha*dz;
            
            dtmp  = dy - (LSCoeff[iNode].r12/LSCoeff[iNode].r11)*dx;
            LSCoeff[iNode].r22  += dtmp*dtmp;
            LSCoeff[iNode].r23  += dtmp*dz;
            LSCoeff[iNode].ro22 += dtmp*dy;
        }

        // Compute r33, ro33
        for (i = crs_IA_Node2Node[iNode]; i < crs_IA_Node2Node[iNode+1]; i++) {
            nid = crs_JA_Node2Node[i];
            dx = coordXYZ[3*nid]     - coordXYZ[3*iNode];
            dy = coordXYZ[3*nid + 1] - coordXYZ[3*iNode + 1];
            dz = coordXYZ[3*nid + 2] - coordXYZ[3*iNode + 2];

            // Unweighted
            if (LSWeight == 0)
                alpha = 1.0;
            else // Inverse Distance
                alpha = 1.0/sqrt(dx*dx + dy*dy + dz*dz);

            dx = alpha*dx;
            dy = alpha*dy;
            dz = alpha*dz;
            
            dtmp  = dz - (LSCoeff[iNode].r13/LSCoeff[iNode].r11)*dx
                    - (LSCoeff[iNode].r23/LSCoeff[iNode].r22)*(dy - (LSCoeff[iNode].r12/LSCoeff[iNode].r11)*dx);
            LSCoeff[iNode].r33  += dtmp*dtmp;
            LSCoeff[iNode].ro33 += dtmp*dz;
        }
    }

    // Set the Flag
    FlagLSCoeff = 1;
}

//------------------------------------------------------------------------------
//! Compute Node Based Least Squares Gradient Coefficients
//------------------------------------------------------------------------------
void Compute_Least_Square_Gradient(int WeightType) {
    int i, ifun, iNode, nid;
    double dx, dy, dz, alpha;
    double Wx, Wy, Wz;

    // Check if any function is added
    if (NFunction <= 0) {
        warn("Compute_Least_Square_Gradient: No Function Set for Gradient Computation");
        return;
    }

    // Check if recompuation of coefficient is required
    if ((FlagLSCoeff == 1) && (LSWeight != WeightType)) {
        LSWeight    = WeightType;
        FlagLSCoeff = -1;
    }
    
    // Check if Least Square Coefficients are computed
    if (FlagLSCoeff != 1)
        Compute_Least_Square_Gradient_Coefficients();
    
    // Loop over the Nodes and Compute the Gradients
    for (iNode = 0; iNode < nNode; iNode++) {
        // Initialize the Gradients to Zero
        for (ifun = 0; ifun < NFunction; ifun++) {
            GradientX[ifun][iNode] = 0.0;
            GradientY[ifun][iNode] = 0.0;
            GradientZ[ifun][iNode] = 0.0;
        }
        
        // Loop over the nodes connected to the "iNode"
        for (i = crs_IA_Node2Node[iNode]; i < crs_IA_Node2Node[iNode+1]; i++) {
            nid = crs_JA_Node2Node[i];
            dx = coordXYZ[3*nid]     - coordXYZ[3*iNode];
            dy = coordXYZ[3*nid + 1] - coordXYZ[3*iNode + 1];
            dz = coordXYZ[3*nid + 2] - coordXYZ[3*iNode + 2];

            // Unweighted
            if (LSWeight == 0)
                alpha = 1.0;
            else // Inverse Distance
                alpha = 1.0/sqrt(dx*dx + dy*dy + dz*dz);

            dx = alpha*dx;
            dy = alpha*dy;
            dz = alpha*dz;

            // Compute the edge weigth: Wz
            Wz = (dz - (LSCoeff[iNode].r13/LSCoeff[iNode].r11)*dx
                    - (LSCoeff[iNode].r23/LSCoeff[iNode].r22)*(dy - (LSCoeff[iNode].r12/LSCoeff[iNode].r11)*dx))/LSCoeff[iNode].ro33;
            // Wy
            Wy = (dy - (LSCoeff[iNode].r12/LSCoeff[iNode].r11)*dx - LSCoeff[iNode].r23*Wz)/LSCoeff[iNode].ro22;
            // Wx
            Wx = (dx - LSCoeff[iNode].r12*Wy - LSCoeff[iNode].r13*Wz)/LSCoeff[iNode].r11;

            // Now Compute the Gradient for a Function contributed by each edge
            for (ifun = 0; ifun < NFunction; ifun++) {
                GradientX[ifun][iNode] += Wx*alpha*(Function[ifun][nid] - Function[ifun][iNode]);
                GradientY[ifun][iNode] += Wy*alpha*(Function[ifun][nid] - Function[ifun][iNode]);
                GradientZ[ifun][iNode] += Wz*alpha*(Function[ifun][nid] - Function[ifun][iNode]);
            }
        }
    }
}

