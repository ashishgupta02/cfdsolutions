/*******************************************************************************
 * File:        MatrixBlock.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _MATRIXBLOCK_H
#define _MATRIXBLOCK_H

#ifdef HAVE_MPI
#include <stdlib.h>

#define MatrixBlockIntSize      4

class MatrixBlock {
public:
    int m, n;
    int istart, jstart, iend, jend;
    double **A;

    // Constructor
    MatrixBlock() {
        m      = -1;
        n      = -1;
        istart = -1;
        jstart = -1;
        iend   = -1;
        A      = NULL;
    }

    // Distructor
    ~MatrixBlock() {
        if (A != NULL) {
            for (int i = 0; i < m; i++)
                delete[] A[i];
            delete[] A;
        }
    }

    // Allocate the Memory to store Matrix
    void Allocate(int dim_m, int dim_n) {
        m = dim_m;
        n = dim_n;
        A = new double*[m];
        for (int i = 0; i < m; i++)
            A[i] = new double[n];
    }
};

#endif /* HAVE_MPI */

#endif /* _MATRIXBLOCK_H */
