/*******************************************************************************
 * File:        ParallelLinearAlgebra.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

// TODO: Do not use: This code has to be made compatible with MC CRS format. 

#include "License.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <math.h>
#ifdef DEBUG
#include <sys/time.h>
#endif
#include "MatrixBlock.h"
#include "LinearAlgebra.h"
#include "ParallelLinearAlgebra.h"
#include "Trim_Utils.h"
#include "List.h"
#include "DoubleList.h"

//------------------------------------------------------------------------------
//! Parallel Sparse GMRES for Linear System Equation Ax = b
//------------------------------------------------------------------------------
void ParallelSparseGMRES(double *CRSMatrix, double *X, double *RHS, int *IA, int *JA, int *IAU, int Dim,
        int NRows, int NCols, int KEnd, int NRestarts, int Precondition, int Rank, int ROOT, int NProc) {
    int i, j, k, istart, iend, NRowsProc, iter;
    double norm, delta, gamma, dtmp1, dtmp2;
    bool converged;
    double *x0    = NULL;
    double *r     = NULL;
    double *v     = NULL;
    double *g     = NULL;
    double *u     = NULL;
    double **H    = NULL;
    double **HS   = NULL;
    double **V    = NULL;
    double *c     = NULL;
    double *s     = NULL;
    double *alpha = NULL;
    double *zk    = NULL;
    double *Mdiag = NULL;
    FILE *fpres   = NULL;

    if (Rank == ROOT) {
        // Open File for Residuals
        if ((fpres = fopen("res.txt", "w")) == NULL)
            error("Unable to Open File to write Residual File - %s", "res.txt");
        else
            fprintf(fpres, "\n");
    }

    NRowsProc = NRows/NProc;
    // Get the Row Range for each process
    istart = Rank * NRowsProc;
    iend   = ((Rank + 1) * NRowsProc) - 1;
    
    g = new double[KEnd + 1];
    r = new double[NRowsProc];
    v = new double[NRows];
    u = new double[NRowsProc];
    H = new double*[KEnd];
    V = new double*[KEnd];
    for (i = 0; i < KEnd; i++) {
        H[i] = NULL;
        V[i] = NULL;
    }
    
    c = new double[KEnd];
    s = new double[KEnd];

    // Get Diagonal Terms for Preconditioning
    if (Precondition) {
        Mdiag = new double[NRowsProc];
        for (i = 0; i < NRowsProc; i++) {
            k = IA[i];
            while (k <= IAU[i] && k < Dim) {
                j = JA[k];
                if (i + istart == j) {
                    Mdiag[i] = CRSMatrix[k];
                    break;
                }
                k++;
            }
        }
    }

    // Xo: Set initial normalized vector
    x0 = new double[NRowsProc];
    norm = sqrt(NRows);
    for (i = 0; i < NRowsProc; i++)
        x0[i] = 1.0 / norm;
    
    converged = false;
    // Number or restarts or converged
    for (iter = 0; iter <= NRestarts && converged == false; iter++) {
        for (i = 0; i < KEnd; i++) {
            c[i] = 0.0;
            s[i] = 0.0;
        }
        V[0] = new double[NRowsProc];

        // A.Xo
        ParallelSparseMatrixVectorMultiple(CRSMatrix, x0, V[0], IA, JA, IAU, NRowsProc, NCols, Rank, NProc);

        // ro = b - A.Xo
        for (i = 0; i < NRowsProc; i++) {
            r[i] = RHS[i] - V[0][i];
             // solve M*r_bar = r
            if (Precondition && fabs(Mdiag[i]) > DBL_ZERO)
                r[i] = r[i] / Mdiag[i];
            V[0][i] = r[i];
        }

        // Vo = ro/||ro||
        dtmp1 = dtmp2 = 0.0;
        for (i = 0; i < NRowsProc; i++)
            dtmp1 += V[0][i]*V[0][i];
        MPI_Allreduce(&dtmp1, &dtmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norm = sqrt(dtmp2);
        if (norm < DBL_ZERO)
            warn("Norm: %E is less then Numerical Zero - Code 1!", norm);
        for (i = 0; i < NRowsProc; i++)
            V[0][i] /= norm;

        // g = ||ro||.e1
        for (i = 1; i < KEnd + 1; i++)
            g[i] = 0.0;
        g[0] = norm;
        
        // k loop
        for (k = 0; k < KEnd; k++) {
            // Allocate memory for next column of H
            H[k] = new double[k + 2];

            // Uk = A.Vk
            ParallelSparseMatrixVectorMultiple(CRSMatrix, V[k], u, IA, JA, IAU, NRowsProc, NCols, Rank, NProc);

            // solve M*u = A*v
            if (Precondition) {
                for (i = 0; i < NRowsProc; i++) {
                    if (Mdiag[i] > DBL_ZERO)
                        u[i] = u[i] / Mdiag[i];
                }
            }

            // Start: Arnoldi Process
            // j loop
            for (j = 0; j <= k; j++) {
                // Hjk = (Vj)'.Uk : Parallel Dot Product
                H[k][j] = 0.0;
                dtmp1 = dtmp2 = 0.0;
                for (i = 0; i < NRowsProc; i++)
                    dtmp1 += u[i] * V[j][i];
                MPI_Allreduce(&dtmp1, &dtmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                H[k][j] = dtmp2;
                
                // Uk = Uk - Hjk.Vj
                for (i = 0; i < NRowsProc; i++)
                    u[i] -= H[k][j] * V[j][i];
            }
            // Hk+1,k = ||Uk||
            dtmp1 = dtmp2 = 0.0;
            for (i = 0; i < NRowsProc; i++)
                dtmp1 += u[i]*u[i];
            MPI_Allreduce(&dtmp1, &dtmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            norm = sqrt(dtmp2);
            if (norm < DBL_ZERO)
                warn("Norm: %E is less then Numerical Zero - Code 2!", norm);
            H[k][k + 1] = norm;
            
            // Vk+1,k = Uk/Hk+1.k
            V[k + 1] = new double[NRowsProc];
            for (i = 0; i < NRowsProc; i++)
                V[k + 1][i] = u[i]/H[k][k + 1];
            // End: Arnoldi Process

            // Start: Givens Rotations
            // Apply previous Givens Rotations to the new column H(:,k)
            for (j = 0; j < k; j++) {
                delta = H[k][j];
                H[k][j] = c[j] * delta + s[j] * H[k][j + 1];
                H[k][j + 1] = -s[j] * delta + c[j] * H[k][j + 1];
            }

            // Apply current Givens rotation to the new column H(:,k)
            gamma = sqrt(H[k][k]*H[k][k] + H[k][k + 1]*H[k][k + 1]);
            c[k] = H[k][k] / gamma;
            s[k] = H[k][k + 1] / gamma;
            H[k][k] = gamma;
            H[k][k + 1] = 0.0;

            // Apply current Givens to the RHS vector g
            delta = g[k];
            g[k] = c[k] * delta + s[k] * g[k + 1];
            g[k + 1] = -s[k] * delta + c[k] * g[k + 1];
            norm = fabs(g[k + 1]);
            // End: Givens Rotations

            if (Rank == ROOT) {
                printf("residual Iterations = %d, k = %d : %17.16e \n", iter, k, norm);
                fprintf(fpres, "%d %17.16e\n", (iter * KEnd) + k, norm);
            }

            // Check for convergence and other stoping conditions
            if (norm < DBL_RES_TOLERANCE || k >= KEnd - 1) {
                if (norm < DBL_RES_TOLERANCE)
                    converged = true;
                if (Rank == ROOT && converged)
                    printf("k iteration converged.\n");

                // Form the Solution
                if (Rank == ROOT)
                    printf("Solving for H*alpha=g, z, X...\n");
                
                alpha = new double[k + 1];
                HS = new double*[k + 1];
                for (i = 0; i < k + 1; i++) {
                    HS[i] = new double[k + 1];
                    for (j = 0; j < k + 1; j++) {
                        if (i > j)
                            HS[i][j] = 0.0;
                        else
                            HS[i][j] = H[j][i];
                    }
                }
                
                Solve_BackSubstitution(HS, alpha, g, k + 1);
                for (i = 0; i < k + 1; i++)
                    delete[] HS[i];
                delete[] HS;

                zk = new double[NRowsProc];
                for (i = 0; i < NRowsProc; i++)
                    zk[i] = 0.0;
                for (j = 0; j < k + 1; j++) {
                    for (i = 0; i < NRowsProc; i++)
                        zk[i] += V[j][i] * alpha[j];
                }
                delete[] alpha;
                
                for (i = 0; i < NRowsProc; i++)
                    x0[i] = x0[i] + zk[i];
                delete[] zk;
                
                // Break the K loop
                break;
            }
        }

        for (i = 0; i < KEnd; i++) {
            if (H[i] != NULL)
                delete[] H[i];
            if (V[i] != NULL)
                delete[] V[i];
            H[i] = NULL;
            V[i] = NULL;
        }
    }
    
    MPI_Gather(x0, NRowsProc, MPI_DOUBLE, X, NRowsProc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (Rank == ROOT)
        fclose(fpres);
    
    delete[] x0;
    delete[] g;
    delete[] r;
    delete[] v;
    delete[] u;
    delete[] c;
    delete[] s;
    delete[] H;
    delete[] V;
    delete[] Mdiag;
}

//------------------------------------------------------------------------------
//! Parallel Sparse Arnodi process for Eigenvalue Computation
//------------------------------------------------------------------------------
void ParallelSparseArnoldi(double *CRSMatrix, int *IA, int *JA, int *IAU, int NRows, int NCols, int KEnd, int Iterations, int Rank, int ROOT, int NProc) {
    double startTime  = 0.0;
    double totalTime0 = 0.0;
    int i, j, k, iter, lastk, NRowsProc;
    double dtmp, dtmp1;
    double *q_k   = NULL;
    double **Q    = NULL;
    double **H    = NULL;
    double **HQR  = NULL;
    FILE *fpeigen = NULL;
    
    // Initialize data
    H = new double*[KEnd];
    Q = new double*[KEnd];
    NRowsProc = NRows/NProc;

    if (Rank == ROOT) {
        // Open File for EigenValue
        if ((fpeigen = fopen("eigenfile.txt", "w")) == NULL)
            error("Unable to Open File to write EigenValue - %s", "eigenfile.txt");
        else
            fprintf(fpeigen, "\n");
    }
    
    // Set the initial normalized search vector
    Q[0] = new double[NRowsProc];
    for (i = 0; i < NRowsProc; i++)
        Q[0][i] = 1.0;
    // Normalize initial search vector
    dtmp = dtmp1 = 0.0;
    for (i = 0; i < NRowsProc; i++)
        dtmp += Q[0][i]*Q[0][i];
    MPI_Allreduce(&dtmp, &dtmp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dtmp = sqrt(dtmp1);
    for (i = 0; i < NRowsProc; i++)
        Q[0][i] /= dtmp;

    // Start the parallel Arnoldi
    lastk = 0;
    q_k = new double[NRows];
    for (iter = 0; iter < Iterations; iter++) {
        for (k = 0; k < KEnd; k++) {
            lastk = k;
            // Allocate memory for next search vector q_k+1/u_k
            Q[k + 1] = new double[NRowsProc];
            // Allocate memory for next column of H
            H[k] = new double[k + 2];

            // Start Parallel Matrix Vector Multiply
            startTime = MPI_Wtime();
            ParallelSparseMatrixVectorMultiple(CRSMatrix, Q[k], Q[k+1], IA, JA, IAU, NRowsProc, NCols, Rank, NProc);
            totalTime0 += MPI_Wtime() - startTime;
            
            for (j = 0; j <= k; j++) {
                //  h_jk = u_k * q_j
                H[k][j] = 0.0;
                dtmp = 0.0;
                for (i = 0; i < NRowsProc; i++)
                    dtmp += Q[k + 1][i] * Q[j][i];
                MPI_Allreduce(&dtmp, &H[k][j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                //  u_k = u_k - h_jk * q_j
                dtmp = H[k][j];
                for (i = 0; i < NRowsProc; i++)
                    Q[k + 1][i] -= dtmp * Q[j][i];
            }
            // h_k+1,k = ||u_k||
            dtmp = 0.0;
            for (i = 0; i < NRowsProc; i++)
                dtmp += Q[k + 1][i] * Q[k + 1][i];
            MPI_Allreduce(&dtmp, &H[k][k + 1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            H[k][k + 1] = sqrt(H[k][k + 1]);

            // Serial Mode - QR Iteration
            if (Rank == ROOT) {
                HQR = new double*[k + 1];
                for (i = 0; i < k + 1; i++) {
                    HQR[i] = new double[k + 1];
                    for (j = 0; j < k + 1; j++) {
                        if (i >= j + 2)
                            HQR[i][j] = 0.0;
                        else
                            HQR[i][j] = H[j][i];
                    }
                }
                
                // QR Iteration
                QR_Iteration(HQR, k + 1, 100);
                // Report Eigenvalues
                double a, b, c, d, discriminant;
                double ew1r, ew1i, ew2r, ew2i, mag;
#ifdef DEBUG
                printf("Schur form: \n");
                for (i = 0; i < k+1; i++) {
                    for (j = 0; j < k+1; j++)
                        printf("%10.5lf ", HQR[i][j]);
                    printf("\n");
                }
#endif
                // Get the Eigenvalues
                if (k == 0)
                    fprintf(fpeigen, "%15.8lf \t %d \n", HQR[0][0], k);
                for (i = 0; i < k; i++) {
                    if (HQR[i + 1][i] < DBL_ZERO) {
                        fprintf(fpeigen, "%15.8lf \t %d \n", HQR[i][i], k);
                        if (i == k-1)
                            fprintf(fpeigen, "%15.8lf \t %d \n", HQR[i + 1][i + 1], k);
                    } else {
                        // Complex Eigenvalues
                        // Solve characteristic polynomial is (a-λ)(d-λ)-bc
                        a = HQR[i][i];
                        b = HQR[i][i + 1];
                        c = HQR[i + 1][i];
                        d = HQR[i + 1][i + 1];
                        ew1r = ew2r = (a + d) / 2.0;
                        discriminant = 4.0*b*c + (a - d)*(a - d);
                        if (discriminant < 0.0) {
                            ew1i = sqrt(fabs(discriminant)) / 2.0;
                            ew2i = -ew1i;
                        } else {
                            ew1i = ew2i = 0.0;
                            ew1r += sqrt(discriminant) / 2.0;
                            ew2r -= sqrt(discriminant) / 2.0;
                        }
                        // Get the magnitude of complex eiqenvalue
                        mag = sqrt(ew1r*ew1r + ew1i*ew1i);
                        fprintf(fpeigen, "%15.8lf \t %d \n", mag, k);
                        mag = sqrt(ew2r*ew2r + ew2i*ew2i);
                        fprintf(fpeigen, "%15.8lf \t %d \n", mag, k);
                        i++;
                    }
                }
                // Free HQR
                if (HQR != NULL) {
                    for (i = 0; i < k + 1; i++)
                        delete[] HQR[i];
                    delete[] HQR;
                    HQR = NULL;
                }
            }
            // if (h_k+1,k == 0.0) then stop
            if (H[k][k + 1] < DBL_ZERO) {
                warn("H_k+1,k = 0.0 on proc %d \n", Rank);
                break;
            }
            // q_k+1 = u_k/h_k+1,k;
            for (i = 0; i < NRowsProc; i++)
                Q[k + 1][i] /= H[k][k + 1];
        }
        for (k = 0; k < lastk; k++) {
            delete[] Q[k + 1];
            delete[] H[k];
        }
    }

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (Rank == ROOT) {
        printf("Elapsed MPI time, matrix-vector multiply = %17.16lf \n", totalTime0);
        fclose(fpeigen);
    }
#endif
    
    delete[] q_k;
    delete[] Q[0];
    delete[] Q;
    delete[] H;
}


//------------------------------------------------------------------------------
//! Parallel Sparse Matrix Vector Multiple
//! y = Ax
// NRows: No of Processor Own Rows
// NCols: Total Number of Columns in Matrix
//------------------------------------------------------------------------------
void ParallelSparseMatrixVectorMultiple(double *A, double *x, double *y, int *IA, int *JA, int *IAU, int NRows, int NCols, int Rank, int NProc) {
    int i, idx, j;

    if (NProc == 1 && NRows == NCols) {
        SparseMatrixVectorMultiply(A, x, y, IA, JA, IAU, NRows, NCols);
        return;
    }

    int istart, iend, owner, nrequests;
    int *nreq_mine   = NULL;
    int *nreq_theirs = NULL;
    struct double_int **req_mine   = NULL;
    struct double_int **req_theirs = NULL;

    // Allocate Memory and Initialize Variables
    istart    = Rank * NRows;
    iend      = ((Rank + 1) * NRows) - 1;
    nrequests = 2*(NProc - 1);
    nreq_mine   = new int[NProc];
    nreq_theirs = new int[NProc];
    req_mine    = new double_int*[NProc];
    req_theirs  = new double_int*[NProc];
    for (i = 0; i < NProc; i++) {
        nreq_mine[i]   = 0;
        nreq_theirs[i] = 0;
        req_mine[i]    = NULL;
        req_theirs[i]  = NULL;
    }

    bool tag[NCols];
    for (i = 0; i < NCols; i++)
        tag[i] = false;
    // Loop over rows
    for (i = 0; i < NRows; i++) {
        // Loop over cols
        for (idx = IA[i]; idx <= IAU[i]; idx++) {
            j = JA[idx];
            owner = j/NRows;
            if (owner != Rank) {
                // tag to make sure we don't 'collect' this column more than once
                if (!tag[j]) {
                    nreq_mine[owner]++;
                    tag[j] = true;
                }
            }
        }
    }

    // allocate each req[p] as new double_int[nreq[p]]
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank && nreq_mine[owner] > 0)
            req_mine[owner] = new double_int[nreq_mine[owner]];
    }
    for (i = 0; i < NCols; i++)
        tag[i] = false;
    // repurpose nreq[p] as index into req[p][nreq[p]]
    for (owner = 0; owner < NProc; owner++)
        nreq_mine[owner] = 0;
    // Loop over the rows
    for (i = 0; i < NRows; i++) {
        // Loop over the cols
        for (idx = IA[i]; idx <= IAU[i]; idx++) {
            j = JA[idx];
            owner = j/NRows;
            if (owner != Rank) {
                // tag to make sure we don't 'collect' this column more than once
                if (!tag[j]) {
                    req_mine[owner][nreq_mine[owner]].i = j;
                    req_mine[owner][nreq_mine[owner]].d = 0.0;
                    nreq_mine[owner]++;
                    tag[j] = true;
                }
            }
        }
    }

    // now nreq[owner] is again the number of values required from owner
    MPI_Request* request = new MPI_Request[nrequests];
    MPI_Status aStatus;
    int aCount;
    int r = 0;

    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank)
            MPI_Isend(req_mine[owner], nreq_mine[owner], MPI_DOUBLE_INT, owner, 0, MPI_COMM_WORLD, &request[r++]);
    }
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank) {
            MPI_Probe(owner, 0, MPI_COMM_WORLD, &aStatus);
            MPI_Get_count(&aStatus, MPI_DOUBLE_INT, &aCount);
            req_theirs[owner] = new double_int[aCount];
            MPI_Irecv(req_theirs[owner], aCount, MPI_DOUBLE_INT, owner, 0, MPI_COMM_WORLD, &request[r++]);
            nreq_theirs[owner] = aCount;
        }
    }
    MPI_Waitall(nrequests, request, MPI_STATUSES_IGNORE);

    // now fill in all req_theirs with our values of x and send it back
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank) {
            for (i = 0; i < nreq_theirs[owner]; i++) {
                j = req_theirs[owner][i].i;
                req_theirs[owner][i].d = x[j - istart];
            }
        }
    }

    r = 0;
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank)
            MPI_Isend(req_theirs[owner], nreq_theirs[owner], MPI_DOUBLE_INT, owner, 0, MPI_COMM_WORLD, &request[r++]);
    }
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank)
            MPI_Irecv(req_mine[owner], nreq_mine[owner], MPI_DOUBLE_INT, owner, 0, MPI_COMM_WORLD, &request[r++]);
    }
    MPI_Waitall(nrequests, request, MPI_STATUSES_IGNORE);

    for (i = 0; i < NProc; i++) {
        if (req_theirs[i] != NULL)
            delete[] req_theirs[i];
    }
    delete[] req_theirs;
    delete[] nreq_theirs;
    delete[] request;

    // Create the tuple (int, double)
    DoubleList Xn;
    for (owner = 0; owner < NProc; owner++) {
        if (owner != Rank) {
            for (i = 0; i < nreq_mine[owner]; i++)
                Xn.Check_List(req_mine[owner][i].i, req_mine[owner][i].d);
        }
    }

    // Use Tuple List to multiple to the matrix-vector multiplication
    double temp;
    for (i = 0; i < NRows; i++) {
        y[i] = 0;
        for (idx = IA[i]; idx <= IAU[i]; idx++) {
            j = JA[idx];
            owner = j/NRows;
            if (owner == Rank)
                temp = x[j - istart];
            else
                temp = Xn.data[Xn.Index(j)];
            y[i] += A[idx] * temp;
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! Parallel Gaussain PA = LU with Block-Cyclic Partition
//------------------------------------------------------------------------------
void ParallelGaussainLUDecompose(double **A, int *P, int *Pinv, int Dim, int ROOT, int Rank, int NProc, int RecursionLevel) {
#ifdef DEBUG
    // For Timing of Algorithm
    double delta_time;
    struct timeval start_time;
    struct timeval end_time;
#endif

    // Temporary Variables
    int i, j;
    int itmp;
    double dtmp;

    // Block Matrix Variables
    int iblock, blocki, blockj, nblocks, blocksize, proc_num, sqrt_np;
    int ranksize, ranksize_n_minus_1, blocks_per_proc, n_per_proc;
    int **RankMatrix    = NULL;
    MatrixBlock *blocklist = NULL;

    // Variable for Distributing and Assembling the Matrix
    double *databuffer = NULL;
    int packbuffer_size = 0;
    char **distbuffer = NULL;
    char **procbuffer = NULL;
    int istart, jstart, iend, jend, position;
    MPI_Request *requests = NULL;
    int nrequests = 0;
    int request = 0;
    struct double_int *dbl_ints = NULL;

    // Variables for LU Decomposition Algorithm
    int col, local_col, local_row;
    List block_vec;
    List receivers;
    int sender;
    double pivot_max;
    int pivot_max_index;
    int pivot_max_proc;
    struct double_int pivot_max_in, pivot_max_out;
    int g_i, g_j, d_i;
    double *m = NULL;
    double *r = NULL;
    double *pivot_gather_dbl = NULL;
    MPI_Status aStatus;
    int columnblock, rowblock;
    int is_receiver;

    // Initialize
    sqrt_np = 0;
    ranksize_n_minus_1 = 0;
    ranksize = 0;
    nblocks = 0;
    blocksize = 0;

    // START: ================== Rank Matrix Computation =======================
    // Check that number of processors has integer square root
    itmp = (int) sqrt((double)NProc);
    if (itmp*itmp == NProc)
        sqrt_np = itmp;
    else
        error("Number of MPI Processors is not a Perfect Square!");

    // Get the Number of Block Matrix
    if (NProc == 1) {
        // Support for Dim%2 == 0 and Recursion Level = 0
        if (Dim%2 == 0 && RecursionLevel == 0) {
            nblocks = (int) pow(4.0, RecursionLevel);
            ranksize = (int) pow(2.0, RecursionLevel);
            ranksize_n_minus_1 = -1;
        } else
            error("Recursion Level Not supported !");
    } else {
        nblocks = (int) pow((double)NProc, RecursionLevel + 1);
        ranksize = (int) pow((double)sqrt_np, RecursionLevel + 1);
        ranksize_n_minus_1 = ranksize / sqrt_np;
    }
    blocks_per_proc = nblocks / NProc;

    // Check if Block Cyclic Partition is Possible and Get Block Size
    dtmp = (double)Dim / sqrt((double)nblocks);
    itmp = (int)dtmp;
    if (fabs(dtmp - itmp) > 0)
        error("Cannot Perform Block Cyclic Partition With Matrix Size %d NProc %d RLevel %d", Dim, NProc, RecursionLevel);
    else
        blocksize = itmp;
    
    n_per_proc = blocks_per_proc * blocksize;

    // Create a Rank Matrix
    RankMatrix = new int*[ranksize];
    for (i = 0; i < ranksize; i++)
        RankMatrix[i] = new int[ranksize];
    
    if (NProc == 1) {
        for (i = 0; i < ranksize; i++)
            for (j = 0; j < ranksize; j++)
                RankMatrix[i][j] = 0;
    } else {
        for (blocki = 0; blocki < ranksize_n_minus_1; blocki++) {
            for (blockj = 0; blockj < ranksize_n_minus_1; blockj++) {
                proc_num = 0;
                for (i = 0; i < sqrt_np; i++) {
                    for (j = 0; j < sqrt_np; j++) {
                        RankMatrix[i + blocki * sqrt_np][j + blockj * sqrt_np] = proc_num;
                        proc_num++;
                    }
                }
            }
        }
    }

#ifdef DEBUG
    if (Rank == 0) {
        printf("RankMatrix:\n");
        for (i = 0; i < ranksize; i++) {
            for (j = 0; j < ranksize; j++) {
                printf("%2d ", RankMatrix[i][j]);
            }
            printf("\n");
        }
        info("No of Blocks: %d", nblocks);
        info("Block Dimensions: %d x %d", blocksize, blocksize);
        info("Rank matrix size: %d x %d", ranksize, ranksize);
    }
#endif
    // END: =================== Rank Matrix Computation ========================

    // Memory set up for Distribution and Assemble of Matrix
    blocklist = new MatrixBlock[blocks_per_proc];
    for (iblock = 0; iblock < blocks_per_proc; iblock++)
        blocklist[iblock].Allocate(blocksize, blocksize);
    
    // Estimate the Distribution Pack Size for Send and Receive
    databuffer = new double[blocksize * blocksize];
    packbuffer_size = 0;
    itmp = 0;
    MPI_Pack_size(MatrixBlockIntSize, MPI_INT, MPI_COMM_WORLD, &itmp);
    packbuffer_size += itmp;
    MPI_Pack_size(blocksize*blocksize, MPI_DOUBLE, MPI_COMM_WORLD, &itmp);
    packbuffer_size += itmp;
    // Memory for Distribution and Assemble of Matrix
    if (Rank == ROOT) {
        distbuffer = new char*[nblocks];
        for (iblock = 0; iblock < nblocks; iblock++)
            distbuffer[iblock] = new char[packbuffer_size];
    }
    // Memory for Each Process to receive Matrix
    procbuffer = new char*[blocks_per_proc];
    for (iblock = 0; iblock < blocks_per_proc; iblock++)
        procbuffer[iblock] = new char[packbuffer_size];
    
    // START: ============== Distribute Matrix From Root Process ===============
    nrequests = blocks_per_proc;
    if (Rank == ROOT)
        nrequests += (NProc * blocks_per_proc);
    requests = new MPI_Request[nrequests];

    // Root Process Sends the Block Data to all Process including itself
    request = 0;
    if (Rank == ROOT) {
        iblock = 0;
        for (blocki = 0; blocki < ranksize; blocki++) {
            for (blockj = 0; blockj < ranksize; blockj++) {
                position = 0;
                // Get the Block Data
                istart = blocki*blocksize;
                jstart = blockj*blocksize;
                iend = ((blocki + 1) * blocksize) - 1;
                jend = ((blockj + 1) * blocksize) - 1;
                for (i = 0; i < blocksize; i++) {
                    for (j = 0; j < blocksize; j++) {
                        databuffer[position] = A[i + istart][j + jstart];
                        position++;
                    }
                }

                // Pack the Block Matrix for Transfer
                position = 0;
                MPI_Pack(&istart,    1, MPI_INT, distbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&jstart,    1, MPI_INT, distbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&iend,      1, MPI_INT, distbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&jend,      1, MPI_INT, distbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
                MPI_Pack(databuffer, blocksize*blocksize, MPI_DOUBLE, distbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
                // Non Blocking Send
                MPI_Isend(distbuffer[iblock], packbuffer_size, MPI_PACKED, RankMatrix[blocki][blockj], 0, MPI_COMM_WORLD, &requests[request]);
                request++;
                iblock++;
            }
        }
    }
    
    // Receive the Matrix Block Data From Root Process -- Root Process also receives the Matrix Block
    for (iblock = 0; iblock < blocks_per_proc; iblock++) {
        MPI_Irecv(procbuffer[iblock], packbuffer_size, MPI_PACKED, ROOT, 0, MPI_COMM_WORLD, &requests[request]);
        request++;
    }
    // Wait Untill all Data is received and will Deallocate "requests" for Non-Blocking Communication
    MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
    // Unpack the Distributed Matrix Data from Buffer
    for (iblock = 0; iblock < blocks_per_proc; iblock++) {
        position = 0;
        MPI_Unpack(procbuffer[iblock], packbuffer_size, &position, &blocklist[iblock].istart, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(procbuffer[iblock], packbuffer_size, &position, &blocklist[iblock].jstart, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(procbuffer[iblock], packbuffer_size, &position, &blocklist[iblock].iend,   1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(procbuffer[iblock], packbuffer_size, &position, &blocklist[iblock].jend,   1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(procbuffer[iblock], packbuffer_size, &position, databuffer, blocksize*blocksize, MPI_DOUBLE, MPI_COMM_WORLD);

        position = 0;
        for (i = 0; i < blocksize; i++) {
            for (j = 0; j < blocksize; j++) {
                blocklist[iblock].A[i][j] = databuffer[position];
                position++;
            }
        }
    }
    // END: =============== Distribute Matrix From Root Process ================

#ifdef DEBUG
    // Synchronize All Process
    MPI_Barrier(MPI_COMM_WORLD);
    // Get the Start Time
    if (Rank == ROOT)
        gettimeofday(&start_time, NULL);
#endif
    
    // =========================================================================
    // START: Gaussian LU Decomposition

    // Initialize P and P-inverse Permutation Vectors
    for (i = 0; i < Dim; i++)
        P[i] = Pinv[i] = i;
    
    // Allocate for multiplier and pivot row lists
    dbl_ints = new double_int[n_per_proc];

    pivot_gather_dbl = new double[NProc];
    m = new double[Dim];
    r = new double[Dim];
    
    // Loop over Columns of Matrix
    for (col = 0; col < Dim; col++) {
        // Step 1: Find the Pivot for the current Column
        for (i = 0; i < NProc; i++)
            pivot_gather_dbl[i] = 0.0;

        // Check the Block List for this column on current process
        block_vec.Reset(0);
        for (iblock = 0; iblock < blocks_per_proc; iblock++)
            if (col >= blocklist[iblock].jstart && col <= blocklist[iblock].jend)
                block_vec.Check_List(iblock);
        
        // Search for Pivot in the selected Blocks having this column
        pivot_max = 0.0;
        pivot_max_index = col;
        pivot_max_proc = -1;
        for (int it = 0; it < block_vec.max ; it++) {
            iblock = block_vec.list[it];
            local_col = col - blocklist[iblock].jstart;
            // Loop over rows in the block
            for (i = 0; i < blocksize; i++) {
                // Get the global row id
                g_i = i + blocklist[iblock].istart;
                // Check if pivot row is greater or equal to the current column
                if (Pinv[g_i] >= col) {
                    if (fabs(blocklist[iblock].A[i][local_col]) > fabs(pivot_max)) {
                        pivot_max = blocklist[iblock].A[i][local_col];
                        pivot_max_index = Pinv[g_i];
                    }
                }
            }
        }
        pivot_max_in.d = fabs(pivot_max);
        pivot_max_in.i = Rank;

        // Reduce: the pivot and the rank
        MPI_Allreduce(&pivot_max_in, &pivot_max_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        pivot_max = pivot_max_out.d;
        pivot_max_proc = pivot_max_out.i;
        
        // Broadcast: form process that owns max pivot its row index
        if (pivot_max_proc >= 0)
            MPI_Bcast(&pivot_max_index, 1, MPI_INT, pivot_max_proc, MPI_COMM_WORLD);

        // Step 2: Update P and Pinv
        if ((pivot_max_index >= 0) && (pivot_max_index < Dim)) {
            itmp = P[col];
            P[col] = P[pivot_max_index];
            P[pivot_max_index] = itmp;

            // Now swap Pinv[P[pivot_max_index]] and Pinv[P[col]]:
            itmp = Pinv[P[col]];
            Pinv[P[col]] = Pinv[P[pivot_max_index]];
            Pinv[P[pivot_max_index]] = itmp;
        } else
            error("Pivot Max Index Out of Bounds: %d", pivot_max_index);

        // Step 3: Compute the mulitpiler for all rows below pivot column in the Block list
        for (i = 0; i < n_per_proc; i++) {
            dbl_ints[i].d = 0.0;
            dbl_ints[i].i = -1;
        }

        // Loop over the Blocks containing the column
        d_i = 0;
        for (int it = 0; it < block_vec.max; it++) {
            iblock = block_vec.list[it];
            local_col = col - blocklist[iblock].jstart;
            // Loop over rows in the block
            for (i = 0; i < blocksize; i++) {
                // Get the global row id
                g_i = i + blocklist[iblock].istart;
                // Check if pivot row is greater to the current column
                if (Pinv[g_i] > col) {
                    blocklist[iblock].A[i][local_col] /= pivot_max;
                    dbl_ints[d_i].d = blocklist[iblock].A[i][local_col];
                    dbl_ints[d_i].i = g_i;
                    d_i++;
                }
            }
        }

        // Step 4: Broadcast multipliers to processes right of pivot column
        receivers.Reset(0);
        blockj = Rank % sqrt_np;
        blocki = (Rank - blockj) / sqrt_np;
        columnblock = (int) floor((float)col / (float)blocksize);
        rowblock    = (int) floor((float)P[col] / (float)blocksize);

        // Step 4a: Identify Current Block is Senders and Receivers of the multipiler
        is_receiver = 0;
        // Sender if RankMatrix[blocki][blockj] is equal to RankMatrix[blocki][columnblock]
        sender = RankMatrix[blocki][columnblock];
        // Check the process is receiver
        for (j = 0; j < sqrt_np; j++) {
            if (RankMatrix[blocki][j] != sender) {
                receivers.Check_List(RankMatrix[blocki][j]);
                if (RankMatrix[blocki][j] == Rank)
                    is_receiver = 1;
            }
        }

        // Step 4b: If Rank == sender, send to all recievers, otherwise receive
        if (Rank == sender) {
            for (int it = 0; it < receivers.max; it++)
                MPI_Send(dbl_ints, d_i, MPI_DOUBLE_INT, receivers.list[it], 0, MPI_COMM_WORLD);
        } else if (is_receiver) {
            MPI_Probe(sender, 0, MPI_COMM_WORLD, &aStatus);
            MPI_Get_count(&aStatus, MPI_DOUBLE_INT, &d_i);
            MPI_Recv(dbl_ints, d_i, MPI_DOUBLE_INT, sender, 0, MPI_COMM_WORLD, &aStatus);
        }

        for (i = 0; i < Dim; i++)
            m[i] = 0.0;

        // Step 4c: Copy multiplier in its "row" in the m vector
        for (i = 0; i < d_i; i++)
            m[dbl_ints[i].i] = dbl_ints[i].d;

        
        // Step 5: Broadcast pivot row to processes below pivot row
        for (i = 0; i < n_per_proc; i++) {
            dbl_ints[i].d = 0.0;
            dbl_ints[i].i = -1;
        }

        // Loop over all blocks on this proc, if block contains pivot row, copy values and indices into dbl_ints
        d_i = 0;
        for (iblock = 0; iblock < blocks_per_proc; iblock++) {
            if (P[col] >= blocklist[iblock].istart && P[col] <= blocklist[iblock].iend) { // block contains pivot row
                local_row = P[col] - blocklist[iblock].istart;
                // Loop over the columns of the block
                for (j = 0; j < blocksize; j++) {
                    g_j = j + blocklist[iblock].jstart;
                    dbl_ints[d_i].d = blocklist[iblock].A[local_row][j];
                    dbl_ints[d_i].i = g_j;
                    d_i++;
                }
            }
        }

        // Step 6: Find the receivers for the pivot row
        receivers.Reset(0);
        is_receiver = 0;
        // which proc owns the pivot row block in this proc's columns
        sender = RankMatrix[rowblock][blockj]; 
        for (i = 0; i < sqrt_np; i++) {
            if (RankMatrix[i][blockj] != sender) {
                receivers.Check_List(RankMatrix[i][blockj]);
                if (RankMatrix[i][blockj] == Rank)
                    is_receiver = 1;
            }
        }

        // Step 6a: if Rank == sender, send to all recievers, otherwise receive
        if (Rank == sender) {
            for (int it = 0; it < receivers.max; it++)
                MPI_Send(dbl_ints, d_i, MPI_DOUBLE_INT, receivers.list[it], 0, MPI_COMM_WORLD);
        } else if (is_receiver) {
            MPI_Probe(sender, 0, MPI_COMM_WORLD, &aStatus);
            MPI_Get_count(&aStatus, MPI_DOUBLE_INT, &d_i);
            MPI_Recv(dbl_ints, d_i, MPI_DOUBLE_INT, sender, 0, MPI_COMM_WORLD, &aStatus);
        }
        for (i = 0; i < Dim; i++)
            r[i] = 0.0;
        for (i = 0; i < d_i; i++)
            r[dbl_ints[i].i] = dbl_ints[i].d;

        // Step 7: A[i][j] = A[i][j] - m[i]*r[j]
        for (iblock = 0; iblock < blocks_per_proc; iblock++) {
            // For each row in the current block
            for (i = blocklist[iblock].istart; i <= blocklist[iblock].iend; i++) {
                if (Pinv[i] > col) {
                    local_row = i - blocklist[iblock].istart;
                    // For each columns of the curren block
                    for (j = 0; j < blocksize; j++) {
                        g_j = j + blocklist[iblock].jstart;
                        if (g_j > col)
                            blocklist[iblock].A[local_row][j] -= m[i] * r[g_j];
                    }
                }
            }
        }

    }
    delete[] dbl_ints;
    delete[] m;
    delete[] r;

    // END: Gaussian LU Decomposition
    // =========================================================================
    
#ifdef DEBUG
    // Synchronize All Process
    MPI_Barrier(MPI_COMM_WORLD);
    // Compute the Total Time
    if (Rank == ROOT) {
        gettimeofday(&end_time, NULL);
        delta_time = (double) (end_time.tv_sec - start_time.tv_sec);
        delta_time += (double) (end_time.tv_usec - start_time.tv_usec) / 1.0e+6;
        info("Elapsed time = %10.5lf", delta_time);
    }
#endif
    
    // START: =============== Asssemble Matrix to Root Process =================
    nrequests = blocks_per_proc;
    if (Rank == ROOT)
        nrequests += (NProc * blocks_per_proc);
    requests = new MPI_Request[nrequests];
    
    // All Process including Root Send The Block Matrix Data to Root Process
    request = 0;
    for (iblock = 0; iblock < blocks_per_proc; iblock++) {
        position = 0;
        // Get the Block Data
        for (i = 0; i < blocksize; i++) {
            for (j = 0; j < blocksize; j++) {
                databuffer[position] = blocklist[iblock].A[i][j];
                position++;
            }
        }

        // Pack the Block Matrix for Transfer
        position = 0;
        MPI_Pack(&blocklist[iblock].istart, 1, MPI_INT, procbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&blocklist[iblock].jstart, 1, MPI_INT, procbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&blocklist[iblock].iend,   1, MPI_INT, procbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&blocklist[iblock].jend,   1, MPI_INT, procbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
        MPI_Pack(databuffer, blocksize*blocksize, MPI_DOUBLE, procbuffer[iblock], packbuffer_size, &position, MPI_COMM_WORLD);
        // Non Blocking Send
        MPI_Isend(procbuffer[iblock], packbuffer_size, MPI_PACKED, ROOT, 0, MPI_COMM_WORLD, &requests[request]);
        request++;
    }

    // Receive the Matrix Block Data From All Process including Root
    if (Rank == ROOT) {
        iblock = 0;
        for (blocki = 0; blocki < ranksize; blocki++) {
            for (blockj = 0; blockj < ranksize; blockj++) {
                MPI_Irecv(distbuffer[iblock], packbuffer_size, MPI_PACKED, RankMatrix[blocki][blockj], 0, MPI_COMM_WORLD, &requests[request]);
                request++;
                iblock++;
            }
        }
    }
    // Wait Untill all Data is received and will Deallocate "requests" for Non-Blocking Communication
    MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
    // Unpack the Distributed Matrix Data from Buffer to Assemble Matrix
    if (Rank == ROOT) {
        for (iblock = 0; iblock < nblocks; iblock++) {
            position = 0;
            MPI_Unpack(distbuffer[iblock], packbuffer_size, &position, &istart,    1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(distbuffer[iblock], packbuffer_size, &position, &jstart,    1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(distbuffer[iblock], packbuffer_size, &position, &iend,      1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(distbuffer[iblock], packbuffer_size, &position, &jend,      1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(distbuffer[iblock], packbuffer_size, &position, databuffer, blocksize*blocksize, MPI_DOUBLE, MPI_COMM_WORLD);
            
            position = 0;
            for (i = istart; i <= iend; i++) {
                for (j = jstart; j <= jend; j++) {
                    A[i][j] = databuffer[position];
                    position++;
                }
            }
        }
    }
    // END: ================ Asssemble Matrix to Root Process ==================
    
    if (Rank == ROOT) {
        for (i = 0; i < nblocks; i++)
            delete[] distbuffer[i];
        delete[] distbuffer;
    }
    for (i = 0; i < blocks_per_proc; i++)
        delete[] procbuffer[i];
    delete[] procbuffer;
    delete[] databuffer;
    
    for (i = 0; i < ranksize; i++)
        delete[] RankMatrix[i];
    delete[] RankMatrix;
    delete[] blocklist;
}

#endif /* HAVE_MPI */
