/*******************************************************************************
 * File:        CommMPI.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <mpi.h>

#include "Utils.h"
#include "Commons.h"
#include "CommMPI.h"
#include "Vector3D.h"

MPI_Datatype MPI_GRAD;
MPI_Datatype MPI_VEC3D;
MPI_Datatype MPI_GHOST;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Initialize(int argc, char *argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Find the number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &NoProc);
    // Find current processor Rank
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    // Set the I/O Rank
    RankIO = 0;

    sendCells.resize(NoProc);
    recvCount.resize(NoProc);

    // Commit custom communication datatypes

    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    int block_lengths[2];
    MPI_Aint displacements[2];

    // MPI_GHOST
    // Declare a dummy instance
    CommMPIGhost dummy;
    displacements[0] = (long) & dummy.globalId - (long) & dummy;
    displacements[1] = (long) & dummy.vars[0] - (long) & dummy;
    block_lengths[0] = 1;
    block_lengths[1] = 6;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_GHOST);
    MPI_Type_commit(&MPI_GHOST);

    // MPI_GRAD
    CommMPIGrad dummy2;
    displacements[0] = (long) & dummy2.globalId - (long) & dummy2;
    displacements[1] = (long) & dummy2.grads[0] - (long) & dummy2;
    block_lengths[0] = 1;
    block_lengths[1] = 15;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_GRAD);
    MPI_Type_commit(&MPI_GRAD);

    // MPI_VEC3D
    CommMPIVec3D dummy4;
    displacements[0] = (long) & dummy4.ids[0] - (long) & dummy4;
    displacements[1] = (long) & dummy4.comp[0] - (long) & dummy4;
    block_lengths[0] = 3;
    block_lengths[1] = 3;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_VEC3D);
    MPI_Type_commit(&MPI_VEC3D);

    return;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Finalize(void) {
    // Free the Infrastructure and Finalize
    MPI_Finalize();
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Handshake(void) {

    int maxGhost = grid.globalCellCount / NoProc * 2;
    int ghosts2receive[NoProc][maxGhost], ghosts2send[NoProc][maxGhost];

    for (int p = 0; p < NoProc; ++p) {
        // Global id's of ghosts to request from each other processors
        // First entry in the array indicates how many of ghosts to be received
        // The rest are global id's
        ghosts2receive[p][0] = 0;
    }

    for (int g = 0; g < grid.ghostCount; ++g) {
        int p = grid.ghost[g].partition;
        ghosts2receive[p][ghosts2receive[p][0] + 1] = grid.ghost[g].globalId;
        ghosts2receive[p][0]++;
    }

    for (int p = 0; p < NoProc; ++p) {
        MPI_Alltoall(ghosts2receive, maxGhost, MPI_INT, ghosts2send, maxGhost,
                MPI_INT, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Transfer data to more efficient containers
    for (int p = 0; p < NoProc; ++p) {
        for (int i = 1; i <= ghosts2send[p][0]; ++i)
            sendCells[p].push_back(ghosts2send[p][i]);
        recvCount[p] = ghosts2receive[p][0];
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Get_Ghost_Centroids(void) {

    for (int p = 0; p < NoProc; ++p) {
        if (Rank != p) {
            CommMPIVec3D sendBuffer[sendCells[p].size()];
            CommMPIVec3D recvBuffer[recvCount[p]];
            int id;
            for (int g = 0; g < (int)sendCells[p].size(); ++g) {
                id = grid.cellGlobal2Local[sendCells[p][g]];
                sendBuffer[g].ids[0] = grid.cell[id].globalId;
                sendBuffer[g].ids[1] = grid.myOffset + id;
                sendBuffer[g].ids[2] = id;
                for (int i = 0; i < 3; ++i)
                    sendBuffer[g].comp[i] = grid.cell[id].centroid[i];
            }

            MPI_Sendrecv(sendBuffer, sendCells[p].size(), MPI_VEC3D, p, 0, 
                    recvBuffer, recvCount[p], MPI_VEC3D, p, MPI_ANY_TAG,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int g = 0; g < recvCount[p]; ++g) {
                id = grid.ghostGlobal2Local[recvBuffer[g].ids[0]];
                grid.ghost[id].matrix_id = recvBuffer[g].ids[1];
                grid.ghost[id].id_in_owner = recvBuffer[g].ids[2];
                for (int i = 0; i < 3; ++i)
                    grid.ghost[id].centroid[i] = recvBuffer[g].comp[i];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Update_Ghost_Primitives(void) {

    for (int p = 0; p < NoProc; ++p) {
        if (Rank != p) {
            CommMPIGhost sendBuffer[sendCells[p].size()];
            CommMPIGhost recvBuffer[recvCount[p]];
            int id;
            for (int g = 0; g < (int)sendCells[p].size(); ++g) {
                id = grid.cellGlobal2Local[sendCells[p][g]];
                sendBuffer[g].globalId = grid.cell[id].globalId;
                sendBuffer[g].vars[0] = grid.cell[id].p;
                sendBuffer[g].vars[1] = grid.cell[id].v[0];
                sendBuffer[g].vars[2] = grid.cell[id].v[1];
                sendBuffer[g].vars[3] = grid.cell[id].v[2];
                sendBuffer[g].vars[4] = grid.cell[id].T;
                sendBuffer[g].vars[5] = grid.cell[id].rho;
            }

            MPI_Sendrecv(sendBuffer, sendCells[p].size(), MPI_GHOST, p, 0,
                    recvBuffer, recvCount[p], MPI_GHOST, p, MPI_ANY_TAG,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int g = 0; g < recvCount[p]; ++g) {
                id = grid.ghostGlobal2Local[recvBuffer[g].globalId];
                if (timeStep == restart + 1) {
                    for (int i = 0; i < 5; ++i) grid.ghost[id].update[i] = 0.;
                } else {
                    grid.ghost[id].update[0] = recvBuffer[g].vars[0] - grid.ghost[id].p;
                    grid.ghost[id].update[1] = recvBuffer[g].vars[1] - grid.ghost[id].v[0];
                    grid.ghost[id].update[2] = recvBuffer[g].vars[2] - grid.ghost[id].v[1];
                    grid.ghost[id].update[3] = recvBuffer[g].vars[3] - grid.ghost[id].v[2];
                    grid.ghost[id].update[4] = recvBuffer[g].vars[4] - grid.ghost[id].T;
                }
                grid.ghost[id].p = recvBuffer[g].vars[0];
                grid.ghost[id].v[0] = recvBuffer[g].vars[1];
                grid.ghost[id].v[1] = recvBuffer[g].vars[2];
                grid.ghost[id].v[2] = recvBuffer[g].vars[3];
                grid.ghost[id].T = recvBuffer[g].vars[4];
                grid.ghost[id].rho = recvBuffer[g].vars[5];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void CommMPI_Update_Ghost_Gradients(void) {

    // Update ghost gradients
    for (int p = 0; p < NoProc; ++p) {
        CommMPIGrad sendBuffer[sendCells[p].size()];
        CommMPIGrad recvBuffer[recvCount[p]];
        int id;
        for (int g = 0; g < (int)sendCells[p].size(); ++g) {
            id = grid.cellGlobal2Local[sendCells[p][g]];
            sendBuffer[g].globalId = grid.cell[id].globalId;
            int count = 0;
            for (int var = 0; var < 5; ++var) {
                for (int comp = 0; comp < 3; ++comp) {
                    sendBuffer[g].grads[count] = grid.cell[id].grad[var][comp];
                    count++;
                }
            }
        }

        MPI_Sendrecv(sendBuffer, sendCells[p].size(), MPI_GRAD, p, 0,
                recvBuffer, recvCount[p], MPI_GRAD, p, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int g = 0; g < recvCount[p]; ++g) {
            id = grid.ghostGlobal2Local[recvBuffer[g].globalId];
            int count = 0;
            for (int var = 0; var < 5; ++var) {
                for (int comp = 0; comp < 3; ++comp) {
                    grid.ghost[id].grad[var].vec[comp] = recvBuffer[g].grads[count];
                    count++;
                }
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
//! Print error message to stderr and exit
//------------------------------------------------------------------------------
void CommMPI_Error(int procid, const char *fmt, ...) {
    if (Rank == procid) {
        va_list args;

        (void) fprintf(stderr, "ERROR  : Rank[%2d ] ", procid);

        va_start(args, fmt);
        (void) vfprintf(stderr, fmt, args);
        va_end(args);

        (void) fprintf(stderr, "\n");

        /* To ensure log files are current */
        (void) fflush(stderr);
    }
    
    CommMPI_Finalize();
    exit(EXIT_FAILURE);
}

//------------------------------------------------------------------------------
//! Print warning message to stderr and continue
//------------------------------------------------------------------------------
void CommMPI_Warn(int procid, const char *fmt, ...) {
    if (Rank == procid) {
        va_list args;

        (void) fprintf(stderr, "WARNING: Rank[%2d ] ", procid);
        va_start(args, fmt);
        (void) vfprintf(stderr, fmt, args);
        va_end(args);

        (void) fprintf(stderr, "\n");
        /* To ensure log files are current */
        (void) fflush(stderr);
    }
}

//------------------------------------------------------------------------------
//! Print information message to stdout and continue
//------------------------------------------------------------------------------
void CommMPI_Info(int procid, const char *fmt, ...) {
    if (Rank == procid) {
        va_list args;

        (void) fprintf(stdout, "INFO   : Rank[%2d ] ", procid);
        va_start(args, fmt);
        (void) vfprintf(stdout, fmt, args);
        va_end(args);

        (void) fprintf(stdout, "\n");
        /* To ensure log files are current */
        (void) fflush(stdout);
    }
}

