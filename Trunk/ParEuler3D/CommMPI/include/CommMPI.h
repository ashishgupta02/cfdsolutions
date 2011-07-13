/*******************************************************************************
 * File:        CommMPI.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _COMMMPI_H
#define	_COMMMPI_H

struct CommMPIGhost {
    int     globalId;
    double  vars[6];
};

struct CommMPIGrad {
    int     globalId;
    double  grads[15];
};

struct CommMPIVec3D {
    int     ids[3]; // contains globalId and matrix_id;
    double  comp[3];
};

void CommMPI_Initialize(int argc, char *argv[]);
void CommMPI_Finalize();
void CommMPI_Handshake(void);
void CommMPI_Get_Ghost_Centroids(void);
void CommMPI_Update_Ghost_Primitives(void);
void CommMPI_Update_Ghost_Gradients(void);
/* Print error message to stderr and exit */
void CommMPI_Error(int procid, const char *fmt, ...);
/* Print warning message to stderr and continue */
void CommMPI_Warn(int procid, const char *fmt, ...);
/* Print information message to stdout and continue */
void CommMPI_Info(int procid, const char *fmt, ...);

#endif	/* _COMMMPI_H */

