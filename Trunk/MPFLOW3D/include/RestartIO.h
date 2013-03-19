/*******************************************************************************
 * File:        RestartIO.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _RESTARTIO_H
#define	_RESTARTIO_H

void Check_Restart(int Iteration);
void Restart_Writer(const char* filename, int verbose);
void Restart_Reader(const char* filename);

#endif	/* _RESTARTIO_H */

