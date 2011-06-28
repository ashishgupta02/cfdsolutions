/*******************************************************************************
 * File:        RestartIO.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef RESTARTIO_H
#define	RESTARTIO_H

void Check_Restart(int Iteration);
void Restart_Writer(const char* filename, int verbose);
void Restart_Reader(const char* filename);

#endif	/* RESTARTIO_H */

