/*******************************************************************************
 * File:        Commons.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _COMMONS_H
#define	_COMMONS_H

#include <vector>
#include <string>

#include "Node3D.h"
#include "Cell3D.h"
#include "Face3D.h"
#include "Ghost3D.h"
#include "Grid3D.h"

using namespace std;

// CommMPI Commons
extern int Rank;
extern int RankIO;
extern int NoProc;
extern vector< vector<int> > sendCells;
extern vector<int> recvCount;

// Grid Database
extern Grid3D grid;

// Options for Boundary Condition Types
#define SYMMETRY    1
#define SLIP        2
#define NOSLIP      3
#define INLET       4
#define OUTLET      5

// Options for Boundary Condition Type variants (kind)
#define NO_REVERSE      1
#define DAMP_REVERSE    2
#define VELOCITY        3
#define MDOT            4

// Face bc type numbering
#define INTERNAL    -1
#define UNASSIGNED  -2
#define GHOST       -3

// MeshIO Commons
extern int      Input_MeshIO_Type;
extern string   Input_Grid_File;
extern string   Input_Solution_File;
extern string   Input_Parameter_File;
extern int      Output_MeshIO_Type;
extern string   Output_Grid_File;
extern string   Output_Solution_File;
extern string   Output_Parameter_File;

// Solver Commons
extern int timeStep;
extern int restart;

// Functions Managing Commons
void Commons_Init();
void Commons_Fini();

#endif	/* _COMMONS_H */

