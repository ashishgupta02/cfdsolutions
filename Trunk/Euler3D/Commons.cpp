/*******************************************************************************
 * File:        Commons.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "Commons.h"

// CommMPI Commons
int                   Rank;
int                   NoProc;
vector< vector<int> > sendCells;
vector<int>           recvCount;

// Grid
Grid3D grid;

// MeshIO Commons
int      Input_MeshIO_Type;
string   Input_Grid_File;
string   Input_Solution_File;
string   Input_Parameter_File;
int      Output_MeshIO_Type;
string   Output_Grid_File;
string   Output_Solution_File;
string   Output_Parameter_File;

// Solver Commons
int timeStep;
int restart;

//------------------------------------------------------------------------------
//! Initializer
//------------------------------------------------------------------------------
void Commons_Init() {
    // MPI Commons
    Rank    = -1;
    NoProc  = 0;

    // MeshIO Commons
    Input_MeshIO_Type       = -1;
    Input_Grid_File.clear();
    Input_Solution_File.clear();
    Input_Parameter_File.clear();
    Output_MeshIO_Type      = -1;
    Output_Grid_File.clear();
    Output_Solution_File.clear();
    Output_Parameter_File.clear();
}

//------------------------------------------------------------------------------
//! Finalize
//------------------------------------------------------------------------------
void Commons_Fini() {
    // MPI Commons
    Rank    = -1;
    NoProc  = 0;

    // MeshIO Commons
    Input_MeshIO_Type       = -1;
    Input_Grid_File.clear();
    Input_Solution_File.clear();
    Input_Parameter_File.clear();
    Output_MeshIO_Type      = -1;
    Output_Grid_File.clear();
    Output_Solution_File.clear();
    Output_Parameter_File.clear();
}

