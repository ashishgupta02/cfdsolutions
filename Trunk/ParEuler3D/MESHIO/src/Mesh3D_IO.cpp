/*******************************************************************************
 * File:        Mesh3D_IO.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "DBCGNS.h"
#include "DBNETCDF.h"
#include "DBUGRID.h"
#include "Mesh3D_IO.h"
#include "Utils.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
Mesh3D_IO::Mesh3D_IO() {
    // Initialize Attributes
    inputDBType   = 0;
    outputDBType  = 0;
    CoreDB        = NULL;
    cgnsIO        = NULL;
    netcdfIO      = NULL;
    ugridIO       = NULL;
}

//------------------------------------------------------------------------------
//! Distructor
//------------------------------------------------------------------------------
Mesh3D_IO::~Mesh3D_IO() {
    // Reset the Pointer of CoreDB
    CoreDB = NULL;
}

//------------------------------------------------------------------------------
//! Set_DB sets the new DB and then frees the original DB 
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_DB(ROOT* object) {
    // To avoid self assignment
    if (CoreDB != object) {
        Reset_DB();
        CoreDB = object;
    }
}

//------------------------------------------------------------------------------
//! Get_DB returns the pointer to the location of CoreDB
//------------------------------------------------------------------------------
ROOT* Mesh3D_IO::Get_DB() const {
    return CoreDB;
}

//------------------------------------------------------------------------------
//! Sets the Input Database Type
//! Type: 0 = CGNS, 1 = NETCDF, 2 = UGRID
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_Input_DataBaseType(int type) {
    if (type >= 0 || type < 2)
        inputDBType = type;
    else
        error("MESH3D_IO: %s", "Invalid Input Database Type");
}

//------------------------------------------------------------------------------
//! Sets the Output Database Type
//! Type: 0 = CGNS, 1 = NETCDF, 2 = UGRID
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_Output_DataBaseType(int type) {
    if (type >= 0 || type < 2)
        outputDBType = type;
    else
        error("MESH3D_IO: %s", "Invalid Output Database Type");
}

//------------------------------------------------------------------------------
//! Sets the Input Grid File Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_InputGrid_Filename(const char* filename) {
    if (filename != NULL)
        inputGFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Input Grid Filename");
}

//------------------------------------------------------------------------------
//! Sets the Input Parameter File Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_InputParam_Filename(const char* filename) {
    if (filename != NULL)
        inputPFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Input Parameter Filename");
}

//------------------------------------------------------------------------------
//! Sets the Input Solution File Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_InputSolution_Filename(const char* filename) {
    if (filename != NULL)
        inputSFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Input Solution Filename");
}

//------------------------------------------------------------------------------
//! Sets the Output Grid File Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_OutputGrid_Filename(const char* filename) {
    if (filename != NULL)
        outputGFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Output Grid Filename");
}

//------------------------------------------------------------------------------
//! Sets the Output Parameter Files Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_OutputParam_Filename(const char* filename) {
    if (filename != NULL)
        outputPFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Output Parameter Filename");
}

//------------------------------------------------------------------------------
//! Sets the Output Solution File Name
//------------------------------------------------------------------------------
void Mesh3D_IO::Set_OutputSolution_Filename(const char* filename) {
    if (filename != NULL)
        outputSFile.append(filename);
    else
        error("MESH3D_IO: %s", "Invalid Output Solution Filename");
}

//------------------------------------------------------------------------------
//! Read Database
//------------------------------------------------------------------------------
void Mesh3D_IO::Read_DB() {
    switch(inputDBType) {
        // CGNS
        case 0:
            cgnsIO = new DBCGNS;
            cgnsIO->Set_InputGrid_Filename(inputGFile.c_str()) ;
            cgnsIO->Set_InputParam_Filename(inputPFile.c_str());
            cgnsIO->Set_InputSolution_Filename(inputSFile.c_str());
            cgnsIO->Read_DB();
            CoreDB = cgnsIO->Get_DB();
            break;
        // NETCDF
        case 1:
            netcdfIO = new DBNETCDF;
            netcdfIO->Set_InputGrid_Filename(inputGFile.c_str()) ;
            netcdfIO->Set_InputParam_Filename(inputPFile.c_str());
            netcdfIO->Set_InputSolution_Filename(inputSFile.c_str());
            netcdfIO->Read_DB();
            CoreDB = netcdfIO->Get_DB();
            break;
        // UGRID
        case 2:
            ugridIO = new DBUGRID;
            ugridIO->Set_InputGrid_Filename(inputGFile.c_str()) ;
            ugridIO->Set_InputParam_Filename(inputPFile.c_str());
            ugridIO->Set_InputSolution_Filename(inputSFile.c_str());
            ugridIO->Read_DB();
            CoreDB = ugridIO->Get_DB();
            break;
    }

    // Check if Reading DB was successful
    if (CoreDB == NULL)
        error("MESH3D_IO: %s", "Reading Database Failed !");
}

//------------------------------------------------------------------------------
//! Write Database
//------------------------------------------------------------------------------
void Mesh3D_IO::Write_DB() {
    if (CoreDB == NULL)
        error("MESH3D_IO: %s", "No Data - Writing Database Failed !");

    switch(outputDBType) {
        // CGNS
        case 0:
            if (cgnsIO == NULL) {
                cgnsIO = new DBCGNS;
                cgnsIO->Set_DB(CoreDB);
            }
            cgnsIO->Set_OutputGrid_Filename(outputGFile.c_str());
            cgnsIO->Set_OutputSolution_Filename(outputSFile.c_str());
            cgnsIO->Set_OutputParam_Filename(outputPFile.c_str());
            // Set to Overwrite Mode
            cgnsIO->Set_OutputGrid_IOMode(3);
            cgnsIO->Set_OutputSolution_IOMode(3);
            cgnsIO->Set_OutputParam_IOMode(3);
            cgnsIO->Write_DB();
            break;
        // NETCDF
        case 1:
            if (netcdfIO == NULL) {
                netcdfIO = new DBNETCDF;
                netcdfIO->Set_DB(CoreDB);
            }
            netcdfIO->Set_OutputGrid_Filename(outputGFile.c_str());
            netcdfIO->Set_OutputSolution_Filename(outputSFile.c_str());
            netcdfIO->Set_OutputParam_Filename(outputPFile.c_str());
            // Set to Overwrite Mode
            netcdfIO->Set_OutputGrid_IOMode(3);
            netcdfIO->Set_OutputSolution_IOMode(3);
            netcdfIO->Set_OutputParam_IOMode(3);
            netcdfIO->Write_DB();
            break;
        // UGRID
        case 2:
            if (ugridIO == NULL) {
                ugridIO = new DBUGRID;
                ugridIO->Set_DB(CoreDB);
            }
            ugridIO->Set_OutputGrid_Filename(outputGFile.c_str());
            ugridIO->Set_OutputSolution_Filename(outputSFile.c_str());
            ugridIO->Set_OutputParam_Filename(outputPFile.c_str());
            // Set to Overwrite Mode
            ugridIO->Set_OutputGrid_IOMode(3);
            ugridIO->Set_OutputSolution_IOMode(3);
            ugridIO->Set_OutputParam_IOMode(3);
            ugridIO->Write_DB();
            break;
    }
}

//------------------------------------------------------------------------------
//! Reset Database
//------------------------------------------------------------------------------
void Mesh3D_IO::Reset_DB() {
    if (CoreDB != NULL) {
        switch(inputDBType) {
            case 0: // CGNS
                if (cgnsIO != NULL) {
                    cgnsIO->Reset();
                    delete cgnsIO;
                    cgnsIO = NULL;
                }
                break;
            case 1: // NETCDF
                if (netcdfIO != NULL) {
                    netcdfIO->Reset();
                    delete netcdfIO;
                    netcdfIO = NULL;
                }
                break;
            case 2: // UGRID
                if (ugridIO != NULL) {
                    ugridIO->Reset();
                    delete ugridIO;
                    ugridIO = NULL;
                }
                break;
        }
    }
    CoreDB = NULL;
}

