/*******************************************************************************
 * File:        Mesh3D_IO.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/


#ifndef _MESH3D_IO_H
#define	_MESH3D_IO_H

#include <string>
#include "corestruct.h"

// Forward Declaration
class DBCGNS;
class DBNETCDF;
class DBUGRID;

using namespace std;

class Mesh3D_IO {
protected:
    ROOT   *CoreDB;
    int    inputDBType;
    int    outputDBType;
    string inputGFile;
    string inputPFile;
    string inputSFile;
    string outputGFile;
    string outputPFile;
    string outputSFile;
private:
    DBCGNS   *cgnsIO;
    DBNETCDF *netcdfIO;
    DBUGRID  *ugridIO;
public:
    //! Constructors and Distructors
    Mesh3D_IO();
    ~Mesh3D_IO();
    //! Operations
    void Set_DB(ROOT* object);
    ROOT* Get_DB() const;
    void Read_DB();
    void Write_DB();
    void Reset_DB();
    void Set_Input_DataBaseType(int type);
    void Set_InputGrid_Filename(const char* filename);
    void Set_InputParam_Filename(const char* filename);
    void Set_InputSolution_Filename(const char* filename);
    void Set_Output_DataBaseType(int type);
    void Set_OutputGrid_Filename(const char* filename);
    void Set_OutputParam_Filename(const char* filename);
    void Set_OutputSolution_Filename(const char* filename);
};

#endif	/* _MESH3D_IO_H */

