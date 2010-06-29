/*******************************************************************************
 * File:        DBNETCDF.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _DBNETCDF_H
#define	_DBNETCDF_H

#include "DB.h"

// Forward Declearation
class DBCGNS;
class DBUGRID;

class DBNETCDF : public DB {
public :
    //! Constructors and Distructors
    DBNETCDF();
    DBNETCDF(const DBNETCDF& other);
    ~DBNETCDF();
    //! Operations
    void Read_DB();
    void Write_DB();
    DBNETCDF& operator=(const DBNETCDF& other);
    int Export_To_CGNS(DBCGNS& other);
    int Import_From_CGNS(const DBCGNS& other);
    int Share_With_CGNS(DBCGNS& other);
    int Export_To_UGRID(DBUGRID& other);
    int Import_From_UGRID(const DBUGRID& other);
    int Share_With_UGRID(DBUGRID& other);
private:
    int Read_NETCDF_DB(int mode);
    int Read_NETCDF_GridFile();
    int Read_NETCDF_SolutionFile();
    int Read_NETCDF_ParamFile();
    int Write_NETCDF_DB(int mode);
    int Write_NETCDF_GridFile(int mode);
    int Write_NETCDF_SolutionFile();
    int Write_NETCDF_ParamFile();
};

#endif	/* _DBNETCDF_H */

