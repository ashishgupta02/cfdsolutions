/*******************************************************************************
 * File:        DBUGRID.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _DBUGRID_H
#define	_DBUGRID_H

#include "DB.h"

// Forward Declearation
class DBCGNS;
class DBNETCDF;

class DBUGRID : public DB {
public :
    //! Constructors and Distructors
    DBUGRID(); //Done
    DBUGRID(const DBUGRID& other); //Done
    ~DBUGRID(); //Done
    //! Operations
    void Read_DB();
    void Write_DB();
    DBUGRID& operator=(const DBUGRID& other); //Done
    int Export_To_CGNS(DBCGNS& other);
    int Import_From_CGNS(const DBCGNS& other);
    int Share_With_CGNS(DBCGNS& other);
    int Export_To_NETCDF(DBNETCDF& other);
    int Import_From_NETCDF(const DBNETCDF& other);
    int Share_With_NETCDF(DBNETCDF& other);
private:
    int Read_UGRID_DB(int mode);
    int Read_UGRID_GridFile();
    int Read_UGRID_SolutionFile();
    int Read_UGRID_ParamFile();
    int Write_UGRID_DB(int mode);
    int Write_UGRID_GridFile(int mode);
    int Write_UGRID_SolutionFile();
    int Write_UGRID_ParamFile();
};

#endif	/* _DBUGRID_H */

