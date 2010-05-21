/*******************************************************************************
 * File:        DBCGNS.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _DBCGNS_H
#define	_DBCGNS_H

#include "DB.h"

// Forward Declearation
class DBNETCDF;

class DBCGNS : public DB {
public :
    //! Constructors and Distructors
    DBCGNS();
    DBCGNS(const DBCGNS& other);
    ~DBCGNS();
    //! Operations
    void Read_DB();
    void Write_DB();
    DBCGNS& operator=(const DBCGNS& other);
    int Export_To_NETCDF(DBNETCDF& other);
    int Import_From_NETCDF(const DBNETCDF& other);
    int Share_With_NETCDF(DBNETCDF& other);
private:
    int Read_CGNS_DB(int mode);
    int Read_CGNS_GridFile();
    int Read_CGNS_SolutionFile();
    int Read_CGNS_ParamFile();
    int Write_CGNS_DB(int mode);
    int Write_CGNS_GridFile(int mode);
    int Write_CGNS_SolutionFile();
    int Write_CGNS_ParamFile();
};

#endif	/* _DBCGNS_H */

