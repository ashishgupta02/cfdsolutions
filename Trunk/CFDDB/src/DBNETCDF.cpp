/*******************************************************************************
 * File:        DBNETCDF.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <cstring>
#include <malloc.h>
#include <stdlib.h>
#include "Utils.h"
#include "DBMANAGER.h"
#include "DBERROR.h"
#include "DBCGNS.h"
#include "DBNETCDF.h"
#include "DBUGRID.h"
// Database Specific Header
#include "corelib.h"
#include "netcdfIO.h"

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
DBNETCDF::DBNETCDF()
{
    // DB Constructor will be automatically called first
}

//------------------------------------------------------------------------------
//! Copy Constructor
//------------------------------------------------------------------------------
DBNETCDF::DBNETCDF(const DBNETCDF& other)
{
    // DB Constructor will be automatically called first
    // Copy the Attributes information
    char *tmp  = NULL;
    size_t len = 0;
    tmp = other.Get_InputGrid_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputGFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputGFile);
        inputGFile = strcpy(inputGFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_InputParam_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputPFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputPFile);
        inputPFile = strcpy(inputPFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_InputSolution_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputSFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputSFile);
        inputSFile = strcpy(inputSFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputGrid_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputGFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputGFile);
        outputGFile = strcpy(outputGFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputParam_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputPFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputPFile);
        outputPFile = strcpy(outputPFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputSolution_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputSFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputSFile);
        outputSFile = strcpy(outputSFile, tmp);
        tmp = NULL;
        len = 0;
    }

    // Copy the containts of other
    if (other.isParent())
        Copy(other.Get_DB());
    else
        Set_DB(other.Get_DB());
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
DBNETCDF::~DBNETCDF()
{
    // All code above only
    // DB Distructor will be called now
}

//------------------------------------------------------------------------------
//! Copy assignment operator
//------------------------------------------------------------------------------
DBNETCDF& DBNETCDF::operator=(const DBNETCDF& other)
{
    // Avoid Self Assignment
    if (&other != this) {
        // Free the memory used for cstrings
        if (inputGFile != NULL)
            free(inputGFile);
        if (inputPFile != NULL)
            free(inputPFile);
        if (inputSFile != NULL)
            free(inputSFile);
        if (outputGFile != NULL)
            free(outputGFile);
        if (outputPFile != NULL)
            free(outputPFile);
        if (outputSFile != NULL)
            free(outputSFile);

        // Initialize the attributes
        inputGFile  = NULL;
        inputPFile  = NULL;
        inputSFile  = NULL;
        outputGFile = NULL;
        outputPFile = NULL;
        outputSFile = NULL;

        // Copy the Attributes information
        char *tmp  = NULL;
        size_t len = 0;
        tmp = other.Get_InputGrid_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputGFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(inputGFile);
            inputGFile = strcpy(inputGFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_InputParam_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputPFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(inputPFile);
            inputPFile = strcpy(inputPFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_InputSolution_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputSFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(inputSFile);
            inputSFile = strcpy(inputSFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputGrid_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputGFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(outputGFile);
            outputGFile = strcpy(outputGFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputParam_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputPFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(outputPFile);
            outputPFile = strcpy(outputPFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputSolution_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputSFile = (char *) malloc((len+1)*sizeof(char));
            str_blank(outputSFile);
            outputSFile = strcpy(outputSFile, tmp);
            tmp = NULL;
            len = 0;
        }
        
        // Check if CoreDB is not resigned back after multiple share
        if (other.Get_DB() == NULL) {
            Set_DB(NULL);
        } else if (CoreDB != other.Get_DB()) {
            if (parent)
                Reset_DB();
            CoreDB = NULL;
            // Copy the containts of other
            if (other.isParent())
                Copy(other.Get_DB());
            else
                Set_DB(other.Get_DB());
        } else {
            if (other.isParent()) {
                Set_DB(NULL);
                Copy(other.Get_DB());
            }
        }
    }
    return *this;
}

//------------------------------------------------------------------------------
//! Read NETCDF Database
//------------------------------------------------------------------------------
void DBNETCDF::Read_DB()
{
    int mode = 0;
    // Check the inputs files names
    if ((inputGFile != NULL) && (file_exists(inputGFile)))
        mode |= 1;
    if ((inputSFile != NULL) && (file_exists(inputSFile)))
        mode |= 2;
    if ((inputPFile != NULL) && (file_exists(inputPFile)))
        mode |= 4;

    // Read only if Grid File is available
    if ((mode & 1) == 1) {
        // Reset the database if it contains data
        if (CoreDB != NULL)
            Set_DB(NULL);
        
        if (!Read_NETCDF_DB(mode)) {
            // Register now with DB Manager
            state  = 1;
            parent = 1;
        }
    }
}

//------------------------------------------------------------------------------
//! Write NETCDF Database
//------------------------------------------------------------------------------
void DBNETCDF::Write_DB()
{
    int mode = 0;
    // Check the output file names
    if ((outputGFile != NULL) && (!file_exists(outputGFile)))
        mode |= 1;
    if ((outputSFile != NULL) && (!file_exists(outputSFile)))
        mode |= 2;
    if ((outputPFile != NULL) && (!file_exists(outputPFile)))
        mode |= 4;

    // Write only if grid the file doesnt exist rest are skiped if exists
    if ((mode & 1) == 1)
        Write_NETCDF_DB(mode);
}

//------------------------------------------------------------------------------
//! Communicate with Other Databases
//! Export to CGNS Class object
//------------------------------------------------------------------------------
int DBNETCDF::Export_To_CGNS(DBCGNS& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Copy the containts to other
    if (CoreDB == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        other.Copy(CoreDB);
    } else if (!other.isParent()) {
        other.Set_DB(NULL);
        other.Copy(CoreDB);
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Import from CGNS Class Object
//------------------------------------------------------------------------------
int DBNETCDF::Import_From_CGNS(const DBCGNS& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Copy the containts of other
    if (other.Get_DB() == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        Copy(other.Get_DB());
    } else if (!parent) {
        Set_DB(NULL);
        Copy(other.Get_DB());
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Share with CGNS Class Object
//------------------------------------------------------------------------------
int DBNETCDF::Share_With_CGNS(DBCGNS& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Share the containts with other
    if (CoreDB == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        other.Set_DB(NULL);
        other.Set_DB(CoreDB);
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Communicate with Other Databases
//! Export to UGRID Class object
//------------------------------------------------------------------------------
int DBNETCDF::Export_To_UGRID(DBUGRID& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Copy the containts to other
    if (CoreDB == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        other.Copy(CoreDB);
    } else if (!other.isParent()) {
        other.Set_DB(NULL);
        other.Copy(CoreDB);
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Import from UGRID Class Object
//------------------------------------------------------------------------------
int DBNETCDF::Import_From_UGRID(const DBUGRID& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Copy the containts of other
    if (other.Get_DB() == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        Copy(other.Get_DB());
    } else if (!parent) {
        Set_DB(NULL);
        Copy(other.Get_DB());
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Share with UGRID Class Object
//------------------------------------------------------------------------------
int DBNETCDF::Share_With_UGRID(DBUGRID& other)
{
    // Check if CoreDB is not resigned back after multiple share
    // Share the containts with other
    if (CoreDB == NULL) {
        return DB_ERROR;
    } else if (CoreDB != other.Get_DB()) {
        other.Set_DB(NULL);
        other.Set_DB(CoreDB);
    }
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the NETCDF File and Get Maximum information available
//------------------------------------------------------------------------------
int DBNETCDF::Read_NETCDF_GridFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the NETCDF File and Append Solution Data to existing Database
//------------------------------------------------------------------------------
int DBNETCDF::Read_NETCDF_SolutionFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the Additional Parameters Describes the Simulation Information
//------------------------------------------------------------------------------
int DBNETCDF::Read_NETCDF_ParamFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the NetCDF DataBase File and Updates the CoreDB
//------------------------------------------------------------------------------
int DBNETCDF::Read_NETCDF_DB(int mode)
{
    // Read Grid File
    if ((mode & 1) == 1) {
        CoreDB = new_root();
        if (Read_NETCDF_GridFile()) {
            return DB_ERROR;
        }
    } else
        return DB_ERROR;

    // Read Solution File
    if ((mode & 2) == 2) {
        if (Read_NETCDF_SolutionFile()) {
            return DB_ERROR;
        }
    }

    // Read Parameter File
    if ((mode & 4) == 4) {
        if (Read_NETCDF_ParamFile()) {
            return DB_ERROR;
        }
    }

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write the Grid and Solution Data into NETCDF File.
//! Solution Data is ignored if Solution Filename is available
//------------------------------------------------------------------------------
int DBNETCDF::Write_NETCDF_GridFile(int mode)
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write only Solution Data in NETCDF File
//------------------------------------------------------------------------------
int DBNETCDF::Write_NETCDF_SolutionFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write the Additional Parameters Describes the Simulation Information
//------------------------------------------------------------------------------
int DBNETCDF::Write_NETCDF_ParamFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Writes the NetCDF DataBase File from CoreDB
//------------------------------------------------------------------------------
int DBNETCDF::Write_NETCDF_DB(int mode)
{
    // Write Grid File
    if (((mode & 1) == 1) && (state == 1)) {
        if (Write_NETCDF_GridFile(mode)) {
            return DB_ERROR;
        }
    } else
        return DB_ERROR;

    // Write Solution File
    if ((mode & 2) == 2) {
        if (Write_NETCDF_SolutionFile()) {
            return DB_ERROR;
        }
    }

    // Write Parameter File
    if ((mode & 4) == 4) {
        if (Write_NETCDF_ParamFile()) {
            return DB_ERROR;
        }
    }
    return DB_OK;
}

