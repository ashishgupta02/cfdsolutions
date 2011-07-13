/*******************************************************************************
 * File:        DBUGRID.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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
#include "ugridIO.h"

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
DBUGRID::DBUGRID()
{
    // DB Constructor will be automatically called first
}

//------------------------------------------------------------------------------
//! Copy Constructor
//------------------------------------------------------------------------------
DBUGRID::DBUGRID(const DBUGRID& other)
{
    // DB Constructor will be automatically called first
    // Copy the Attributes information
    char *tmp = NULL;
    size_t len = 0;
    tmp = other.Get_InputGrid_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputGFile = (char *) malloc((len + 1) * sizeof (char));
        str_blank(inputGFile);
        inputGFile = strcpy(inputGFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_InputParam_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputPFile = (char *) malloc((len + 1) * sizeof (char));
        str_blank(inputPFile);
        inputPFile = strcpy(inputPFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_InputSolution_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        inputSFile = (char *) malloc((len + 1) * sizeof (char));
        str_blank(inputSFile);
        inputSFile = strcpy(inputSFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputGrid_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputGFile = (char *) malloc((len + 1) * sizeof (char));
        str_blank(outputGFile);
        outputGFile = strcpy(outputGFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputParam_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputPFile = (char *) malloc((len + 1) * sizeof (char));
        str_blank(outputPFile);
        outputPFile = strcpy(outputPFile, tmp);
        tmp = NULL;
        len = 0;
    }
    tmp = other.Get_OutputSolution_Filename();
    if (tmp != NULL) {
        len = strlen(tmp);
        outputSFile = (char *) malloc((len + 1) * sizeof (char));
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
DBUGRID::~DBUGRID()
{
    // All code above only
    // DB Distructor will be called now
}

//------------------------------------------------------------------------------
//! Copy assignment operator
//------------------------------------------------------------------------------
DBUGRID& DBUGRID::operator=(const DBUGRID& other)
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
        inputGFile = NULL;
        inputPFile = NULL;
        inputSFile = NULL;
        outputGFile = NULL;
        outputPFile = NULL;
        outputSFile = NULL;

        // Copy the Attributes information
        char *tmp = NULL;
        size_t len = 0;
        tmp = other.Get_InputGrid_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputGFile = (char *) malloc((len + 1) * sizeof (char));
            str_blank(inputGFile);
            inputGFile = strcpy(inputGFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_InputParam_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputPFile = (char *) malloc((len + 1) * sizeof (char));
            str_blank(inputPFile);
            inputPFile = strcpy(inputPFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_InputSolution_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            inputSFile = (char *) malloc((len + 1) * sizeof (char));
            str_blank(inputSFile);
            inputSFile = strcpy(inputSFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputGrid_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputGFile = (char *) malloc((len + 1) * sizeof (char));
            str_blank(outputGFile);
            outputGFile = strcpy(outputGFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputParam_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputPFile = (char *) malloc((len + 1) * sizeof (char));
            str_blank(outputPFile);
            outputPFile = strcpy(outputPFile, tmp);
            tmp = NULL;
            len = 0;
        }
        tmp = other.Get_OutputSolution_Filename();
        if (tmp != NULL) {
            len = strlen(tmp);
            outputSFile = (char *) malloc((len + 1) * sizeof (char));
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
//! Read UGRID Database
//------------------------------------------------------------------------------
void DBUGRID::Read_DB()
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

        if (!Read_UGRID_DB(mode)) {
            // Register now with DB Manager
            state  = 1;
            parent = 1;
        }
    }
}

//------------------------------------------------------------------------------
//! Write UGRID Database
//------------------------------------------------------------------------------
void DBUGRID::Write_DB()
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
        Write_UGRID_DB(mode);
}

//------------------------------------------------------------------------------
//! Communicate with Other Databases
//! Export to CGNS Class Object
//------------------------------------------------------------------------------
int DBUGRID::Export_To_CGNS(DBCGNS& other)
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
int DBUGRID::Import_From_CGNS(const DBCGNS& other)
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
int DBUGRID::Share_With_CGNS(DBCGNS& other)
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
//! Export to NetCDF Class object
//------------------------------------------------------------------------------
int DBUGRID::Export_To_NETCDF(DBNETCDF& other)
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
//! Import from NetCDF Class Object
//------------------------------------------------------------------------------
int DBUGRID::Import_From_NETCDF(const DBNETCDF& other)
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
//! Share with NetCDF Class Object
//------------------------------------------------------------------------------
int DBUGRID::Share_With_NETCDF(DBNETCDF& other)
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
//! Read the UGRID File and Get Maximum information available
//------------------------------------------------------------------------------
int DBUGRID::Read_UGRID_GridFile()
{
    // Initialize the UGRIDIO Library
    UGRIDIO_INIT();

    // Open UGRID Grid File to Read in Mode Read only
    int mode = 1;
    if(open_ugrid(inputGFile, mode))
        return DB_ERROR;

    // Read the UGRID Database and Update the data structue
    // Get the Number of Bases
    CoreDB->nbases = 1;

    CoreDB->bases = new_base(CoreDB->nbases);
    if (CoreDB->bases == NULL)
        return DB_ERROR;
    
    BASE *qbase = NULL;
    qbase = read_ugrid_grid();
    if (qbase != NULL) {
        // Get the value of current base
        CoreDB->bases[0].celldim   = qbase->celldim;
        CoreDB->bases[0].phydim    = qbase->phydim;
        CoreDB->bases[0].baseclass = qbase->baseclass;
        CoreDB->bases[0].ndesc     = qbase->ndesc;
        CoreDB->bases[0].desc      = qbase->desc;
        CoreDB->bases[0].nzones    = qbase->nzones;
        CoreDB->bases[0].zones     = qbase->zones;
        for (int i = 0; i < 5; i++)
            CoreDB->bases[0].baseunits[i] = qbase->baseunits[i];
        str_blank(CoreDB->bases[0].name);
        strcpy(CoreDB->bases[0].name, qbase->name);

        // Now Reset the data and free the qbase
        qbase->zones = NULL;
        qbase->desc  = NULL;
        del_base(qbase);
        qbase = NULL;
    }

    // Close UGRID Database
    if (close_ugrid())
        return DB_ERROR;

    // Finalize the UGRIDIO Library
    UGRIDIO_FINI();

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the UGRID File and Append Solution Data to existing Database
//------------------------------------------------------------------------------
int DBUGRID::Read_UGRID_SolutionFile()
{
    // Initialize the UGRIDIO Library
    UGRIDIO_INIT();

    // Open UGRID Solution File to Read in Mode Read only
    int mode = 1;
    if(open_ugrid(inputSFile, mode))
        return DB_ERROR;

    // Read the UGRID Database and Update the data structue
    if (CoreDB->nbases > 0) {
        int pass = 0;
        BASE *qbase = NULL;
        int nb = 1;

        qbase = read_ugrid_solution();
        if (qbase != NULL) {
            // Get the value of current base
            pass = 0;
            if (CoreDB->bases[nb-1].celldim == qbase->celldim)
                pass++;
            if (CoreDB->bases[nb-1].phydim  == qbase->phydim)
                pass++;
            if (pass == 2) {
                if (CoreDB->bases[nb-1].nzones != qbase->nzones) {
                    del_base(qbase);
                    qbase = NULL;
                    return DB_ERROR;
                }
                // Update the Solution if already exists add new solution
                for (int nz = 1; nz <= CoreDB->bases[nb-1].nzones; nz++) {
                    if (CoreDB->bases[nb-1].zones[nz-1].nsols == 0) {
                        CoreDB->bases[nb-1].zones[nz-1].nsols = qbase->zones[nz-1].nsols;
                        CoreDB->bases[nb-1].zones[nz-1].sols  = qbase->zones[nz-1].sols;
                    } else {
                        SOLUTION *sol = NULL;
                        int nsol = 0;
                        nsol = CoreDB->bases[nb-1].zones[nz-1].nsols;
                        sol  = CoreDB->bases[nb-1].zones[nz-1].sols;

                        // Allocate Additional Resources
                        CoreDB->bases[nb-1].zones[nz-1].nsols = nsol + qbase->zones[nz-1].nsols;
                        CoreDB->bases[nb-1].zones[nz-1].sols = NULL;
                        CoreDB->bases[nb-1].zones[nz-1].sols = new_solution(CoreDB->bases[nb-1].zones[nz-1].nsols);
                        if (CoreDB->bases[nb-1].zones[nz-1].sols == NULL)
                            return DB_ERROR;

                        // Get Back the Orignal Data
                        for (int ns = 1; ns <= nsol; ns++) {
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].id        = sol[ns-1].id;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].location  = sol[ns-1].location;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].size      = sol[ns-1].size;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].dataclass = sol[ns-1].dataclass;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].nflds     = sol[ns-1].nflds;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].flds      = sol[ns-1].flds;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].ndesc     = sol[ns-1].ndesc;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].desc      = sol[ns-1].desc;
                            str_blank(CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].name);
                            strcpy(CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].name, sol[ns-1].name);
                            for (int i = 0; i < 3; i++)
                                for (int j = 0; j < 2; j++)
                                    CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].rind[i][j] = sol[ns-1].rind[i][j];
                            for (int i = 0; i < 5; i++)
                                CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].units[i] = sol[ns-1].units[i];

                            // Set the values to avoid freeing when delete
                            sol[ns-1].nflds = 0;
                            sol[ns-1].flds  = NULL;
                            sol[ns-1].ndesc = 0;
                            sol[ns-1].desc  = NULL;
                        }

                        // Now Add the additional Data
                        for (int ns = nsol+1; ns <= CoreDB->bases[nb-1].zones[nz-1].nsols; ns++) {
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].id        = ns;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].location  = qbase->zones[nz-1].sols[ns-nsol-1].location;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].size      = qbase->zones[nz-1].sols[ns-nsol-1].size;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].dataclass = qbase->zones[nz-1].sols[ns-nsol-1].dataclass;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].nflds     = qbase->zones[nz-1].sols[ns-nsol-1].nflds;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].flds      = qbase->zones[nz-1].sols[ns-nsol-1].flds;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].ndesc     = qbase->zones[nz-1].sols[ns-nsol-1].ndesc;
                            CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].desc      = qbase->zones[nz-1].sols[ns-nsol-1].desc;
                            // Check if names are not same
                            if (strcmp(qbase->zones[nz-1].sols[ns-nsol-1].name, CoreDB->bases[nb-1].zones[nz-1].sols[ns-2].name)) {
                                str_blank(CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].name);
                                strcpy(CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].name, qbase->zones[nz-1].sols[ns-nsol-1].name);
                            }
                            for (int i = 0; i < 3; i++)
                                for (int j = 0; j < 2; j++)
                                    CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].rind[i][j] = qbase->zones[nz-1].sols[ns-nsol-1].rind[i][j];
                            for (int i = 0; i < 5; i++)
                                CoreDB->bases[nb-1].zones[nz-1].sols[ns-1].units[i] = qbase->zones[nz-1].sols[ns-nsol-1].units[i];

                            // Set the values to avoid freeing when delete
                            qbase->zones[nz-1].sols[ns-nsol-1].nflds = 0;
                            qbase->zones[nz-1].sols[ns-nsol-1].flds  = NULL;
                            qbase->zones[nz-1].sols[ns-nsol-1].ndesc = 0;
                            qbase->zones[nz-1].sols[ns-nsol-1].desc  = NULL;
                        }

                        // Free the Solution
                        del_solution(sol);
                        sol = NULL;
                    }
                    qbase->zones[nz-1].nsols = 0;
                    qbase->zones[nz-1].sols  = NULL;
                }
            }

            // Now Reset the data and free the qbase
            del_base(qbase);
            qbase = NULL;
        }
    }

    // Close UGRID Database
    if (close_ugrid())
        return DB_ERROR;

    // Finalize the UGRIDIO Library
    UGRIDIO_FINI();

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the Additional Parameters Describes the Simulation Information
//------------------------------------------------------------------------------
int DBUGRID::Read_UGRID_ParamFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Read the UGRID DataBase File and Updates the CoreDB
//------------------------------------------------------------------------------
int DBUGRID::Read_UGRID_DB(int mode)
{
    // Read Grid File
    if ((mode & 1) == 1) {
        CoreDB = new_root();
        if (Read_UGRID_GridFile()) {
            return DB_ERROR;
        }
    } else
        return DB_ERROR;

    // Read Solution File
    if ((mode & 2) == 2) {
        if (Read_UGRID_SolutionFile()) {
            return DB_ERROR;
        }
    }

    // Read Parameter File
    if ((mode & 4) == 4) {
        if (Read_UGRID_ParamFile()) {
            return DB_ERROR;
        }
    }

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write the Grid and Solution Data into UGRID File.
//! Solution Data is ignored if Solution Filename is available
//------------------------------------------------------------------------------
int DBUGRID::Write_UGRID_GridFile(int mode)
{
    // Initialize the UGRIDIO Library
    UGRIDIO_INIT();

    // Open UGRID Grid File
    if(open_ugrid(outputGFile, outputGIOMode))
        return DB_ERROR;

    // Write the Database if Data is available
    if (CoreDB->nbases > 0) {
        ZONE *z = NULL;
        // Temporary Make Solution Parameter to Zero so avoid solution write
        if((mode &2) == 2) {
            z = new_zone(CoreDB->bases[0].nzones);
            if (z != NULL) {
                for (int i = 0; i < CoreDB->bases[0].nzones; i++) {
                    z[i].nsols = CoreDB->bases[0].zones[i].nsols;
                    CoreDB->bases[0].zones[i].nsols = 0;
                }
            }
        }

        // Now Write Base to File
        if (write_ugrid_grid(&CoreDB->bases[0]))
            return DB_ERROR;

        // Restore back the Solution Parameter
        if((mode &2) == 2) {
            if (z != NULL) {
                for (int i = 0; i < CoreDB->bases[0].nzones; i++) {
                    CoreDB->bases[0].zones[i].nsols = z[i].nsols;
                    z[i].nsols = 0;
                }
                del_zone(z);
                z = NULL;
            }
        }
    }

    // Close UGRID Database
    if (close_ugrid())
        return DB_ERROR;

    // Finalize the UGRIDIO Library
    UGRIDIO_FINI();

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write only Solution Data in UGRID File
//------------------------------------------------------------------------------
int DBUGRID::Write_UGRID_SolutionFile()
{
    // Initialize the UGRIDIO Library
    UGRIDIO_INIT();

    // Open UGRID Solution File
    if(open_ugrid(outputSFile, outputSIOMode))
        return DB_ERROR;
    
    // Write the Database if Data is available
    if (CoreDB->nbases > 0) {
        ZONE *z = NULL;
        // Temporary Make Parameter so only solution is writeen
        z = new_zone(CoreDB->bases[0].nzones);
        if (z != NULL) {
            for (int i = 0; i < CoreDB->bases[0].nzones; i++) {
                z[i].nverts = CoreDB->bases[0].zones[i].nverts;
                z[i].nesets = CoreDB->bases[0].zones[i].nesets;
                z[i].nints  = CoreDB->bases[0].zones[i].nints;
                z[i].nconns = CoreDB->bases[0].zones[i].nconns;
                z[i].nbocos = CoreDB->bases[0].zones[i].nbocos;
                z[i].ndesc  = CoreDB->bases[0].zones[i].ndesc;

                CoreDB->bases[0].zones[i].nverts = 0;
                CoreDB->bases[0].zones[i].esets  = 0;
                CoreDB->bases[0].zones[i].nints  = 0;
                CoreDB->bases[0].zones[i].nconns = 0;
                CoreDB->bases[0].zones[i].nbocos = 0;
                CoreDB->bases[0].zones[i].ndesc  = 0;
            }
        }

        // Now Write Solution Base to File
        if (write_ugrid_solution(&CoreDB->bases[0]))
            return DB_ERROR;

        // Restore back the Parameter
        if (z != NULL) {
            for (int i = 0; i < CoreDB->bases[0].nzones; i++) {
                CoreDB->bases[0].zones[i].nverts = z[i].nverts;
                CoreDB->bases[0].zones[i].nesets = z[i].nesets;
                CoreDB->bases[0].zones[i].nints  = z[i].nints;
                CoreDB->bases[0].zones[i].nconns = z[i].nconns;
                CoreDB->bases[0].zones[i].nbocos = z[i].nbocos;
                CoreDB->bases[0].zones[i].ndesc  = z[i].ndesc;

                z[i].nverts = 0;
                z[i].nesets = 0;
                z[i].nints  = 0;
                z[i].nconns = 0;
                z[i].nbocos = 0;
                z[i].ndesc  = 0;
            }
            del_zone(z);
            z = NULL;
        }
    }

    // Close UGRID Database
    if (close_ugrid())
        return DB_ERROR;

    // Finalize the UGRIDIO Library
    UGRIDIO_FINI();

    return DB_OK;
}

//------------------------------------------------------------------------------
//! Write the Additional Parameters Describes the Simulation Information
//------------------------------------------------------------------------------
int DBUGRID::Write_UGRID_ParamFile()
{
    return DB_OK;
}

//------------------------------------------------------------------------------
//! Writes the UGRID DataBase File from CoreDB
//------------------------------------------------------------------------------
int DBUGRID::Write_UGRID_DB(int mode)
{
    // Write Grid File
    if (((mode & 1) == 1) && (state == 1)) {
        if (Write_UGRID_GridFile(mode)) {
            return DB_ERROR;
        }
    } else
        return DB_ERROR;

    // Write Solution File
    if ((mode & 2) == 2) {
        if (Write_UGRID_SolutionFile()) {
            return DB_ERROR;
        }
    }

    // Write Parameter File
    if ((mode & 4) == 4) {
        if (Write_UGRID_ParamFile()) {
            return DB_ERROR;
        }
    }
    return DB_OK;
}

