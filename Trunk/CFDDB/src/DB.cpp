/*******************************************************************************
 * File:        DB.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <cstring>
#include <malloc.h>
#include "Utils.h"
#include "corelib.h"
#include "DB.h"

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
DB::DB() {
    // Initialize Attributes
    CoreDB        = NULL;
    inputGFile    = NULL;
    inputPFile    = NULL;
    inputSFile    = NULL;
    outputGFile   = NULL;
    outputPFile   = NULL;
    outputSFile   = NULL;
    parent        = 0;
    state         = 0;
    outputGIOMode = 3;
    outputPIOMode = 3;
    outputSIOMode = 3;
}

//------------------------------------------------------------------------------
//! Copy Constructor: It is deep copy
//------------------------------------------------------------------------------
DB::DB(const DB& other) {
    // Initialize Attributes
    CoreDB      = NULL;
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

    outputGIOMode = other.Get_OuputGrid_IOMode();
    outputPIOMode = other.Get_OuputParam_IOMode();
    outputSIOMode = other.Get_OutputSolution_IOMode();
    
    // Copy the containts of other
    if (other.isParent())
        Copy(other.Get_DB());
    else
        Set_DB(other.Get_DB());
}

//------------------------------------------------------------------------------
//! Distructor
//------------------------------------------------------------------------------
DB::~DB() {
    // Free the memory used for CoreDB
    if (parent)
        Reset_DB();
    
    CoreDB = NULL;
    
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
}

//------------------------------------------------------------------------------
//! = assignment operator
//------------------------------------------------------------------------------
DB& DB::operator=(const DB& other) {
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

        outputGIOMode = other.Get_OuputGrid_IOMode();
        outputPIOMode = other.Get_OuputParam_IOMode();
        outputSIOMode = other.Get_OutputSolution_IOMode();

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
//! Set_DB sets the new DB and then frees the original DB if parent
//! Demoted to child if CoreDB are not same object else no change
//------------------------------------------------------------------------------
void DB::Set_DB(ROOT* object) {
    // To avoid self assignment
    if (CoreDB != object) {
        Reset_DB();
        CoreDB = object;
        if (object != NULL) {
            state  = 1;
            parent = 0;
        } else {
            state  = 0;
            parent = 0;
        }
    }
}

//------------------------------------------------------------------------------
//! Get_DB returns the pointer to the location of CoreDB
//! Else NULL if no data available
//------------------------------------------------------------------------------
ROOT* DB::Get_DB() const {
    if (state)
        return CoreDB;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Shows the Data is Child = 0 or Parent = 1
//------------------------------------------------------------------------------
unsigned int DB::isParent() const {
        return parent;
}

//------------------------------------------------------------------------------
//! Shows the Data state of Object No Data = 0 or Data = 1
//------------------------------------------------------------------------------
unsigned int DB::getState() const{
    return state;
}

//------------------------------------------------------------------------------
//! Calls the internal Reset DB function to free the database and resets
//------------------------------------------------------------------------------
void DB::Reset() {
    // Free the memory used for CoreDB
    if(parent)
        Reset_DB();
    
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
    CoreDB      = NULL;
    inputGFile  = NULL;
    inputPFile  = NULL;
    inputSFile  = NULL;
    outputGFile = NULL;
    outputPFile = NULL;
    outputSFile = NULL;
    parent      = 0;
    state       = 0;
}

//------------------------------------------------------------------------------
//! Sets the Input Grid File Name
//------------------------------------------------------------------------------
void DB::Set_InputGrid_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        inputGFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputGFile);
        inputGFile = strcpy(inputGFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Sets the Input Parameter File Name
//------------------------------------------------------------------------------
void DB::Set_InputParam_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        inputPFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputPFile);
        inputPFile = strcpy(inputPFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Sets the Input Solution File Name
//------------------------------------------------------------------------------
void DB::Set_InputSolution_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        inputSFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(inputSFile);
        inputSFile = strcpy(inputSFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Sets the Output Grid File Name
//------------------------------------------------------------------------------
void DB::Set_OutputGrid_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        outputGFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputGFile);
        outputGFile = strcpy(outputGFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Sets the Output Parameter Files Name
//------------------------------------------------------------------------------
void DB::Set_OutputParam_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        outputPFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputPFile);
        outputPFile = strcpy(outputPFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Sets the Output Solution File Name
//------------------------------------------------------------------------------
void DB::Set_OutputSolution_Filename(const char* filename) {
    if (filename != NULL) {
        size_t len = strlen(filename);
        outputSFile = (char *) malloc((len+1)*sizeof(char));
        str_blank(outputSFile);
        outputSFile = strcpy(outputSFile, filename);
    }
}

//------------------------------------------------------------------------------
//! Retrives the Input Grid File Name
//------------------------------------------------------------------------------
char* DB::Get_InputGrid_Filename() const {
    if (inputGFile != NULL)
        return inputGFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Retrives the Input Parameter File Name
//------------------------------------------------------------------------------
char* DB::Get_InputParam_Filename() const {
    if (inputPFile != NULL)
        return inputPFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Retrives the Input Solution File Name
//------------------------------------------------------------------------------
char* DB::Get_InputSolution_Filename() const {
    if (inputSFile != NULL)
        return inputSFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Retrives the Output Grid File Name
//------------------------------------------------------------------------------
char* DB::Get_OutputGrid_Filename() const {
    if (outputGFile != NULL)
        return outputGFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Retrives the Output Paramter File Name
//------------------------------------------------------------------------------
char* DB::Get_OutputParam_Filename() const {
    if (outputPFile != NULL)
        return outputPFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Retrives the Output Solution File Name
//------------------------------------------------------------------------------
char* DB::Get_OutputSolution_Filename() const {
    if (outputSFile != NULL)
        return outputSFile;
    else
        return NULL;
}

//------------------------------------------------------------------------------
//! Set the Input-Output Mode for Grid File
//! mode = 1 => Read Only
//! mode = 2 => Modify
//! mode = 3 => Write Only
//------------------------------------------------------------------------------
void DB::Set_OutputGrid_IOMode(int mode) {
    if (mode < 0)
        outputGIOMode = 1;
    else
        outputGIOMode = mode;
}

//------------------------------------------------------------------------------
//! Set the Input-Output Mode for Parameter File
//! mode = 1 => Read Only
//! mode = 2 => Modify
//! mode = 3 => Write Only
//------------------------------------------------------------------------------
void DB::Set_OutputParam_IOMode(int mode) {
    if (mode < 0)
        outputPIOMode = 1;
    else
        outputPIOMode = mode;
}

//------------------------------------------------------------------------------
//! Set the Input-Output Mode for Solution File
//! mode = 1 => Read Only
//! mode = 2 => Modify
//! mode = 3 => Write Only
//------------------------------------------------------------------------------
void DB::Set_OutputSolution_IOMode(int mode) {
    if (mode < 0)
        outputSIOMode = 1;
    else
        outputSIOMode = mode;
}

//------------------------------------------------------------------------------
//! Get the Input-Output Mode for Grid File
//! 1 => Read Only
//! 2 => Modify
//! 3 => Write Only
//------------------------------------------------------------------------------
int DB::Get_OuputGrid_IOMode() const {
    return outputGIOMode;
}

//------------------------------------------------------------------------------
//! Get the Input-Output Mode for Parameter File
//! 1 => Read Only
//! 2 => Modify
//! 3 => Write Only
//------------------------------------------------------------------------------
int DB::Get_OuputParam_IOMode() const {
    return outputPIOMode;
}

//------------------------------------------------------------------------------
//! Get the Input-Output Mode for Solution File
//! 1 => Read Only
//! 2 => Modify
//! 3 => Write Only
//------------------------------------------------------------------------------
int DB::Get_OutputSolution_IOMode() const {
    return outputSIOMode;
}

//------------------------------------------------------------------------------
//! Copies the bit-by-bit CoreDB from other DB
//------------------------------------------------------------------------------
void DB::Copy(const ROOT* object) {
    int desc_len = 0;
    
    // To avoid self assignment
    if ((object != NULL) && (CoreDB != object)) {
        // Clean the database
        Reset_DB();
        // ==================================
        // Create the Database
        // ==================================
        CoreDB = new_root();
        if (object->nbases > 0) {
            // ==================================
            // Create the Bases
            // ==================================
            BASE *b_o = object->bases;
            if (b_o != NULL) {
                CoreDB->nbases = object->nbases;
                CoreDB->bases = new_base(CoreDB->nbases);
                BASE *b_n = CoreDB->bases;
                for (int nb = 0; nb < CoreDB->nbases; nb++, b_o++, b_n++) {
                    // Copy the base attributes
                    b_n->id        = b_o->id;
                    b_n->celldim   = b_o->celldim;
                    b_n->phydim    = b_o->phydim;
                    b_n->baseclass = b_o->baseclass;
                    for (int i = 0; i < 5; i++)
                        b_n->baseunits[i] = b_o->baseunits[i];
                    strcpy(b_n->name, b_o->name);
                    if (b_o->nzones > 0) {
                        // ==================================
                        // Create the Zones
                        // ==================================
                        ZONE *z_o = b_o->zones;
                        if (z_o != NULL) {
                            b_n->nzones = b_o->nzones;
                            b_n->zones  = new_zone(b_o->nzones);
                            ZONE *z_n = b_n->zones;
                            for (int nz = 0; b_o->nzones; nz++, z_o++, z_n++) {
                                // Copy the zone attributes
                                z_n->id = z_o->id;
                                strcpy(z_n->name, z_o->name);
                                z_n->type = z_o->type;
                                z_n->idim = z_o->idim;
                                for (int i = 0; i < 3; i++)
                                    z_n->dim[i] = z_o->dim[i];
                                for (int i = 0; i < 5; i++)
                                    z_n->units[i] = z_o->units[i];
                                z_n->dataclass = z_o->dataclass;
                                z_n->datatype  = z_o->datatype;
                                z_n->vertflags = z_o->vertflags;
                                if (z_o->nverts > 0) {
                                    // ==================================
                                    // Create the Vertex
                                    // ==================================
                                    VERTEX *ver_o = z_o->verts;
                                    if (ver_o != NULL) {
                                        z_n->nverts = z_o->nverts;
                                        z_n->verts  = new_vertex(z_o->nverts);
                                        VERTEX *ver_n = z_n->verts;
                                        for (int nv = 0; nv < z_o->nverts; nv++, ver_o++, ver_n++) {
                                            ver_n->id = ver_o->id;
                                            ver_n->x  = ver_o->x;
                                            ver_n->y  = ver_o->y;
                                            ver_n->z  = ver_o->z;
                                            ver_n->w  = ver_o->w;
                                        }
                                    }
                                }
                                if (z_o->nesets > 0) {
                                    // ==================================
                                    // Create the Element Set
                                    // ==================================
                                    ELEMSET *eset_o = z_o->esets;
                                    if (eset_o != NULL) {
                                        z_n->nesets = z_o->nesets;
                                        z_n->esets  = new_elemset(z_o->nesets);
                                        ELEMSET *eset_n = z_n->esets;
                                        for (int ne = 0; ne < z_o->nesets; ne++, eset_o++, eset_n++) {
                                            eset_n->id = eset_o->id;
                                            strcpy(eset_n->name, eset_o->name);
                                            eset_n->type   = eset_o->type;
                                            eset_n->start  = eset_o->start;
                                            eset_n->end    = eset_o->end;
                                            eset_n->nbndry = eset_o->nbndry;
                                            // Copy the Connectivity
                                            if (eset_o->conn != NULL) {
                                                eset_n->csize  = eset_o->csize;
                                                eset_n->conn = (int *) malloc(eset_o->csize*sizeof(int));
                                                if (eset_n->conn != NULL) {
                                                    for (int i = 0; i < eset_o->csize; i++)
                                                        eset_n->conn[i] = eset_o->conn[i];
                                                }
                                            }
                                            // Copy the Parent Connectivity
                                            if (eset_o->parent != NULL) {
                                                eset_n->psize  = eset_o->psize;
                                                eset_n->parent = (int *) malloc(eset_o->psize*sizeof(int));
                                                if (eset_n->parent != NULL) {
                                                    for (int i = 0; i < eset_o->psize; i++)
                                                        eset_n->parent[i] = eset_o->parent[i];
                                                }
                                            }
                                        }
                                    }
                                }
                                if (z_o->nints > 0) {
                                    // ==================================
                                    // Create the Interface
                                    // ==================================
                                    INTERFACE *ints_o = z_o->ints;
                                    if (ints_o != NULL) {
                                        z_n->nints = z_o->nints;
                                        z_n->ints  = new_interface(z_o->nints);
                                        INTERFACE *ints_n = z_n->ints;
                                        for (int ni = 0; ni < z_o->nints; ni++, ints_o++, ints_n++) {
                                            ints_n->id = ints_o->id;
                                            strcpy(ints_n->name, ints_o->name);
                                            for (int i = 0; i < 3; i++) {
                                                ints_n->transform[i] = ints_o->transform[i];
                                                for (int j = 0; j < 2; j++) {
                                                    ints_n->range[i][j] = ints_o->range[i][j];
                                                    ints_n->d_range[i][j] = ints_o->d_range[i][j];
                                                }
                                            }
                                            strcpy(ints_n->d_name, ints_o->d_name);
                                            ints_n->d_zone = ints_o->d_zone;
                                        }
                                    }
                                }
                                if (z_o->nconns > 0) {
                                    // ==================================
                                    // Create Zone Grid Connectivity
                                    // ==================================
                                    CONNECT *conn_o =  z_o->conns;
                                    if (conn_o != NULL) {
                                        z_n->nconns = z_o->nconns;
                                        z_n->conns  = new_connect(z_o->nconns);
                                        CONNECT *conn_n = z_n->conns;
                                        for (int nc = 0; nc < z_o->nconns; nc++, conn_o++, conn_n++) {
                                            conn_n->id       = conn_o->id;
                                            strcpy(conn_n->name, conn_o->name);
                                            conn_n->type     = conn_o->type;
                                            conn_n->location = conn_o->location;
                                            conn_n->ptype    = conn_o->ptype;
                                            conn_n->pnts = (int *) calloc(conn_o->npnts*z_o->idim, sizeof(int));
                                            if (conn_n->pnts != NULL) {
                                                conn_n->npnts    = conn_o->npnts;
                                                for (int i = 0; i < (conn_o->npnts*z_o->idim); i++)
                                                    conn_n->pnts[i] = conn_o->pnts[i];
                                            }
                                            strcpy(conn_n->d_name, conn_o->d_name);
                                            conn_n->d_ztype  = conn_o->d_ztype;
                                            conn_n->d_ptype  = conn_o->d_ptype;
                                            conn_n->d_pnts = (int *) calloc(conn_o->d_npnts*z_o->idim, sizeof(int));
                                            if (conn_n->d_pnts != NULL) {
                                                conn_n->d_npnts  = conn_o->d_npnts;
                                                for (int i = 0; i < (conn_o->d_npnts*z_o->idim); i++)
                                                    conn_n->d_pnts[i] = conn_o->d_pnts[i];
                                            }
                                            conn_n->d_zone   = conn_o->d_zone;
                                        }
                                    }
                                }
                                if (z_o->nbocos > 0) {
                                    // ==================================
                                    // Create Boundary Conditions
                                    // ==================================
                                    BOCO *boco_o = z_o->bocos;
                                    if (boco_o != NULL) {
                                        z_n->nbocos  = z_o->nbocos;
                                        z_n->bocos   = new_boco(z_o->nbocos);
                                        BOCO *boco_n = z_n->bocos;
                                        for (int nb = 0; nb < z_o->nbocos; nb++, boco_o++, boco_n++) {
                                            boco_n->id = boco_o->id;
                                            strcpy(boco_n->name, boco_o->name);
                                            boco_n->type  = boco_o->type;
                                            boco_n->ptype = boco_o->ptype;
                                            boco_n->pnts  =(int *) calloc(boco_o->npnts*z_o->idim, sizeof(int));
                                            if (boco_n->pnts != NULL) {
                                                boco_n->npnts = boco_o->npnts;
                                                for (int i = 0; i < (boco_o->npnts*z_o->idim); i++)
                                                    boco_n->pnts[i] = boco_o->pnts[i];
                                            }
                                            for (int i = 0; i < 3; i++)
                                                boco_n->n_index[i] = boco_o->n_index[i];
                                            boco_n->n_type = boco_o->n_type;
                                            boco_n->n_list = (double *) calloc(boco_o->n_cnt, sizeof(double));
                                            if (boco_n->n_list != NULL) {
                                                boco_n->n_cnt  = boco_o->n_cnt;
                                                for (int i = 0; i < boco_o->n_cnt; i++)
                                                    boco_n->n_list[i] = boco_o->n_list[i];
                                            }
                                        }
                                    }
                                }
                                if (z_o->nsols > 0) {
                                    // ==================================
                                    // Create Solution
                                    // ==================================
                                    SOLUTION *sol_o = z_o->sols;
                                    if (sol_o != NULL) {
                                        z_n->nsols = z_o->nsols;
                                        z_n->sols  = new_solution(z_o->nsols);
                                        SOLUTION *sol_n = z_n->sols;
                                        for (int ns = 0; ns < z_o->nsols; ns++, sol_o++, sol_n++) {
                                            sol_n->id = sol_o->id;
                                            strcpy(sol_n->name, sol_o->name);
                                            sol_n->location = sol_o->location;
                                            for (int i = 0; i < 3; i++)
                                                for (int j = 0; j < 2; j++)
                                                    sol_n->rind[i][j] = sol_o->rind[i][j];
                                            sol_n->size = sol_o->size;
                                            for (int i = 0; i < 5; i++)
                                                sol_n->units[i] = sol_o->units[i];
                                            sol_n->dataclass = sol_o->dataclass;
                                            if (sol_o->nflds > 0) {
                                                // ==================================
                                                // Create Solution Fields
                                                // ==================================
                                                FIELD *fld_o = sol_o->flds;
                                                if (fld_o != NULL) {
                                                    sol_n->nflds = sol_o->nflds;
                                                    sol_n->flds  = new_field(sol_o->nflds, 0);
                                                    FIELD *fld_n = sol_n->flds;
                                                    for (int fd = 0; fd < sol_o->nflds; fd++, fld_o++, fld_n++) {
                                                        fld_n->id = fld_o->id;
                                                        strcpy(fld_n->name, fld_o->name);
                                                        fld_n->datatype = fld_o->datatype;
                                                        for (int i = 0; i < 5; i++) {
                                                            fld_n->units[i]    = fld_o->units[i];
                                                            fld_n->exponent[i] = fld_o->exponent[i];
                                                        }
                                                        for (int i = 0; i < 2; i++)
                                                            fld_n->dataconv[i] = fld_o->dataconv[i];
                                                        fld_n->dataclass = fld_o->dataclass;
                                                        fld_n->convtype  = fld_o->convtype;
                                                        fld_n->exptype   = fld_o->exptype;
                                                        if (fld_o->data != NULL) {
                                                            fld_n->data = (double *) malloc(sol_o->size*sizeof(double));
                                                            if (fld_n->data != NULL)
                                                                for (int i = 0; i < sol_o->size; i++)
                                                                    fld_n->data[i] = fld_o->data[i];
                                                        }
                                                    }
                                                }
                                            }
                                            if (sol_o->ndesc > 0) {
                                                // ==================================
                                                // Create the Solution Descriptors
                                                // ==================================
                                                DESC *des_o = sol_o->desc;
                                                if (des_o != NULL) {
                                                    sol_n->ndesc = sol_o->ndesc;
                                                    sol_n->desc  = new_desc(sol_o->ndesc);
                                                    DESC *des_n  = sol_n->desc;
                                                    for (int nd = 0; nd < sol_o->ndesc; nd++, des_o++, des_n++) {
                                                        des_n->id = des_o->id;
                                                        strcpy(des_n->name, des_o->name);
                                                        if (des_o->desc != NULL) {
                                                            desc_len = strlen(des_o->desc);
                                                            des_n->desc = (char *) malloc(desc_len*sizeof(char));
                                                            strcpy(des_n->desc, des_o->name);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if (z_o->ndesc > 0) {
                                    // ==================================
                                    // Create the Zone Descriptors
                                    // ==================================
                                    DESC *des_o = z_o->desc;
                                    if (des_o != NULL) {
                                        z_n->ndesc = z_o->ndesc;
                                        z_n->desc  = new_desc(z_o->ndesc);
                                        DESC *des_n = z_n->desc;
                                        for (int nd = 0; nd < z_o->ndesc; nd++, des_o++, des_n++) {
                                            des_n->id = des_o->id;
                                            strcpy(des_n->name, des_o->name);
                                            if (des_o->desc != NULL) {
                                                desc_len = strlen(des_o->desc);
                                                des_n->desc = (char *) malloc(desc_len*sizeof(char));
                                                strcpy(des_n->desc, des_o->desc);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (b_o->ndesc > 0) {
                        // ==================================
                        // Create the Base Descriptors
                        // ==================================
                        DESC *des_o = b_o->desc;
                        if (des_o != NULL) {
                            b_n->ndesc = b_o->ndesc;
                            b_n->desc  = new_desc(b_o->ndesc);
                            DESC *des_n = b_n->desc;
                            for (int nd = 0; nd < b_o->ndesc; nd++, des_o++, des_n++) {
                                des_n->id = des_o->id;
                                strcpy(des_n->name, des_o->name);
                                if (des_o->desc != NULL) {
                                    desc_len = strlen(des_o->desc);
                                    des_n->desc = (char *) malloc(desc_len*sizeof(char));
                                    strcpy(des_n->desc, des_o->desc);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (object->ndesc > 0) {
            // ==================================
            // Create the ROOT Descriptors
            // ==================================
            DESC *des_o = object->desc;
            if (des_o != NULL) {
                CoreDB->ndesc  = object->ndesc;
                CoreDB->desc = new_desc(CoreDB->ndesc);
                DESC *des_n = CoreDB->desc;
                for (int nd = 0; nd < object->ndesc; nd++, des_o++, des_n++) {
                    des_n->id = des_o->id;
                    strcpy(des_n->name, des_o->name);
                    if (des_o->desc != NULL) {
                        desc_len = strlen(des_o->desc);
                        des_n->desc = (char *) malloc(desc_len*sizeof(char));
                        strcpy(des_n->desc, des_o->desc);
                    }
                }
            }
        }
        // Change state and parent
        state = 1;
        parent = 1;
    }
}

//------------------------------------------------------------------------------
//! Frees the memory used if it is parent else simple reset of pointer if child
//------------------------------------------------------------------------------
void DB::Reset_DB() {
    if ((CoreDB != NULL) && (parent == 1)) {
        // ==================================
        // Check for Bases
        // ==================================
        BASE *b = CoreDB->bases;
        if (b != NULL) {
            for (int nb = 0; nb < CoreDB->nbases; nb++, b++) {
                // ==================================
                // Check for Zones
                // ==================================
                ZONE *z = b->zones;
                if (z != NULL) {
                    for (int nz = 0; nz < b->nzones; nz++, z++) {
                        // ==================================
                        // Check for vertex
                        // ==================================
                        if (z->verts != NULL)
                            free(z->verts);
                        // ==================================
                        // Check for Element Set
                        // ==================================
                        ELEMSET *eset = z->esets;
                        if (eset != NULL) {
                            for (int ne = 0; ne < z->nesets; ne++, eset++) {
                                // Check for Connectivity
                                if (eset->conn != NULL)
                                    free(eset->conn);
                                // Check for Parents
                                if (eset->parent != NULL)
                                    free(eset->parent);
                            }
                            // Finally Free Element Set
                            free(z->esets);
                        }
                        // ==================================
                        // Check for Interface
                        // ==================================
                        if (z->ints != NULL)
                            free(z->ints);
                        // ==================================
                        // Check for Zone Grid Connectivity
                        // ==================================
                        CONNECT *conn = z->conns;
                        if (conn != NULL) {
                            for (int nc = 0; nc < z->nconns; nc++, conn++) {
                                // Check for Connectivity Points arrays
                                if (conn->pnts != NULL)
                                    free(conn->pnts);
                                // Check for Donor Connectivity Points arrays
                                if (conn->d_pnts != NULL)
                                    free(conn->d_pnts);
                            }
                            // Finally Free Zone Grid Connectivity
                            free(z->conns);
                        }
                        // ==================================
                        // Check for Boundary Conditions
                        // ==================================
                        BOCO *boco = z->bocos;
                        if (boco != NULL) {
                            for (int nbc = 0; nbc < z->nbocos; nbc++, boco++) {
                                // Check for Boundary Points Arrays
                                if (boco->pnts != NULL)
                                    free(boco->pnts);
                                // Check for Boundary Normal Arrays
                                if (boco->n_list != NULL)
                                    free(boco->n_list);
                            }
                            // Finally Free Boundary Condition
                            free(z->bocos);
                        }
                        // ==================================
                        // Check for Solutions Fields
                        // ==================================
                        SOLUTION *sol = z->sols;
                        if (sol != NULL) {
                            for (int ns = 0; ns < z->nsols; ns++, sol++) {
                                // Check for Solution Field
                                FIELD *fld = sol->flds;
                                if (fld != NULL) {
                                    for (int nf = 0; nf < sol->nflds; nf++, fld++) {
                                        if (fld->data != NULL)
                                            free(fld->data);
                                    }
                                    // Finally Free Field
                                    free(sol->flds);
                                }
                                // Check for Solution Descriptor
                                DESC *des = sol->desc;
                                if (des != NULL) {
                                    for (int nd = 0; nd < sol->ndesc; nd++, des++) {
                                        if (des->desc != NULL)
                                            free(des->desc);
                                    }
                                    // Finally Free the Solution Descriptor
                                    free(sol->desc);
                                }
                            }
                            // Finally Free Solution Fields
                            free(z->sols);
                        }
                        // ==================================
                        // Check for Zone Discriptors
                        // ==================================
                        DESC *des = z->desc;
                        if (des != NULL) {
                            for (int nd = 0; nd < z->ndesc; nd++, des++) {
                                if (des->desc != NULL)
                                    free(des->desc);
                            }
                            // Finally Free the Zone Descriptors
                            free(z->desc);
                        }
                    }
                    // Finally Free the Zones
                    free(b->zones);
                }
                // ==================================
                // Check for Base Descriptors
                // ==================================
                DESC *des = b->desc;
                if (des != NULL) {
                    for (int nd = 0; nd < b->ndesc; nd++, des++) {
                        if (des->desc != NULL)
                            free(des->desc);
                    }
                    // Finally Free the Base Descriptors
                    free(b->desc);
                }
            }
            // Finally Free Bases
            free(CoreDB->bases);
        }
        // ==================================
        // Check for ROOT Descriptors
        // ==================================
        DESC *des = CoreDB->desc;
        if (des != NULL) {
            for (int nd = 0; nd < CoreDB->ndesc; nd++, des++) {
                if (des->desc !=NULL)
                    free(des->desc);
            }
            // Finally Free the CoreDB Descriptors
            free(CoreDB->desc);
        }
        // Finally Free CoreDB
        free(CoreDB);
    }
    // Initialize the CoreDB
    CoreDB = NULL;
    state = 0;
    parent = 0;
}
