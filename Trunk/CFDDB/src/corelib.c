/*******************************************************************************
 * File:        corelib.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>

/* Application Header Files */
#include "corelib.h"

/*---------- CORELIB_FATAL ---------------------------------------------
 * exit with error message
 *---------------------------------------------------------------------*/

static void CORELIB_FATAL(char *procname, char *errmsg) {
    char *msg = errmsg;
    
    if (NULL == msg) {
        msg = "unknown error";
    }
    fflush(stdout);
    if (NULL == procname || !*procname)
        fprintf(stderr, "%s\n", msg);
    else
        fprintf(stderr, "%s:%s\n", procname, msg);
    
    exit(CORELIB_ERROR);
}

/*---------- new_root -------------------------------------------------
 * create new root(s)
 *---------------------------------------------------------------------*/

ROOT *new_root(void) {
    ROOT *newroot = NULL;

    newroot = (ROOT *) calloc(1, sizeof(ROOT));
    if (NULL == newroot)
        CORELIB_FATAL("new_root", "calloc failed for new root");

    newroot->nbases = 0;
    newroot->bases  = NULL;
    newroot->ndesc  = 0;
    newroot->desc   = NULL;
    
    return newroot;
}

/*---------- new_base -------------------------------------------------
 * create new base(s)
 *---------------------------------------------------------------------*/

BASE *new_base(int count) {
    int n, i;
    BASE *newbase = NULL;
    
    newbase = (BASE *) calloc(count, sizeof(BASE));
    if (NULL == newbase)
        CORELIB_FATAL("new_base", "calloc failed for new base");
    for (n = 0; n < count; n++) {
        newbase[n].id = n + 1;
        sprintf(newbase[n].name, "Base%d", n + 1);
        /* Physical and Cell Dim */
        newbase[n].celldim = 0;
        newbase[n].phydim = 0;
        newbase[n].baseclass = 0;
        for (i = 0; i < 5; i++)
            newbase[n].baseunits[i] = 0;
        newbase[n].nzones = 0;
        newbase[n].zones = NULL;
        newbase[n].ndesc = 0;
        newbase[n].desc = NULL;
    }
    return newbase;
}

/*---------- new_zone -------------------------------------------------
 * create new zone(s)
 *---------------------------------------------------------------------*/

ZONE *new_zone(int count) {
    int n;
    ZONE *newzone = NULL;
    
    newzone = (ZONE *) calloc(count, sizeof(ZONE));
    if (NULL == newzone)
        CORELIB_FATAL("new_zone", "calloc failed for new zones");
    for (n = 0; n < count; n++) {
        newzone[n].id = n + 1;
        sprintf(newzone[n].name, "Zone%d", n + 1);
        /* 2 = Structure */
        newzone[n].type = 2;
        newzone[n].vertflags = 7;
        /* 4 = RealDouble */
        newzone[n].datatype = 4;
    }
    return newzone;
}

/*---------- new_vertex -----------------------------------------------
 * create coordinate array for a zone
 *---------------------------------------------------------------------*/

VERTEX *new_vertex(int nverts) {
    int n;
    VERTEX *verts = NULL;
    
    verts = (VERTEX *) calloc(nverts, sizeof(VERTEX));
    if (NULL == verts)
        CORELIB_FATAL("new_vertex", "calloc failed for new vertex array");
    for (n = 0; n < nverts; n++) {
        verts[n].id = n + 1;
        verts[n].w = 1.0;
    }
    return verts;
}

/*---------- new_elemset ----------------------------------------------
 * create element set array
 *---------------------------------------------------------------------*/

ELEMSET *new_elemset(int nesets) {
    int n;
    ELEMSET *esets = NULL;
    
    esets = (ELEMSET *) calloc(nesets, sizeof(ELEMSET));
    if (NULL == esets)
        CORELIB_FATAL("new_elemset", "calloc failed for new element set array");
    for (n = 0; n < nesets; n++) {
        esets[n].id     = n + 1;
        sprintf(esets[n].name, "ElemSet%d", n + 1);
        esets[n].type   = 0;
        esets[n].start  = 0;
        esets[n].end    = 0;
        esets[n].nbndry = 0;
        esets[n].csize  = 0;
        esets[n].psize  = 0;
        esets[n].conn   = NULL;
        esets[n].parent = NULL;

    }
    return esets;
}

/*---------- new_interface --------------------------------------------
 * create grid 1to1 interface array for a zone
 *---------------------------------------------------------------------*/

INTERFACE *new_interface(int nints) {
    int n;
    INTERFACE *ints = NULL;
    
    ints = (INTERFACE *) calloc(nints, sizeof(INTERFACE));
    if (NULL == ints)
        CORELIB_FATAL("new_interface", "calloc failed for new interface array");
    for (n = 0; n < nints; n++) {
        ints[n].id = n + 1;
        sprintf(ints[n].name, "Interface%d", n + 1);
    }
    return ints;
}

/*---------- new_connect ----------------------------------------------
 * create grid connectivity array for a zone
 *---------------------------------------------------------------------*/

CONNECT *new_connect(int nconns) {
    int n;
    CONNECT *conns = NULL;
    
    conns = (CONNECT *) calloc(nconns, sizeof(CONNECT));
    if (NULL == conns)
        CORELIB_FATAL("new_connect", "calloc failed for new connectivity array");
    for (n = 0; n < nconns; n++) {
        conns[n].id = n + 1;
        sprintf(conns[n].name, "Connectivity%d", n + 1);
    }
    return conns;
}

/*---------- new_boco -------------------------------------------------
 * create boundary condition array
 *---------------------------------------------------------------------*/

BOCO *new_boco(int nbocos) {
    int n;
    BOCO *bocos = NULL;
    
    bocos = (BOCO *) calloc(nbocos, sizeof(BOCO));
    if (NULL == bocos)
        CORELIB_FATAL("new_boco", "calloc failed for new boco array");
    for (n = 0; n < nbocos; n++) {
        bocos[n].id = n + 1;
        sprintf(bocos[n].name, "Boco%d", n + 1);
    }
    return bocos;
}

/*---------- new_solution ---------------------------------------------
 * create solution array for a zone
 *---------------------------------------------------------------------*/

SOLUTION *new_solution(int nsols) {
    int n;
    SOLUTION *sols = NULL;
    
    sols = (SOLUTION *) calloc(nsols, sizeof(SOLUTION));
    if (NULL == sols)
        CORELIB_FATAL("new_solution", "calloc failed for new solution array");
    for (n = 0; n < nsols; n++) {
        sols[n].id = n + 1;
        sprintf(sols[n].name, "FlowSolution%d", n + 1);
    }
    return sols;
}

/*---------- new_field ------------------------------------------------
 * create solution variable array for a zone
 *---------------------------------------------------------------------*/

FIELD *new_field(int nflds, int size) {
    int n;
    FIELD *flds = NULL;
    
    flds = (FIELD *) calloc(nflds, sizeof(FIELD));
    if (NULL == flds)
        CORELIB_FATAL("new_field", "calloc failed for new field array");
    for (n = 0; n < nflds; n++) {
        flds[n].id = n + 1;
        sprintf(flds[n].name, "Field%d", n+1);
        /* 4 = RealDouble */
        flds[n].datatype = 4;
        if (size > 0) {
            flds[n].data = (double *) calloc(size, sizeof(double));
            if (NULL == flds[n].data)
                CORELIB_FATAL("new_field", "calloc failed for field data array");
        }
    }
    return flds;
}

/*---------- new_desc -------------------------------------------------
 * create descriptor array
 *---------------------------------------------------------------------*/

DESC *new_desc(int ndesc) {
    int n;
    DESC *desc = NULL;
    
    desc = (DESC *) calloc(ndesc, sizeof(DESC));
    if (NULL == desc)
        CORELIB_FATAL("new_desc", "calloc failed for new descriptor array");
    for (n = 0; n < ndesc; n++) {
        desc[n].id = n + 1;
        sprintf(desc[n].name, "Descriptor%d", n + 1);
    }
    return desc;
}

/*---------- del_root -------------------------------------------------
 * Delete the memory associated in ROOT
 *---------------------------------------------------------------------*/

int del_root(ROOT *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_base -------------------------------------------------
 * Delete the memory associated in BASE
 *---------------------------------------------------------------------*/

int del_base(BASE *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;
    
    return CORELIB_OK;
}

/*---------- del_zone -------------------------------------------------
 * Delete the memory associated in ZONE
 *---------------------------------------------------------------------*/
int del_zone(ZONE *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_vertex -----------------------------------------------
 * Delete the memory associated in VERTEX
 *---------------------------------------------------------------------*/
int del_vertex(VERTEX *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_elemset ----------------------------------------------
 * Delete the memory associated in ELEMSET
 *---------------------------------------------------------------------*/
int del_elemset(ELEMSET *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_interface --------------------------------------------
 * Delete the memory associated in INTERFACE
 *---------------------------------------------------------------------*/
int del_interface(INTERFACE *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_connect ----------------------------------------------
 * Delete the memory associated in CONNECT
 *---------------------------------------------------------------------*/
int del_connect(CONNECT *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_boco -------------------------------------------------
 * Delete the memory associated in BOCO
 *---------------------------------------------------------------------*/
int del_boco(BOCO *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_desc -------------------------------------------------
 * Delete the memory associated in DESC
 *---------------------------------------------------------------------*/
int del_desc(DESC *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_solution ----------------------------------------------
 * Delete the memory associated in SOLUTION
 *---------------------------------------------------------------------*/
int del_solution(SOLUTION *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}

/*---------- del_field ------------------------------------------------
 * Delete the memory associated in FIELD
 *---------------------------------------------------------------------*/
int del_field(FIELD *POINTER) {
    if (POINTER == NULL) return CORELIB_OK;

    return CORELIB_OK;
}
