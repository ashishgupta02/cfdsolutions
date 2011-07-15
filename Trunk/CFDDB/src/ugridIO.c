/*******************************************************************************
 * File:        ugridIO.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Application Header Files */
#include "Utils.h"
#include "corelib.h"
#include "ugridIO.h"

/* Static Variables */
static FILE *ugridfn      = NULL;

/*---------- UGRIDIO_INIT ----------------------------------------------
 * Initialize the UGRID Data Structure
 *---------------------------------------------------------------------*/
void UGRIDIO_INIT(void) {
    ugridfn = NULL;
}

/*---------- UGRIDIO_FINI ----------------------------------------------
 * Finalize the UGRID Data Structure
 *---------------------------------------------------------------------*/
void UGRIDIO_FINI(void) {
    ugridfn = NULL;
}

/*---------- open_ugrid ------------------------------------------------
 * open a UGRID file
 *  mode:	1 = Read Only
 *		2 = Write Only
 *		3 = Create New Database if non exist
 *---------------------------------------------------------------------*/
int open_ugrid(const char *ugridfile, int mode) {
    
    switch (mode) {
        case 1:
            ugridfn = fopen(ugridfile, "r");
            break;
        case 2:
            ugridfn = fopen(ugridfile, "w");
            break;
        case 3:
            ugridfn = fopen(ugridfile, "w+");
            break;
    }
    
    if (ugridfn == NULL)
        error("open_ugrid: %s", "Unable to open/create file");
    
    return OK;
}

/*---------- close_ugrid -----------------------------------------------
 * Close a UGRID file
 *---------------------------------------------------------------------*/
int close_ugrid() {
    if (ugridfn != NULL) {
        if (fclose(ugridfn) == EOF)
            error("close_ugrid: %s", "Unable to close file");
        ugridfn = NULL;
    }
    return OK;
}

/*---------- read_ugrid_grid -------------------------------------------
 * Read the ugrid file and Get the required Grid Data
 * File Format
 * nnodes, ntri, nquad, ntet, npyramid, nprism, nhex
 * Coordinates # 3 doubles per node
 * Surface connectivity for Triangles # 3 integers per triangle
 * Surface connectivity for Quads # 4 integers per quad
 * Surface tags for Triangles # 1 integer (tag) per triangle
 * Surface tags for Quads # 1 integer (tag) per quad
 * Volume Connectivity for Tets # 4 integers per tet
 * Volume connectivity for Pyramids # 5 integers per pyramid
 * Volume connectivity for Prisms # 6 integers per prism
 * Volume connectivity for Hex # 8 integers per hex
 *---------------------------------------------------------------------*/
BASE* read_ugrid_grid(void) {
    int i, j, k, idum, ne, ne_start, ne_end, nb = 1;
    BASE *qbase = NULL;
    int nnode, ntri, nquad, ntet, npyra, nprism, nhex, nsurf, sflag;
    int pyra_translate[5] = {0, 3, 4, 1, 2};
    int pyra_tmp[5]       = {0, 0, 0, 0, 0};
    int *c2n_tri   = NULL;
    int *c2n_quad  = NULL;
    int *surf_tri  = NULL;
    int *surf_quad = NULL;
    int *nsur_tri  = NULL;
    int *nsur_quad = NULL;

    /* Initialize */
    nnode = ntri = nquad = ntet = npyra = nprism = nhex = nsurf = sflag = 0;
    
    if (ugridfn != NULL) {
        /* Create Base and Assign Default Values */
        qbase = new_base(nb);
        qbase->id = nb;
        /* Set to handle only 3D Volume Grids */
        qbase->celldim = 3;
        qbase->phydim  = 3;
        /* Set all to default values SI Units*/
        qbase->baseunits[0] = 2; /* Kilogram */
        qbase->baseunits[1] = 2; /* Meter */
        qbase->baseunits[2] = 2; /* Second */
        qbase->baseunits[3] = 2; /* Kelvin */
        qbase->baseunits[4] = 2; /* Degree */
        qbase->baseclass = 2; /* Dimensional */
        
        /* Create a Descriptor for Base*/
        qbase->ndesc = 1;
        qbase->desc = new_desc(qbase->ndesc);
        strcpy(qbase->desc[0].name, "Information ");
        qbase->desc[0].desc = (char *) malloc(50*sizeof(char));
        strcpy(qbase->desc[0].desc, "Converted From Simcenter UGrid ASCII Format \0");
        
        /* Read the Data from File and Create Zone*/
        qbase->nzones = 1;
        qbase->zones = NULL;
        qbase->zones = new_zone(qbase->nzones);
        qbase->zones[0].type = 3; /* Unstructured */
        qbase->zones[0].vertflags = 7; /* For 3D Vector */
        qbase->zones[0].idim   = 1;
        qbase->zones[0].dim[0] = 0; /* NVertex3D */
        qbase->zones[0].dim[1] = 0; /* NCell3D */
        qbase->zones[0].dim[2] = 0; /* NBoundVertex */
        qbase->zones[0].dataclass = qbase->baseclass;
        for (i = 0; i < 5; i++)
            qbase->zones[0].units[i] = qbase->baseunits[i];
        
        /* Read the Grid Attributes */
        idum = fscanf(ugridfn, "%d %d %d %d %d %d %d\n",
                &nnode, &ntri, &nquad, &ntet, &npyra, &nprism, &nhex);
        qbase->zones[0].dim[0] = nnode;
        qbase->zones[0].dim[1] = ntet + npyra + nprism + nhex;
        qbase->zones[0].nverts = nnode;
        /* Allocate Memory to create Verts */
        qbase->zones[0].verts = new_vertex(nnode);
        for (i = 0; i < nnode; i++) {
            idum = fscanf(ugridfn, "%lf %lf %lf", &(qbase->zones[0].verts[i].x),
                    &(qbase->zones[0].verts[i].y), &(qbase->zones[0].verts[i].z));
        }

        /* Determine the Number of Element Sets */
        qbase->zones[0].nesets = 0;
        /* Read Triangle connectivity */
        if (ntri > 0) {
            c2n_tri = (int *) calloc(3*ntri, sizeof(int));
            if (c2n_tri == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 1");
            for (i = 0; i < ntri; i++)
                idum = fscanf(ugridfn, "%d %d %d",
                        &c2n_tri[3*i + 0], &c2n_tri[3*i + 1], &c2n_tri[3*i + 2]);
        }
        /* Read Quadrilateral Connectivity */
        if (nquad > 0) {
            c2n_quad = (int *) calloc(4*nquad, sizeof(int));
            if (c2n_quad == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 2");
            for (i = 0; i < nquad; i++)
                idum = fscanf(ugridfn, "%d %d %d %d", &c2n_quad[4*i + 0], &c2n_quad[4*i + 1],
                        &c2n_quad[4*i + 2], &c2n_quad[4*i + 3]);
        }
        /* Read Triangle Surface Tag */
        if (ntri > 0) {
            surf_tri = (int *) calloc(ntri, sizeof(int));
            if (surf_tri == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 3");
            for (i = 0; i < ntri; i++) {
                idum = fscanf(ugridfn, "%d\n", &surf_tri[i]);
                nsurf = MAX(nsurf, surf_tri[i]);
            }
        }
        /* Read Quadrilateral Surface Tag */
        if (nquad > 0) {
            surf_quad = (int *) calloc(nquad, sizeof(int));
            if (surf_quad == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 4");
            for (i = 0; i < nquad; i++) {
                idum = fscanf(ugridfn, "%d\n", &surf_quad[i]);
                nsurf = MAX(nsurf, surf_quad[i]);
            }
        }
        qbase->zones[0].nesets += nsurf;
        
        /* Check for Tetrahera */
        if (ntet > 0)
            qbase->zones[0].nesets++;
        /* Check for Pyramids */
        if (npyra > 0)
            qbase->zones[0].nesets++;
        /* Check for Prism */
        if (nprism > 0)
            qbase->zones[0].nesets++;
        /* Check for Hexahyderon */
        if (nhex > 0)
            qbase->zones[0].nesets++;

        /* Finally Create the Elemsets */
        qbase->zones[0].esets = new_elemset(qbase->zones[0].nesets);
        ne = 0;
        ne_start = 0;
        ne_end   = 0;
        /* Read Tetrahydera Connectivity */
        if (ntet > 0) {
            qbase->zones[0].esets[ne].id = ne+1;
            /* Set the Name */
            strcpy(qbase->zones[0].esets[ne].name, "TETRAHEDRON ");
            /* Set the Type */
            qbase->zones[0].esets[ne].type = 10;
            /* Make this unsorted boundary */
            qbase->zones[0].esets[ne].nbndry = 0;
            /* Set for non boundary Connectivity */
            qbase->zones[0].esets[ne].parent = NULL;
            /* Set the Start and End Point */
            ne_start = ne_end + 1;
            ne_end   = ne_start + ntet - 1;
            qbase->zones[0].esets[ne].start = ne_start;
            qbase->zones[0].esets[ne].end   = ne_end;
            /* Allocate memory to store connectivity */
            qbase->zones[0].esets[ne].csize = 4*ntet;
            qbase->zones[0].esets[ne].conn = (int *) calloc(4*ntet, sizeof(int));
            if (qbase->zones[0].esets[ne].conn == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 5");
            for (i = 0; i < ntet; i++)
                idum = fscanf(ugridfn, "%d %d %d %d",
                        &qbase->zones[0].esets[ne].conn[4*i + 0],
                        &qbase->zones[0].esets[ne].conn[4*i + 1],
                        &qbase->zones[0].esets[ne].conn[4*i + 2],
                        &qbase->zones[0].esets[ne].conn[4*i + 3]);
            ne++;
        }
        /* Read Pyramids Connectivity */
        if (npyra > 0) {
            qbase->zones[0].esets[ne].id = ne+1;
            /* Set the Name */
            strcpy(qbase->zones[0].esets[ne].name, "PYRAMIDS    ");
            /* Set the Type */
            qbase->zones[0].esets[ne].type = 12;
            /* Make this unsorted boundary */
            qbase->zones[0].esets[ne].nbndry = 0;
            /* Set for non boundary Connectivity */
            qbase->zones[0].esets[ne].parent = NULL;
            /* Set the Start and End Point */
            ne_start = ne_end + 1;
            ne_end   = ne_start + npyra - 1;
            qbase->zones[0].esets[ne].start = ne_start;
            qbase->zones[0].esets[ne].end   = ne_end;
            /* Allocate memory to store connectivity */
            qbase->zones[0].esets[ne].csize = 5*npyra;
            qbase->zones[0].esets[ne].conn = (int *) calloc(5*npyra, sizeof(int));
            if (qbase->zones[0].esets[ne].conn == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 5");
            for (i = 0; i < npyra; i++) {
                idum = fscanf(ugridfn, "%d %d %d %d %d",
                        &pyra_tmp[0], &pyra_tmp[1], &pyra_tmp[2], &pyra_tmp[3], &pyra_tmp[4]);
                /* Convert the node ordering to SIDS convention */
                for (j = 0; j < 5; j++)
                    qbase->zones[0].esets[ne].conn[5*i + pyra_translate[j]] = pyra_tmp[j];
            }
            ne++;
        }
        /* Read Prism Connectivity */
        if (nprism > 0) {
            qbase->zones[0].esets[ne].id = ne+1;
            /* Set the Name */
            strcpy(qbase->zones[0].esets[ne].name, "PRISM       ");
            /* Set the Type */
            qbase->zones[0].esets[ne].type = 14;
            /* Make this unsorted boundary */
            qbase->zones[0].esets[ne].nbndry = 0;
            /* Set for non boundary Connectivity */
            qbase->zones[0].esets[ne].parent = NULL;
            /* Set the Start and End Point */
            ne_start = ne_end + 1;
            ne_end   = ne_start + nprism - 1;
            qbase->zones[0].esets[ne].start = ne_start;
            qbase->zones[0].esets[ne].end   = ne_end;
            /* Allocate memory to store connectivity */
            qbase->zones[0].esets[ne].csize = 6*nprism;
            qbase->zones[0].esets[ne].conn = (int *) calloc(6*nprism, sizeof(int));
            if (qbase->zones[0].esets[ne].conn == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 7");
            for (i = 0; i < nprism; i++)
                idum = fscanf(ugridfn, "%d %d %d %d %d %d",
                        &qbase->zones[0].esets[ne].conn[6*i + 0],
                        &qbase->zones[0].esets[ne].conn[6*i + 1],
                        &qbase->zones[0].esets[ne].conn[6*i + 2],
                        &qbase->zones[0].esets[ne].conn[6*i + 3],
                        &qbase->zones[0].esets[ne].conn[6*i + 4],
                        &qbase->zones[0].esets[ne].conn[6*i + 5]);
            ne++;
        }
        /* Read Hexahyderon Connectivity */
        if (nhex > 0) {
            qbase->zones[0].esets[ne].id = ne+1;
            /* Set the Name */
            strcpy(qbase->zones[0].esets[ne].name, "HEXAHEDRON  ");
            /* Set the Type */
            qbase->zones[0].esets[ne].type = 17;
            /* Make this unsorted boundary */
            qbase->zones[0].esets[ne].nbndry = 0;
            /* Set for non boundary Connectivity */
            qbase->zones[0].esets[ne].parent = NULL;
            /* Set the Start and End Point */
            ne_start = ne_end + 1;
            ne_end   = ne_start + nhex - 1;
            qbase->zones[0].esets[ne].start = ne_start;
            qbase->zones[0].esets[ne].end   = ne_end;
            /* Allocate memory to store connectivity */
            qbase->zones[0].esets[ne].csize = 8*nhex;
            qbase->zones[0].esets[ne].conn = (int *) calloc(8*nhex, sizeof(int));
            if (qbase->zones[0].esets[ne].conn == NULL)
                error("read_ugrid_grid: %s", "Memory Allocation Failed - 8");
            for (i = 0; i < nhex; i++)
                idum = fscanf(ugridfn, "%d %d %d %d %d %d %d %d",
                        &qbase->zones[0].esets[ne].conn[8*i + 0],
                        &qbase->zones[0].esets[ne].conn[8*i + 1],
                        &qbase->zones[0].esets[ne].conn[8*i + 2],
                        &qbase->zones[0].esets[ne].conn[8*i + 3],
                        &qbase->zones[0].esets[ne].conn[8*i + 4],
                        &qbase->zones[0].esets[ne].conn[8*i + 5],
                        &qbase->zones[0].esets[ne].conn[8*i + 6],
                        &qbase->zones[0].esets[ne].conn[8*i + 7]);
            ne++;
        }
        /* Separate triangles and quadrilaterals to make surface group */
        nsur_tri = (int *) calloc(nsurf, sizeof(int));
        if (nsur_tri == NULL)
            error("read_ugrid_grid: %s", "Memory Allocation Failed - 9");
        nsur_quad = (int *) calloc(nsurf, sizeof(int));
        if (nsur_quad == NULL)
            error("read_ugrid_grid: %s", "Memory Allocation Failed - 10");
        for (i = 0; i < nsurf; i++)
            nsur_tri[i] = nsur_quad[i] = 0;
        if (ntri > 0)
            for (i = 0; i < ntri; i++)
                nsur_tri[surf_tri[i]-1]++;
        if (nquad > 0)
            for (i = 0; i < nquad; i++)
                nsur_quad[surf_quad[i]-1]++;
        /* Algorithm to separate Triangle, Quadrilateral and MIXED */
        for (i = 0; i < nsurf; i++) {
            /* Triangle Surface */
            if ((nsur_tri[i] > 0) && (nsur_quad[i] == 0)) {
                qbase->zones[0].esets[ne].id = ne + 1;
                /* Set the Type */
                qbase->zones[0].esets[ne].type = 5;
                /* Make this unsorted boundary */
                qbase->zones[0].esets[ne].nbndry = 0;
                /* Set for non boundary Connectivity */
                qbase->zones[0].esets[ne].parent = NULL;
                /* Set the Start and End Point */
                ne_start = ne_end + 1;
                ne_end = ne_start + nsur_tri[i] - 1;
                qbase->zones[0].esets[ne].start = ne_start;
                qbase->zones[0].esets[ne].end = ne_end;
                /* Allocate memory to store connectivity */
                qbase->zones[0].esets[ne].csize = 3 * nsur_tri[i];
                qbase->zones[0].esets[ne].conn = (int *) calloc(3 * nsur_tri[i], sizeof(int));
                if (qbase->zones[0].esets[ne].conn == NULL)
                    error("read_ugrid_grid: %s", "Memory Allocation Failed - 11");
                k = 0;
                for (j = 0; j < ntri; j++) {
                    if (surf_tri[j]-1 != i)
                        continue;
                    qbase->zones[0].esets[ne].conn[3*k + 0] = c2n_tri[3*j + 0];
                    qbase->zones[0].esets[ne].conn[3*k + 1] = c2n_tri[3*j + 1];
                    qbase->zones[0].esets[ne].conn[3*k + 2] = c2n_tri[3*j + 2];
                    k++;
                }
            }
            /* Quadrilateral Surface */
            if ((nsur_quad[i] > 0) && (nsur_tri[i] == 0)) {
                qbase->zones[0].esets[ne].id = ne + 1;
                /* Set the Type */
                qbase->zones[0].esets[ne].type = 7;
                /* Make this unsorted boundary */
                qbase->zones[0].esets[ne].nbndry = 0;
                /* Set for non boundary Connectivity */
                qbase->zones[0].esets[ne].parent = NULL;
                /* Set the Start and End Point */
                ne_start = ne_end + 1;
                ne_end = ne_start + nsur_quad[i] - 1;
                qbase->zones[0].esets[ne].start = ne_start;
                qbase->zones[0].esets[ne].end = ne_end;
                /* Allocate memory to store connectivity */
                qbase->zones[0].esets[ne].csize = 4 * nsur_quad[i];
                qbase->zones[0].esets[ne].conn = (int *) calloc(4 * nsur_quad[i], sizeof(int));
                if (qbase->zones[0].esets[ne].conn == NULL)
                    error("read_ugrid_grid: %s", "Memory Allocation Failed - 12");
                k = 0;
                for (j = 0; j < nquad; j++) {
                    if (surf_quad[j]-1 != i)
                        continue;
                    qbase->zones[0].esets[ne].conn[4*k + 0] = c2n_quad[4*j + 0];
                    qbase->zones[0].esets[ne].conn[4*k + 1] = c2n_quad[4*j + 1];
                    qbase->zones[0].esets[ne].conn[4*k + 2] = c2n_quad[4*j + 2];
                    qbase->zones[0].esets[ne].conn[4*k + 3] = c2n_quad[4*j + 3];
                    k++;
                }
            }
            /* MIXED Element Surface */
            if ((nsur_tri[i] > 0) && (nsur_quad[i] > 0)) {
                qbase->zones[0].esets[ne].id = ne + 1;
                /* Set the Type */
                qbase->zones[0].esets[ne].type = 20;
                /* Make this unsorted boundary */
                qbase->zones[0].esets[ne].nbndry = 0;
                /* Set for non boundary Connectivity */
                qbase->zones[0].esets[ne].parent = NULL;
                /* Set the Start and End Point */
                ne_start = ne_end + 1;
                ne_end = ne_start + nsur_tri[i] + nsur_quad[i] - 1;
                qbase->zones[0].esets[ne].start = ne_start;
                qbase->zones[0].esets[ne].end = ne_end;
                /* Allocate memory to store connectivity */
                qbase->zones[0].esets[ne].csize = 4*nsur_tri[i] + 5*nsur_quad[i];
                qbase->zones[0].esets[ne].conn = (int *) calloc(4*nsur_tri[i] + 5*nsur_quad[i], sizeof(int));
                if (qbase->zones[0].esets[ne].conn == NULL)
                    error("read_ugrid_grid: %s", "Memory Allocation Failed - 13");
                /* First Triangle */
                k = 0;
                for (j = 0; j < ntri; j++) {
                    if (surf_tri[j]-1 != i)
                        continue;
                    qbase->zones[0].esets[ne].conn[4*k + 0] = 5;
                    qbase->zones[0].esets[ne].conn[4*k + 1] = c2n_tri[3*j + 0];
                    qbase->zones[0].esets[ne].conn[4*k + 2] = c2n_tri[3*j + 1];
                    qbase->zones[0].esets[ne].conn[4*k + 3] = c2n_tri[3*j + 2];
                    k++;
                }
                /* Now Quadrilaterals */
                k = 0;
                for (j = 0; j < nquad; j++) {
                    if (surf_quad[j]-1 != i)
                        continue;
                    qbase->zones[0].esets[ne].conn[4*nsur_tri[i] + 5*k + 0] = 7;
                    qbase->zones[0].esets[ne].conn[4*nsur_tri[i] + 5*k + 1] = c2n_quad[4*j + 0];
                    qbase->zones[0].esets[ne].conn[4*nsur_tri[i] + 5*k + 2] = c2n_quad[4*j + 1];
                    qbase->zones[0].esets[ne].conn[4*nsur_tri[i] + 5*k + 3] = c2n_quad[4*j + 2];
                    qbase->zones[0].esets[ne].conn[4*nsur_tri[i] + 5*k + 4] = c2n_quad[4*j + 3];
                    k++;
                }
            }
            ne++;
        }

        /* Create the Interface Data Structure */
        /* 1-to-1 Connectivity Not Valid for UnStructured Grid */
        qbase->zones[0].nints = 0;
        qbase->zones[0].ints  = NULL;

        /* Create the Connectivity Data Structure */
        /* Get Connectivity Between Zones - Only One Zone Supported Setting to NULL */
        qbase->zones[0].nconns = 0;
        qbase->zones[0].conns  = NULL;

        /* Create the Boundary Connectivity Data Structure */
        qbase->zones[0].nbocos = nsurf;
        qbase->zones[0].bocos  = new_boco(qbase->zones[0].nbocos);
        ne  = qbase->zones[0].nesets - nsurf;
        for (i = 0; i < qbase->zones[0].nbocos; i++) {
            qbase->zones[0].bocos[i].id = i+1;
            /* Set the boundary name */
            sprintf(qbase->zones[0].bocos[i].name, "Boundary %d", i+1);
            /* Set the Boundary Type */
            qbase->zones[0].bocos[i].type = 1; /* BCTypeUserDefined */
            /* Set for No Normal Information */
            qbase->zones[0].bocos[i].n_index[0] = 0;
            qbase->zones[0].bocos[i].n_index[1] = 0;
            qbase->zones[0].bocos[i].n_index[2] = 0;
            qbase->zones[0].bocos[i].n_cnt      = 0;
            qbase->zones[0].bocos[i].n_type     = 0; /* DataTypeNull */
            qbase->zones[0].bocos[i].n_list     = NULL;
            /* Set the ptyoe as ElementRange */
            qbase->zones[0].bocos[i].ptype = 6; /* ElementRange */
            qbase->zones[0].bocos[i].npnts = 2;
            qbase->zones[0].bocos[i].pnts = (int *) malloc(2*sizeof(int));
            qbase->zones[0].bocos[i].pnts[0] = qbase->zones[0].esets[ne].start;
            qbase->zones[0].bocos[i].pnts[1] = qbase->zones[0].esets[ne].end;
            ne++;
        }
    }

    /* Create Solution Data Structure */
    qbase->zones[0].nsols = 0;
    qbase->zones[0].sols  = NULL;

    /* Create Discriptor */
    qbase->zones[0].ndesc = 0;
    qbase->zones[0].desc  = NULL;

    /* Free the Used Resources */
    if (c2n_tri != NULL)
        free(c2n_tri);
    if (c2n_quad != NULL)
        free(c2n_quad);
    if (surf_tri != NULL)
        free(surf_tri);
    if (surf_quad != NULL)
        free(surf_quad);
    if (nsur_tri != NULL)
        free(nsur_tri);
    if (nsur_quad != NULL)
        free(nsur_quad);
    
    return qbase;
}

/*---------- read_ugrid_solution ---------------------------------------
 * Read the ugrid file and Get the required Solution data
 *---------------------------------------------------------------------*/
BASE* read_ugrid_solution(void) {
    BASE *qbase = NULL;

    return qbase;
}

/*---------- write_ugrid_grid ------------------------------------------
 * Write the ugrid file and Set the required Grid and Solution data
 *---------------------------------------------------------------------*/
int write_ugrid_grid(BASE *base) {
    return OK;
}

/*---------- write_ugrid_solution --------------------------------------
 * Write the ugrid file and Set the required Solution data
 *---------------------------------------------------------------------*/
int write_ugrid_solution(BASE *base) {
    return OK;
}

