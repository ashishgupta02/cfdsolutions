/*******************************************************************************
 * File:        cgnsIO.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cgnslib.h>

/* Application Header Files */
#include "Utils.h"
#include "corelib.h"
#include "cgnsIO.h"

/* Static Variables */
static int cgnsfn        = 0;
static char basename[33] = "";
static int baseclass     = 0;
static int baseunits[5]  = {0, 0, 0, 0, 0};
static int cgnsbase      = 0;
static int CellDim       = 0;
static int PhyDim        = 0;
static int nZones        = 0;
static ZONE *Zones;
static int nDescs        = 0;
static DESC *Descs;

/*---------- CGNSIO_FATAL ---------------------------------------------
 * exit with error message
 *---------------------------------------------------------------------*/
void CGNSIO_FATAL(char *procname, char *errmsg) {
    char *msg = errmsg;
    
    if (NULL == msg) {
        msg = (char *)cg_get_error();
        if (NULL == msg || !*msg)
            msg = "unknown error";
    }
    fflush(stdout);
    if (NULL == procname || !*procname)
        fprintf(stderr, "%s\n", msg);
    else
        fprintf(stderr, "%s:%s\n", procname, msg);
    if (cgnsfn) cg_close(cgnsfn);
    exit(CG_ERROR);
}

/*---------- CGNSIO_INIT ----------------------------------------------
 * Initialise the CGNS Data Structure
 *---------------------------------------------------------------------*/
void CGNSIO_INIT(void) {
    short i = 0;
    
    /* Initialise the Static Variables */
    nZones  = 0;
    CellDim = 0;
    PhyDim  = 0;
    Zones   = NULL;
    nDescs  = 0;
    Descs   = NULL;
    baseclass = 0;
    for (i = 0; i < 5; i++) baseunits[i] = 0;
    cgnsbase = 0;
    cgnsfn = 0;
    str_blank(basename);
}

/*---------- CGNSIO_FINI ----------------------------------------------
 * Finalize the CGNS Data Structure
 *---------------------------------------------------------------------*/
void CGNSIO_FINI(void) {
    short i = 0;

    /* Finalize the Static Variables */
    nZones  = 0;
    CellDim = 0;
    PhyDim  = 0;
    Zones   = NULL;
    nDescs  = 0;
    Descs   = NULL;
    baseclass = 0;
    for (i = 0; i < 5; i++) baseunits[i] = 0;
    cgnsbase = 0;
    cgnsfn = 0;
    str_blank(basename);
}

/*---------- open_cgns ------------------------------------------------
 * open a CGNS file
 *  mode:	1 = Read Only
 *		2 = Write Only
 *		3 = Create New Database if non exist
 *---------------------------------------------------------------------*/
int open_cgns(const char *cgnsfile, int mode) {
    int ier = 0;
    
    switch (mode) {
        case 1:
            ier = cg_open(cgnsfile, MODE_READ, &cgnsfn);
            break;
        case 2:
            ier = cg_open(cgnsfile, MODE_MODIFY, &cgnsfn);
            break;
        case 3:
            ier = cg_open(cgnsfile, MODE_WRITE, &cgnsfn);
            break;
    }
    
    if (ier) CGNSIO_FATAL("open_cgns", NULL);
    
    return CG_OK;
}

/*---------- close_cgns -----------------------------------------------
 * Close a CGNS file
 *---------------------------------------------------------------------*/
int close_cgns() {
    if (cgnsfn != 0) {
        if (cg_close(cgnsfn))
            CGNSIO_FATAL("close_cgns", NULL);
        cgnsfn = 0;
    }
    return CG_OK;
}

/*---------- get_cgns_nbases -------------------------------------------
 * Get the Number of Bases available in CGNS File
 *---------------------------------------------------------------------*/
int get_cgns_nbases(void) {
    int nbases = 0;
    if (cgnsfn != 0) {
        if(cg_nbases(cgnsfn, &nbases))
            CGNSIO_FATAL("get_cgns_nbases", NULL);
    }
    return nbases;
}

/*---------- get_cgns_base ---------------------------------------------
 * Read the CGNS file and Get the required Base data
 *---------------------------------------------------------------------*/
BASE* get_cgns_base(int ibase) {
    int i, nb = 1;
    BASE *qbase = NULL;
    if (cgnsfn != 0) {
        /* Init the static values */
        cgnsbase = 0;
        nZones  = 0;
        CellDim = 0;
        PhyDim  = 0;
        Zones   = NULL;
        nDescs  = 0;
        Descs   = NULL;
        baseclass = 0;
        for (i = 0; i < 5; i++) baseunits[i] = 0;
        str_blank(basename);
        /* Read the CGNS Base data from File */
        qbase = new_base(nb);
        qbase->id = ibase;
        read_base(ibase);
        qbase->baseclass = baseclass;
        qbase->celldim   = CellDim;
        qbase->phydim    = PhyDim;
        qbase->nzones    = nZones;
        qbase->zones     = Zones;
        qbase->ndesc     = nDescs;
        qbase->desc      = Descs;
        for (i = 0; i < 5; i++)
            qbase->baseunits[i] = baseunits[i];
        strcpy(qbase->name, basename);
    }
    return qbase;
}

/*---------- set_cgns_base ---------------------------------------------
 * Set Base pointed data to Static variables
 *---------------------------------------------------------------------*/
int set_cgns_base(BASE *base) {
    int i;
    if (cgnsfn != 0) {
        /* Init the static values */
        cgnsbase = 0;
        nZones  = 0;
        CellDim = 0;
        PhyDim  = 0;
        Zones   = NULL;
        nDescs  = 0;
        Descs   = NULL;
        baseclass = 0;
        for (i = 0; i < 5; i++) baseunits[i] = 0;
        str_blank(basename);
        /* Set the values to static variable pointed by base */
        if (base != NULL) {
            cgnsbase  = base->id;
            baseclass = base->baseclass;
            CellDim   = base->celldim;
            PhyDim    = base->phydim;
            nZones    = base->nzones;
            Zones     = base->zones;
            nDescs    = base->ndesc;
            Descs     = base->desc;
            for (i = 0; i < 5; i++)
                baseunits[i] = base->baseunits[i];
            strcpy(basename, base->name);
            /* Write the data into CGNS File */
            write_base();
        } else
            return CG_ERROR;
    }
    return CG_OK;
}

/*---------- read_base -------------------------------------------------
 * read the CGNS Base from file
 *---------------------------------------------------------------------*/
void read_base(int ibase){
    int nz, ns, nd, celldim;
    char buff[33];
    DataClass_t tmpdata = DataClassNull;

    cgnsbase = ibase;
    if (cg_goto(cgnsfn, cgnsbase, "end"))
        CGNSIO_FATAL("read_base", NULL);

    read_units(baseunits);
    
    if (cg_dataclass_read(&tmpdata))
        baseclass = 0;
    else
        baseclass = tmpdata;
    
    if (cg_base_read(cgnsfn, cgnsbase, buff, &celldim, &nd))
        CGNSIO_FATAL("read_base", NULL);
    CellDim = celldim;
    PhyDim  = nd;
    strcpy(basename, buff);

    /* Read Zones and Zone Data */
    read_zones();
    for (nz = 1; nz <= nZones; nz++) {
        read_zone_grid(nz);
        read_zone_element(nz);
        read_zone_interface(nz);
        read_zone_connect(nz);
        read_zone_boco(nz);
        read_zone_solution(nz);
        for (ns = 1; ns <= Zones[nz-1].nsols; ns++)
            read_solution_field(nz, ns, 0);
    }

    /* Read the Descriptors of the Base */
    if (cg_goto(cgnsfn, cgnsbase, "end"))
        CGNSIO_FATAL("read_base", NULL);
    if (cg_ndescriptors(&nd))
        CGNSIO_FATAL("read_base", NULL);
    nDescs = nd;
    if (nd) {
        nDescs = nd;
        Descs = new_desc(nd);
        for (nd = 0; nd < nDescs; nd++) {
            if (cg_descriptor_read(nd+1, Descs[nd].name, &Descs[nd].desc))
                CGNSIO_FATAL("read_base", NULL);
        }
   }
}

/*---------- read_zones -----------------------------------------------
 * read zone information from CGNS file
 *---------------------------------------------------------------------*/
int read_zones(void){
    int n, nz, nd, celldim;
    int sizes[9];
    ZoneType_t zonetype;
    char buff[33];
    
    if (cg_goto(cgnsfn, cgnsbase, "end"))
        CGNSIO_FATAL("read_zones", NULL);

    celldim = CellDim;
    
    if (cg_nzones(cgnsfn, cgnsbase, &nZones))
        CGNSIO_FATAL("read_zones", NULL);
    Zones = new_zone(nZones);
    
    /* read the zone information */
    for (nz = 0; nz < nZones; nz++) {
        if (cg_zone_read(cgnsfn, cgnsbase, nz+1, buff, sizes) ||
                cg_zone_type(cgnsfn, cgnsbase, nz+1, &zonetype))
            CGNSIO_FATAL("read_zones", NULL);
        if (zonetype != Structured && zonetype != Unstructured)
            CGNSIO_FATAL("read_zones", "invalid zone type");
        Zones[nz].id = nz + 1;
        strcpy(Zones[nz].name, buff);
        Zones[nz].type = zonetype;
        Zones[nz].idim = zonetype == Structured ? celldim : 1;
        for (n = 0; n < 3; n++)
            Zones[nz].dim[n] = sizes[n];
        
        /* get units */
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz+1, "end"))
            CGNSIO_FATAL("read_zones", NULL);
        if (!read_units(Zones[nz].units)) {
            for (n = 0; n < 5; n++)
                Zones[nz].units[n] = baseunits[n];
        }
        if (cg_dataclass_read((DataClass_t *)&Zones[nz].dataclass))
            Zones[nz].dataclass = baseclass;
        
        /* get descriptors */
        if (cg_ndescriptors(&nd))
            CGNSIO_FATAL("cg_ndescriptors", NULL);
        if (nd) {
            Zones[nz].ndesc = nd;
            Zones[nz].desc = new_desc(nd);
            for (n = 0; n < nd; n++) {
                if (cg_descriptor_read(n+1, Zones[nz].desc[n].name,
                        &Zones[nz].desc[n].desc))
                    CGNSIO_FATAL("cg_descriptor_read", NULL);
            }
        }
        
        /* get zone counts */
        if (cg_nsections(cgnsfn, cgnsbase, nz+1, &Zones[nz].nesets) ||
                cg_n1to1(cgnsfn, cgnsbase, nz+1, &Zones[nz].nints) ||
                cg_nconns(cgnsfn, cgnsbase, nz+1, &Zones[nz].nconns) ||
                cg_nbocos(cgnsfn, cgnsbase, nz+1, &Zones[nz].nbocos) ||
                cg_nsols(cgnsfn, cgnsbase, nz+1, &Zones[nz].nsols))
            CGNSIO_FATAL("read_zones", NULL);
    }
    return nZones;
}

/*---------- read_zone_grid -------------------------------------------
 * read zone grid coordinates
 *---------------------------------------------------------------------*/
int read_zone_grid(int nz) {
    int n, nverts, nc, ncoords, rng[2][3];
    DataType_t datatype;
    char buff[33];
    double *xyz;
    ZONE *z = &Zones[nz-1];
    
    if (z->type == Structured) {
        nverts = z->dim[0] * z->dim[1] * z->dim[2];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = z->dim[n];
        }
    }
    else {
        nverts = z->dim[0];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = nverts;
        }
    }
    
    xyz = (double *) malloc(nverts * sizeof(double));
    if (NULL == xyz)
        CGNSIO_FATAL("read_zone_grid", "malloc failed for coordinate working array");
    z->vertflags = 0;
    z->nverts = nverts;
    z->verts = new_vertex(nverts);
    
    
    /* read the nodes */
    
    if (cg_ncoords(cgnsfn, cgnsbase, nz, &ncoords))
        CGNSIO_FATAL("read_zone_grid", NULL);
    for (nc = 1; nc <= ncoords; nc++) {
        if (cg_coord_info(cgnsfn, cgnsbase, nz, nc, &datatype, buff) ||
                cg_coord_read(cgnsfn, cgnsbase, nz, buff, RealDouble,
                rng[0], rng[1], xyz))
            CGNSIO_FATAL("read_zone_grid", NULL);
        if (!strcmp(buff, "CoordinateX")) {
            z->vertflags |= 1;
            z->datatype = datatype;
            for (n = 0; n < nverts; n++)
                z->verts[n].x = xyz[n];
            continue;
        }
        if (!strcmp(buff, "CoordinateY")) {
            z->vertflags |= 2;
            for (n = 0; n < nverts; n++)
                z->verts[n].y = xyz[n];
            continue;
        }
        if (!strcmp(buff, "CoordinateZ")) {
            z->vertflags |= 4;
            for (n = 0; n < nverts; n++)
                z->verts[n].z = xyz[n];
        }
    }
    free(xyz);
    return nverts;
}

/*---------- read_zone_element ----------------------------------------
 * read zone element sets and elements
 *---------------------------------------------------------------------*/
int read_zone_element(int nz) {
    int ns, size, iparent;
    ZONE *z = &Zones[nz-1];
    ELEMSET *eset;
    
    if (cg_nsections(cgnsfn, cgnsbase, nz, &z->nesets))
        CGNSIO_FATAL("read_zone_element", NULL);
    if (z->nesets) {
        z->esets = eset = new_elemset(z->nesets);
        for (ns = 1; ns <= z->nesets; ns++, eset++) {
            if (cg_section_read(cgnsfn, cgnsbase, nz, ns,
                    eset->name, (ElementType_t *)&eset->type,
                    &eset->start, &eset->end, &eset->nbndry, &iparent) ||
                    cg_ElementDataSize(cgnsfn, cgnsbase, nz, ns, &size))
                CGNSIO_FATAL("read_zone_element", NULL);
            eset->csize = size;
            eset->conn  = (int *) malloc(size * sizeof(int));
            if (NULL == eset->conn)
                CGNSIO_FATAL("read_zone_element",
                        "malloc failed for element connectivity");
            if (iparent) {
                size = 4 * (eset->end - eset->start + 1);
                eset->psize  = size;
                eset->parent = (int *) malloc(size * sizeof(int));
                if (NULL == eset->conn)
                    CGNSIO_FATAL("read_zone_element", "malloc failed for parent data");
            }
            if (cg_elements_read(cgnsfn, cgnsbase, nz, ns,
                    eset->conn, eset->parent))
                CGNSIO_FATAL("read_zone_element", NULL);
        }
    }
    return z->nesets;
}

/*---------- read_zone_interface --------------------------------------
 * read zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/
int read_zone_interface(int nz) {
    int i, j, n, ni, range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz-1];
    INTERFACE *ints;
    
    if (cg_n1to1(cgnsfn, cgnsbase, nz, &z->nints))
        CGNSIO_FATAL("read_zone_interface", NULL);
    if (z->nints) {
        z->ints = ints = new_interface(z->nints);
        for (ni = 1; ni <= z->nints; ni++, ints++) {
            if (cg_1to1_read(cgnsfn, cgnsbase, nz, ni,
                    ints->name, ints->d_name, (int *)range,
                    (int *)d_range, (int *)ints->transform))
                CGNSIO_FATAL("read_zone_interface", NULL);
            for (j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++) {
                    ints->range[i][j] = range[j][i];
                    ints->d_range[i][j] = d_range[j][i];
                }
            }
            for (n = 0; n < nZones; n++) {
                if (!strcmp(Zones[n].name, ints->d_name)) {
                    ints->d_zone = n + 1;
                    break;
                }
            }
        }
    }
    return z->nints;
}

/*---------- read_zone_connect ----------------------------------------
 * read zone connectivities
 *---------------------------------------------------------------------*/
int read_zone_connect(int nz) {
    int n, nc, npnts;
    GridLocation_t location;
    GridConnectivityType_t type;
    PointSetType_t ptype, d_ptype;
    ZoneType_t d_ztype;
    DataType_t d_datatype;
    ZONE *z = &Zones[nz-1];
    CONNECT *conns;
    
    if (cg_nconns(cgnsfn, cgnsbase, nz, &z->nconns))
        CGNSIO_FATAL("read_zone_connect", NULL);
    if (z->nconns) {
        z->conns = conns = new_connect(z->nconns);
        for (nc = 1; nc <= z->nconns; nc++, conns++) {
            if (cg_conn_info(cgnsfn, cgnsbase, nz, nc,
                    conns->name, &location, &type, &ptype,
                    &conns->npnts, conns->d_name, &d_ztype,
                    &d_ptype, &d_datatype, &conns->d_npnts))
                CGNSIO_FATAL("read_zone_connect", NULL);
            conns->location = location;
            conns->type = type;
            conns->ptype = ptype;
            conns->d_ztype = d_ztype;
            conns->d_ptype = d_ptype;
            npnts = conns->npnts * z->idim;
            conns->pnts = (int *) calloc(npnts, sizeof(int));
            npnts = conns->d_npnts * z->idim;
            conns->d_pnts = (int *) calloc(npnts, sizeof(int));
            if (NULL == conns->pnts || NULL == conns->d_pnts)
                CGNSIO_FATAL("read_zone_connect",
                        "malloc failed for connectivity point arrays");
            if (cg_conn_read(cgnsfn, cgnsbase, nz, nc,
                    conns->pnts, Integer, conns->d_pnts))
                CGNSIO_FATAL("read_zone_connect", NULL);
            for (n = 0; n < nZones; n++) {
                if (!strcmp(Zones[n].name, conns->d_name)) {
                    conns->d_zone = n + 1;
                    break;
                }
            }
        }
    }
    return z->nconns;
}

/*---------- read_zone_boco  ------------------------------------------
 * read zone boundary conditions
 *---------------------------------------------------------------------*/
int read_zone_boco(int nz) {
    int nb, npnts, ndatasets;
    BCType_t bctype;
    PointSetType_t ptype;
    DataType_t datatype;
    ZONE *z = &Zones[nz-1];
    BOCO *bocos;
    
    if (cg_nbocos(cgnsfn, cgnsbase, nz, &z->nbocos))
        CGNSIO_FATAL("read_zone_boco", NULL);
    if (z->nbocos) {
        z->bocos = bocos = new_boco(z->nbocos);
        for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
            if (cg_boco_info(cgnsfn, cgnsbase, nz, nb, bocos->name,
                    &bctype, &ptype, &bocos->npnts, bocos->n_index,
                    &bocos->n_cnt, &datatype, &ndatasets))
                CGNSIO_FATAL("read_zone_boco", NULL);
            bocos->type = bctype;
            bocos->ptype = ptype;
            bocos->n_type = datatype;
            npnts = bocos->npnts * z->idim;
            bocos->pnts = (int *) calloc(npnts, sizeof(int));
            if (NULL == bocos->pnts)
                CGNSIO_FATAL("read_zone_boco",
                        "calloc failed for boco point arrays");
            if (bocos->n_cnt) {
                bocos->n_list = (double *) calloc(bocos->n_cnt, sizeof(double));
                if (NULL == bocos->n_list)
                    CGNSIO_FATAL("read_zone_boco",
                            "calloc failed for boco normal list");
            }
            if (cg_boco_read(cgnsfn, cgnsbase, nz, nb,
                    bocos->pnts, bocos->n_list))
                CGNSIO_FATAL("read_zone_boco", NULL);
        }
    }
    return z->nbocos;
}

/*---------- read_zone_solution ---------------------------------------
 * read zone solution
 *---------------------------------------------------------------------*/
int read_zone_solution(int nz) {
    int i, j, ns, nd;
    DataType_t datatype;
    GridLocation_t location;
    ZONE *z = &Zones[nz-1];
    SOLUTION *sols;
    
    if (cg_nsols(cgnsfn, cgnsbase, nz, &z->nsols))
        CGNSIO_FATAL("read_zone_solution", NULL);
    if (z->nsols) {
        z->sols = sols = new_solution(z->nsols);
        for (ns = 1; ns <= z->nsols; ns++, sols++) {
            if (cg_sol_info(cgnsfn, cgnsbase, nz, ns,
                    sols->name, &location))
                CGNSIO_FATAL("read_zone_solution", NULL);
            sols->location = location;
            if (z->type == Structured) {
                if (sols->location == Vertex) {
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 2; j++)
                            sols->rind[i][j] = 0;
                    sols->size = z->dim[0] * z->dim[1] * z->dim[2];
                }
                else if (sols->location == CellCenter) {
                    if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                            "FlowSolution_t", ns, "end"))
                        CGNSIO_FATAL("read_zone_solution", NULL);
                    if (cg_rind_read((int *)sols->rind)) {
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 2; j++)
                                sols->rind[i][j] = 0;
                    }
                    sols->size = 1;
                    for (i = 0; i < 3; i++) {
                        sols->size *= (z->dim[i] - 1 +
                                sols->rind[i][0] + sols->rind[i][1]);
                    }
                }
                else
                    CGNSIO_FATAL("read_zone_solution",
                            "solution location not Vertex or CellCenter");
            }
            else {
                sols->size = sols->location == Vertex ? z->dim[0] : z->dim[1];
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 2; j++)
                        sols->rind[i][j] = 0;
            }
            if (cg_nfields(cgnsfn, cgnsbase, nz, ns, &sols->nflds))
                CGNSIO_FATAL("read_zone_solution", NULL);
            if (sols->nflds) {
                sols->flds = new_field(sols->nflds, 0);
                for (i = 0; i < sols->nflds; i++) {
                    if (cg_field_info(cgnsfn, cgnsbase, nz, ns, i+1,
                            &datatype, sols->flds[i].name))
                        CGNSIO_FATAL("read_zone_solution", NULL);
                }
            }
            
            /* get units */
            
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                    "FlowSolution_t", ns, "end"))
                CGNSIO_FATAL("read_zone_solution", NULL);
            if (!read_units(sols->units)) {
                for (i = 0; i < 5; i++)
                    sols->units[i] = z->units[i];
            }
            if (cg_dataclass_read((DataClass_t *)&sols->dataclass))
                sols->dataclass = z->dataclass;
            
            /* get descriptors */
            
            if (cg_ndescriptors(&nd))
                CGNSIO_FATAL("cg_ndescriptors", NULL);
            if (nd) {
                sols->ndesc = nd;
                sols->desc = new_desc(nd);
                for (i = 0; i < nd; i++) {
                    if (cg_descriptor_read(i+1, sols->desc[i].name,
                            &sols->desc[i].desc))
                        CGNSIO_FATAL("cg_descriptor_read", NULL);
                }
            }
        }
    }
    return z->nsols;
}

/*---------- read_solution_field --------------------------------------
 * read solution field data for a zone
 *---------------------------------------------------------------------*/
int read_solution_field(int nz, int ns, int nf) {
    int n, is, ie, min[3], max[3];
    DataType_t datatype;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];
    FIELD *f;
    
    if (z->type == Structured) {
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = z->dim[n];
        }
        if (s->location == CellCenter) {
            for (n = 0; n < 3; n++)
                max[n] += s->rind[n][0] + s->rind[n][1] - 1;
        }
    }
    else {
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = s->size;
        }
    }
    if (nf) {
        is = ie = nf;
    }
    else {
        is = 1;
        ie = s->nflds;
    }
    f = &s->flds[is-1];
    for (nf = is; nf <= ie; nf++, f++) {
        if (cg_field_info(cgnsfn, cgnsbase, nz, ns, nf,
                &datatype, f->name))
            CGNSIO_FATAL("read_solution_field", NULL);
        f->id = nf;
        f->datatype = datatype;
        f->data = (double *) malloc(s->size * sizeof(double));
        if (NULL == f->data)
            CGNSIO_FATAL("read_solution_field",
                    "malloc failed for solution field data");
        if (cg_field_read(cgnsfn, cgnsbase, nz, ns, f->name,
                RealDouble, min, max, f->data))
            CGNSIO_FATAL("read_solution_field", NULL);
        
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                "FlowSolution_t", ns, "DataArray_t", nf, "end"))
            CGNSIO_FATAL("read_solution_field", NULL);
        if (!read_units(f->units)) {
            for (n = 0; n < 5; n++)
                f->units[n] = s->units[n];
        }
        
        /* read data class, conversion and exponents */
        
        if (cg_dataclass_read((DataClass_t *)&f->dataclass))
            f->dataclass = s->dataclass;
        
        if (cg_conversion_info(&datatype))
            f->dataconv[0] = 1.0;
        else {
            f->convtype = datatype;
            if (datatype == RealSingle) {
                float conv[2];
                if (cg_conversion_read(conv))
                    CGNSIO_FATAL("read_solution_field", NULL);
                for (n = 0; n < 2; n++)
                    f->dataconv[n] = conv[n];
            }
            else if (datatype == RealDouble) {
                if (cg_conversion_read(f->dataconv))
                    CGNSIO_FATAL("read_solution_field", NULL);
            }
            else
                CGNSIO_FATAL("cg_conversion_info", "invalid data type");
        }
        
        if (!cg_exponents_info(&datatype)) {
            f->exptype = datatype;
            if (datatype == RealSingle) {
                float exp[5];
                if (cg_exponents_read(exp))
                    CGNSIO_FATAL("read_solution_field", NULL);
                for (n = 0; n < 5; n++)
                    f->exponent[n] = exp[n];
            }
            else if (datatype == RealDouble) {
                if (cg_exponents_read(f->exponent))
                    CGNSIO_FATAL("read_solution_field", NULL);
            }
            else
                CGNSIO_FATAL("cg_exponents_info", "invalid data type");
        }
    }
    return s->size;
}

/*---------- read_units -----------------------------------------------
 * read unit specifications
 *---------------------------------------------------------------------*/
int read_units(int units[5]) {
    int n;
    MassUnits_t mass;
    LengthUnits_t length;
    TimeUnits_t time;
    TemperatureUnits_t temp;
    AngleUnits_t angle;
    
    if (cg_units_read(&mass, &length, &time, &temp, &angle)) {
        for (n = 0; n < 5; n++)
            units[n] = 0;
        return 0;
    }
    units[0] = mass;
    units[1] = length;
    units[2] = time;
    units[3] = temp;
    units[4] = angle;
    return 1;
}

/*---------- write_base -----------------------------------------------
 * write the Base Data structure to CGNS file
 *---------------------------------------------------------------------*/
void write_base(void) {
    int n, nz, ns;

    /* Create the new base in the CGNS File */
    if (cg_base_write(cgnsfn, basename, CellDim, PhyDim, &cgnsbase))
        CGNSIO_FATAL("write_base", NULL);

    for (n = 0; n < 5; n++) {
        if (baseunits[n]) {
            if (cg_goto(cgnsfn, cgnsbase, "end") ||
                    cg_units_write((MassUnits_t)baseunits[0],
                    (LengthUnits_t)baseunits[1],
                    (TimeUnits_t)baseunits[2],
                    (TemperatureUnits_t)baseunits[3],
                    (AngleUnits_t)baseunits[4]))
                CGNSIO_FATAL("write_base", NULL);
            break;
        }
    }
    if (baseclass) {
        if (cg_goto(cgnsfn, cgnsbase, "end") ||
                cg_dataclass_write((DataClass_t)baseclass))
            CGNSIO_FATAL("write_base", NULL);
    }

    /* Write Zone and Zone Data to File */
    write_zones();
    for (nz = 1; nz <= nZones; nz++) {
        write_zone_grid(nz);
        write_zone_element(nz);
        write_zone_interface(nz);
        write_zone_connect(nz);
        write_zone_boco(nz);
        for (ns = 1; ns <= Zones[nz-1].nsols; ns++) {
            write_zone_solution(nz, ns);
            write_solution_field(nz, ns, 0);
        }
    }

    /* Write the Base Discriptors */
    if (nDescs) {
        if (cg_goto(cgnsfn, cgnsbase, "end"))
                CGNSIO_FATAL("write_base", NULL);
        for (n = 0; n < nDescs; n++) {
            if (cg_descriptor_write(Descs[n].name, Descs[n].desc))
                CGNSIO_FATAL("write_base", NULL);
        }
    }
}

/*---------- write_zones ----------------------------------------------
 * write zone information to CGNS file
 *---------------------------------------------------------------------*/
void write_zones(void) {
    int n, nz, sizes[3][3];
    ZONE *z = Zones;
    
    /* write the zone information */
    for (nz = 0; nz < nZones; nz++, z++) {
        if (!z->id) continue;
        for (n = 0; n < 3; n++) {
            sizes[0][n] = z->dim[n];
            sizes[1][n] = z->dim[n] - 1;
            sizes[2][n] = 0;
        }
        if (cg_zone_write(cgnsfn, cgnsbase, z->name,
                (int *)sizes, (ZoneType_t)z->type, &z->id))
            CGNSIO_FATAL("write_zones", NULL);
        
        for (n = 0; n < 5; n++) {
            if (z->units[n] && z->units[n] != baseunits[n]) {
                if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                        cg_units_write((MassUnits_t)z->units[0],
                        (LengthUnits_t)z->units[1],
                        (TimeUnits_t)z->units[2],
                        (TemperatureUnits_t)z->units[3],
                        (AngleUnits_t)z->units[4]))
                    CGNSIO_FATAL("write_zones", NULL);
                break;
            }
        }
        
        if (z->dataclass && z->dataclass != baseclass) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                    cg_dataclass_write((DataClass_t)z->dataclass))
                CGNSIO_FATAL("write_zones", NULL);
        }
        
        if (z->ndesc) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end"))
                CGNSIO_FATAL("write_zones", NULL);
            for (n = 0; n < z->ndesc; n++) {
                if (cg_descriptor_write(z->desc[n].name, z->desc[n].desc))
                    CGNSIO_FATAL("cg_descriptor_write", NULL);
            }
        }
    }
}

/*---------- write_zone_grid ------------------------------------------
 * write zone grid coordinates
 *---------------------------------------------------------------------*/
void write_zone_grid(int nz) {
    int n, nc;
    ZONE *z = &Zones[nz-1];
    
    if (z->verts == NULL || (z->nverts < 1) || (z->vertflags & 7) == 0) return;
    
    if (z->datatype == RealSingle) {
        float *xyz = (float *) malloc(z->nverts * sizeof(float));
        if (NULL == xyz)
            CGNSIO_FATAL("write_zone_grid",
                    "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].x;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateX", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].y;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateY", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float)z->verts[n].z;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateZ", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        free(xyz);
    }
    else {
        double *xyz = (double *) malloc(z->nverts * sizeof(double));
        if (NULL == xyz)
            CGNSIO_FATAL("write_zone_grid",
                    "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].x;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateX", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].y;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateY", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].z;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateZ", xyz, &nc)) CGNSIO_FATAL("write_zone_grid", NULL);
        }
        free(xyz);
    }
}

/*---------- write_zone_element ---------------------------------------
 * write zone element sets and elements
 *---------------------------------------------------------------------*/
void write_zone_element(int nz) {
    int ns;
    ZONE *z = &Zones[nz-1];
    ELEMSET *eset = z->esets;
    
    if (eset == NULL || (z->nesets < 1) || z->type == Structured) return;
    
    for (ns = 1; ns <= z->nesets; ns++, eset++) {
        if (eset->id) {
            if (cg_section_write(cgnsfn, cgnsbase, z->id,
                    eset->name, (ElementType_t)eset->type,
                    eset->start, eset->end, eset->nbndry,
                    eset->conn, &eset->id))
                CGNSIO_FATAL("write_zone_element", NULL);
            if (eset->parent != NULL &&
                    cg_parent_data_write(cgnsfn, cgnsbase, z->id,
                    eset->id, eset->parent))
                CGNSIO_FATAL("write_zone_element", NULL);
        }
    }
}

/*---------- write_zone_interface -------------------------------------
 * write zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/
void write_zone_interface(int nz) {
    int i, j, ni, range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz-1];
    INTERFACE *ints = z->ints;

    if (z->nints < 1) return;
    if (ints == NULL) return;
    
    for (ni = 1; ni <= z->nints; ni++, ints++) {
        if (ints->id) {
            for (j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++) {
                    range[j][i] = ints->range[i][j];
                    d_range[j][i] = ints->d_range[i][j];
                }
            }
            if (cg_1to1_write(cgnsfn, cgnsbase, z->id,
                    ints->name, ints->d_name, (int *)range,
                    (int *)d_range, ints->transform, &ints->id))
                CGNSIO_FATAL("write_zone_interface", NULL);
        }
    }
}

/*---------- write_zone_connect ---------------------------------------
 * write zone connectivities
 *---------------------------------------------------------------------*/
void write_zone_connect(int nz) {
    int nc;
    ZONE *z = &Zones[nz-1];
    CONNECT *conns = z->conns;

    if (z->nconns < 1) return;
    if (conns == NULL) return;
    
    for (nc = 1; nc <= z->nconns; nc++, conns++) {
        if (conns->id &&
                cg_conn_write(cgnsfn, cgnsbase, z->id,
                conns->name, (GridLocation_t)conns->location,
                (GridConnectivityType_t)conns->type,
                (PointSetType_t)conns->ptype, conns->npnts, conns->pnts,
                conns->d_name, (ZoneType_t)conns->d_ztype,
                (PointSetType_t)conns->d_ptype, Integer, conns->d_npnts,
                conns->d_pnts, &conns->id))
            CGNSIO_FATAL("write_zone_connect", NULL);
    }
}

/*---------- write_zone_bocos -----------------------------------------
 * write zone boundary conditions
 *---------------------------------------------------------------------*/
void write_zone_boco(int nz) {
    int nb;
    ZONE *z = &Zones[nz-1];
    BOCO *bocos = z->bocos;

    if (z->nbocos < 1) return;
    if (bocos == NULL) return;
    
    for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
        if ((bocos->id) &&
                (cg_boco_write(cgnsfn, cgnsbase, z->id,
                bocos->name, (BCType_t)bocos->type,
                (PointSetType_t)bocos->ptype, bocos->npnts,
                bocos->pnts, &bocos->id) ||
                cg_boco_normal_write(cgnsfn, cgnsbase, z->id, bocos->id,
                bocos->n_index, bocos->n_cnt, (DataType_t)bocos->n_type,
                bocos->n_list)))
            CGNSIO_FATAL("write_zone_boco", NULL);
    }
}

/*---------- write_zone_solution --------------------------------------
 * write zone solution
 *---------------------------------------------------------------------*/
void write_zone_solution(int nz, int ns) {
    int n;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s;
    
    if (z->sols == NULL || ns < 1 || ns > z->nsols) return;
    s = &z->sols[ns-1];
    
    if (cg_sol_write(cgnsfn, cgnsbase, z->id, s->name,
            (GridLocation_t)s->location, &s->id))
        CGNSIO_FATAL("write_zone_solution", NULL);
    if (z->type == Structured && s->location == CellCenter) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end") ||
                cg_rind_write((int *)s->rind))
            CGNSIO_FATAL("write_zone_solution", NULL);
    }
    
    for (n = 0; n < 5; n++) {
        if (s->units[n] && s->units[n] != z->units[n]) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "end") ||
                    cg_units_write((MassUnits_t)s->units[0],
                    (LengthUnits_t)s->units[1],
                    (TimeUnits_t)s->units[2],
                    (TemperatureUnits_t)s->units[3],
                    (AngleUnits_t)s->units[4]))
                CGNSIO_FATAL("write_zone_solution", NULL);
            break;
        }
    }
    
    if (s->dataclass && s->dataclass != z->dataclass) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end") ||
                cg_dataclass_write((DataClass_t)s->dataclass))
            CGNSIO_FATAL("write_zone_solution", NULL);
    }
    
    if (s->ndesc) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end"))
            CGNSIO_FATAL("write_zone_solution", NULL);
        for (n = 0; n < s->ndesc; n++) {
            if (cg_descriptor_write(s->desc[n].name, s->desc[n].desc))
                CGNSIO_FATAL("cg_descriptor_write", NULL);
        }
    }
}

/*---------- write_solution_field -------------------------------------
 * write solution field data for a zone
 *---------------------------------------------------------------------*/
void write_solution_field(int nz, int ns, int nf) {
    int n, is, ie;
    float *data = NULL;
    ZONE *z = &Zones[nz-1];
    SOLUTION *s = &z->sols[ns-1];
    FIELD *f;
    
    if (nf) {
        is = ie = nf;
    }
    else {
        is = 1;
        ie = s->nflds;
    }
    f = &s->flds[is-1];
    
    for (nf = is; nf <= ie; nf++, f++) {
        if (f->data == NULL) continue;
        if (f->datatype == RealSingle) {
            if (data == NULL) {
                data = (float *) malloc(s->size * sizeof(float));
                if (NULL == data)
                    CGNSIO_FATAL("write_solution_field",
                            "malloc failed for working array");
            }
            for (n = 0; n < s->size; n++)
                data[n] = (float)f->data[n];
            if (cg_field_write(cgnsfn, cgnsbase, z->id, s->id,
                    RealSingle, f->name, data, &f->id))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        else {
            if (cg_field_write(cgnsfn, cgnsbase, z->id, s->id,
                    RealDouble, f->name, f->data, &f->id))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        for (n = 0; n < 5; n++) {
            if (f->units[n] && f->units[n] != s->units[n]) {
                if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                        "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                        cg_units_write((MassUnits_t)f->units[0],
                        (LengthUnits_t)f->units[1],
                        (TimeUnits_t)f->units[2],
                        (TemperatureUnits_t)f->units[3],
                        (AngleUnits_t)f->units[4]))
                    CGNSIO_FATAL("write_solution_field", NULL);
                break;
            }
        }
        if (f->dataclass && f->dataclass != s->dataclass) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_dataclass_write((DataClass_t)f->dataclass))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        if (f->convtype == RealSingle) {
            float conv[2];
            for (n = 0; n < 2; n++)
                conv[n] = (float)f->dataconv[n];
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_conversion_write(RealSingle, conv))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        else if (f->convtype == RealDouble) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_conversion_write(RealDouble, f->dataconv))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        else {}
        if (f->exptype == RealSingle) {
            float exp[5];
            for (n = 0; n < 5; n++)
                exp[n] = (float)f->dataconv[n];
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_exponents_write(RealSingle, exp))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        else if (f->exptype == RealDouble) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_exponents_write(RealDouble, f->exponent))
                CGNSIO_FATAL("write_solution_field", NULL);
        }
        else {}
    }
    if (data != NULL) free(data);
}
