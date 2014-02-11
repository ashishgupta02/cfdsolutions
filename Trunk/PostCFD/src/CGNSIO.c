/*[
 * Copyright 2006   Ashish Gupta
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
 * with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

/*
 * File		CGNSIO.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "CGNSIO.h"

#ifdef _WIN32
#define PATH_DELIM ';'
#include <direct.h>
#include <io.h>
#else
#define PATH_DELIM ':'
#include <unistd.h>
#endif

float Version = 0;
int nZones = 0;
ZONE *Zones;
int baseclass = 0;
int baseunits[5] = {0, 0, 0, 0, 0};
int cgnsbase = 1;
int cgnsfn = 0;
int element_node_counts[] = {
    0, 0, /* ElementTypeNull, ElementTypeUserDefined */
    1, 2, 3, /* NODE, BAR_2, BAR_3 */
    3, 6, /* TRI_3, TRI_6 */
    4, 8, 9, /* QUAD_4, QUAD_8, QUAD_9 */
    4, 10, /* TETRA_4, TETRA_10 */
    5, 14, /* PYRA_5, PYRA_14 */
    6, 15, 18, /* PENTA_6, PENTA_15, PENTA_18 */
    8, 20, 27, /* HEXA_8, HEXA_20, HEXA_27 */
    0, 0 /* MIXED, NGON_n */
};
/* To get temp file name */
static char cgnstemp[16] = "";

/*---------- file_exists ----------------------------------------------
 * check if a file exists
 *---------------------------------------------------------------------*/
int file_exists(const char *file) {
    struct stat st;

    if (access(file, 0) || stat(file, &st) ||
            S_IFREG != (st.st_mode & S_IFMT)) return 0;
    return 1;
}

/*---------- FATAL ----------------------------------------------------
 * exit with error message
 *---------------------------------------------------------------------*/
void FATAL(char *procname, char *errmsg) {
    char *msg = errmsg;

    if (NULL == msg) {
        msg = (char *) cg_get_error();
        if (NULL == msg || !*msg)
            msg = "unknown error";
    }
    fflush(stdout);
    if (NULL == procname || !*procname)
        fprintf(stderr, "%s\n", msg);
    else
        fprintf(stderr, "%s:%s\n", procname, msg);
    if (cgnsfn) cg_close(cgnsfn);
    exit(1);
}

/*---------- temporary_file -------------------------------------------
 * create a temporary file
 *---------------------------------------------------------------------*/

char *temporary_file(void) {
    char temp[16];

    strcpy(temp, "cgnsXXXXXX");

#ifdef _WIN32
    if (mktemp(temp) == NULL)
        FATAL("temporary_file", "failed to create temporary filename");

#else
    if (mkstemp(temp) == -1)
        FATAL("temporary_file", "failed to create temporary filename");

#endif
    /*
            if (mktemp (temp) == NULL)
                    FATAL ("temporary_file", "failed to create temporary filename");
     */
    strcpy(cgnstemp, temp);

    return cgnstemp;
}

/*---------- copy_file ------------------------------------------------
 * make a copy of a file
 *---------------------------------------------------------------------*/
void copy_file(char *oldfile, char *newfile) {
    int c;
    FILE *oldfp, *newfp;

    if (NULL == (oldfp = fopen(oldfile, "rb")))
        FATAL("copy_file", "error opening input file for reading");
    if (NULL == (newfp = fopen(newfile, "w+b"))) {
        fclose(oldfp);
        FATAL("copy_file", "error opening output file for writing");
    }
    while (EOF != (c = getc(oldfp)))
        putc(c, newfp);
    fclose(oldfp);
    fclose(newfp);
}

/*---------- new_zone -------------------------------------------------
 * create new zone(s)
 *---------------------------------------------------------------------*/
ZONE *new_zone(int count) {
    int n;
    ZONE *z;

    z = (ZONE *) calloc(count, sizeof (ZONE));
    if (NULL == z)
        FATAL("new_zone", "calloc failed for new zones");
    for (n = 0; n < count; n++) {
        z[n].id = n + 1;
        sprintf(z[n].name, "Zone%d", n + 1);
        z[n].type = Structured;
        z[n].vertflags = 7;
        z[n].datatype = RealDouble;
    }
    return z;
}

/*---------- new_vertex -----------------------------------------------
 * create coordinate array for a zone
 *---------------------------------------------------------------------*/
VERTEX *new_vertex(int nverts) {
    int n;
    VERTEX *verts;

    verts = (VERTEX *) calloc(nverts, sizeof (VERTEX));
    if (NULL == verts)
        FATAL("new_vertex", "calloc failed for new vertex array");
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
    ELEMSET *esets;

    esets = (ELEMSET *) calloc(nesets, sizeof (ELEMSET));
    if (NULL == esets)
        FATAL("new_elemset", "calloc failed for new element set array");
    for (n = 0; n < nesets; n++) {
        esets[n].id = n + 1;
        sprintf(esets[n].name, "ElemSet%d", n + 1);
    }
    return esets;
}

/*---------- new_interface --------------------------------------------
 * create grid 1to1 interface array for a zone
 *---------------------------------------------------------------------*/
INTERFACE *new_interface(int nints) {
    int n;
    INTERFACE *ints;

    ints = (INTERFACE *) calloc(nints, sizeof (INTERFACE));
    if (NULL == ints)
        FATAL("new_interface", "calloc failed for new interface array");
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
    CONNECT *conns;

    conns = (CONNECT *) calloc(nconns, sizeof (CONNECT));
    if (NULL == conns)
        FATAL("new_connect", "calloc failed for new connectivity array");
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
    BOCO *bocos;

    bocos = (BOCO *) calloc(nbocos, sizeof (BOCO));
    if (NULL == bocos)
        FATAL("new_boco", "calloc failed for new boco array");
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
    SOLUTION *sols;

    sols = (SOLUTION *) calloc(nsols, sizeof (SOLUTION));
    if (NULL == sols)
        FATAL("new_solution", "calloc failed for new solution array");
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
    FIELD *flds;

    flds = (FIELD *) calloc(nflds, sizeof (FIELD));
    if (NULL == flds)
        FATAL("new_field", "calloc failed for new field array");
    for (n = 0; n < nflds; n++) {
        flds[n].id = n + 1;
        sprintf(flds[n].name, "Field%d", n + 1);
        flds[n].datatype = RealDouble;
        if (size > 0) {
            flds[n].data = (double *) calloc(size, sizeof (double));
            if (NULL == flds[n].data)
                FATAL("new_field", "calloc failed for field data array");
        }
    }
    return flds;
}

/*---------- new_desc -------------------------------------------------
 * create descriptor array
 *---------------------------------------------------------------------*/
DESC *new_desc(int ndesc) {
    int n;
    DESC *desc;

    desc = (DESC *) calloc(ndesc, sizeof (DESC));
    if (NULL == desc)
        FATAL("new_desc", "calloc failed for new descriptor array");
    for (n = 0; n < ndesc; n++) {
        desc[n].id = n + 1;
        sprintf(desc[n].name, "Descriptor%d", n + 1);
    }
    return desc;
}

/*---------- vertex_index ---------------------------------------------
 * get index in vertex array for structured grid
 *---------------------------------------------------------------------*/
int vertex_index(ZONE *z, int i, int j, int k) {
    return i - 1 + z->dim[0] * ((j - 1) + z->dim[1] * (k - 1));
}

/*---------- cell_index -----------------------------------------------
 * read index in cell array for structured grid
 *---------------------------------------------------------------------*/
int cell_index(ZONE *z, int i, int j, int k) {
    return i - 1 + (z->dim[0] - 1) * ((j - 1) + (z->dim[1] - 1) * (k - 1));
}

/*---------- solution_index -------------------------------------------
 * read index in solution array for structured grid
 *---------------------------------------------------------------------*/
int solution_index(ZONE *z, SOLUTION *s, int i, int j, int k) {
    int ni, nj;

    ni = z->dim[0] - 1 + s->rind[0][0] + s->rind[0][1];
    nj = z->dim[1] - 1 + s->rind[1][0] + s->rind[1][1];
    i += s->rind[0][0];
    j += s->rind[1][0];
    k += s->rind[2][0];
    return i - 1 + ni * ((j - 1) + nj * (k - 1));
}

/*---------- open_cgns ------------------------------------------------
 * open a CGNS file
 *-------------------------------------------------------------------*/
int open_cgns(const char *cgnsfile, int mode) {
    int nbases = 0;

    if (cgnsfn) {
        cg_close(cgnsfn);
        cgnsfn = 0;
    }
    switch (mode) {
        case 0:
            if (cg_open(cgnsfile, MODE_READ, &cgnsfn) ||
                    cg_nbases(cgnsfn, &nbases))
                FATAL("open_cgns", NULL);
            return nbases;
            break;
        case 1:
            if (cg_open(cgnsfile, MODE_MODIFY, &cgnsfn) ||
                    cg_nbases(cgnsfn, &nbases))
                FATAL("open_cgns", NULL);
            return nbases;
            break;
        case 2:
            if (cg_open(cgnsfile, MODE_WRITE, &cgnsfn))
                FATAL("open_cgns", NULL);
            return 0;
            break;
    }

    return 0;
}

/*---------- find_base ------------------------------------------------
 * find base id from base name
 *---------------------------------------------------------------------*/
int find_base(char *basename) {
    int nbases, nb, idum;
    char buff[33];

    if (cg_nbases(cgnsfn, &nbases))
        FATAL("find_base", NULL);
    for (nb = 1; nb <= nbases; nb++) {
        if (cg_base_read(cgnsfn, nb, buff, &idum, &idum))
            FATAL("find_base", NULL);
        if (!strcmp(buff, basename))
            return nb;
    }
    return 0;
}

/*---------- read_cgns ------------------------------------------------
 * read the CGNS file
 *---------------------------------------------------------------------*/
void read_cgns(void) {
    int nz, ns;

    read_zones();
    for (nz = 1; nz <= nZones; nz++) {
        read_zone_grid(nz);
        read_zone_element(nz);
        read_zone_interface(nz);
        read_zone_connect(nz);
        read_zone_boco(nz);
        read_zone_solution(nz);
        for (ns = 1; ns <= Zones[nz - 1].nsols; ns++)
            read_solution_field(nz, ns, 0);
    }
}

/*---------- read_zones -----------------------------------------------
 * read zone information from CGNS file
 *---------------------------------------------------------------------*/
int read_zones(void) {
    int n, nz, nd, celldim, phydim, sizes[9];
    ZoneType_t zonetype = ZoneTypeNull;
    char buff[33];

    if (cg_goto(cgnsfn, cgnsbase, "end"))
        FATAL("read_zones", NULL);

    read_units(baseunits);

    if (cg_dataclass_read((DataClass_t *) & baseclass))
        baseclass = 0;

    if (cg_base_read(cgnsfn, cgnsbase, buff, &celldim, &phydim) ||
            cg_nzones(cgnsfn, cgnsbase, &nZones))
        FATAL("read_zones", NULL);

    Zones = new_zone(nZones);

    /* read the zone information */
    for (nz = 0; nz < nZones; nz++) {

        if (cg_zone_read(cgnsfn, cgnsbase, nz + 1, buff, sizes) ||
                cg_zone_type(cgnsfn, cgnsbase, nz + 1, &zonetype))
            FATAL("read_zones", NULL);

        if (zonetype != Structured && zonetype != Unstructured)
            FATAL("read_zones", "invalid zone type");

        Zones[nz].id = nz + 1;
        strcpy(Zones[nz].name, buff);
        Zones[nz].type = zonetype;
        Zones[nz].idim = zonetype == Structured ? celldim : 1;

        switch (celldim) {
            case 1:
                /* Unstructured Grid */
                for (n = 0; n < 3; n++)
                    Zones[nz].dim[n] = sizes[n];
                break;
            case 2:
                /* Surface Grid */
                for (n = 0; n < 2; n++)
                    Zones[nz].dim[n] = sizes[n];
                Zones[nz].dim[2] = 1;
                break;
            case 3:
                /* Volume Grid */
                for (n = 0; n < 3; n++)
                    Zones[nz].dim[n] = sizes[n];
                break;
        }
        /*		
        for (n = 0; n < 3; n++)
                Zones[nz].dim[n] = sizes[n];
         */
        /* get units */
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz + 1, "end"))
            FATAL("read_zones", NULL);

        if (!read_units(Zones[nz].units)) {
            for (n = 0; n < 5; n++)
                Zones[nz].units[n] = baseunits[n];
        }

        if (cg_dataclass_read((DataClass_t *) & Zones[nz].dataclass))
            Zones[nz].dataclass = baseclass;

        /* get descriptors */
        if (cg_ndescriptors(&nd))
            FATAL("cg_ndescriptors", NULL);
        if (nd) {
            Zones[nz].ndesc = nd;
            Zones[nz].desc = new_desc(nd);
            for (n = 0; n < nd; n++) {
                if (cg_descriptor_read(n + 1, Zones[nz].desc[n].name,
                        &Zones[nz].desc[n].desc))
                    FATAL("cg_descriptor_read", NULL);
            }
        }
        /* get zone counts */
        if (cg_nsections(cgnsfn, cgnsbase, nz + 1, &Zones[nz].nesets) ||
                cg_n1to1(cgnsfn, cgnsbase, nz + 1, &Zones[nz].nints) ||
                cg_nconns(cgnsfn, cgnsbase, nz + 1, &Zones[nz].nconns) ||
                cg_nbocos(cgnsfn, cgnsbase, nz + 1, &Zones[nz].nbocos) ||
                cg_nsols(cgnsfn, cgnsbase, nz + 1, &Zones[nz].nsols))
            FATAL("read_zones", NULL);
    }
    return nZones;
}

/*---------- read_zone_data -------------------------------------------
 * read all zone data
 *---------------------------------------------------------------------*/
void read_zone_data(int nz) {
    int ns, nsols;

    read_zone_grid(nz);
    read_zone_element(nz);
    read_zone_interface(nz);
    read_zone_connect(nz);
    read_zone_boco(nz);
    nsols = read_zone_solution(nz);
    for (ns = 1; ns <= nsols; ns++)
        read_solution_field(nz, ns, 0);
}

/*---------- read_zone_grid -------------------------------------------
 * read zone grid coordinates
 *---------------------------------------------------------------------*/
int read_zone_grid(int nz) {
    int n, nverts, nc, ncoords, rng[2][3];
    DataType_t datatype;
    char buff[33];
    double *xyz;
    ZONE *z = &Zones[nz - 1];

    if (z->type == Structured) {
        /* Structured Surface and Volume Grid */
        nverts = z->dim[0] * z->dim[1] * z->dim[2];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = z->dim[n];
        }
    } else {
        /* Unstructured Grid */
        nverts = z->dim[0];
        for (n = 0; n < 3; n++) {
            rng[0][n] = 1;
            rng[1][n] = nverts;
        }
    }

    xyz = (double *) malloc(nverts * sizeof (double));
    if (NULL == xyz)
        FATAL("read_zone_grid", "malloc failed for coordinate working array");

    z->vertflags = 0;
    z->nverts = nverts;
    z->verts = new_vertex(nverts);

    /* read the nodes */
    if (cg_ncoords(cgnsfn, cgnsbase, nz, &ncoords))
        FATAL("read_zone_grid", NULL);
    for (nc = 1; nc <= ncoords; nc++) {
        if (cg_coord_info(cgnsfn, cgnsbase, nz, nc, &datatype, buff) ||
                cg_coord_read(cgnsfn, cgnsbase, nz, buff, RealDouble,
                rng[0], rng[1], xyz))
            FATAL("read_zone_grid", NULL);
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
    ZONE *z = &Zones[nz - 1];
    ELEMSET *eset;

    if (cg_nsections(cgnsfn, cgnsbase, nz, &z->nesets))
        FATAL("read_zone_element", NULL);
    if (z->nesets) {
        z->esets = eset = new_elemset(z->nesets);
        for (ns = 1; ns <= z->nesets; ns++, eset++) {
            if (cg_section_read(cgnsfn, cgnsbase, nz, ns,
                    eset->name, (ElementType_t *) & eset->type,
                    &eset->start, &eset->end, &eset->nbndry, &iparent) ||
                    cg_ElementDataSize(cgnsfn, cgnsbase, nz, ns, &size))
                FATAL("read_zone_element", NULL);
            eset->conn = (int *) malloc(size * sizeof (int));
            if (NULL == eset->conn)
                FATAL("read_zone_element",
                    "malloc failed for element connectivity");
            if (iparent) {
                size = 4 * (eset->end - eset->start + 1);
                eset->parent = (int *) malloc(size * sizeof (int));
                if (NULL == eset->conn)
                    FATAL("read_zone_element", "malloc failed for parent data");
            }
            if (cg_elements_read(cgnsfn, cgnsbase, nz, ns,
                    eset->conn, eset->parent))
                FATAL("read_zone_element", NULL);
        }
    }
    return z->nesets;
}

/*---------- structured_elements --------------------------------------
 * build elements from structured zone
 *---------------------------------------------------------------------*/
int structured_elements(int nz) {
    int i, j, k, n, nelems = 0;
    ZONE *z = &Zones[nz - 1];
    ELEMSET *eset;

    if (z->type != Structured)
        return 0;

    switch (z->idim) {
        case 2:
            /* Structured Surface Grid */
            FATAL("structured_elements", "Surface Grid not supported");
            break;
        case 3:
            /* Structured Volume Grid */
            nelems = (z->dim[0] - 1) * (z->dim[1] - 1) * (z->dim[2] - 1);
            break;
    }
    if (nelems == 0)
        return 0;

    z->nesets = 1;
    z->esets = eset = new_elemset(1);
    strcpy(eset->name, "StructuredGridElements");
    eset->type = HEXA_8;
    eset->start = 1;
    eset->end = nelems;
    eset->nbndry = 0;
    eset->conn = (int *) malloc(8 * nelems * sizeof (int));
    if (NULL == eset->conn)
        FATAL("structured_elements", "malloc failed for element connectivity");
    eset->parent = NULL;
    for (n = 0, k = 1; k < z->dim[2]; k++) {
        for (j = 1; j < z->dim[1]; j++) {
            for (i = 1; i < z->dim[0]; i++) {
                eset->conn[n++] = vertex_index(z, i, j, k) + 1;
                eset->conn[n++] = vertex_index(z, i + 1, j, k) + 1;
                eset->conn[n++] = vertex_index(z, i + 1, j + 1, k) + 1;
                eset->conn[n++] = vertex_index(z, i, j + 1, k) + 1;
                eset->conn[n++] = vertex_index(z, i, j, k + 1) + 1;
                eset->conn[n++] = vertex_index(z, i + 1, j, k + 1) + 1;
                eset->conn[n++] = vertex_index(z, i + 1, j + 1, k + 1) + 1;
                eset->conn[n++] = vertex_index(z, i, j + 1, k + 1) + 1;
            }
        }
    }
    return 1;
}

/*---------- read_zone_interface --------------------------------------
 * read zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/
int read_zone_interface(int nz) {
    int i, j, n, ni, range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz - 1];
    INTERFACE *ints;

    if (cg_n1to1(cgnsfn, cgnsbase, nz, &z->nints))
        FATAL("read_zone_interface", NULL);
    if (z->nints) {
        z->ints = ints = new_interface(z->nints);
        for (ni = 1; ni <= z->nints; ni++, ints++) {
            if (cg_1to1_read(cgnsfn, cgnsbase, nz, ni,
                    ints->name, ints->d_name, (int *) range,
                    (int *) d_range, (int *) ints->transform))
                FATAL("read_zone_interface", NULL);
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
    ZONE *z = &Zones[nz - 1];
    CONNECT *conns;

    if (cg_nconns(cgnsfn, cgnsbase, nz, &z->nconns))
        FATAL("read_zone_connect", NULL);
    if (z->nconns) {
        z->conns = conns = new_connect(z->nconns);
        for (nc = 1; nc <= z->nconns; nc++, conns++) {
            if (cg_conn_info(cgnsfn, cgnsbase, nz, nc,
                    conns->name, &location, &type, &ptype,
                    &conns->npnts, conns->d_name, &d_ztype,
                    &d_ptype, &d_datatype, &conns->d_npnts))
                FATAL("read_zone_connect", NULL);
            conns->location = location;
            conns->type = type;
            conns->ptype = ptype;
            conns->d_ztype = d_ztype;
            conns->d_ptype = d_ptype;
            npnts = conns->npnts * z->idim;
            conns->pnts = (int *) calloc(npnts, sizeof (int));
            npnts = conns->d_npnts * z->idim;
            conns->d_pnts = (int *) calloc(npnts, sizeof (int));
            if (NULL == conns->pnts || NULL == conns->d_pnts)
                FATAL("read_zone_connect",
                    "malloc failed for connectivity point arrays");
            if (cg_conn_read(cgnsfn, cgnsbase, nz, nc,
                    conns->pnts, Integer, conns->d_pnts))
                FATAL("read_zone_connect", NULL);
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
    ZONE *z = &Zones[nz - 1];
    BOCO *bocos;

    if (cg_nbocos(cgnsfn, cgnsbase, nz, &z->nbocos))
        FATAL("read_zone_boco", NULL);
    if (z->nbocos) {
        z->bocos = bocos = new_boco(z->nbocos);
        for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
            if (cg_boco_info(cgnsfn, cgnsbase, nz, nb, bocos->name,
                    &bctype, &ptype, &bocos->npnts, bocos->n_index,
                    &bocos->n_cnt, &datatype, &ndatasets))
                FATAL("read_zone_boco", NULL);
            bocos->type = bctype;
            bocos->ptype = ptype;
            bocos->n_type = datatype;
            npnts = bocos->npnts * z->idim;
            bocos->pnts = (int *) calloc(npnts, sizeof (int));
            if (NULL == bocos->pnts)
                FATAL("read_zone_boco",
                    "calloc failed for boco point arrays");
            if (bocos->n_cnt) {
                bocos->n_list = (double *) calloc(bocos->n_cnt, sizeof (double));
                if (NULL == bocos->n_list)
                    FATAL("read_zone_boco",
                        "calloc failed for boco normal list");
            }
            if (cg_boco_read(cgnsfn, cgnsbase, nz, nb,
                    bocos->pnts, bocos->n_list))
                FATAL("read_zone_boco", NULL);
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
    ZONE *z = &Zones[nz - 1];
    SOLUTION *sols;

    if (cg_nsols(cgnsfn, cgnsbase, nz, &z->nsols))
        FATAL("read_zone_solution", NULL);

    if (z->nsols) {
        z->sols = sols = new_solution(z->nsols);

        for (ns = 1; ns <= z->nsols; ns++, sols++) {
            if (cg_sol_info(cgnsfn, cgnsbase, nz, ns,
                    sols->name, &location))
                FATAL("read_zone_solution", NULL);

            sols->location = location;

            if (z->type == Structured) {
                if (sols->location == Vertex) {
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 2; j++)
                            sols->rind[i][j] = 0;
                    sols->size = z->dim[0] * z->dim[1] * z->dim[2];
                } else if (sols->location == CellCenter) {
                    if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                            "FlowSolution_t", ns, "end"))
                        FATAL("read_zone_solution", NULL);

                    if (cg_rind_read((int *) sols->rind)) {
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 2; j++)
                                sols->rind[i][j] = 0;
                    }

                    sols->size = 1;

                    switch (z->idim) {
                        case 2:
                            /* Structured Surface Grid */
                            for (i = 0; i < 2; i++) {
                                sols->size *= (z->dim[i] - 1 +
                                        sols->rind[i][0] + sols->rind[i][1]);
                            }
                            break;
                        case 3:
                            /* Structured Volume Grid */
                            for (i = 0; i < 3; i++) {
                                sols->size *= (z->dim[i] - 1 +
                                        sols->rind[i][0] + sols->rind[i][1]);
                            }
                            break;
                    }
                } else
                    FATAL("read_zone_solution",
                        "solution location not Vertex or CellCenter");
            }                /* For Unstructured Grid: Only Vertex and CellCenter Support */
            else {
                sols->size = sols->location == Vertex ? z->dim[0] : z->dim[1];
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 2; j++)
                        sols->rind[i][j] = 0;
            }

            if (cg_nfields(cgnsfn, cgnsbase, nz, ns, &sols->nflds))
                FATAL("read_zone_solution", NULL);

            if (sols->nflds) {
                sols->flds = new_field(sols->nflds, 0);
                for (i = 0; i < sols->nflds; i++) {
                    if (cg_field_info(cgnsfn, cgnsbase, nz, ns, i + 1,
                            &datatype, sols->flds[i].name))
                        FATAL("read_zone_solution", NULL);
                }
            }

            /* get units */
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                    "FlowSolution_t", ns, "end"))
                FATAL("read_zone_solution", NULL);

            if (!read_units(sols->units)) {
                for (i = 0; i < 5; i++)
                    sols->units[i] = z->units[i];
            }

            if (cg_dataclass_read((DataClass_t *) & sols->dataclass))
                sols->dataclass = z->dataclass;

            /* get descriptors */
            if (cg_ndescriptors(&nd))
                FATAL("cg_ndescriptors", NULL);

            if (nd) {
                sols->ndesc = nd;
                sols->desc = new_desc(nd);

                for (i = 0; i < nd; i++) {
                    if (cg_descriptor_read(i + 1, sols->desc[i].name,
                            &sols->desc[i].desc))
                        FATAL("cg_descriptor_read", NULL);
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
    ZONE *z = &Zones[nz - 1];
    SOLUTION *s = &z->sols[ns - 1];
    FIELD *f;

    if (z->type == Structured) {
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = z->dim[n];
        }
        if (s->location == CellCenter) {
            switch (z->idim) {
                case 2:
                    /* Structured Surface Grid */
                    for (n = 0; n < 2; n++)
                        max[n] += s->rind[n][0] + s->rind[n][1] - 1;
                    break;
                case 3:
                    /* Structured Volume Grid */
                    for (n = 0; n < 3; n++)
                        max[n] += s->rind[n][0] + s->rind[n][1] - 1;
                    break;
            }
        }
    } else {
        /* Unstructured Grid */
        for (n = 0; n < 3; n++) {
            min[n] = 1;
            max[n] = s->size;
        }
    }

    if (nf) {
        is = ie = nf;
    } else {
        is = 1;
        ie = s->nflds;
    }

    f = &s->flds[is - 1];

    for (nf = is; nf <= ie; nf++, f++) {
        if (cg_field_info(cgnsfn, cgnsbase, nz, ns, nf,
                &datatype, f->name))
            FATAL("read_solution_field", NULL);

        f->id = nf;
        f->datatype = datatype;
        f->data = (double *) malloc(s->size * sizeof (double));

        if (NULL == f->data)
            FATAL("read_solution_field",
                "malloc failed for solution field data");

        if (cg_field_read(cgnsfn, cgnsbase, nz, ns, f->name,
                RealDouble, min, max, f->data))
            FATAL("read_solution_field", NULL);

        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", nz,
                "FlowSolution_t", ns, "DataArray_t", nf, "end"))
            FATAL("read_solution_field", NULL);

        if (!read_units(f->units)) {
            for (n = 0; n < 5; n++)
                f->units[n] = s->units[n];
        }

        /* read data class, conversion and exponents */
        if (cg_dataclass_read((DataClass_t *) & f->dataclass))
            f->dataclass = s->dataclass;

        if (cg_conversion_info(&datatype))
            f->dataconv[0] = 1.0;
        else {
            f->convtype = datatype;
            if (datatype == RealSingle) {
                float conv[2];
                if (cg_conversion_read(conv))
                    FATAL("read_solution_field", NULL);
                for (n = 0; n < 2; n++)
                    f->dataconv[n] = conv[n];
            } else if (datatype == RealDouble) {
                if (cg_conversion_read(f->dataconv))
                    FATAL("read_solution_field", NULL);
            } else
                FATAL("cg_conversion_info", "invalid data type");
        }

        if (!cg_exponents_info(&datatype)) {
            f->exptype = datatype;
            if (datatype == RealSingle) {
                float exp[5];
                if (cg_exponents_read(exp))
                    FATAL("read_solution_field", NULL);
                for (n = 0; n < 5; n++)
                    f->exponent[n] = exp[n];
            } else if (datatype == RealDouble) {
                if (cg_exponents_read(f->exponent))
                    FATAL("read_solution_field", NULL);
            } else
                FATAL("cg_exponents_info", "invalid data type");
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

/*---------- write_cgns -----------------------------------------------
 * write the CGNS file
 *---------------------------------------------------------------------*/
void write_cgns(void) {
    int nz, ns;

    write_zones();
    for (nz = 1; nz <= nZones; nz++) {
        write_zone_grid(nz);
        write_zone_element(nz);
        write_zone_interface(nz);
        write_zone_connect(nz);
        write_zone_boco(nz);
        for (ns = 1; ns <= Zones[nz - 1].nsols; ns++) {
            write_zone_solution(nz, ns);
            write_solution_field(nz, ns, 0);
        }
    }
}

/*---------- write_zones ----------------------------------------------
 * write zone information to CGNS file
 *---------------------------------------------------------------------*/
void write_zones(void) {
    int n, nz, sizes[3][3];
    ZONE *z = Zones;

    for (n = 0; n < 5; n++) {
        /* Pre-intialized baseunits[5]={0,0,0,0,0} */
        if (baseunits[n]) {
            if (cg_goto(cgnsfn, cgnsbase, "end") ||
                    cg_units_write((MassUnits_t) baseunits[0],
                    (LengthUnits_t) baseunits[1],
                    (TimeUnits_t) baseunits[2],
                    (TemperatureUnits_t) baseunits[3],
                    (AngleUnits_t) baseunits[4]))
                FATAL("write_zones", NULL);
            break;
        }
    }
    /* Pre-initialized baseclass = 0 */
    if (baseclass) {
        if (cg_goto(cgnsfn, cgnsbase, "end") ||
                cg_dataclass_write((DataClass_t) baseclass))
            FATAL("write_zones", NULL);
    }
    /* write the zone information */
    for (nz = 0; nz < nZones; nz++, z++) {
        if (!z->id) continue;
        for (n = 0; n < 3; n++) {
            sizes[0][n] = z->dim[n];
            sizes[1][n] = z->dim[n] - 1;
            sizes[2][n] = 0;
        }
        if (cg_zone_write(cgnsfn, cgnsbase, z->name,
                (int *) sizes, (ZoneType_t) z->type, &z->id))
            FATAL("write_zones", NULL);
        for (n = 0; n < 5; n++) {
            if (z->units[n] && z->units[n] != baseunits[n]) {
                if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                        cg_units_write((MassUnits_t) z->units[0],
                        (LengthUnits_t) z->units[1],
                        (TimeUnits_t) z->units[2],
                        (TemperatureUnits_t) z->units[3],
                        (AngleUnits_t) z->units[4]))
                    FATAL("write_zones", NULL);
                break;
            }
        }
        if (z->dataclass && z->dataclass != baseclass) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                    cg_dataclass_write((DataClass_t) z->dataclass))
                FATAL("write_zones", NULL);
        }
        if (z->ndesc) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end"))
                FATAL("write_zones", NULL);
            for (n = 0; n < z->ndesc; n++) {
                if (cg_descriptor_write(z->desc[n].name, z->desc[n].desc))
                    FATAL("cg_descriptor_write", NULL);
            }
        }
    }
}

/*---------- write_zone_data ------------------------------------------
 * write all zone data
 *---------------------------------------------------------------------*/
void write_zone_data(int nz) {
    int n, ns, sizes[3][3];
    ZONE *z = &Zones[nz - 1];

    /* write the zone information */
    for (n = 0; n < 3; n++) {
        sizes[0][n] = z->dim[n];
        sizes[1][n] = z->dim[n] - 1;
        sizes[2][n] = 0;
    }
    if (cg_zone_write(cgnsfn, cgnsbase, z->name,
            (int *) sizes, (ZoneType_t) z->type, &z->id))
        FATAL("write_zone_data", NULL);
    for (n = 0; n < 5; n++) {
        if (z->units[n] && z->units[n] != baseunits[n]) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id, "end") ||
                    cg_units_write((MassUnits_t) z->units[0],
                    (LengthUnits_t) z->units[1],
                    (TimeUnits_t) z->units[2],
                    (TemperatureUnits_t) z->units[3],
                    (AngleUnits_t) z->units[4]))
                FATAL("write_zone_data", NULL);
            break;
        }
    }
    write_zone_grid(nz);
    write_zone_element(nz);
    write_zone_interface(nz);
    write_zone_connect(nz);
    write_zone_boco(nz);
    for (ns = 1; ns <= Zones[nz - 1].nsols; ns++) {
        write_zone_solution(nz, ns);
        write_solution_field(nz, ns, 0);
    }
}

/*---------- write_zone_grid ------------------------------------------
 * write zone grid coordinates
 *---------------------------------------------------------------------*/
void write_zone_grid(int nz) {
    int n, nc;
    ZONE *z = &Zones[nz - 1];

    if (z->verts == NULL || (z->vertflags & 7) == 0)
        return;
    if (z->datatype == RealSingle) {
        float *xyz = (float *) malloc(z->nverts * sizeof (float));
        if (NULL == xyz)
            FATAL("write_zone_grid",
                "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float) z->verts[n].x;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateX", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float) z->verts[n].y;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateY", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = (float) z->verts[n].z;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealSingle,
                    "CoordinateZ", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        free(xyz);
    } else {
        double *xyz = (double *) malloc(z->nverts * sizeof (double));
        if (NULL == xyz)
            FATAL("write_zone_grid",
                "malloc failed for coordinate working array");
        if ((z->vertflags & 1) == 1) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].x;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateX", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 2) == 2) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].y;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateY", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        if ((z->vertflags & 4) == 4) {
            for (n = 0; n < z->nverts; n++)
                xyz[n] = z->verts[n].z;
            if (cg_coord_write(cgnsfn, cgnsbase, z->id, RealDouble,
                    "CoordinateZ", xyz, &nc)) FATAL("write_zone_grid", NULL);
        }
        free(xyz);
    }
}

/*---------- write_zone_element ---------------------------------------
 * write zone element sets and elements
 *---------------------------------------------------------------------*/
void write_zone_element(int nz) {
    int ns;
    ZONE *z = &Zones[nz - 1];
    ELEMSET *eset = z->esets;

    if (eset == NULL || z->type == Structured)
        return;
    for (ns = 1; ns <= z->nesets; ns++, eset++) {
        if (eset->id) {
            if (cg_section_write(cgnsfn, cgnsbase, z->id,
                    eset->name, (ElementType_t) eset->type,
                    eset->start, eset->end, eset->nbndry,
                    eset->conn, &eset->id))
                FATAL("write_zone_element", NULL);
            if (eset->parent != NULL &&
                    cg_parent_data_write(cgnsfn, cgnsbase, z->id,
                    eset->id, eset->parent))
                FATAL("write_zone_element", NULL);
        }
    }
}

/*---------- write_zone_interface -------------------------------------
 * write zone 1 to 1 interfaces
 *---------------------------------------------------------------------*/
void write_zone_interface(int nz) {
    int i, j, ni, range[2][3], d_range[2][3];
    ZONE *z = &Zones[nz - 1];
    INTERFACE *ints = z->ints;

    if (ints == NULL)
        return;
    for (ni = 1; ni <= z->nints; ni++, ints++) {
        if (ints->id) {
            for (j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++) {
                    range[j][i] = ints->range[i][j];
                    d_range[j][i] = ints->d_range[i][j];
                }
            }
            if (cg_1to1_write(cgnsfn, cgnsbase, z->id,
                    ints->name, ints->d_name, (int *) range,
                    (int *) d_range, ints->transform, &ints->id))
                FATAL("write_zone_interface", NULL);
        }
    }
}

/*---------- write_zone_connect ---------------------------------------
 * write zone connectivities
 *---------------------------------------------------------------------*/
void write_zone_connect(int nz) {
    int nc;
    ZONE *z = &Zones[nz - 1];
    CONNECT *conns = z->conns;

    if (conns == NULL)
        return;
    for (nc = 1; nc <= z->nconns; nc++, conns++) {
        if (conns->id &&
                cg_conn_write(cgnsfn, cgnsbase, z->id,
                conns->name, (GridLocation_t) conns->location,
                (GridConnectivityType_t) conns->type,
                (PointSetType_t) conns->ptype, conns->npnts, conns->pnts,
                conns->d_name, (ZoneType_t) conns->d_ztype,
                (PointSetType_t) conns->d_ptype, Integer, conns->d_npnts,
                conns->d_pnts, &conns->id))
            FATAL("write_zone_connect", NULL);
    }
}

/*---------- write_zone_bocos -----------------------------------------
 * write zone boundary conditions
 *---------------------------------------------------------------------*/
void write_zone_boco(int nz) {
    int nb;
    ZONE *z = &Zones[nz - 1];
    BOCO *bocos = z->bocos;

    if (bocos == NULL)
        return;
    for (nb = 1; nb <= z->nbocos; nb++, bocos++) {
        if ((bocos->id) &&
                (cg_boco_write(cgnsfn, cgnsbase, z->id,
                bocos->name, (BCType_t) bocos->type,
                (PointSetType_t) bocos->ptype, bocos->npnts,
                bocos->pnts, &bocos->id) ||
                cg_boco_normal_write(cgnsfn, cgnsbase, z->id, bocos->id,
                bocos->n_index, bocos->n_cnt, (DataType_t) bocos->n_type,
                bocos->n_list)))
            FATAL("write_zone_boco", NULL);
    }
}

/*---------- write_zone_solution --------------------------------------
 * write zone solution
 *---------------------------------------------------------------------*/
void write_zone_solution(int nz, int ns) {
    int n;
    ZONE *z = &Zones[nz - 1];
    SOLUTION *s;

    if (z->sols == NULL || ns < 1 || ns > z->nsols)
        return;
    s = &z->sols[ns - 1];
    if (cg_sol_write(cgnsfn, cgnsbase, z->id, s->name,
            (GridLocation_t) s->location, &s->id))
        FATAL("write_zone_solution", NULL);
    if (z->type == Structured && s->location == CellCenter) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end") ||
                cg_rind_write((int *) s->rind))
            FATAL("write_zone_solution", NULL);
    }
    for (n = 0; n < 5; n++) {
        if (s->units[n] && s->units[n] != z->units[n]) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "end") ||
                    cg_units_write((MassUnits_t) s->units[0],
                    (LengthUnits_t) s->units[1],
                    (TimeUnits_t) s->units[2],
                    (TemperatureUnits_t) s->units[3],
                    (AngleUnits_t) s->units[4]))
                FATAL("write_zone_solution", NULL);
            break;
        }
    }
    if (s->dataclass && s->dataclass != z->dataclass) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end") ||
                cg_dataclass_write((DataClass_t) s->dataclass))
            FATAL("write_zone_solution", NULL);
    }
    if (s->ndesc) {
        if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                "FlowSolution_t", s->id, "end"))
            FATAL("write_zone_solution", NULL);
        for (n = 0; n < s->ndesc; n++) {
            if (cg_descriptor_write(s->desc[n].name, s->desc[n].desc))
                FATAL("cg_descriptor_write", NULL);
        }
    }
}

/*---------- write_solution_field -------------------------------------
 * write solution field data for a zone
 *---------------------------------------------------------------------*/
void write_solution_field(int nz, int ns, int nf) {
    int n, is, ie;
    float *data = NULL;
    ZONE *z = &Zones[nz - 1];
    SOLUTION *s = &z->sols[ns - 1];
    FIELD *f;

    if (nf) {
        is = ie = nf;
    } else {
        is = 1;
        ie = s->nflds;
    }
    f = &s->flds[is - 1];
    for (nf = is; nf <= ie; nf++, f++) {
        if (f->data == NULL) continue;
        if (f->datatype == RealSingle) {
            if (data == NULL) {
                data = (float *) malloc(s->size * sizeof (float));
                if (NULL == data)
                    FATAL("write_solution_field",
                        "malloc failed for working array");
            }
            for (n = 0; n < s->size; n++)
                data[n] = (float) f->data[n];
            if (cg_field_write(cgnsfn, cgnsbase, z->id, s->id,
                    RealSingle, f->name, data, &f->id))
                FATAL("write_solution_field", NULL);
        } else {
            if (cg_field_write(cgnsfn, cgnsbase, z->id, s->id,
                    RealDouble, f->name, f->data, &f->id))
                FATAL("write_solution_field", NULL);
        }
        for (n = 0; n < 5; n++) {
            if (f->units[n] && f->units[n] != s->units[n]) {
                if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                        "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                        cg_units_write((MassUnits_t) f->units[0],
                        (LengthUnits_t) f->units[1],
                        (TimeUnits_t) f->units[2],
                        (TemperatureUnits_t) f->units[3],
                        (AngleUnits_t) f->units[4]))
                    FATAL("write_solution_field", NULL);
                break;
            }
        }
        if (f->dataclass && f->dataclass != s->dataclass) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_dataclass_write((DataClass_t) f->dataclass))
                FATAL("write_solution_field", NULL);
        }
        if (f->convtype == RealSingle) {
            float conv[2];
            for (n = 0; n < 2; n++)
                conv[n] = (float) f->dataconv[n];
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_conversion_write(RealSingle, conv))
                FATAL("write_solution_field", NULL);
        } else if (f->convtype == RealDouble) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_conversion_write(RealDouble, f->dataconv))
                FATAL("write_solution_field", NULL);
        } else {
        }
        if (f->exptype == RealSingle) {
            float exp[5];
            for (n = 0; n < 5; n++)
                exp[n] = (float) f->dataconv[n];
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_exponents_write(RealSingle, exp))
                FATAL("write_solution_field", NULL);
        } else if (f->exptype == RealDouble) {
            if (cg_goto(cgnsfn, cgnsbase, "Zone_t", z->id,
                    "FlowSolution_t", s->id, "DataArray_t", f->id, "end") ||
                    cg_exponents_write(RealDouble, f->exponent))
                FATAL("write_solution_field", NULL);
        } else {
        }
    }
    if (data != NULL) free(data);
}

/*---------- print_interface ------------------------------------------
 * Prints Interface Information
 *---------------------------------------------------------------------*/
void print_interface(ZONE *zone) {
    int ni;
    INTERFACE *ints;

    read_zone_interface(zone->id);
    ints = zone->ints;
    for (ni = 0; ni < zone->nints; ni++) {
        printf("  interface %d - %s\n", ints[ni].id, ints[ni].name);
        printf("    range       = %d %d %d %d %d %d\n",
                ints[ni].range[0][0], ints[ni].range[0][1],
                ints[ni].range[1][0], ints[ni].range[1][1],
                ints[ni].range[2][0], ints[ni].range[2][1]);
        printf("    donor name  = %s\n", ints[ni].d_name);
        printf("    donor range = %d %d %d %d %d %d\n",
                ints[ni].d_range[0][0], ints[ni].d_range[0][1],
                ints[ni].d_range[1][0], ints[ni].d_range[1][1],
                ints[ni].d_range[2][0], ints[ni].d_range[2][1]);
        printf("    transform   = %d %d %d\n", ints[ni].transform[0],
                ints[ni].transform[1], ints[ni].transform[2]);
        printf("    donor zone  = %d\n", ints[ni].d_zone);
    }
}

/*---------- print_connect --------------------------------------------
 * Prints Connectivity Information
 *---------------------------------------------------------------------*/
void print_connect(ZONE *zone) {
    int nc;
    CONNECT *conns;

    read_zone_connect(zone->id);
    conns = zone->conns;
    for (nc = 0; nc < zone->nconns; nc++) {
        printf("  connectivity %d - %s\n", conns[nc].id, conns[nc].name);
        printf("    type          = %d\n", conns[nc].type);
        printf("    location      = %d\n", conns[nc].location);
        printf("    pt type       = %d\n", conns[nc].ptype);
        printf("    points        = %d\n", conns[nc].npnts);
        printf("    donor name    = %s\n", conns[nc].d_name);
        printf("    donor pt type = %d\n", conns[nc].d_ptype);
        printf("    donor points  = %d\n", conns[nc].d_npnts);
        printf("    donor zone    = %d\n", conns[nc].d_zone);
    }
}

/*---------- print_solution -------------------------------------------
 * Prints Solution Information
 *---------------------------------------------------------------------*/
void print_solution(ZONE *zone) {
    int ns, nf;
    SOLUTION *sols;

    read_zone_solution(zone->id);
    sols = zone->sols;
    for (ns = 0; ns < zone->nsols; ns++) {
        printf("Solution No = %d\n", sols[ns].id);
        printf("Solution Name = %s\n", sols[ns].name);
        printf("Location = %d\n", sols[ns].location);
        printf("Rind Data = %d %d %d %d %d %d\n",
                sols[ns].rind[0][0], sols[ns].rind[0][1],
                sols[ns].rind[1][0], sols[ns].rind[1][1],
                sols[ns].rind[2][0], sols[ns].rind[2][1]);
        printf("Size = %d\n", sols[ns].size);
        printf("No of Fields = %d\n", sols[ns].nflds);
        for (nf = 0; nf < sols[ns].nflds; nf++)
            printf("\t%s\n", sols[ns].flds[nf].name);
    }
}

/*---------- print_GridLocation ---------------------------------------
 * Prints Grid Location
 *---------------------------------------------------------------------*/
void print_GridLocation(int location) {
    switch (location) {
        case 0:
            printf("Grid Location = GridLocationNuLL\n");
            break;
        case 1:
            printf("Grid Location = GridLocationUserDefined\n");
            break;
        case 2:
            printf("Grid Location = Vertex\n");
            break;
        case 3:
            printf("Grid Location = CellCenter\n");
            break;
        case 4:
            printf("Grid Location = FaceCenter\n");
            break;
        case 5:
            printf("Grid Location = IFaceCenter\n");
            break;
        case 6:
            printf("Grid Location = JFaceCenter\n");
            break;
        case 7:
            printf("Grid Location = KFaceCenter\n");
            break;
        case 8:
            printf("Grid Location = EdgeCenter\n");
            break;
    }
}

/*---------- print_BCDataType -----------------------------------------
 * Prints Boundary Condition Data Type
 *---------------------------------------------------------------------*/
void print_BCDataType(int type) {
    switch (type) {
        case 0:
            printf("Boundary Condition Data Type = BCDataTypeNull\n");
            break;
        case 1:
            printf("Boundary Condition Data Type = BCDataTypeUserDefined\n");
            break;
        case 2:
            printf("Boundary Condition Data Type = Dirichlet\n");
            break;
        case 3:
            printf("Boundary Condition Data Type = Neumann\n");
            break;
    }
}

/*---------- print_GridConnectivityType -------------------------------
 * Prints Grid Connectivity Type
 *---------------------------------------------------------------------*/
void print_GridConnectivityType(int type) {
    switch (type) {
        case 0:
            printf("Grid Connectivity Type = GridConnectivityTypeNull\n");
            break;
        case 1:
            printf("Grid Connectivity Type = GridConnectivityTypeUserDefined\n");
            break;
        case 2:
            printf("Grid Connectivity Type = Overset\n");
            break;
        case 3:
            printf("Grid Connectivity Type = Abutting\n");
            break;
        case 4:
            printf("Grid Connectivity Type = Abutting1to1\n");
            break;
    }
}

/*---------- print_PointSetType ---------------------------------------
 * Prints Point Set Type
 *---------------------------------------------------------------------*/
void print_PointSetType(int type) {
    switch (type) {
        case 0:
            printf("Point Set Type = PointSetTypeNull\n");
            break;
        case 1:
            printf("Point Set Type = PointSetTypeUserDefined\n");
            break;
        case 2:
            printf("Point Set Type = PointList\n");
            break;
        case 3:
            printf("Point Set Type = PointListDonor\n");
            break;
        case 4:
            printf("Point Set Type = PointRange\n");
            break;
        case 5:
            printf("Point Set Type = PointRangeDonor\n");
            break;
        case 6:
            printf("Point Set Type = ElementRange\n");
            break;
        case 7:
            printf("Point Set Type = ElementList\n");
            break;
        case 8:
            printf("Point Set Type = CellListDonor\n");
            break;
    }
}

/*---------- print_BCType ---------------------------------------------
 * Prints Boundary Condition Type
 *---------------------------------------------------------------------*/
void print_BCType(int type) {
    switch (type) {
        case 0:
            printf("Boundary Condition Type = BCTypeNull\n");
            break;
        case 1:
            printf("Boundary Condition Type = BCTypeUserDefined\n");
            break;
        case 2:
            printf("Boundary Condition Type = BCAxisymmetricWedge\n");
            break;
        case 3:
            printf("Boundary Condition Type = BCDegenerateLine\n");
            break;
        case 4:
            printf("Boundary Condition Type = BCDegeneratePoint\n");
            break;
        case 5:
            printf("Boundary Condition Type = BCDirichlet\n");
            break;
        case 6:
            printf("Boundary Condition Type = BCExtrapolate\n");
            break;
        case 7:
            printf("Boundary Condition Type = BCFarfield\n");
            break;
        case 8:
            printf("Boundary Condition Type = BCGeneral\n");
            break;
        case 9:
            printf("Boundary Condition Type = BCInflow\n");
            break;
        case 10:
            printf("Boundary Condition Type = BCInflowSubsonic\n");
            break;
        case 11:
            printf("Boundary Condition Type = BCInflowSupersonic\n");
            break;
        case 12:
            printf("Boundary Condition Type = BCNeumann\n");
            break;
        case 13:
            printf("Boundary Condition Type = BCOutflow\n");
            break;
        case 14:
            printf("Boundary Condition Type = BCOutflowSubsonic\n");
            break;
        case 15:
            printf("Boundary Condition Type = BCOutflowSupersonic\n");
            break;
        case 16:
            printf("Boundary Condition Type = BCSymmetryPlane\n");
            break;
        case 17:
            printf("Boundary Condition Type = BCSymmetryPolar\n");
            break;
        case 18:
            printf("Boundary Condition Type = BCTunnelInflow\n");
            break;
        case 19:
            printf("Boundary Condition Type = BCTunnelOutflow\n");
            break;
        case 20:
            printf("Boundary Condition Type = BCWall\n");
            break;
        case 21:
            printf("Boundary Condition Type = BCWallInviscid\n");
            break;
        case 22:
            printf("Boundary Condition Type = BCWallViscous\n");
            break;
        case 23:
            printf("Boundary Condition Type = BCWallViscousHeatFlux\n");
            break;
        case 24:
            printf("Boundary Condition Type = BCWallViscousIsothermal\n");
            break;
        case 25:
            printf("Boundary Condition Type = FamilySpecified\n");
            break;
    }
}

/*---------- print_DataType -------------------------------------------
 * Prints Data Type
 *---------------------------------------------------------------------*/
void print_DataType(int type) {
    switch (type) {
        case 0:
            printf("Data Type = DataTypeNull\n");
            break;
        case 1:
            printf("Data Type = DataTypeUserDefined\n");
            break;
        case 2:
            printf("Data Type = Integer\n");
            break;
        case 3:
            printf("Data Type = RealSingle\n");
            break;
        case 4:
            printf("Data Type = RealDouble\n");
            break;
        case 5:
            printf("Data Type = Character\n");
            break;
    }
}

/*---------- print_ZoneType -------------------------------------------
 * Prints Zone Type
 *---------------------------------------------------------------------*/
void print_ZoneType(int type) {
    switch (type) {
        case 0:
            printf("Zone Type = ZoneTypeNull\n");
            break;
        case 1:
            printf("Zone Type = ZoneTypeUserDefined\n");
            break;
        case 2:
            printf("Zone Type = Structured\n");
            break;
        case 3:
            printf("Zone Type = Unstructured\n");
            break;
    }
}

