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
 * File		Read_CGNS.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "CGNSIO.h"
#include "Read_CGNS.h"

/*---------- Read ----------------------------------------------
 * Read CGNS File
 * -------------------------------------------------------------*/
int Read(const char *file) {
    int nbases, celldim, phydim;
    char basename[33];
    ZONE *z;
    BOCO *bocos;
    int nz, nb;

    /* Checks for existance of file */
    if (!file_exists(file))
        FATAL(NULL, "File does not exist");

    /* Read CGNS file */
    printf("Reading CGNS file from %s\n", file);
    fflush(stdout);

    /* Open CGNS File */
    nbases = open_cgns(file, 0);
    if (!nbases)
        FATAL(NULL, "No bases in CGNS file");

    cg_version(cgnsfn, &Version);
    printf("File version = %lf\n", Version);

    printf("No of Bases in file = %d\n", nbases);

    if (nbases == 1)
        cgnsbase = 1;
    else {
        printf("Give base No to browse : ");
        scanf("%d", &cgnsbase);

        if (cgnsbase < 1 || cgnsbase > nbases)
            FATAL(NULL, "Invailed base index");
    }

    if (cg_base_read(cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL(NULL, NULL);

    printf("Using base %d - %s\n", cgnsbase, basename);
    printf("Cell dimension     = %d\n", celldim);
    printf("Physical dimension = %d\n", phydim);

    read_cgns();

    printf("No of zones nodes = %d\n", nZones);

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        printf("\nZone No = %d\nZone Name = %s\n", z->id, z->name);

        print_ZoneType(z->type);
        switch (z->type) {
            case 2:
                printf("Dimensions = %d x %d x %d\n",
                        z->dim[0], z->dim[1], z->dim[2]);
                break;
            case 3:
                printf("Dimensions = %d\n", z->dim[0]);
                break;
        }

        if (z->nbocos) {
            printf("No of boundary conditions = %d\n", z->nbocos);
            for (bocos = z->bocos, nb = 1; nb <= z->nbocos; nb++, bocos++) {
                printf("Boundry condition name = %s\n", bocos->name);

                print_BCType(bocos->type);
                print_PointSetType(bocos->ptype);
                print_BCDataType(bocos->n_type);
            }
        }

        if (z->nsols) {
            printf("No of solutions node = %d\n", z->nsols);
            print_solution(z);
        }
    }

    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);

    return 0;
}

/*---------- Info ----------------------------------------------
 * Information CGNS File
 * -------------------------------------------------------------*/
int Info(const char *file) {
    int nz, nbases, celldim, phydim;
    char basename[33];
    ZONE *z;

    if (!file_exists(file))
        FATAL(NULL, "CGNSfile does not exist or is not a file");

    /* read CGNS file */
    printf("reading CGNS file from %s\n", file);
    fflush(stdout);
    nbases = open_cgns(file, 0);
    if (!nbases) FATAL(NULL, "no bases in CGNS file");

    cg_version(cgnsfn, &Version);
    printf("File version = %lf\n", Version);
    printf("No of Bases in file = %d\n", nbases);

    if (nbases == 1)
        cgnsbase = 1;
    else {
        printf("Give base No to browse");
        scanf("%d", &cgnsbase);

        if (cgnsbase < 1 || cgnsbase > nbases)
            FATAL(NULL, "Invailed base index");
    }

    if (cg_base_read(cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL(NULL, NULL);

    printf("using base %d - %s\n", cgnsbase, basename);
    printf("cell dimension     = %d\n", celldim);
    printf("physical dimension = %d\n", phydim);

    read_zones();

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        printf("\nzone %d - %s\n", z->id, z->name);
        printf("type      = %d\n", z->type);
        printf("dimension = %d x %d x %d\n", z->dim[0],
                z->dim[1], z->dim[2]);
        printf("1to1      = %d\n", z->nints);
        print_interface(z);

        printf("connects  = %d\n", z->nconns);
        print_connect(z);

        printf("solutions = %d\n", z->nsols);
        print_solution(z);
    }

    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);
    exit(0);
}

