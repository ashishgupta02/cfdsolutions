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
 * File		CGNSPatch.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "CGNSIO.h"
#include "CGNSPatch.h"

/*--------------------------------------------------------------*/
int CGNSPatch_Connectivity(const char *file1, const char *file2) {
    int nbases;
    int i, nz, e;
    ZONE *z;
    ELEMSET *eset;

    /* Checks for existance of file */
    if (!file_exists(file1))
        FATAL(NULL, "File 1 does not exist");

    if (!file_exists(file1))
        FATAL(NULL, "File 2 does not exist");

    /* Open CGNS file 1 Read*/
    nbases = open_cgns(file1, 0);
    if (!nbases)
        FATAL(NULL, "No bases in CGNS file 1");

    if (nbases == 1)
        cgnsbase = 1;
    else {
        printf("Give base No to Post-Processed : ");
        scanf("%d", &cgnsbase);

        if (cgnsbase < 1 || cgnsbase > nbases)
            FATAL(NULL, "Invailed base index");
    }

    /* Read CGNS file 1 */
    read_cgns();

    /* Close cgns file 1 Interface */
    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);

    printf("Hi \n");

    /* Open CGNS file 2 Modify*/
    nbases = open_cgns(file2, 1);
    if (!nbases)
        FATAL(NULL, "No bases in CGNS file 2");

    /* Apply Connectivity Patch */
    for (nz = 1; nz <= nZones; nz++) {
        z = &Zones[nz - 1];
        /* Below will work properly only when nsections = 1 */
        for (eset = z->esets, e = 1; e <= z->nesets; e++, eset++) {
            for (i = 0; i < 3 * eset->end; i++)
                eset->conn[i] = eset->conn[i] + 1;
        }
        write_zone_element(nz);
    }

    /* Close cgns file 2 Interface */
    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);

    printf("CGNS File Connectivity Sucessfully Patched\n");

    return 0;
}

/*--------------------------------------------------------------*/
int CGNSPatch(const char *file1, const char *file2) {
    int nbases;
    int i, nz, bc;
    ZONE *z;
    BOCO *boco;

    /* Checks for existance of file */
    if (!file_exists(file1))
        FATAL(NULL, "File 1 does not exist");

    if (!file_exists(file1))
        FATAL(NULL, "File 2 does not exist");

    /* Open CGNS file 1 Read*/
    nbases = open_cgns(file1, 0);
    if (!nbases)
        FATAL(NULL, "No bases in CGNS file 1");

    if (nbases == 1)
        cgnsbase = 1;
    else {
        printf("Give base No to Post-Processed : ");
        scanf("%d", &cgnsbase);

        if (cgnsbase < 1 || cgnsbase > nbases)
            FATAL(NULL, "Invailed base index");
    }

    /* Read CGNS file 1 */
    read_cgns();

    /* Close cgns file 1 Interface */
    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);

    /* Open CGNS file 2 Modify*/
    nbases = open_cgns(file2, 1);
    if (!nbases)
        FATAL(NULL, "No bases in CGNS file 2");

    /* Read Boundary Condition */
    for (nz = 1; nz <= nZones; nz++) {
        z = &Zones[nz - 1];
        for (boco = z->bocos, bc = 1; bc <= z->nbocos; bc++, boco++) {
            for (i = 0; i < boco->npnts; i++)
                boco->pnts[i] = boco->pnts[i] - 1;
        }
        write_zone_boco(nz);
    }

    /* Close cgns file 2 Interface */
    if (cg_close(cgnsfn))
        FATAL(NULL, NULL);

    printf("CGNS File Sucessfully Patched\n");

    return 0;
}
