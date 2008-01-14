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
 * File		PostProcessing_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "CGNSIO.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPInitializeSolution_2D.h"
#include "FluidProperties.h"
#include "ReferenceQuantities.h"
#include "PostAnalysis_2D.h"
#include "PPSolutionGradients_2D.h"
#include "Graphics_2D.h"

/*---------- PostProcesssing --------------------------------------
 * Post-Analysis of CFD Data
 * Developed only for CFD-Tutor with single base - single zone 2D
 * -------------------------------------------------------------*/
int PostProcessing (const char *file) {
	int visual;
	int nbases, ncoords, celldim, phydim;
	char basename[33];
	ZONE *z;
	int i, j, nz;
	SOLUTION *s;

	/* Checks for existance of file */
	if (!file_exists(file))
		FATAL (NULL, "File does not exist");
	
	/* Read CGNS file */
	printf ("Reading CGNS file from %s\n", file);
	fflush (stdout);
	
	/* Open CGNS File */
	nbases = open_cgns (file, 0);
	if (!nbases)
		FATAL (NULL, "No bases in CGNS file");

	cg_version (cgnsfn, &Version);
	printf("File version = %lf\n", Version);

	printf("No of Bases = %d\n", nbases);

	if (nbases == 1)
		cgnsbase = 1;
	else {
		printf ("Give base No to Post-Processed : ");
		scanf ("%d", &cgnsbase);

		if (cgnsbase < 1 || cgnsbase > nbases)
			FATAL (NULL, "Invailed base index");
	}

	if (cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim) ||
		cg_nzones (cgnsfn, cgnsbase, &nZones))
		FATAL(NULL, NULL);
	
	if (celldim != 2 || phydim != 2)
		FATAL (NULL, "Not A 2D Grid");
	
	if (nZones > 1) 
		FATAL (NULL, "Not A Single Zone");
	

	printf ("Base Number = %d\n", cgnsbase);
	printf ("Base Name = %s\n", basename);
	printf ("Cell Dimension = %d\n", celldim);
	printf ("Physical Dimension = %d\n", phydim);

	read_cgns ();
	
	printf ("No of Zones = %d\n", nZones);
	
	/* Print Zone Information */
	for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
		printf ("Zone No = %d\n", z->id);
		printf ("Zone Name = %s\n", z->name);
		if (cg_ncoords (cgnsfn, cgnsbase, nz, &ncoords))
			FATAL ("Post-Processing ncoords", NULL);
		if(ncoords != 2)
			FATAL (NULL, "No 3D Support Now"); 

		print_ZoneType (z->type);
		switch (z->type) {
			case 2:
				/* For Structured Grid */
				printf ("Dimensions = %d x %d x %d\n",z->dim[0], z->dim[1], z->dim[2]);
				break;
			case 3:
				/* For Unstructured Grid */
				printf ("No of Nodes = %d\n", z->dim[0]);
				printf ("No of Cells = %d\n", z->dim[1]);
				break;
			default:
				/* Unknown Type Grids */
				FATAL(NULL,"Unknown Grid Type");
		}
	}
	
	/* Print Available Solution Fields */
	for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
		if (z->nsols) {
			printf ("No of Solutions Nodes = %d\n", z->nsols);
			for (s = z->sols, i = 1; i <= z->nsols; i++, s++) {
				printf("Solution Node = %d\n", i);
				
				if(s->nflds == 0)
					FATAL(NULL,"No Solution Fields: Exiting");
					
				printf("\tNo of Fields = %d\n", s->nflds);
				for (j = 0; j < s->nflds; j++) {
					printf("\tField = %s\n", s->flds[j].name);
				}
			}
		}
		else
			FATAL(NULL,"No Solution Node Available: Exiting");
	}
		
	/* Initialize Data Structure Only Once */
	for (z = Zones, nz = 1; nz <= nZones; nz++, z++)
		InitializeGrid_2D(z);
	

	/* Do Post-Processing Operation Below */
	/*****************/
	
	/* Post Analysis Function Calls */
	for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
		InitializeSolutionCGNS_2D(z);
		InitializeSolution_2D();
		GetSolutionList_2D();
			
		/* Visual Mode Option */
		printf("Goto Visual Mode (0/1):");
		scanf("%d", &visual);
		
		if (visual == 1)
			Graphics_2D();
		else {
			GetFluidProperties();
			GetReferenceQuantities();
			PostAnalysis_2D();
		}	
	}
	
	/*****************/
	/* Do Post-Processing Operation Above */
		
	/* Close cgns Interface */
	if (cg_close (cgnsfn))
		FATAL (NULL, NULL);

	return 0;
}
