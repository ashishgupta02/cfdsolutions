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
 * File		Main.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "Read_CGNS.h"
#include "CGNSPatch.h"
#include "PostProcessing_2D.h"

/* Command line options */

static char options[] = "acirput";

static char *usgmsg[] = {
                         "usage  : Post-Processing [options]",
                         "options:",
                         "	-a		= Append to CGNS file",
                         "	-c		= Convert to CGNS formate",
                         "	-i		= CGNS Information",
                         "	-r		= Read a CGNS file",
                         "	-p		= Post-Processing",
						 "	-u		= Patch CGNS File",
						 "	-t		= Patch CGNS Connectivity", 
                         NULL
                 };

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static void usage (char **usgmsg) {
	int n;

	for (n = 0; NULL != usgmsg[n]; n++)
		fprintf (stderr, "%s\n", usgmsg[n]);

	exit (1);
}

/*---------- arguments ---------------------------------------------------
 * Get option letter from argument vector or terminates on error
 *----------------------------------------------------------------------*/

static int arguments (int argc, char **argv) {
	unsigned int n;

	if (argc < 2) {
		printf (" Too few arguments\n");
		usage (usgmsg);
	}

	if (argc < 3) {
		printf ("No input file\n");
		exit (1);
	}
	if ((argc == 3) || (argc == 4)) {
		for (n = 0; n < strlen(options); n++)
			if (argv[1][1] == options[n])
				return n;
		MSG("Invalid option");
		usage (usgmsg);
	}
	
	return 0;
}

/*-----------------------------------------------------------------------------
 * Main Function
 * ---------------------------------------------------------------------------*/

int main (int argc, char *argv[]) {

	int err, sel;
	char File[100], File1[100]; 
	
	sel = arguments (argc, argv);
	strcpy (File,argv[2]);
	
	/* Selecting the module */
	switch(sel){
	case 0:
		/* Append CGNS File */
		MSG("main: Not Yet Implemented");
		break;
	case 1:
		/* Write New CGNS File */
		/*
		if((err = convert(File)) == 1) {
			FATAL("Convert","Unknown error");
			exit(1);
		}
		*/
		MSG("main: Not Yet Implemented");
		break;
	case 2:
		/* Cgns File Information */
		err = Info(File);
		if(err == 1) {
			MSG("main: Info CGNS");
			exit(1);
		}
		break;
	case 3:
		/* Open CGNS File In Read Mode */
		if((err = Read(File)) == 1) {
			MSG("main: Read CGNS");
			exit(1);
		}
		break;
	case 4:
		/* Calls the Post-Processing Module */
		if((err = PostProcessing(File)) == 1) {
			MSG("main: Post-Processing CGNS");
			exit(1);
		}
		break;
	case 5:
		/* CGNS Patch */
		strcpy(File1, argv[3]);
		if((err = CGNSPatch(File, File1)) == 1) {
			MSG("main: Patch Failed");
			exit(1);
		}
		break;
	case 6:
		/* CGNS Connectivity Patch */
		strcpy(File1, argv[3]);
		if((err = CGNSPatch_Connectivity(File, File1)) == 1) {
			MSG("main: Patch Connectivity Failed");
			exit(1);
		}
		break;
	}
	return 0;
}
