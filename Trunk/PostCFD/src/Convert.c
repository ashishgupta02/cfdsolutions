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
 * File		Convert.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include <cgnslib.h>
#include "CGNSIO.h"

/* This function opens a input file and convert it into CGNS formate */
int convert(const char *file){
	FILE *fp;
	char str[10], name[33];
	double *xcoord = NULL, *ycoord = NULL, *zcoord = NULL;
	int i, type, *curve, *range;
	
	/* Checks for existance of file */
	if (!file_exists( file )){
		FATAL(NULL, "File does not exist"); 
	}
	
	strcpy(name,file);
	strcat(name,".cgns");
	printf("Output File: %s\n", name);
	
	/* Open the input file */
	if ((fp = fopen(file,"r")) == NULL){
		printf("Cannot open input file. \n");
		exit(1);
	}
	
	fscanf(fp, "%d%s\n", &type, str);
	
	switch(type){
		case 31:
			printf("Data of Type : %s\n", str);
			curve = (int *)malloc( sizeof(int));
			range = (int *)malloc( sizeof(int));
			fscanf(fp, "%d%d\n", &curve[0], &range[0]);
			xcoord = (double *)malloc( range[0]*sizeof(double));
			ycoord = (double *)malloc( range[0]*sizeof(double));
			zcoord = (double *)malloc( range[0]*sizeof(double));
			for(i = 0; i < range[0]; i++){
				fscanf(fp,"%lf%lf%lf", &xcoord[i], &ycoord[i], &zcoord[i]);
				printf("%lf\t%lf\t%lf\n", xcoord[i], ycoord[i], zcoord[i]);
			}	
			break;
		case 51:
			printf("Data of Type : %s\n", str);
			curve = (int *)malloc( 2*sizeof(int));
			range = (int *)malloc( 2*sizeof(int));
			fscanf(fp, "%d%d%d%d\n", &curve[0], &curve[1], &range[0], &range[1]);
			xcoord = (double *)malloc( range[0]*range[1]*sizeof(double));
			ycoord = (double *)malloc( range[0]*range[1]*sizeof(double));
			zcoord = (double *)malloc( range[0]*range[1]*sizeof(double));
			for(i = 0; i < range[0]*range[1]; i++){
				fscanf(fp,"%lf%lf%lf", &xcoord[i], &ycoord[i], &zcoord[i]);
				printf("%lf\t%lf\t%lf\n", xcoord[i], ycoord[i], zcoord[i]);
			}
			break;
		default:
			printf("Not Implemented\n");
			break;
	}
	
	/* Close the input file */
	fclose(fp);
	
	/* Creat a New CGNS file */
	open_cgns (name, 2);
	/* Write in CGNS file */
	write_cgns();
	if(cg_close(cgnsfn))
		FATAL(NULL, NULL);
	
	/* Releasing the allocated memory */
	if(xcoord != NULL)
		free(xcoord);
	if(ycoord != NULL)
		free(ycoord);
	if(zcoord != NULL)
		free(zcoord);

	return 0;
}

