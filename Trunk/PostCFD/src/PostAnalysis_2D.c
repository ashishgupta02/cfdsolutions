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
 * File		PostAnalysis_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "PPInitializeSolution_2D.h"
#include "PostAnalysis_2D.h"
#include "PPScalarFlowVariables_2D.h"
#include "PPVectorFlowVariables_2D.h"
#include "PPTensorFlowVariables_2D.h"

/* Post-Processing options */

/* static char options[] = "123"; */

static char *PPOptions[] = {
	"Post-Processing [options]",
	"options:",
	"	1	= Scalar Variables",
	"	2	= Vector Variables",
	"	3	= Tensor Variables",
	NULL
};

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static int usage (char **usgmsg) {
	int n;

	for (n = 0; NULL != usgmsg[n]; n++)
		fprintf (stderr, "%s\n", usgmsg[n]);

	return 0;
}

/*---------------------------------------------------------------*/
int PostAnalysis_2D (void) {
	int Option;
loop1:		
	/* Getting Post-Processing Option */
	usage(PPOptions);
	printf("Give Post-Processing Module Option: ");
	scanf("%d", &Option);
	if(Option < 1 || Option > 3) goto loop1;

	GetSolutionList_2D();
	
	switch (Option) {
		case 1:
			/* Calculate Scalar Variables */
			ScalarFlowVariableOptions_2D();
			break;
		case 2:
			/* Calculate Vector Variables */
			//VectorFlowVariableOptions_2D();
			break;
		case 3:
			/* Calculate Tensor Variables */
			TensorFlowVariableOptions_2D();
			break;
		default:
			/* Invalid Option */
			FATAL(NULL,"Invalid Post-Processing Option");
	}
	
	printf("Perform Another Post-Processing Operation (1/0):");
	scanf("%d", &Option);
	
	if (Option == 1)
		goto loop1;
	
	return 0;
}
