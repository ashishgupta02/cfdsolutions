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
 * File		ColorCode_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "ColorCode_2D.h"

/* Color Encode */
static int ColorRange = 256;
double Red[256];
double Green[256];
double Blue[256];
static int colorflag = 0;

/*--------------------------------------------------------------*/
/* Color Varies From:
	Blue -> Blue + Green
	Blue + Green -> Green + Red
	Green + Red -> Red
*/
void ColorEncode_2D (void) {
	int i, tmp;
	
	if (colorflag == 0 ) {
		tmp = (ColorRange - 1)/3;
		/* Blue -> Blue + Green */
		for (i = 0; i <= tmp; i++) {
			Red[i]   = 0.0;
			Green[i] = ((double)i)/((double)tmp);
			Blue[i]  = 1.0;
		}
		
		/* Blue + Green -> Green + Red */
		for (i = 0; i <= tmp; i++) {
			Red[tmp+i]   = ((double)i)/((double)tmp);
			Green[tmp+i] = 1.0;
			Blue[tmp+i]  = ((double)(tmp-i))/((double)tmp);
		}
		
		/* Green + Red -> Red */
		for(i = 0; i <= tmp; i++) {
			Red[2*tmp+i]   = 1;
			Green[2*tmp+i] = ((double)(tmp-i))/((double)tmp);
			Blue[2*tmp+i]  = 0;
		}
		colorflag = 1;
	}
}
