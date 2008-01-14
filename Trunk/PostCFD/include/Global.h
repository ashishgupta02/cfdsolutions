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
 * File		Global.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

#ifndef __GLOBAL_H__
#define __GLOBAL_H__
/* This contains the global header information */

#undef WIN32
#undef LINUX
#undef MacOS
#undef UNIX

#define WIN32
//#define LINUX
//#define MacOS
//#define UNIX

#ifdef WIN32
//#include <windows.h>
#else
#include <unistd.h>
#endif

#ifdef MacOS
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <ctype.h>
#include <signal.h>
#include <assert.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#endif
