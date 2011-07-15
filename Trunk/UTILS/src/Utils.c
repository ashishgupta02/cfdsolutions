/*******************************************************************************
 * File:        Utils.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#define PATH_DELIM ';'
#include <direct.h>
#include <io.h>
#else
#define PATH_DELIM ':'
#include <unistd.h>
#endif

/* Application Header Files */
#include "Utils.h"
