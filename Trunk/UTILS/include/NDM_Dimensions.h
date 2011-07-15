/*******************************************************************************
 * File:        NDM_Dimensions.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _NDM_DIMENSIONS_H
#define _NDM_DIMENSIONS_H

/*==============================================================================
| NOTE: the following defines sorted in alphabetic order
==============================================================================*/
/*------------------------------------------------------------------------------
| Define the length of keynames, string parameter and internal buffers
------------------------------------------------------------------------------*/
#define NDM_BUF_SZ          128
#define NDM_MAX_KEYNAME_LEN 124
#define NDM_MAX_STR_PARAMETER_LEN NDM_MAX_KEYNAME_LEN

/*------------------------------------------------------------------------------
| Definition to avoid the use of FILENAME_MAX which is possibly too small
------------------------------------------------------------------------------*/
#define NDM_MAX_FILENAME_LEN 1024

/*------------------------------------------------------------------------------
| Maximum markers allowed for each boundary part
| Note: we encountered the 1st grid with more than 1024 marker
------------------------------------------------------------------------------*/
#define NDM_MAX_MARKERS 4096

/*------------------------------------------------------------------------------
| Define maximum length of marker names
------------------------------------------------------------------------------*/
#define NDM_MAX_MARKERNAME_LEN 50

/*------------------------------------------------------------------------------
| Maximum number of points per volume element
------------------------------------------------------------------------------*/
#define NDM_MAX_POINTS_PER_VOLELEMENT 8

/*------------------------------------------------------------------------------
| Maximum argument leng for communication (kaps)
------------------------------------------------------------------------------*/
#define NDM_MAX_ARG_LEN 255

#endif /* _NDM_DIMENSIONS_H */

