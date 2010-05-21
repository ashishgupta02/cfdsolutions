/*******************************************************************************
 * File:        corelib.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _CORELIB_H
#define	_CORELIB_H

#include "corestruct.h"

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
extern "C" {
#endif
    
    #define CORELIB_OK      0
    #define CORELIB_ERROR   1

    /* Create new root */
    ROOT *new_root(void);
    /* Create new base(s) */
    BASE *new_base(int count);
    /* Create new zone(s) */
    ZONE *new_zone(int count);
    /* Create coordinate array for a zone */
    VERTEX *new_vertex(int nverts);
    /* Create element set array */
    ELEMSET *new_elemset(int nsets);
    /* Create grid 1to1 interface array for a zone */
    INTERFACE *new_interface(int nints);
    /* Create grid connectivity array for a zone */
    CONNECT *new_connect(int nconns);
    /* Create boundary condition array */
    BOCO *new_boco(int nbocos);
    /* Create descriptor array */
    DESC *new_desc(int ndesc);
    /* Create solution array for a zone */
    SOLUTION *new_solution(int nsols);
    /* Create solution variable array for a zone */
    FIELD *new_field(int nfields, int size);

    /* Delete the memory associated in ROOT */
    int del_root(ROOT *POINTER);
    /* Delete the memory associated in BASE */
    int del_base(BASE *POINTER);
    /* Delete the memory associated in ZONE */
    int del_zone(ZONE *POINTER);
    /* Delete the memory associated in VERTEX */
    int del_vertex(VERTEX *POINTER);
    /* Delete the memory associated in ELEMSET */
    int del_elemset(ELEMSET *POINTER);
    /* Delete the memory associated in INTERFACE */
    int del_interface(INTERFACE *POINTER);
    /* Delete the memory associated in CONNECT */
    int del_connect(CONNECT *POINTER);
    /* Delete the memory associated in BOCO */
    int del_boco(BOCO *POINTER);
    /* Delete the memory associated in DESC */
    int del_desc(DESC *POINTER);
    /* Delete the memory associated in SOLUTION */
    int del_solution(SOLUTION *POINTER);
    /* Delete the memory associated in FIELD */
    int del_field(FIELD *POINTER);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
}
#endif

#endif	/* _CORELIB_H */

