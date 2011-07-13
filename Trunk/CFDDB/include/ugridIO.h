/*******************************************************************************
 * File:        ugridIO.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _UGRIDIO_H
#define	_UGRIDIO_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
extern "C" {
#endif

    #include "corestruct.h"

    /* Define Function Libraries */
    /* Initialise the UGRID Data Structure */
    void UGRIDIO_INIT(void);
    /* Finalize the UGRID Data Structure */
    void UGRIDIO_FINI(void);
    /*************************************************/
    /* UGRID file I/O Methods */
    int open_ugrid(const char *ugridfile, int mode);
    int close_ugrid();
    BASE* read_ugrid_grid(void);
    BASE* read_ugrid_solution(void);
    int write_ugrid_grid(BASE *base);
    int write_ugrid_solution(BASE *base);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
}
#endif

#endif	/* _UGRIDIO_H */

