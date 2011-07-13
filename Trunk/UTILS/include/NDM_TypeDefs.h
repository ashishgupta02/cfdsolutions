/*******************************************************************************
 * File:        NDM_TypeDefs.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _NDM_TYPEDEFS_H
#define	_NDM_TYPEDEFS_H

#ifdef	__cplusplus
extern "C" {
#endif

/*******************************************************************************
* we define own types in order to allow to switch from int to short
* or from double to float if needed or if it is of advantage
* depending on the architecture used for compilation
*
* NOTE: NDMIndex is used for index-list containing point-numbers only!
*       face-numbers or element-numbers can be up to one magnitued
*       larger, thus we restrict 'NDMIndex' to point-indicees as
*       in tetrapnt, hexapnt, etc...
*       This will allow to reduce memory for large index-lists
*       without the risk to exceed the range e.g. of short for
*       large grids
*
* The macro NDM_USE_DOUBLE has to be defined if NDMDouble is set to double
* and should be undefined for float
*******************************************************************************/
#define NDMIndex  int
#define NDMDouble double

#define NDM_USE_DOUBLE

#ifdef	__cplusplus
}
#endif

#endif	/* _NDM_TYPEDEFS_H */

