/*******************************************************************************
 * File:        Check_Scanf.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _CHECK_SCANF_H
#define _CHECK_SCANF_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
* an alternative to the standard scanf library function:
*  -  this function use fgets to read the line and vsscanf to convert it
*  -  the given input line has to match the format list
*  -  the routine does not stop until the requred input is given
*******************************************************************************/
int check_scanf(char *fmt, ...);


/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _CHECK_SCANF_H */
