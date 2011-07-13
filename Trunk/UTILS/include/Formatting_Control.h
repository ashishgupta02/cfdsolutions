/*******************************************************************************
 * File:        Formatting_Control.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _FORMATTING_CONTROL_H
#define _FORMATTING_CONTROL_H

#include "NDM_TypeDefs.h"

/*******************************************************************************
* Compute the word length required for displaying all signiificant figures
* of a number
*******************************************************************************/
void determine_significant_digits(NDMDouble val);

/*******************************************************************************
*
*******************************************************************************/
int *return_significant_digits_for_time(void);

/*******************************************************************************
*
*******************************************************************************/
void ndm_sprintf(char *c, const char *format, ...);

#endif /** _FORMATTING_CONTROL_H **/
