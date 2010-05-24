/*******************************************************************************
 * File:        Machine_Parameters.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MACHINE_PARAMETERS_H
#define _MACHINE_PARAMETERS_H

#include "NDM_TypeDefs.h"

typedef struct machine_parameters {
    int ibeta;
    int it;
    int irnd;
    int ngrd;
    int machep;
    int negep;
    int iexp;
    int minexp;
    int maxexp;
    NDMDouble eps;
    NDMDouble epsneg;
    NDMDouble xmin;
    NDMDouble xmax;
} MachineParameters;

/*******************************************************************************
 *
 ******************************************************************************/
MachineParameters *return_machine_parameters(void);

/*******************************************************************************
 *
 ******************************************************************************/
void print_machine_parameters(void);

/*******************************************************************************
 *
 ******************************************************************************/
void show_machine_parameters(void);

#endif /** _MACHINE_PARAMETERS_H **/
