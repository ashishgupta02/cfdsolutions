/*******************************************************************************
 * File:        Output_Level.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _OUTPUT_LEVEL_H
#define _OUTPUT_LEVEL_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
* substitute replaced functions by macros for backward compability
* till they are renamed everywhere in the code ....
*******************************************************************************/
#include "Utils.h"

#define ndm_message ndm_print2
#define screen_output(a) ndm_print_output(a)
#define debug_level1(a)  have_debug_level(a)

/*******************************************************************************
*define output-level for stdout/stderr control
*******************************************************************************/

/*******************************************************************************
* Default output level
*******************************************************************************/
#define DEFAULT_OUTPUT_LEVEL              5
#define DEFAULT_PRINTF_LEVEL              5

/*******************************************************************************
* Ouptut levels for specific information
*******************************************************************************/
/* PRINT_ALWAYS for error-messages only, to allow switch off all other output */
#define PRINT_ALWAYS                      0
#define PRINT_MINIMUM_LEVEL               1

#define PRINT_DEFAULT_LEVEL               5

#define PRINT_MONITORING_LEVEL            2
#define PRINT_BOUNDARY_MAPPING_LEVEL      3
#define PRINT_PARAMETER_LEVEL             3
#define PRINT_CHECK_GRID_LOW_LEVEL        3
#define PRINT_BOUNDARY_TREATMENT_LEVEL    5

#define PRINT_GRID_CONTENTS_LEVEL         6
#define PRINT_TRANSITION_MAPPING_LEVEL    6

#define PRINT_CHECK_GRID_HIGH_LEVEL       7
#define PRINT_DUALGRID_PREP_PARA          7
#define PRINT_ADAPTATION_LEVEL            7

#define PRINT_PARTITIONING_LEVEL         10
#define PRINT_BANDWIDTH_LEVEL            10
#define PRINT_BUILDLAYERS_LEVEL          10

#define PRINT_GRID_STATISTICS_LEVEL      15
#define PRINT_TRANSITION                 15
#define PRINT_PERIODIC_BDRY_INFO_LEVEL   15
#define PRINT_ENGINE_BDRY_INFO_LEVEL     15
#define PRINT_HEAT_EXCHANGER_BDRY_INFO_LEVEL     15
#define PRINT_ACTUATOR_BDRY_INFO_LEVEL   15
#define PRINT_ADAP_LAYER_INFO_LEVEL      15
#define PRINT_ADAP_YPLUS_INFO_LEVEL      15

#define PRINT_FACE_TYPE_LEVEL            20
#define PRINT_NEAR_POINT_LEVEL           20
#define PRINT_WALL_DIST_LEVEL            20
#define PRINT_COLORING_LEVEL             20
#define PRINT_ADAPTATION_INFO_LEVEL      20

#define PRINT_FNGB_LEVEL                 20

#define PRINT_SUBGRIDS_LEVEL             25
#define PRINT_CHECK_SUBGRIDS_LEVEL       25

#define PRINT_FURTHER_ADAP_INFO_LEVEL    30

#define PRINT_PARALLEL_EXTENSIONS        40

#define PRINT_TIMING_OUTPUT              50

#define PRINT_MAXIMUM_LEVEL              50
#define DEBUG_LEVEL_1                    51

#define PRINT_FURTHER_GRID_STATISTICS    60


/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _OUTPUT_LEVEL_H */

