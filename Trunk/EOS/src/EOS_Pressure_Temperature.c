/*******************************************************************************
 * File:        EOS_Pressure_Temperature.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS_Internal.h"

// ------------------ Pressure and Temperature Formulation ---------------------
// This functions are valid Only in this Regions:
// 1) Super Critical State
// 2) Subcooled Compressed Liquid
// 3) Super Heated Vapor
// -----------------------------------------------------------------------------

