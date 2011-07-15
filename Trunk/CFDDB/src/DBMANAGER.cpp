/*******************************************************************************
 * File:        DBMANAGER.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <iostream>
#include "DBMANAGER.h"

//------------------------------------------------------------------------------
//! Define the static DBMANAGER Pointer
//------------------------------------------------------------------------------
DBMANAGER* DBMANAGER::inst_ = NULL;

//------------------------------------------------------------------------------
//! Private Constructor
//------------------------------------------------------------------------------
DBMANAGER::DBMANAGER() {
    // Initialize the Attributes
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
DBMANAGER::~DBMANAGER() {
    if (inst_ != NULL) {
        delete(inst_);
    }
    // Free the Memory Used by Attributes
}

//------------------------------------------------------------------------------
//! Get the Only Single Instance of DBMANAGER
//------------------------------------------------------------------------------
DBMANAGER* DBMANAGER::getInstance() {
    if (inst_ == NULL) {
        inst_ = new DBMANAGER();
    }
    return(inst_);
}

