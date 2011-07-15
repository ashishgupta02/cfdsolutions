/*******************************************************************************
 * File:        DBERROR.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _DBERROR_H
#define	_DBERROR_H

class DBError {
public:
    static DBError* getInstance();
    // Member Functions to Register Errors
    
    // Distructor
    ~DBError();
protected:
    // Attributes to Register Error Mesages
    
private:
    // The one, Single Instance
    static DBError* inst_;
    // Private Constructors
    DBError();
};

#endif	/* _DBERROR_H */

