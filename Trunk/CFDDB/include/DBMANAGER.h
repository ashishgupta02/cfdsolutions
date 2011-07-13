/*******************************************************************************
 * File:        DBMANAGER.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _DBMANAGER_H
#define	_DBMANAGER_H

class DBMANAGER {
public:
    static DBMANAGER* getInstance();
    // Member Functions to Register and Deregister Instances of DB Derived Classes
    ~DBMANAGER();
protected:
    // Attributes to Register Instances of DB Derived Classes
    
private:
    // The one, Single Instance
    static DBMANAGER* inst_;
    // Private Constructors
    DBMANAGER();
};

#endif	/* _DBMANAGER_H */

