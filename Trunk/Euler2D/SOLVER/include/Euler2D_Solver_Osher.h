/* 
 * File:   Euler2D_Solver_Osher.h
 * Author: Ashish Gupta
 *
 * Created on March 21, 2010, 11:04 PM
 */

#ifndef _EULER2D_SOLVER_OSHER_H
#define	_EULER2D_SOLVER_OSHER_H

#include "Euler2D_Mesh.h"

class Euler2D_Solver_Osher: virtual public Euler2D_Mesh {
protected:
    // Constants
    double Gamma;
    // Reference Conditions
    double Ref_Mach;
    double Ref_Alpha;
    double Ref_Pressure;
public:
    Euler2D_Solver_Osher();
    virtual ~Euler2D_Solver_Osher();
    void Solve();
protected:
    // Initialize the Solution
    void Initialize_Solution();
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_OSHER_H */

