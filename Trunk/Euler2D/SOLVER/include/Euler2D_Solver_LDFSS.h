/* 
 * File:   Euler2D_Solver_LDFSS.h
 * Author: Ashish Gupta
 *
 * Created on March 21, 2010, 11:03 PM
 */

#ifndef _EULER2D_SOLVER_LDFSS_H
#define	_EULER2D_SOLVER_LDFSS_H

#include "Euler2D_Mesh.h"

class Euler2D_Solver_LDFSS: virtual public Euler2D_Mesh {
protected:
    // Constants
    double Gamma;
    // Reference Conditions
    double Ref_Mach;
    double Ref_Alpha;
    double Ref_Pressure;
public:
    Euler2D_Solver_LDFSS();
    virtual ~Euler2D_Solver_LDFSS();
    void Solve();
protected:
    // Initialize the Solution
    void Initialize_Solution();
private:
    void Init();
};

#endif	/* _EULER2D_SOLVER_LDFSS_H */

