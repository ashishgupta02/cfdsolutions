/* 
 * File:   Bowyer_Watson.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#ifndef _BOWYER_WATSON_H
#define	_BOWYER_WATSON_H

    #include "List.h"
    
    int Bowyer_Watson_Triangulation(int nn, int tdim, double x[], double y[],
        int &ntri, int tri[][3], int cell2cell[][3], List* node2cell[]);

#endif	/* _BOWYER_WATSON_H */

