/*******************************************************************************
 * File:        Grid_Utils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _GRID_UTILS_H
#define	_GRID_UTILS_H

    #include "List.h"

    int Point_In_Circle(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

    int Point_In_Triangle_Mesh(int seed, double xt, double yt, double x[], double y[], int tri[][3], int nbr[][3]);

    void Get_Domain_Extent(int nn, double x[], double y[],
        double *Xmin, double *Ymin, double *Xmax, double *Ymax);
    
    void Create_BoundingBox_Domain(double x[], double y[],
        double Xmin, double Ymin, double Xmax, double Ymax);
    
    void Delete_BoundingBox_Domain(int ntri, int tri[][3]);

    void Create_Node2Cell_Connectivity(int nn, int ntri, int tri[][3], List* node2cell[]);

    void Create_Cell2Cell_Connectivity(int ntri, int tri[][3], int cell2cell[][3], List* node2cell[]);

    void Update_Cell2Cell_Connectivity(int tri[][3], int cell2cell[][3],
        List* node2cell[], List* new_tri);

    void Find_Neighbour_Triangles_Violating_CircleTest(int tid, int pid, double x[], double y[],
        int tri[][3], int (*cell2cell)[3], List *violated_tri);

    int Boundary_Edge_Recovery(int nn, int nb, int nbs[], int ***bs,
        double x[], double y[], int ntri, int tri[][3],
        int cell2cell[][3], List* node2cell[]);

    int Flood_Fill_Algorithm(int nb, int nbs[], int ***bs,
        int &ntri, int tri[][3], int cell2cell[][3], List* node2cell[]);

#endif	/* _GRID_UTILS_H */

