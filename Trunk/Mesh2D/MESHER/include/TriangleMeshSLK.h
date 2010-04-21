/* 
 * File:   TriangleMeshSLK.h
 * Author: Ashish Gupta
 *
 * Created on April 20, 2010, 5:24 PM
 */

#ifndef _TRIANGLEMESHSLK_H
#define	_TRIANGLEMESHSLK_H

#include "Point2D.h"
#include "Linked_List.h"

class TriangleMeshSLK {
public:
    TriangleMeshSLK();
    ~TriangleMeshSLK();
    
    int trimesh(int npt, int tdim, int nb, int *nbs, int ***bs, double x[], double y[], int tri[][3]);
    int trimesh(int npt, int tdim, int nb, int *nbs, int ***bs, Point2D pt[], int tri[][3]);
protected:
    void compress_list(Linked_List *del, Linked_List **Lhash, Linked_List *mknbr, int tri[][3], int nbr[][3], int &nt);
    void nbr_search(Point2D pt, Point2D p[], int tri[][3], int nbr[][3], Linked_List *del, int m);
    void make_nbrs(Linked_List *mknbr, int nn, int nt, int tri[][3], int nbr[][3], Linked_List **Lhash);
    int circle_test(Point2D t1, Point2D t2, Point2D t3, Point2D t);
    int search(int t, Point2D pt, Point2D p[], int tri[][3], int nbr[][3]);
};

#endif	/* _TRIANGLEMESHSLK_H */

