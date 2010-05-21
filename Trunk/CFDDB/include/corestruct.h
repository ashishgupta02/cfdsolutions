/*******************************************************************************
 * File:        corestruct.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _CORESTRUCT_H
#define	_CORESTRUCT_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
extern "C" {
#endif
    
    /* Structure Declarations */
    /* Structure for Descriptor */
    typedef struct _DESC {
        int id;
        char name[33];
        char *desc;
    } DESC;
    
    /* Structure for Vertex */
    typedef struct _VERTEX {
        int id;
        double x, y, z, w;
    } VERTEX;
    
    /* Structure for Element Set */
    typedef struct _ELEMSET {
        int id;
        char name[33];
        int type;
        int start;
        int end;
        int nbndry;
        int csize;
        int *conn;
        int psize;
        int *parent;
    } ELEMSET;
    
    /* Structure for Interface */
    typedef struct _INTERFACE {
        int id;
        char name[33];
        int range[3][2];
        char d_name[33];
        int d_range[3][2];
        int transform[3];
        int d_zone;
    } INTERFACE;
    
    /* Structure for Zone Grid Connectivity */
    typedef struct _CONNECT {
        int id;
        char name[33];
        int type;
        int location;
        int ptype;
        int npnts;
        int *pnts;
        char d_name[33];
        int d_ztype;
        int d_ptype;
        int d_npnts;
        int *d_pnts;
        int d_zone;
    } CONNECT;
    
    /* Structure for Boundary Condition */
    typedef struct _BOCO {
        int id;
        char name[33];
        int type;
        int ptype;
        int npnts;
        int *pnts;
        int n_index[3];
        int n_cnt;
        int n_type;
        double *n_list;
    } BOCO;
    
    /* Structure for Field */
    typedef struct _FIELD {
        int id;
        char name[33];
        int datatype;
        int units[5];
        int dataclass;
        int convtype;
        double dataconv[2];
        int exptype;
        double exponent[5];
        double *data;
    } FIELD;
    
    /* Structure for Solution */
    typedef struct _SOLUTION {
        int id;
        char name[33];
        int location;
        int rind[3][2];
        int size;
        int units[5];
        int dataclass;
        int nflds;
        FIELD *flds;
        int ndesc;
        DESC *desc;
    } SOLUTION;
    
    /* Structure for Zone */
    typedef struct _ZONE {
        int id;
        char name[33];
        int type;
        int idim;
        int dim[3];
        int units[5];
        int dataclass;
        int datatype;
        int vertflags;
        int nverts;
        VERTEX *verts;
        int nesets;
        ELEMSET *esets;
        int nints;
        INTERFACE *ints;
        int nconns;
        CONNECT *conns;
        int nbocos;
        BOCO *bocos;
        int nsols;
        SOLUTION *sols;
        int ndesc;
        DESC *desc;
    } ZONE;
    
    /* Struture for Base */
    typedef struct _BASE {
        int id;
        int celldim;
        int phydim;
        int baseclass;
        int baseunits[5];
        char name[33];
        int nzones;
        ZONE *zones;
        int ndesc;
        DESC *desc;
    } BASE;
    
    /* Struture for CFDDB */
    typedef struct _ROOT {
        int nbases;
        BASE *bases;
        int ndesc;
        DESC *desc;
    } ROOT;
    
/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef	__cplusplus
}
#endif

#endif	/* _CORESTRUCT_H */

