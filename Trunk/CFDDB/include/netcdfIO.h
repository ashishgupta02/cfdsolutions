/*******************************************************************************
 * File:        netcdfIO.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _NETCDFIO_H
#define	_NETCDFIO_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
    
    /* Structure Declarations */
    
    typedef union _NCPOINTER {
        char *cp;
        short *sp;
        int *ip;
        float *fp;
        double *dp;
    } NCPOINTER;
    
    /* Structure to store value of Attribute */
    typedef struct _VATTRIBUTE {
        unsigned int atype;
        int length;
        char *name;
        NCPOINTER av;
    } VATTRIBUTE;
    
    /* Structure to manage multiple Attribute */
    typedef struct _ATTRIBUTE {
        int nattrib;
        VATTRIBUTE *pattrib;
    } ATTRIBUTE;
    
    /* Structure to store value of Dimension with Attributes */
    typedef struct _VDIMENSION {
        int value;
        char *name;
    } VDIMENSION;
    
    /* Structure to manage multiple Dimension */
    typedef struct _DIMENSION {
        int ndim;
        VDIMENSION *pdim;
    } DIMENSION;
    
    /* Structure to store variable of anytype with Dimensions and Attribute */
    typedef struct  _VVARIABLE {
        unsigned int vtype;
        unsigned int state;
        char *name;
        NCPOINTER vv;
        DIMENSION dim;
        ATTRIBUTE attrib;
    } VVARIABLE;
    
    /* Structure to manage multiple Variable */
    typedef struct _VARIABLE {
        int nvar;
        VVARIABLE *pvar;
    } VARIABLE;
    
    /* Global Variables */
    extern DIMENSION gDim;
    extern ATTRIBUTE gAttrib;
    extern VARIABLE gVar;
    extern int netcdffn;
    
    /* Define Function Libraries */
    /* Exit with error message */
    void NETCDFIO_FATAL(char *procname, char *errmsg);
    /* Initialise the NetCDF Data Structure */
    int NETCDFIO_INIT(void);
    /* Resets the NetCDF Data Structure safely */
    int NETCDFIO_RESET(void);
    /* Print NetCDF File Containts Info */
    int NETCDFIO_INFO(void);
    /* Open a NETCDF file */
    int open_netcdf(const char *netcdffile, int mode);
    
    /* Read the NETCDF file */
    int read_netcdf(void);
    /* Read global Dimension information from NETCDF file */
    int read_global_dimension(void);
    /* Read global Attribute information from NETCDF file */
    int read_global_attribute(void);
    /* Read global Variable information from NETCDF file */
    int read_global_variable(void);
    /* Get global Variable from NETCDF file
   only one variable is maintained to avoid memory overhead */
    int get_global_variable(int varid);
    
    /* Write the NETCDF file */
    void write_netcdf(void);
    /* Write global Dimension information to NETCDF file */
    void write_global_dimension(void);
    /* Write global Attribute information to NETCDF file */
    void write_global_attribute(void);
    /* Write global Variable information to NETCDF file */
    void write_global_variable(void);
    
    /* Dictionary Operations */
    int inquire_var(char *name, int *status);
    int inquire_dim(char *name, int *status);
    
/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef __cplusplus
}
#endif

#endif	/* _NETCDFIO_H */

