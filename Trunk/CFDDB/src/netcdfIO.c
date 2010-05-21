/*******************************************************************************
 * File:        netcdfIO.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

/* Custom header files */
#include "Utils.h"
#include "netcdfIO.h"

/* Defination of Global Variable */
DIMENSION gDim;
ATTRIBUTE gAttrib;
VARIABLE gVar;
int netcdffn = 0;

/* Local Variables */
static int errstatus = NC_NOERR;
static char netcdftemp[16] = "";

/* Static Dictionary for TAU File */
/* Global Dimension Name Dictionary */
static char *global_dimension_name[] = {
    "no_of_points",
    "no_of_elements",
    "no_of_tetraeders",
    "points_per_tetraeder",
    "no_of_pyramids",
    "points_per_pyramid",
    "no_of_prisms",
    "points_per_prism",
    "no_of_hexaeders",
    "points_per_hexaeder",
    "no_of_surfacetriangles",
    "points_per_surfacetriangle",
    "no_of_surfacequadrilaterials",
    "points_per_surfacequadrilateral",
    "no_of_surfaceelements",
    "no_of_markers",
    "no_of_hanging_faces",
    "points_per_hanging_face",
    "no_of_child_faces",
    "child_face_data",
    NULL
};

/* Global Variables Name Dictionary */
static char *global_variable_name[] = {
    "points_of_tetraeders",
    "points_of_pyramids",
    "points_of_prisms",
    "points_of_hexaeders",
    "points_of_surfacetriangles",
    "points_of_surfacequadrilaterals",
    "boundarymarker_of_surfaces",
    "boundarypanel_of_surfaces",
    "points_xc",
    "points_yc",
    "points_zc",
    "marker",
    "centaur_key",
    NULL
};

/* Global Attributes Name Dictionary */
static char *global_attribute_name[] = {
    "type",
    "marker_xx",
    "deformation_timestep",
    "deformation_alpha",
    "deformation_dalpha",
    "deformation_h",
    "deformation_dh",
    "deformation_step",
    "adoptation_level",
    "connectivity_check",
    "volume_check",
    "point_check",
    "pyra_quadrilat_check",
    "pyra_vertical_check",
    "test_lib_check",
    "periodic_pairs_check",
    "is_balanced",
    "preserve_piles",
    "nownpoints",
    "naddpoints",
    "nextpoints",
    "nparpoints",
    "output_type",
    "parameters",
    "first_residual",
    NULL
};

#define NofDimensionsNames		20
#define NofVariableNames		13
#define NofAttributeNames		25

static int *valid_dimsp;
static int *valid_attrp;
static int *valid_varp;

/*---------- NETCDFIO_FATAL --------------------------------------------
 * Exit with error message
 *---------------------------------------------------------------------*/

void NETCDFIO_FATAL(char *procname, char *errmsg) {
    char *msg = errmsg;
    
    if (NULL == msg) {
        if (errstatus != NC_NOERR) {
            msg = (char *)nc_strerror(errstatus);
            if (NULL == msg || !*msg)
                msg = "unknown error";
        }
    }
    fflush(stdout);
    if (NULL == procname || !*procname)
        fprintf(stderr, "%s\n", msg);
    else
        fprintf(stderr, "%s:%s\n", procname, msg);
    if (netcdffn) nc_close(netcdffn);
    if (netcdftemp[0]) unlink_name(netcdftemp);
    exit(1);
}

/*---------- open_netcdf -----------------------------------------------
 * open a NETCDF file
 * 	mode:	1 = Read Only
 *		2 = Write
 *		3 = Create New Database if non exist
 *---------------------------------------------------------------------*/

int open_netcdf(const char *netcdffile, int mode) {
    
    switch (mode) {
        case 1:
            errstatus = nc_open(netcdffile, NC_NOWRITE, &netcdffn);
            break;
        case 2:
            errstatus = nc_open(netcdffile, NC_WRITE, &netcdffn);
            break;
        case 3:
            errstatus = nc_create(netcdffile, NC_NOCLOBBER, &netcdffn);
            break;
    }
    
    if (errstatus) NETCDFIO_FATAL("open_netcdf", NULL);
    
    return NC_NOERR;
}

/*---------- reset_tau_dictionary -------------------------------------
 * Reset TAU Grid/Solution File Dictionary Map
 *---------------------------------------------------------------------*/

static void reset_tau_dictionary(void) {
    if (valid_dimsp != NULL)
        free(valid_dimsp);
    if (valid_attrp != NULL)
        free(valid_attrp);
    if (valid_varp != NULL)
        free(valid_varp);
}

/*---------- init_tau_dictionary ---------------------------------------
 * Build TAU Grid/Solution File Dictionary Map
 *---------------------------------------------------------------------*/

static int init_tau_dictionary(void) {
    int i, j;
    int ndims, nvars, ngatts, unlimdimid;
    size_t length;
    char recname [NC_MAX_NAME+1];
    
    /* Dictionary Initialise */
    valid_dimsp = NULL;
    valid_attrp = NULL;
    valid_varp  = NULL;
    
    /* Initialise */
    ndims = 0;
    nvars = 0;
    ngatts = 0;
    unlimdimid = 0;
    
    /* Inquire NetCDF file for containts of the file */
    errstatus = nc_inq(netcdffn, &ndims, &nvars, &ngatts, &unlimdimid);
    if (errstatus)
        NETCDFIO_FATAL("init_tau_dictionary", NULL);
    
    /* Check of possible errors */
    if (ndims > NofDimensionsNames)
        NETCDFIO_FATAL("init_tau_dictionary", "Error 1");
    if (nvars > NofVariableNames)
        NETCDFIO_FATAL("init_tau_dictionary", "Error 2");
    if (ngatts > NofAttributeNames)
        warn("init_tau_dictionary %s", "More attributes are found then documented");
    
    /* Allocate the Valid Array lookup */
    valid_dimsp = malloc(NofDimensionsNames*sizeof(int));
    valid_varp = malloc(NofVariableNames*sizeof(int));
    valid_attrp = malloc(NofAttributeNames*sizeof(int));
    
    /* Check for malloc errors */
    if (valid_dimsp == NULL || valid_varp == NULL || valid_attrp == NULL)
        NETCDFIO_FATAL("init_tau_dictionary", "Error 3");
    
    /* Initialise the new arrays to zero */
    for (i = 0; i < NofDimensionsNames; i++)
        valid_dimsp[i] = 0;
    
    for (i = 0; i < NofVariableNames; i++)
        valid_varp[i] = 0;
    
    for (i = 0; i < NofAttributeNames; i++)
        valid_attrp[i] = 0;
    
    /* Get the valid Dimensions Names */
    for (i = 0; i < ndims; i++) {
        errstatus = nc_inq_dim(netcdffn, i, recname, &length);
        if (errstatus)
            NETCDFIO_FATAL("init_tau_dictionary", NULL);
        
        for (j = 0; j < NofDimensionsNames; j++) {
            errstatus = strcmp(global_dimension_name[j], recname);
            if (errstatus == 0) {
                valid_dimsp[j] = 1;
                break;
            }
        }
    }
    
    /* Get the valid Variables Names */
    for (i = 0; i < nvars; i++) {
        errstatus = nc_inq_varname(netcdffn, i, recname);
        if (errstatus)
            NETCDFIO_FATAL("init_tau_dictionary", NULL);
        
        for (j = 0; j < NofVariableNames; j++) {
            errstatus = strcmp(global_variable_name[j], recname);
            if (errstatus == 0) {
                valid_varp[j] = 1;
                break;
            }
        }
    }
    
    /* Get the valid Global Attributes Names */
    for (i = 0; i < ngatts; i++) {
        errstatus = nc_inq_attname(netcdffn, NC_GLOBAL, i, recname);
        if (errstatus)
            NETCDFIO_FATAL("init_tau_dictionary", NULL);
        
        for (j = 0; j < NofAttributeNames; j++) {
            errstatus = strcmp(global_attribute_name[j], recname);
            if (errstatus == 0) {
                valid_attrp[j] = 1;
                break;
            }
        }
    }
    
    return NC_NOERR;
}

/*---------- inquire_var ----------------------------------------------
 * Checks if Variable Exists in NetCDF Database
 *---------------------------------------------------------------------*/

int inquire_var(char *name, int *status) {
    int i, out;
    
    *status = 0;
    if (valid_varp == NULL)
        return 1;
    
    out = 1;
    for (i = 0; i < NofVariableNames; i++) {
        if (valid_varp[i]) {
            out = strcmp(name, global_variable_name[i]);
            if (out == 0) {
                *status = 1;
                break;
            }
        }
    }
    
    return NC_NOERR;
}

/*---------- inquire_dim ----------------------------------------------
 * Checks if Dimensions Exists in NetCDF Database
 *---------------------------------------------------------------------*/

int inquire_dim(char *name, int *status) {
    int i, out;
    
    *status = 0;
    if (valid_dimsp == NULL)
        return 1;
    
    out = 1;
    for (i = 0; i < NofDimensionsNames; i++) {
        if (valid_dimsp[i]) {
            out = strcmp(name, global_dimension_name[i]);
            if (out == 0) {
                *status = 1;
                break;
            }
        }
    }
    
    return NC_NOERR;
}

/*---------- read_netcdf ----------------------------------------------
 * Read the NETCDF file
 *---------------------------------------------------------------------*/

int read_netcdf(void) {
    /* Initialise TAU Dictionary */
    errstatus = init_tau_dictionary();
    if (errstatus)
        NETCDFIO_FATAL("read_netcdf", NULL);
    /* Get the Global Dimensions */
    errstatus = read_global_dimension();
    if (errstatus)
        NETCDFIO_FATAL("read_netcdf", NULL);
    /* Get the Global Attributes */
    errstatus = read_global_attribute();
    if (errstatus)
        NETCDFIO_FATAL("read_netcdf", NULL);
    /* Get the Global Variables */
    errstatus = read_global_variable();
    if (errstatus)
        NETCDFIO_FATAL("read_netcdf", NULL);
    return NC_NOERR;
}

/*---------- delete_global_dimension -----------------------------------
 * Delete global Dimension information from memory
 *---------------------------------------------------------------------*/

static void delete_global_dimension(void) {
    int i;
    
    if ((gDim.ndim > 0) && (gDim.pdim != NULL)) {
        /* Free the internal memory of Data Structure */
        for (i = 0; i < gDim.ndim; i++) {  
            if (gDim.pdim[i].name != NULL)
                free(gDim.pdim[i].name);
        }
        free(gDim.pdim);
        /* Reset to NULL */
        gDim.pdim = NULL;
        gDim.ndim = 0;
    }
}

/*---------- read_global_dimension -------------------------------------
 * Read global Dimension information from NETCDF file
 *---------------------------------------------------------------------*/

int read_global_dimension(void) {
    int i, ndims;
    size_t length;
    char recname [NC_MAX_NAME+1];
    
    /* Initialise */
    ndims = 0;
    
    /* Inquire for Number of Global Dimensions */
    errstatus = nc_inq_ndims(netcdffn, &ndims);
    if (errstatus)
        NETCDFIO_FATAL("read_global_dimension", NULL);
    
    /* Update in Date Structure */
    gDim.ndim = ndims;
    if (ndims == 0) {
        gDim.pdim = NULL;
        return NC_NOERR;
    }
    /* Allocate the memory to store the values */
    gDim.pdim = malloc(ndims*sizeof(VDIMENSION));
    if (gDim.pdim == NULL)
        NETCDFIO_FATAL("read_global_dimension", "Error 1");
    
    /* Get the value for all Dimensions */
    for (i = 0; i<ndims; i++) {
        str_blank(recname);
        errstatus = nc_inq_dim(netcdffn, i, recname, &length);
        if (errstatus)
            NETCDFIO_FATAL("read_global_dimension", NULL);
        gDim.pdim[i].value = length;
        gDim.pdim[i].name = malloc((NC_MAX_NAME+1)*sizeof(char));
        str_blank(gDim.pdim[i].name);
        strcpy(gDim.pdim[i].name, recname);
    }
    
    return NC_NOERR;
}

/*---------- delete_global_attribute -----------------------------------
 * Delete global Attribute information from memory
 *---------------------------------------------------------------------*/

static void delete_global_attribute(void) {
    int i;
    
    if ((gAttrib.nattrib > 0) && (gAttrib.pattrib != NULL)) {
        for (i = 0; i < gAttrib.nattrib; i++) {
            /* Free the memory for name */
            free(gAttrib.pattrib[i].name);
            /* Free memory used to store attributes values */
            switch (gAttrib.pattrib[i].atype) {
                case NC_CHAR:
                    if (gAttrib.pattrib[i].av.cp != NULL)
                        free(gAttrib.pattrib[i].av.cp);
                    break;
                case NC_SHORT:
                    if (gAttrib.pattrib[i].av.sp != NULL)
                        free(gAttrib.pattrib[i].av.sp);
                    break;
                case NC_INT:
                    if (gAttrib.pattrib[i].av.ip != NULL)
                        free(gAttrib.pattrib[i].av.ip);
                    break;
                case NC_FLOAT:
                    if (gAttrib.pattrib[i].av.fp != NULL)
                        free(gAttrib.pattrib[i].av.fp);
                    break;
                case NC_DOUBLE:
                    if (gAttrib.pattrib[i].av.dp != NULL)
                        free(gAttrib.pattrib[i].av.dp);
                    break;
            }
        }
        /* Finally free memory with pattrib */
        free(gAttrib.pattrib);
        gAttrib.pattrib = NULL;
        gAttrib.nattrib = 0;
    }
}

/*---------- read_global_attribute -------------------------------------
 * Read global Attribute information from NETCDF file
 *---------------------------------------------------------------------*/

int read_global_attribute(void) {
    int i, natts;
    nc_type atype;
    size_t length;
    char name[NC_MAX_NAME +1];
    
    /* Initialise */
    natts = 0;
    length = 0;
    atype = NC_NAT;
    str_blank(name);
    
    /* Inquire for Number of Global Attributes */
    errstatus = nc_inq_natts(netcdffn, &natts);
    if (errstatus)
        NETCDFIO_FATAL("read_global_attribute", NULL);
    /* Update in Date Structure */
    gAttrib.nattrib = natts;
    if (natts == 0) {
        gAttrib.pattrib = NULL;
        return NC_NOERR;
    }
    /* Allocate the memory to store the values */
    gAttrib.pattrib = malloc(natts*sizeof(VATTRIBUTE));
    if (gAttrib.pattrib == NULL)
        NETCDFIO_FATAL("read_global_attribute", "Error 1");
    /* Get the value for all Attributes */
    for (i = 0; i < natts; i++) {
        errstatus = nc_inq_attname(netcdffn, NC_GLOBAL, i, name);
        if (errstatus)
            NETCDFIO_FATAL("read_global_attribute", NULL);
        
        /* Get the type of the attribute */
        errstatus = nc_inq_atttype(netcdffn, NC_GLOBAL, name, &atype);
        if (errstatus)
            NETCDFIO_FATAL("read_global_attribute", NULL);
        gAttrib.pattrib[i].atype = atype;
        
        /* Get the length of the Attribute */
        errstatus = nc_inq_attlen(netcdffn, NC_GLOBAL, name, &length);
        if (errstatus)
            NETCDFIO_FATAL("read_global_attribute", NULL);
        gAttrib.pattrib[i].length = length;
        
        /* Get the value of the Attribute */
        switch (atype) {
            case NC_CHAR:
                /* Allocate the Space */
                gAttrib.pattrib[i].av.cp = malloc((length+1)*sizeof(char));
                errstatus = nc_get_att_text(netcdffn, NC_GLOBAL, name,
                        gAttrib.pattrib[i].av.cp);
                if (gAttrib.pattrib[i].av.cp != NULL)
                    gAttrib.pattrib[i].av.cp[length] = '\0';
                break;
            case NC_SHORT:
                gAttrib.pattrib[i].av.sp = malloc(length*sizeof(short));
                errstatus = nc_get_att_short(netcdffn, NC_GLOBAL, name,
                        gAttrib.pattrib[i].av.sp);
                break;
            case NC_INT:
                gAttrib.pattrib[i].av.ip = malloc(length*sizeof(int));
                errstatus = nc_get_att_int(netcdffn, NC_GLOBAL, name,
                        gAttrib.pattrib[i].av.ip);
                break;
            case NC_FLOAT:
                gAttrib.pattrib[i].av.fp = malloc(length*sizeof(float));
                errstatus = nc_get_att_float(netcdffn, NC_GLOBAL, name,
                        gAttrib.pattrib[i].av.fp);
                break;
            case NC_DOUBLE:
                gAttrib.pattrib[i].av.dp = malloc(length*sizeof(double));
                errstatus = nc_get_att_double(netcdffn, NC_GLOBAL, name,
                        gAttrib.pattrib[i].av.dp);
                break;
            default:
                NETCDFIO_FATAL("read_global_attribute", "Error 2");
        }
        
        if (errstatus)
            NETCDFIO_FATAL("read_global_attribute", NULL);
        
        /* Copy the name */
        length = strlen(name);
        gAttrib.pattrib[i].name = malloc((length+1)*sizeof(char));
        str_blank(gAttrib.pattrib[i].name);
        strcpy(gAttrib.pattrib[i].name, name);
        
        /* Reset the Value for next use */
        length = 0;
        atype = NC_NAT;
        str_blank(name);
    }
    
    return NC_NOERR;
}

/*---------- delete_global_variable -----------------------------------
 * Delete global Variable information from memory
 *---------------------------------------------------------------------*/

static void delete_global_variable(void) {
    int i, j;
    
    /* If there are any Variable */
    if ((gVar.nvar > 0) && (gVar.pvar != NULL)) {
        for (i = 0; i < gVar.nvar; i++) {
            /* Free memory used to store name */
            if (gVar.pvar[i].name != NULL)
                free(gVar.pvar[i].name);
            /* Check State and Free memory */
            if (gVar.pvar[i].state) {
                switch (gVar.pvar[i].vtype) {
                    case NC_CHAR:
                        if (gVar.pvar[i].vv.cp != NULL)
                            free(gVar.pvar[i].vv.cp);
                        break;
                    case NC_SHORT:
                        if (gVar.pvar[i].vv.sp != NULL)
                            free(gVar.pvar[i].vv.sp);
                        break;
                    case NC_INT:
                        if (gVar.pvar[i].vv.ip != NULL)
                            free(gVar.pvar[i].vv.ip);
                        break;
                    case NC_FLOAT:
                        if (gVar.pvar[i].vv.fp != NULL)
                            free(gVar.pvar[i].vv.fp);
                        break;
                    case NC_DOUBLE:
                        if (gVar.pvar[i].vv.dp != NULL)
                            free(gVar.pvar[i].vv.dp);
                        break;
                }
            }
            
            /* Free the memory used by dimensions */
            if ((gVar.pvar[i].dim.ndim > 0) && (gVar.pvar[i].dim.pdim != NULL)) {
                for (j = 0; j < gVar.pvar[i].dim.ndim; j++)
                    if (gVar.pvar[i].dim.pdim[j].name != NULL)
                        free(gVar.pvar[i].dim.pdim[j].name);
                
                if (gVar.pvar[i].dim.pdim != NULL)
                    free(gVar.pvar[i].dim.pdim);
            }
            
            /* Free the memory used by attributes */
            if ((gVar.pvar[i].attrib.nattrib > 0) && (gVar.pvar[i].attrib.pattrib != NULL)) {
                for (j = 0; j < gVar.pvar[i].attrib.nattrib; j++) {
                    /* Free the memory for name */
                    if (gVar.pvar[i].attrib.pattrib[j].name != NULL)
                        free(gVar.pvar[i].attrib.pattrib[j].name);
                    /* Free memory used to store attributes values */
                    switch (gVar.pvar[i].attrib.pattrib[j].atype) {
                        case NC_CHAR:
                            if (gVar.pvar[i].attrib.pattrib[j].av.cp != NULL)
                                free(gVar.pvar[i].attrib.pattrib[j].av.cp);
                            break;
                        case NC_SHORT:
                            if (gVar.pvar[i].attrib.pattrib[j].av.sp != NULL)
                                free(gVar.pvar[i].attrib.pattrib[j].av.sp);
                            break;
                        case NC_INT:
                            if (gVar.pvar[i].attrib.pattrib[j].av.ip != NULL)
                                free(gVar.pvar[i].attrib.pattrib[j].av.ip);
                            break;
                        case NC_FLOAT:
                            if (gVar.pvar[i].attrib.pattrib[j].av.fp != NULL)
                                free(gVar.pvar[i].attrib.pattrib[j].av.fp);
                            break;
                        case NC_DOUBLE:
                            if (gVar.pvar[i].attrib.pattrib[j].av.dp != NULL)
                                free(gVar.pvar[i].attrib.pattrib[j].av.dp);
                            break;
                    }
                }
                /* Free memory with pattrib */
                if (gVar.pvar[i].attrib.pattrib != NULL)
                    free(gVar.pvar[i].attrib.pattrib);
            }
        }
        /* Free the root */
        if (gVar.pvar != NULL)
            free(gVar.pvar);
        gVar.pvar = NULL;
        gVar.nvar = 0;
    }
}

/*---------- get_global_variable --------------------------------------
 * Get global Variable from NETCDF file only one variable is maintained
 * to avoid memory overhead
 *---------------------------------------------------------------------*/

int get_global_variable(int varid) {
    int i, ndim, memsize;
    
    /* Check the validity of varid */
    if ((varid < 0) || (varid > (gVar.nvar-1)))
        NETCDFIO_FATAL("get_global_variable", "Error 1");
    
    /* Check if any variable is allocated and free */
    for (i = 0; i < gVar.nvar; i++) {
        if (gVar.pvar[i].state) {
            switch (gVar.pvar[i].vtype) {
                case NC_CHAR:
                    free(gVar.pvar[i].vv.cp);
                    gVar.pvar[i].vv.cp = NULL;
                    break;
                case NC_SHORT:
                    free(gVar.pvar[i].vv.sp);
                    gVar.pvar[i].vv.sp = NULL;
                    break;
                case NC_INT:
                    free(gVar.pvar[i].vv.ip);
                    gVar.pvar[i].vv.ip = NULL;
                    break;
                case NC_FLOAT:
                    free(gVar.pvar[i].vv.fp);
                    gVar.pvar[i].vv.fp = NULL;
                    break;
                case NC_DOUBLE:
                    free(gVar.pvar[i].vv.dp);
                    gVar.pvar[i].vv.dp = NULL;
                    break;
            }
            gVar.pvar[i].state = 0;
        }
    }
    
    /* Get the value of Number of Dimensions */
    ndim = gVar.pvar[varid].dim.ndim;
    memsize = 1;
    for (i = 0; i < ndim; i++) {
        memsize = gVar.pvar[varid].dim.pdim[i].value * memsize;
    }
    /* Now Get the requested variable */
    switch (gVar.pvar[varid].vtype) {
        case NC_CHAR:
            /* Allocate the memory for the variable */
            gVar.pvar[varid].vv.cp = malloc(memsize*sizeof(char));
            if (gVar.pvar[varid].vv.cp == NULL)
                NETCDFIO_FATAL("rget_global_variable", "Error 2");
            
            errstatus = nc_get_var_text(netcdffn, varid, gVar.pvar[varid].vv.cp);
            break;
        case NC_SHORT:
            /* Allocate the memory for the variable */
            gVar.pvar[varid].vv.sp = malloc(memsize*sizeof(short));
            if (gVar.pvar[varid].vv.sp == NULL)
                NETCDFIO_FATAL("rget_global_variable", "Error 3");
            
            errstatus = nc_get_var_short(netcdffn, varid, gVar.pvar[varid].vv.sp);
            break;
        case NC_INT:
            /* Allocate the memory for the variable */
            gVar.pvar[varid].vv.ip = malloc(memsize*sizeof(int));
            if (gVar.pvar[varid].vv.ip == NULL)
                NETCDFIO_FATAL("rget_global_variable", "Error 4");
            
            errstatus = nc_get_var_int(netcdffn, varid, gVar.pvar[varid].vv.ip);
            break;
        case NC_FLOAT:
            /* Allocate the memory for the variable */
            gVar.pvar[varid].vv.fp = malloc(memsize*sizeof(float));
            if (gVar.pvar[varid].vv.fp == NULL)
                NETCDFIO_FATAL("rget_global_variable", "Error 5");
            
            errstatus = nc_get_var_float(netcdffn, varid, gVar.pvar[varid].vv.fp);
            break;
        case NC_DOUBLE:
            /* Allocate the memory for the variable */
            gVar.pvar[varid].vv.dp = malloc(memsize*sizeof(double));
            if (gVar.pvar[varid].vv.dp == NULL)
                NETCDFIO_FATAL("rget_global_variable", "Error 6");
            
            errstatus = nc_get_var_double(netcdffn, varid, gVar.pvar[varid].vv.dp);
            break;
        default:
            NETCDFIO_FATAL("rget_global_variable", "Data Type Not Supported");
    }
    if (errstatus)
        NETCDFIO_FATAL("get_global_variable", NULL);
    
    gVar.pvar[varid].state = 1;
    
    return NC_NOERR;
}

/*---------- read_global_variable --------------------------------------
 * Read global Variable information from NETCDF file
 *---------------------------------------------------------------------*/

int read_global_variable(void) {
    int i, j, nvars, ndims, natts;
    int *dimids;
    size_t alength;
    nc_type vtype, atype;
    char name[NC_MAX_NAME+1];
    char *aname;
    
    /* Initialise */
    nvars = 0;
    ndims = 0;
    natts = 0;
    dimids = NULL;
    alength = 0;
    vtype = NC_NAT;
    atype = NC_NAT;
    aname = NULL;
    str_blank(name);
    
    /* Get the number of variables */
    errstatus = nc_inq_nvars(netcdffn, &nvars);
    if (errstatus)
        NETCDFIO_FATAL("read_global_variable", NULL);
    
    /* Allocate the memory to store the variable information */
    gVar.nvar = nvars;
    if (nvars == 0) {
        gVar.pvar = NULL;
        return NC_NOERR;
    }
    /* Allocate the memory to store the values */
    gVar.pvar = malloc(nvars*sizeof(VVARIABLE));
    if (gVar.pvar == NULL)
        NETCDFIO_FATAL("read_global_variable", "Error 1");
    
    /* Get the value for all Variables */
    for (i = 0; i < nvars; i++) {
        /* Get the name of the variable */
        errstatus = nc_inq_varname(netcdffn, i, name);
        if (errstatus)
            NETCDFIO_FATAL("read_global_variable", NULL);
        gVar.pvar[i].name = malloc((NC_MAX_NAME+1)*sizeof(char));
        str_blank(gVar.pvar[i].name);
        strcpy(gVar.pvar[i].name, name);
        /* Set the State of Variable to inactive */
        gVar.pvar[i].state = 0;
        
        /* Get the type of the variable */
        errstatus = nc_inq_vartype(netcdffn, i, &vtype);
        if (errstatus)
            NETCDFIO_FATAL("read_global_variable", NULL);
        gVar.pvar[i].vtype = vtype;
        
        /* Initialise the pointers */
        switch (vtype) {
            case NC_CHAR:
                gVar.pvar[i].vv.cp = NULL;
                break;
            case NC_SHORT:
                gVar.pvar[i].vv.sp = NULL;
                break;
            case NC_INT:
                gVar.pvar[i].vv.ip = NULL;
                break;
            case NC_FLOAT:
                gVar.pvar[i].vv.fp = NULL;
                break;
            case NC_DOUBLE:
                gVar.pvar[i].vv.dp = NULL;
                break;
            default:
                NETCDFIO_FATAL("read_global_variable", "Error 2");
        }
        
        /* Get the number of dimensions associated with variable */
        errstatus = nc_inq_varndims(netcdffn, i, &ndims);
        if (errstatus)
            NETCDFIO_FATAL("read_global_variable", NULL);
        gVar.pvar[i].dim.ndim = ndims;
        gVar.pvar[i].dim.pdim = NULL;
        
        /* Get the Dimensions Ids of dimensions for perticular variable */
        if (ndims > 0) {
            dimids = malloc(ndims * sizeof(int));
            if (dimids == NULL)
                NETCDFIO_FATAL("read_global_variable", "Error 3");
            errstatus = nc_inq_vardimid(netcdffn, i, dimids);
            if (errstatus)
                NETCDFIO_FATAL("read_global_variable", NULL);
            /* Copy the information from Global Dimensions */
            gVar.pvar[i].dim.pdim = malloc(ndims*sizeof(VDIMENSION));
            for (j = 0; j < ndims; j++) {
                gVar.pvar[i].dim.pdim[j].value = gDim.pdim[dimids[j]].value;
                gVar.pvar[i].dim.pdim[j].name = malloc((NC_MAX_NAME+1)*sizeof(char));
                str_blank(gVar.pvar[i].dim.pdim[j].name);
                strcpy(gVar.pvar[i].dim.pdim[j].name, gDim.pdim[dimids[j]].name);
            }
        }
        
        /* Get the number of attributes associated with variable */
        errstatus = nc_inq_varnatts(netcdffn, i, &natts);
        if (errstatus)
            NETCDFIO_FATAL("read_global_variable", NULL);
        gVar.pvar[i].attrib.nattrib = natts;
        gVar.pvar[i].attrib.pattrib = NULL;
        
        /* Get the value of Attribute for each variable */
        if (natts > 0) {
            gVar.pvar[i].attrib.pattrib = malloc(natts*sizeof(VATTRIBUTE));
            
            /* Get the value for all Attributes */
            for (j = 0; j < natts; j++) {
                errstatus = nc_inq_attname(netcdffn, i, j, aname);
                if (errstatus)
                    NETCDFIO_FATAL("read_global_variable", NULL);
                
                /* Get the type of the attribute */
                errstatus = nc_inq_atttype(netcdffn, i, aname, &atype);
                if (errstatus)
                    NETCDFIO_FATAL("read_global_variable", NULL);
                gVar.pvar[i].attrib.pattrib[j].atype = atype;
                
                /* Get the length of the Attribute */
                errstatus = nc_inq_attlen(netcdffn, i, aname, &alength);
                if (errstatus)
                    NETCDFIO_FATAL("read_global_variable", NULL);
                gVar.pvar[i].attrib.pattrib[j].length = alength;
                
                /* Get the value of the Attribute */
                switch (atype) {
                    case NC_CHAR:
                        /* Allocate the Space */
                        gVar.pvar[i].attrib.pattrib[j].av.cp =
                                malloc((alength+1)*sizeof(char));
                        errstatus = nc_get_att_text(netcdffn, i, aname,
                                gVar.pvar[i].attrib.pattrib[j].av.cp);
                        if (gVar.pvar[i].attrib.pattrib[j].av.cp != NULL)
                            gVar.pvar[i].attrib.pattrib[j].av.cp[alength] = '\0';
                        break;
                    case NC_SHORT:
                        gVar.pvar[i].attrib.pattrib[j].av.sp =
                                malloc(alength*sizeof(short));
                        errstatus = nc_get_att_short(netcdffn, i, aname,
                                gVar.pvar[i].attrib.pattrib[j].av.sp);
                        break;
                    case NC_INT:
                        gVar.pvar[i].attrib.pattrib[j].av.ip =
                                malloc(alength*sizeof(int));
                        errstatus = nc_get_att_int(netcdffn, i, aname,
                                gVar.pvar[i].attrib.pattrib[j].av.ip);
                        break;
                    case NC_FLOAT:
                        gVar.pvar[i].attrib.pattrib[j].av.fp =
                                malloc(alength*sizeof(float));
                        errstatus = nc_get_att_float(netcdffn, i, aname,
                                gVar.pvar[i].attrib.pattrib[j].av.fp);
                        break;
                    case NC_DOUBLE:
                        gVar.pvar[i].attrib.pattrib[j].av.dp =
                                malloc(alength*sizeof(double));
                        errstatus = nc_get_att_double(netcdffn, i, aname,
                                gVar.pvar[i].attrib.pattrib[j].av.dp);
                        break;
                    default:
                        NETCDFIO_FATAL("read_global_variable", "Error 4");
                }
                
                if (errstatus)
                    NETCDFIO_FATAL("read_global_variable", NULL);
                
                /* Copy the name */
                alength = strlen(aname);
                gVar.pvar[i].attrib.pattrib[j].name = malloc((alength+1)*sizeof(char));
                str_blank(gVar.pvar[i].attrib.pattrib[j].name);
                strcpy(gVar.pvar[i].attrib.pattrib[j].name, aname);
                
                /* Reset the Value for next use */
                alength = 0;
                atype = NC_NAT;
                free(aname);
                aname = NULL;
            }
        }
        
        /* Reset the temporary variables */
        ndims = 0;
        natts = 0;
        free(dimids);
        vtype = NC_NAT;
        str_blank(name);
    }
    nvars = 0;
    
    return NC_NOERR;
}

/*---------- NETCDFIO_INIT ---------------------------------------------
 * Initialise the NetCDF Data Structure
 *---------------------------------------------------------------------*/

int NETCDFIO_INIT(void) {
    
    /* Initialise the global variables */
    gDim.ndim = 0;
    gDim.pdim = NULL;
    gAttrib.nattrib = 0;
    gAttrib.pattrib = NULL;
    gVar.nvar = 0;
    gVar.pvar = NULL;
    netcdffn = 0;
    str_blank(netcdftemp);
    errstatus = NC_NOERR;
    
    return NC_NOERR;
}

/*---------- NETCDFIO_RESET --------------------------------------------
 * Resets the NetCDF Data Structure safely
 *---------------------------------------------------------------------*/

int NETCDFIO_RESET(void) {
    
    reset_tau_dictionary();
    delete_global_dimension();
    delete_global_attribute();
    delete_global_variable();
    
    errstatus = nc_close(netcdffn);
    if (errstatus)
        NETCDFIO_FATAL("NETCDFIO_RESET", NULL);
    
    return NC_NOERR;
}

/*---------- NETCDFIO_INFO- --------------------------------------------
 * Print NetCDF File Containts Info
 *---------------------------------------------------------------------*/

int NETCDFIO_INFO(void) {
    int i;
    
    /* Print the Global Dimension Information */
    if (gDim.ndim > 0) {
        fprintf(stdout, "No of Global Dimensions = %d\n", gDim.ndim+1);
        for (i = 0; i < gDim.ndim; i++)
            fprintf(stdout, "\t%s = %d\n",
                    gDim.pdim[i].name, gDim.pdim[i].value);
    }
    
    /* Print the Global Attribute */
    if (gAttrib.nattrib > 0) {
        fprintf(stdout, "No of Global Attributes = %d \n", gAttrib.nattrib);
        for (i = 0; i < gAttrib.nattrib; i++) {
            switch (gAttrib.pattrib[i].atype) {
                case NC_CHAR:
                    fprintf(stdout, "\t%s = %s\n",
                            gAttrib.pattrib[i].name, gAttrib.pattrib[i].av.cp);
                    break;
                case NC_SHORT:
                    break;
                case NC_INT:
                    break;
                case NC_FLOAT:
                    break;
                case NC_DOUBLE:
                    break;
            }
        }
    }
    
    /* Print the Global Variable */
    if (gVar.nvar > 0) {
        fprintf(stdout, "No of Global Variables = %d \n", gVar.nvar);
        for (i = 0; i < gVar.nvar; i++)
            fprintf(stdout, "\t%s\n", gVar.pvar[i].name);
    }
    
    return NC_NOERR;
}

/*---------- write_global_dimension ------------------------------------
 * Write global Dimension information to NETCDF file
 *---------------------------------------------------------------------*/

void write_global_dimension(void) {
    int i, dimid;
    
    /* Check if the availablity */
    if ((gDim.ndim < 1) || (gDim.pdim == NULL)) return;
    
    /* Write the Dimensions */
    for (i = 0; i < gDim.ndim; i++) {
        dimid = 0;
        errstatus = nc_def_dim(netcdffn, gDim.pdim[i].name, gDim.pdim[i].value, &dimid);
        if (errstatus)
            NETCDFIO_FATAL("write_global_dimension", NULL);
    }
}

/*---------- write_global_attribute ------------------------------------
 * Write global Attribute information to NETCDF file
 *---------------------------------------------------------------------*/

void write_global_attribute(void) {
    int i, atype, attribid;
    
    /* Check if the availablity */
    if ((gAttrib.nattrib < 1) || (gAttrib.pattrib == NULL)) return;
    
    /* Write the Attributes */
    for (i = 0; i < gAttrib.nattrib; i++) {
        attribid = 0;
        atype = gAttrib.pattrib[i].atype;
        /* Get the value of the Attribute */
        switch (atype) {
            case NC_CHAR:
                errstatus = nc_put_att_text(netcdffn, NC_GLOBAL, gAttrib.pattrib[i].name,
                        gAttrib.pattrib[i].length, gAttrib.pattrib[i].av.cp);
                break;
            case NC_SHORT:
                errstatus = nc_put_att_short(netcdffn, NC_GLOBAL, gAttrib.pattrib[i].name,
                        atype, gAttrib.pattrib[i].length, gAttrib.pattrib[i].av.sp);
                break;
            case NC_INT:
                errstatus = nc_put_att_int(netcdffn, NC_GLOBAL, gAttrib.pattrib[i].name,
                        atype, gAttrib.pattrib[i].length, gAttrib.pattrib[i].av.ip);
                break;
            case NC_FLOAT:
                errstatus = nc_put_att_float(netcdffn, NC_GLOBAL, gAttrib.pattrib[i].name,
                        atype, gAttrib.pattrib[i].length, gAttrib.pattrib[i].av.fp);
                break;
            case NC_DOUBLE:
                errstatus = nc_put_att_double(netcdffn, NC_GLOBAL, gAttrib.pattrib[i].name,
                        atype, gAttrib.pattrib[i].length, gAttrib.pattrib[i].av.dp);
                break;
            default:
                NETCDFIO_FATAL("write_global_attribute", "Error 1");
        }
        
        if (errstatus)
            NETCDFIO_FATAL("write_global_attribute", NULL);
    }
}

/*---------- write_global_variable ------------------------------------
 * Write global Variable information to NETCDF file
 *---------------------------------------------------------------------*/

void write_global_variable(void) {
    int i, j, varid, vtype, id;
    int *dimids;
    
    /* Initialize */
    dimids = NULL;
    
    /* Check if the availablity */
    if ((gVar.nvar < 1) || (gVar.pvar == NULL)) return;
    
    for (i = 0; i < gVar.nvar; i++) {
        varid = 0;
        /* Get the dimensions IDs */
        if (dimids != NULL) free(dimids);
        dimids = (int *) malloc(gVar.pvar[i].dim.ndim*sizeof(int));
        for (j = 0; j < gVar.pvar[i].dim.ndim; j++) {
            id = 0;
            errstatus = nc_inq_dimid(netcdffn, gVar.pvar[i].dim.pdim[j].name, &id);
            if (errstatus)
                NETCDFIO_FATAL("write_global_variable", NULL);
            dimids[j] = id;
        }
        
        errstatus = nc_def_var(netcdffn, gVar.pvar[i].name, gVar.pvar[i].vtype,
                gVar.pvar[i].dim.ndim, dimids, &varid);
        if (errstatus)
            NETCDFIO_FATAL("write_global_variable", NULL);
        
        /* End the Define mode */
        errstatus = nc_enddef(netcdffn);
        if (errstatus)
            NETCDFIO_FATAL("write_global_variable", NULL);
        
        /* Now write the data into variable */
        vtype = gVar.pvar[i].vtype;
        switch (vtype) {
            case NC_CHAR:
                errstatus = nc_put_var_text(netcdffn, varid, gVar.pvar[i].vv.cp);
                break;
            case NC_SHORT:
                errstatus = nc_put_var_short(netcdffn, varid, gVar.pvar[i].vv.sp);
                break;
            case NC_INT:
                errstatus = nc_put_var_int(netcdffn, varid, gVar.pvar[i].vv.ip);
                break;
            case NC_FLOAT:
                errstatus = nc_put_var_float(netcdffn, varid, gVar.pvar[i].vv.fp);
                break;
            case NC_DOUBLE:
                errstatus = nc_put_var_double(netcdffn, varid, gVar.pvar[i].vv.dp);
                break;
            default:
                NETCDFIO_FATAL("read_global_variable", "Error 1");
        }
        if (errstatus)
            NETCDFIO_FATAL("write_global_variable", NULL);
        
        /* Start the define mode */
        errstatus = nc_redef(netcdffn);
        if (errstatus)
            NETCDFIO_FATAL("write_global_variable", NULL);
    }
    
    /* Free memory */
    if (dimids != NULL) free(dimids);
}

/*---------- write_netcdf ----------------------------------------------
 * Write the NETCDF file
 *---------------------------------------------------------------------*/

void write_netcdf(void) {
    /* Write Dimensions */
    write_global_dimension();
    /* Write Attributes */
    write_global_attribute();
    /* Write Variables */
    write_global_variable();
}
