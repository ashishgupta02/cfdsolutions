/*******************************************************************************
 * File:        cgnsIO.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _CGNSIO_H
#define	_CGNSIO_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
    
    #include "corestruct.h"
    
    /* Define Function Libraries */
    /* Exit with error message */
    void CGNSIO_FATAL(char *procname, char *errmsg);
    /* Initialise the CGNS Data Structure */
    void CGNSIO_INIT(void);
    /* Finalize the CGNS Data Structure */
    void CGNSIO_FINI(void);
    /*************************************************/
    /* Open a CGNS file */
    int open_cgns(const char *cgnsfile, int mode);
    int get_cgns_nbases(void);
    BASE* get_cgns_base(int ibase);
    int set_cgns_base(BASE *base);
    int close_cgns();
    /*************************************************/
    /* Read the CGNS file */
    /* Read Base information from CGNS file */
    void read_base(int ibase);
    /* Read zone information from CGNS file */
    int read_zones(void);
    /* Read zone grid coordinates */
    int read_zone_grid(int izone);
    /* Read zone element sets and elements */
    int read_zone_element(int izone);
    /* Read zone 1 to 1 interfaces */
    int read_zone_interface(int izone);
    /* Read zone connectivities */
    int read_zone_connect(int izone);
    /* Read zone boundary conditions */
    int read_zone_boco(int izone);
    /* Read zone solution */
    int read_zone_solution(int izone);
    /* Read solution field data for a zone */
    int read_solution_field(int izone, int isol, int ifld);
    /* Read unit specifications */
    int read_units(int units[5]);
    /*************************************************/
    /* Write the CGNS file */
    /* Write Base information to CGNS file */
    void write_base(void);
    /* Write zone information to CGNS file */
    void write_zones(void);
    /* Write zone grid coordinates */
    void write_zone_grid(int izone);
    /* Write zone element sets and elements */
    void write_zone_element(int izone);
    /* Write zone 1 to 1 interfaces */
    void write_zone_interface(int izone);
    /* Write zone connectivities */
    void write_zone_connect(int izone);
    /* Write zone boundary conditions */
    void write_zone_boco(int izone);
    /* Write zone solution */
    void write_zone_solution(int izone, int isol);
    /* Write solution field data for a zone */
    void write_solution_field(int izone, int isol, int ifld);
    /*************************************************/

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#ifdef __cplusplus
}
#endif

#endif	/* _CGNSIO_H */
