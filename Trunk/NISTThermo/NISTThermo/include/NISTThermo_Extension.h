/*******************************************************************************
 * File:        NISTThermo_Extension.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _NISTTHERMO_EXTENSION_H
#define	_NISTTHERMO_EXTENSION_H

/* Include All Original NIST Reference Properties Interface */
#include "NISTThermo.h"

#define NIST_EX_THERM2_DIM      29

#ifndef WIN32
#if defined (UNDERSCORE)

// Extended Temperature Density Based Functions
#define TDSSOUND    tdssound_
#define TDDPDT2     tddpdt2_
#define TDDPDTD     tddpdtd_
#define TDDHDTCD    tddhdtcd_
#define TDDHDTCP    tddhdtcp_
#define TDDHDDCT    tddhddct_
#define TDDHDDCP    tddhddcp_
#define TDDHDPCT    tddhdpct_
#define TDDHDPCD    tddhdpcd_
// Extended Multiphase Temperature Density Based Functions
#define TDETDFLSH   tdetdflsh_
#define TDETDFLSH2  tdetdflsh2_
#define TDESSOUND   tdessound_
#define TDEPRESS    tdepress_
#define TDEENTRO    tdeentro_
#define TDEENTHAL   tdeenthal_
#define TDEENERGY   tdeenergy_
#define TDECVCP     tdecvcp_
#define TDEDPDD     tdedpdd_
#define TDEDPDT     tdedpdt_
#define TDEDDDT     tdedddt_
#define TDEDPDD2    tdedpdd2_
#define TDEDPDT2    tdedpdt2_
#define TDEDPDTD    tdedpdtd_
#define TDEDHDTCD   tdedhdtcd_
#define TDEDHDTCP   tdedhdtcp_
#define TDEDHDDCT   tdedhddct_
#define TDEDHDDCP   tdedhddcp_
#define TDEDHDPCT   tdedhdpct_
#define TDEDHDPCD   tdedhdpcd_

// Extended Pressure Density Based Functions
#define PDSSOUND    pdssound_
// Extended Multiphase Pressure Density Based Functions
#define PDEPDFLSH   pdepdflsh_
#define PDEPDFLSH2  pdepdflsh2_
#define PDESSOUND   pdessound_
#define PDETEMP     pdetemp_
#define PDEENTRO    pdeentro_
#define PDEENTHAL   pdeenthal_
#define PDEENERGY   pdeenergy_
#define PDECVCP     pdecvcp_
#define PDEDPDD     pdedpdd_
#define PDEDPDT     pdedpdt_
#define PDEDDDT     pdedddt_
#define PDEDPDD2    pdedpdd2_
#define PDEDPDT2    pdedpdt2_
#define PDEDPDTD    pdedpdtd_
#define PDEDHDTCD   pdedhdtcd_
#define PDEDHDTCP   pdedhdtcp_
#define PDEDHDDCT   pdedhddct_
#define PDEDHDDCP   pdedhddcp_
#define PDEDHDPCT   pdedhdpct_
#define PDEDHDPCD   pdedhdpcd_

#else

// Extended Temperature Density Based Functions
#define TDSSOUND    tdssound
#define TDDPDT2     tddpdt2
#define TDDPDTD     tddpdtd
#define TDDHDTCD    tddhdtcd
#define TDDHDTCP    tddhdtcp
#define TDDHDDCT    tddhddct
#define TDDHDDCP    tddhddcp
#define TDDHDPCT    tddhdpct
#define TDDHDPCD    tddhdpcd
// Extended Multiphase Temperature Density Based Functions
#define TDETDFLSH   tdetdflsh
#define TDETDFLSH2  tdetdflsh2
#define TDESSOUND   tdessound
#define TDEPRESS    tdepress
#define TDEENTRO    tdeentro
#define TDEENTHAL   tdeenthal
#define TDEENERGY   tdeenergy
#define TDECVCP     tdecvcp
#define TDEDPDD     tdedpdd
#define TDEDPDT     tdedpdt
#define TDEDDDT     tdedddt
#define TDEDPDD2    tdedpdd2
#define TDEDPDT2    tdedpdt2
#define TDEDPDTD    tdedpdtd
#define TDEDHDTCD   tdedhdtcd
#define TDEDHDTCP   tdedhdtcp
#define TDEDHDDCT   tdedhddct
#define TDEDHDDCP   tdedhddcp
#define TDEDHDPCT   tdedhdpct
#define TDEDHDPCD   tdedhdpcd

// Extended Pressure Density Based Functions
#define PDSSOUND    pdssound
// Extended Multiphase Pressure Density Based Functions
#define PDEPDFLSH   pdepdflsh
#define PDEPDFLSH2  pdepdflsh2
#define PDESSOUND   pdessound
#define PDETEMP     pdetemp
#define PDEENTRO    pdeentro
#define PDEENTHAL   pdeenthal
#define PDEENERGY   pdeenergy
#define PDECVCP     pdecvcp
#define PDEDPDD     pdedpdd
#define PDEDPDT     pdedpdt
#define PDEDDDT     pdedddt
#define PDEDPDD2    pdedpdd2
#define PDEDPDT2    pdedpdt2
#define PDEDPDTD    pdedpdtd
#define PDEDHDTCD   pdedhdtcd
#define PDEDHDTCP   pdedhdtcp
#define PDEDHDDCT   pdedhddct
#define PDEDHDDCP   pdedhddcp
#define PDEDHDPCT   pdedhdpct
#define PDEDHDPCD   pdedhdpcd

#endif /* UNDERSCORE */
#endif

#ifdef	__cplusplus
extern "C" {
#endif

// Definitions of the Refprop extension types
PREFIX TDETDFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TDETDFLSH2 __ProtoGlarp__((double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PDEPDFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PDEPDFLSH2 __ProtoGlarp__((double *,double *,double *,double *,double *,double *,int *,char*,int ));

#ifdef	__cplusplus
}
#endif

#endif	/* _NISTTHERMO_EXTENSION_H */

