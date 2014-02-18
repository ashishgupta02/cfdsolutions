/*******************************************************************************
 * File:        NISTThermo.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _NISTTHERMO_H
#define	_NISTTHERMO_H

/* Defines */
#define NISTTHERMO_REF_STR_LEN      3
#define NISTTHERMO_GEN_STR_LEN      255
#define NISTTHERMO_FILE_STR_LEN     5100 /* 255x20 */
#define NISTTHERMO_MAX_COMPONENT    20
#define NISTTHERMO_MAX_PARAMETER    72
#define NISTTHERMO_MAX_COEFFICIENT  50

/* For windows applications */
#ifdef WIN32
#ifdef __MINGW32__
#undef WIN32
#endif
#endif

#if defined(WIN32)
#ifdef ABSOFT
#undef WIN32
#define PREFIX extern void __cdecl
#undef  UNDERSCORE
#else
#ifdef IFORT
#define PREFIX extern void
#else
#define PREFIX extern void __stdcall
#endif /* IFORT */
#endif /* ABSOFT */
#endif

/* For non windows applications:
  (i)  Change case of routine names
  (ii) Add underscores where appropriate */
#if defined(mips) || defined(__mips)  || defined(__MACH__) || defined(__alpha) || defined(__linux) || defined(__sparc) || defined (__MINGW32__)
/* SGI & DECStation; SGI;  MAC OS X; DEC Alpha; Linux; SUN */
#define PREFIX extern void
#define UNDERSCORE
#else 
/* HP 9000/700; IBM RS6000 */
#if defined(__hpux) || defined(_IBMR2)
#define PREFIX extern void
#undef  UNDERSCORE
#endif 
#endif

#ifndef WIN32
#if defined (UNDERSCORE)
#define ABFL1   abfl1_
#define ABFL2   abfl2_
#define ACTVY   actvy_
#define AG      ag_
#define CCRIT   ccrit_
#define CP0     cp0_
#define CRITP   critp_
#define CSATK   csatk_
#define CV2PK   cv2pk_
#define CVCPK   cvcpk_
#define CVCP    cvcp_
#define DBDT    dbdt_
#define DBFL1   dbfl1_
#define DBFL2   dbfl2_
#define DDDP    dddp_
#define DDDT    dddt_
#define DEFLSH  deflsh_
#define DHD1    dhd1_
#define DHFLSH  dhflsh_
#define DIELEC  dielec_
#define DOTFILL dotfill_
#define DPDD2   dpdd2_
#define DPDDK   dpddk_
#define DPDD    dpdd_
#define DPDTK   dpdtk_
#define DPDT    dpdt_
#define DPTSATK dptsatk_
#define DSFLSH  dsflsh_
#define ENTHAL  enthal_
#define ENTRO   entro_
#define ESFLSH  esflsh_
#define FGCTY   fgcty_
#define FPV     fpv_
#define GERG04  gerg04_
#define GETFIJ  getfij_
#define GETKTV  getktv_
#define GIBBS   gibbs_
#define HSFLSH  hsflsh_
#define INFO    info_
#define LIMITK  limitk_
#define LIMITS  limits_
#define LIMITX  limitx_
#define MELTP   meltp_
#define MELTT   meltt_
#define MLTH2O  mlth2o_
#define NAME    name_
#define PDFL1   pdfl1_
#define PDFLSH  pdflsh_
#define PEFLSH  peflsh_
#define PHFL1   phfl1_
#define PHFLSH  phflsh_
#define PQFLSH  pqflsh_
#define PREOS   preos_
#define PRESS   press_
#define PSFL1   psfl1_
#define PSFLSH  psflsh_
#define PUREFLD purefld_
#define QMASS   qmass_
#define QMOLE   qmole_
#define SATD    satd_
#define SATE    sate_
#define SATH    sath_
#define SATP    satp_
#define SATS    sats_
#define SATT    satt_
#define SETAGA  setaga_
#define SETKTV  setktv_
#define SETMIX  setmix_
#define SETMOD  setmod_
#define SETREF  setref_
#define SETUP   setup_
#define SETPATH setpath_
#define SPECGR  specgr_
#define SUBLP   sublp_
#define SUBLT   sublt_
#define SURFT   surft_
#define SURTEN  surten_
#define TDFLSH  tdflsh_
#define TEFLSH  teflsh_
#define THERM0  therm0_
#define THERM2  therm2_
#define THERM3  therm3_
#define THERM   therm_
#define THFLSH  thflsh_
#define TPFLSH  tpflsh_
#define TPRHO   tprho_
#define TQFLSH  tqflsh_
#define TRNPRP  trnprp_
#define TSFLSH  tsflsh_
#define VIRB    virb_
#define VIRC    virc_
#define WMOL    wmol_
#define XMASS   xmass_
#define XMOLE   xmole_

#else
#define ABFL1   abfl1
#define ABFL2   abfl2
#define ACTVY   actvy
#define AG      ag
#define CCRIT   ccrit
#define CP0     cp0
#define CRITP   critp
#define CSATK   csatk
#define CV2PK   cv2pk
#define CVCPK   cvcpk
#define CVCP    cvcp
#define DBDT    dbdt
#define DBFL1   dbfl1
#define DBFL2   dbfl2
#define DDDP    dddp
#define DDDT    dddt
#define DEFLSH  deflsh
#define DHD1    dhd1
#define DHFLSH  dhflsh
#define DIELEC  dielec
#define DOTFILL dotfill
#define DPDD2   dpdd2
#define DPDDK   dpddk
#define DPDD    dpdd
#define DPDTK   dpdtk
#define DPDT    dpdt
#define DPTSATK dptsatk
#define DSFLSH  dsflsh
#define ENTHAL  enthal
#define ENTRO   entro
#define ESFLSH  esflsh
#define FGCTY   fgcty
#define FPV     fpv
#define GERG04  gerg04
#define GETFIJ  getfij
#define GETKTV  getktv
#define GIBBS   gibbs
#define HSFLSH  hsflsh
#define INFO    info
#define LIMITK  limitk
#define LIMITS  limits
#define LIMITX  limitx
#define MELTP   meltp
#define MELTT   meltt
#define MLTH2O  mlth2o
#define NAME    name
#define PDFL1   pdfl1
#define PDFLSH  pdflsh
#define PEFLSH  peflsh
#define PHFL1   phfl1
#define PHFLSH  phflsh
#define PQFLSH  pqflsh
#define PREOS   preos
#define PRESS   press
#define PSFL1   psfl1
#define PSFLSH  psflsh
#define PUREFLD purefld
#define QMASS   qmass
#define QMOLE   qmole
#define SATD    satd
#define SATE    sate
#define SATH    sath
#define SATP    satp
#define SATS    sats
#define SATT    satt
#define SETAGA  setaga
#define SETKTV  setktv
#define SETMIX  setmix
#define SETMOD  setmod
#define SETREF  setref
#define SETUP   setup
#define SETPATH setpath
#define SPECGR  specgr
#define SUBLP   sublp
#define SUBLT   sublt
#define SURFT   surft
#define SURTEN  surten
#define TDFLSH  tdflsh
#define TEFLSH  teflsh
#define THERM0  therm0
#define THERM2  therm2
#define THERM3  therm3
#define THERM   therm
#define THFLSH  thflsh
#define TPFLSH  tpflsh
#define TPRHO   tprho
#define TQFLSH  tqflsh
#define TRNPRP  trnprp
#define TSFLSH  tsflsh
#define VIRB    virb
#define VIRC    virc
#define WMOL    wmol
#define XMASS   xmass
#define XMOLE   xmole

#endif /* UNDERSCORE */
#endif

#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef	__cplusplus
extern "C" {
#endif

// Definitions of the Refprop types
#if defined(WIN32) && !defined(IFORT)
PREFIX GETFIJ  __ProtoGlarp__((char*,int,double *,char*,int,char*,int ));
PREFIX GETKTV  __ProtoGlarp__((int *,int *,char*,int,double *,char*,int,char*,int,char*,int,char*,int ));
PREFIX LIMITK  __ProtoGlarp__((char*,int,int *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX LIMITS  __ProtoGlarp__((char*,int,double *,double *,double *,double *,double * ));
PREFIX LIMITX  __ProtoGlarp__((char*,int,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SETKTV  __ProtoGlarp__((int *,int *,char*,int,double *,char*,int,int *,char*,int ));
PREFIX SETMIX  __ProtoGlarp__((char*,int,char*,int,char*,int,int *,char*,int,double *,int *,char*,int ));
PREFIX SETMOD  __ProtoGlarp__((int *,char*,int,char*,int,char*,int,int *,char*,int ));
PREFIX SETREF  __ProtoGlarp__((char*,int,int *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SETUP   __ProtoGlarp__((int *,char*,int,char*,int,char*,int,int *,char*,int ));
#else
PREFIX GETFIJ  __ProtoGlarp__((char*,double *,char*,char*,int ,int ,int ));
PREFIX GETKTV  __ProtoGlarp__((int *,int *,char*,double *,char*,char*,char*,char*,int ,int ,int ,int ,int ));
PREFIX LIMITK  __ProtoGlarp__((char*,int *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ,int ));
PREFIX LIMITS  __ProtoGlarp__((char*,double *,double *,double *,double *,double *,int ));
PREFIX LIMITX  __ProtoGlarp__((char*,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ,int ));
PREFIX SETKTV  __ProtoGlarp__((int *,int *,char*,double *,char*,int *,char*,int ,int ,int ));
PREFIX SETMIX  __ProtoGlarp__((char*,char*,char*,int *,char*,double *,int *,char*,int ,int ,int ,int ,int ));
PREFIX SETMOD  __ProtoGlarp__((int *,char*,char*,char*,int *,char*,int ,int ,int ,int ));
PREFIX SETREF  __ProtoGlarp__((char*,int *,double *,double *,double *,double *,double *,int *,char*,int ,int ));
PREFIX SETUP   __ProtoGlarp__((int *,char*,char*,char*,int *,char*,int ,int ,int ,int ));
#endif

PREFIX ABFL1   __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX ABFL2   __ProtoGlarp__((double *,double *,double *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX ACTVY   __ProtoGlarp__((double *,double *,double *,double *));
PREFIX AG      __ProtoGlarp__((double *,double *,double *,double *,double *));
PREFIX CCRIT   __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX CP0     __ProtoGlarp__((double *,double *,double *));
PREFIX CRITP   __ProtoGlarp__((double *,double *,double *,double *,int *,char*,int ));
PREFIX CSATK   __ProtoGlarp__((int *,double *,int *,double *,double *,double *,int *,char*,int ));
PREFIX CV2PK   __ProtoGlarp__((int *,double *,double *,double *,double *,int *,char*,int ));
PREFIX CVCPK   __ProtoGlarp__((int *,double *,double *,double *,double *));
PREFIX CVCP    __ProtoGlarp__((double *,double *,double *,double *,double *));
PREFIX DBDT    __ProtoGlarp__((double *,double *,double *));
PREFIX DBFL1   __ProtoGlarp__((double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX DBFL2   __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX DDDP    __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DDDT    __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DEFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX DHD1    __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX DHFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX DIELEC  __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DOTFILL __ProtoGlarp__((int *,double *,double *,double *,int *,char*,int ));
PREFIX DPDD2   __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DPDDK   __ProtoGlarp__((int *,double *,double *,double *));
PREFIX DPDD    __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DPDTK   __ProtoGlarp__((int *,double *,double *,double *));
PREFIX DPDT    __ProtoGlarp__((double *,double *,double *,double *));
PREFIX DPTSATK __ProtoGlarp__((int *,double *,int *,double *,double *,double *,double *,int *,char*,int ));
PREFIX DSFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX ENTHAL  __ProtoGlarp__((double *,double *,double *,double *));
PREFIX ENTRO   __ProtoGlarp__((double *,double *,double *,double *));
PREFIX ESFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX FGCTY   __ProtoGlarp__((double *,double *,double *,double *));
PREFIX FPV     __ProtoGlarp__((double *,double *,double *,double *,double *));
PREFIX GERG04  __ProtoGlarp__((int *,int *,int *,char*,int ));
PREFIX GIBBS   __ProtoGlarp__((double *,double *,double *,double *,double *));
PREFIX HSFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX INFO    __ProtoGlarp__((int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX MELTP   __ProtoGlarp__((double *,double *,double *,int *,char*,int ));
PREFIX MELTT   __ProtoGlarp__((double *,double *,double *,int *,char*,int ));
PREFIX MLTH2O  __ProtoGlarp__((double *,double *,double *));
PREFIX NAME    __ProtoGlarp__((int *,char*,char*,char*,int ,int ,int ));
PREFIX PDFL1   __ProtoGlarp__((double *,double *,double *,double *,int *,char*,int ));
PREFIX PDFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PEFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PHFL1   __ProtoGlarp__((double *,double *,double *,int *,double *,double *,int *,char*,int ));
PREFIX PHFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PQFLSH  __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PREOS   __ProtoGlarp__((int *));
PREFIX PRESS   __ProtoGlarp__((double *,double *,double *,double *));
PREFIX PSFL1   __ProtoGlarp__((double *,double *,double *,int *,double *,double *,int *,char*,int ));
PREFIX PSFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX PUREFLD __ProtoGlarp__((int *));
PREFIX QMASS   __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX QMOLE   __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SATD    __ProtoGlarp__((double *,double *,int *,int *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SATE    __ProtoGlarp__((double *,double *,int *,int *,int *,double *,double *,double *,int *,double *,double *,double *,int *,char*,int ));
PREFIX SATH    __ProtoGlarp__((double *,double *,int *,int *,int *,double *,double *,double *,int *,double *,double *,double *,int *,char*,int ));
PREFIX SATP    __ProtoGlarp__((double *,double *,int *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SATS    __ProtoGlarp__((double *,double *,int *,int *,int *,double *,double *,double *,int *,double *,double *,double *,int *,double *,double *,double *,int *,char*,int ));
PREFIX SATT    __ProtoGlarp__((double *,double *,int *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX SETAGA  __ProtoGlarp__((int *,char*,int ));
PREFIX SETPATH __ProtoGlarp__((char*,int ));
PREFIX SPECGR  __ProtoGlarp__((double *,double *,double *,double *));
PREFIX SUBLP   __ProtoGlarp__((double *,double *,double *,int *,char*,int ));
PREFIX SUBLT   __ProtoGlarp__((double *,double *,double *,int *,char*,int ));
PREFIX SURFT   __ProtoGlarp__((double *,double *,double *,double *,int *,char*,int ));
PREFIX SURTEN  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TDFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TEFLSH  __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX THERM0  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX THERM2  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX THERM3  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX THERM   __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *));
PREFIX THFLSH  __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TPFLSH  __ProtoGlarp__((double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TPRHO   __ProtoGlarp__((double *,double *,double *,int *,int *,double *,int *,char*,int ));
PREFIX TQFLSH  __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TRNPRP  __ProtoGlarp__((double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX TSFLSH  __ProtoGlarp__((double *,double *,double *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,char*,int ));
PREFIX VIRB    __ProtoGlarp__((double *,double *,double *));
PREFIX VIRC    __ProtoGlarp__((double *,double *,double *));
PREFIX WMOL    __ProtoGlarp__((double *,double *));
PREFIX XMASS   __ProtoGlarp__((double *,double *,double *));
PREFIX XMOLE   __ProtoGlarp__((double *,double *,double *));

#ifdef	__cplusplus
}
#endif

#endif	/* _NISTTHERMO_H */

