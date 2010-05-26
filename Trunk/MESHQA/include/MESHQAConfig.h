/*******************************************************************************
 * File:        MESHQAConfig.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MESHQACONFIG_H
#define	_MESHQACONFIG_H


#if defined(WIN32)
/* ====================== WIN32 ==================== */
#define MESHQA_TARGETPLATFORM        "WIN32"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1

#elif defined(LINUX_IA32)
/* ====================== Linux IA32 ==================== */
#define MESHQA_TARGETPLATFORM        "LINUX_IA32"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1
#define MESHQA_HAVE_EXECINFO         1
#ifndef _GNU_SOURCE
#define _GNU_SOURCE                  1
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS            64
#endif

#elif defined(LINUX_AMD64)
/* ====================== Linux AMD64 ==================== */
#define MESHQA_TARGETPLATFORM        "LINUX_AMD64"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1
#define MESHQA_HAVE_EXECINFO         1
#ifndef _GNU_SOURCE
#define _GNU_SOURCE                  1
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS            64
#endif

#elif defined(LINUX_IA64)

/* ====================== Linux IA64 ==================== */
#define MESHQA_TARGETPLATFORM        "LINUX_IA64"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1
#define MESHQA_HAVE_EXECINFO         1 /* may not be true */

#elif defined(AIX)
/* ====================== AIX 5.3 ==================== */
#define MESHQA_TARGETPLATFORM        "AIX"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1

#elif defined(DARWIN)
/* ====================== DARWIN ==================== */
#define MESHQA_TARGETPLATFORM        "DARWIN"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1
#define MESHQA_HAVE_EXECINFO         1 /* only available starting with 10.5 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE                  1
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS            64
#endif

#elif defined(HPUX_64)
/* ====================== HP ==================== */
#define MESHQA_TARGETPLATFORM        "HPUX_64"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1

#elif defined(IRIX)
/* ====================== IRIX ==================== */
#define MESHQA_TARGETPLATFORM        "IRIX"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1

#elif defined(SUNOS)
/* ====================== SUNOS =================== */
#define MESHQA_TARGETPLATFORM        "SUNOS"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1

#elif defined(OSF1_ALPHA)
/* ====================== DEC ALPHA =================== */
/* NOTE: we use the GNU compiler, not aCC! */
#define MESHQA_TARGETPLATFORM        "OSF1_ALPHA"
#define MESHQA_HAVE_STD_STREAM       1
#define MESHQA_HAVE_STD_STL          1
#define MESHQA_HAVE_SNPRINTF         1
#ifndef _GNU_SOURCE
#define _GNU_SOURCE                  1
#endif

#else
#error ---> This is an error. Some system must be defined at this stage.
#endif


#ifndef WIN32

#ifndef __stdcall
#define __stdcall
#endif

#ifndef __cdecl
#define __cdecl
#endif

#endif


#ifdef __cplusplus

#ifdef MESHQA_HAVE_STD_STREAM
#define MESHQA_IS_OPEN(x)   x.is_open()
#define MESHQA_IOS_FLAGS(f) std::ios::fmtflags(f)
#else
#define MESHQA_IS_OPEN(x)   x
#define MESHQA_IOS_FLAGS(f) f
#endif

#ifdef MESHQA_HAVE_STD_STREAM
#define MESHQA_STREAM(x) std::x
#define MESHQA_ENDL      std::endl
#define MESHQA_FLUSH     std::flush
#else
#define MESHQA_STREAM(x) x
#define MESHQA_ENDL      endl
#define MESHQA_FLUSH     flush
#endif

#ifdef MESHQA_HAVE_STD_STL
#define MESHQA_STL(x) std::x
#else
#define MESHQA_STL(x) x
#endif

#endif  /* __cplusplus */

#endif	/* _MESHQACONFIG_H */

