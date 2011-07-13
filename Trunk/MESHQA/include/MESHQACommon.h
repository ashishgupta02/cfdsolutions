/*******************************************************************************
 * File:        MESHQACommon.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MESHQACOMMON_H
#define	_MESHQACOMMON_H

#include "MESHQAConfig.h"

/*-----------------------------------------------------------------------------
//
//  Common includes
*/
#ifndef MESHQA_NO_AUTO_INCLUDES

/* --- standard C-includes --- */
/*
  Pre-include the usually needed system header files to avoid the same #ifdef-switches
  throughout most of the own header files.
*/
#include <sys/types.h>
#include <sys/stat.h>
#ifndef WIN32
#include <sys/times.h>
#endif

#include <fcntl.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h> /* <--- for exit() */
#include <errno.h>
#include <time.h>

#ifdef WIN32
#include <tchar.h>
#include <io.h>
#include <direct.h>
#include <process.h>
#else
#include <unistd.h>
#endif

#ifdef MESHQA_HAVE_PTHREAD
#include <pthread.h>
#endif

/* --- large file I/O --- */
#ifdef __USE_LARGEFILE64
#define MESHQA_ADD_FLAG_LARGEFILE(flags) (flags |= O_LARGEFILE)
#else
#define MESHQA_ADD_FLAG_LARGEFILE(flags)
#endif

/* --- math --- */
#include <math.h>

#if defined(IRIX) || defined(SUNOS)
#include <ieeefp.h>
#endif

/* --- C++ streams --- */
#ifdef MESHQA_HAVE_STD_STREAM
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#ifdef IRIX
#include <string.h>
#else
#include <cstring>
#endif
#else
#include <iostream.h>
#include <strstream.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>
#endif


/* --- STL includes --- */
#ifdef MESHQA_HAVE_STD_STL
#include <algorithm>
#include <list>
#include <vector>
#include <set>
#include <map>
#else
#include <algorithm.h>
#include <list.h>
#include <vector.h>
#include <set.h>
#include <map.h>
#endif

#include <string>

#endif // MESHQA_NO_AUTO_INCLUDES

/*-----------------------------------------------------------------------------
//
//  Macros
*/

/*!
  This is a macro which may be used in case of fatal errors.
*/
#if defined(MESHQA_ERROR)
#error ---> This should not be defined here!
#endif

#define MESHQA_ERROR " *** ERROR *** ERROR *** ERROR *** ERROR *** ERROR *** ERROR *** "


/*!
  This is a macro for checking if a pointer is zero or not.
*/
#if defined(MESHQA_ZERO_PTR)
#error ---> This should not be defined here!
#endif

#define MESHQA_ZERO_PTR    (static_cast<const char*>(0))


/*!
  These are macros for the de-allocation of dynamic memory segments.
*/
#if defined(MESHQA_DELETE) || defined(MESHQA_DELETE_ARRAY) || defined(MESHQA_FREE)
#error ---> This should not be defined here!
#endif

#define MESHQA_DELETE(ptr)       {if((ptr) != 0){ delete   (ptr); (ptr)=0;}}
#define MESHQA_DELETE_ARRAY(ptr) {if((ptr) != 0){ delete[] (ptr); (ptr)=0;}}

#define MESHQA_FREE(ptr)         {if((ptr) != NULL){ free (ptr); (ptr)=NULL;}}

/* Note: MESHQA_DELETE(ptr), etc. can be safely used multiple times without problems,
         de-allocation of memory will only performed once. */


/*!
  This is a macro for returning a boolean value.
*/
#if defined(MESHQA_BOOL)
#error ---> This should not be defined here!
#endif

#define MESHQA_BOOL(a)    ((a) ? true : false)


/*!
  This is a macro for returning a 'boolean' string.
*/
#if defined(MESHQA_BOOL2STRING)
#error ---> This should not be defined here!
#endif

#define MESHQA_BOOL2STRING(a)    ((a) != 0 ? "yes" : "no")


/*!
  This is a macro for casting a variable to an integer number.
  It is intendend to be used for char (byte) variables only.
*/
#if defined(MESHQA_BYTE2INT)
#error ---> This should not be defined here!
#endif

#define MESHQA_BYTE2INT(a)       (static_cast<int>(a))


/*!
  This is a macro for computing the square of a value.
*/
#if defined(MESHQA_SQUARED)
#error ---> This should not be defined here!
#endif

#define MESHQA_SQUARED(a)      ((a) * (a))


/*!
  This is a macro for checking whether a value is finite or not.
*/
#if defined(MESHQA_IS_FINITE)
#error ---> This should not be defined here!
#endif

#ifdef WIN32
# define MESHQA_IS_FINITE(a)      (_finite(a))
#else
# ifdef _HPUX_SOURCE
#  ifdef __GNUC__
#   define MESHQA_IS_FINITE(a)      (1)
#  else
#   define MESHQA_IS_FINITE(a)      (isfinite(a))
#  endif /* __GCC__ */
# else
#  define MESHQA_IS_FINITE(a)      (finite(a))
# endif /* _HPUX_SOURCE */
#endif


/*!
  This is a macro for checking whether a value is nan or not.
*/
#if defined(MESHQA_IS_NAN)
#error ---> This should not be defined here!
#endif

#ifdef WIN32
#define MESHQA_IS_NAN(a)      (_isnan(a))
#else
#define MESHQA_IS_NAN(a)      (isnan(a))
#endif


/*!
  This macro defines the path separator depending on the system.
*/
#if defined(MESHQA_PATH_SEPARATOR)
#error ---> This should not be defined here!
#endif

#ifdef WIN32
#define MESHQA_PATH_SEPARATOR '\\'
#else
#define MESHQA_PATH_SEPARATOR '/'
#endif


/*!
  This macro should be used in case of a fatal error.
*/
#if defined(MESHQAErrorMacro)
#error ---> This should not be defined here!
#endif

#define MESHQAErrorMacro(x) { char _msg[4096]; sprintf(_msg, " in %s in line %d:", __FILE__, __LINE__); exit(-1); }

/*!
  These macros may be used to explicitly define physical dimensions in loops, etc.
*/
#if defined(MESHQA_1D) || defined(MESHQA_2D) || defined(MESHQA_3D)
#error ---> This should not be defined here!
#endif
#if defined(MESHQA_4D) || defined(MESHQA_5D) || defined(MESHQA_6D)
#error ---> This should not be defined here!
#endif
#if defined(MESHQA_7D) || defined(MESHQA_8D) || defined(MESHQA_9D)
#error ---> This should not be defined here!
#endif

#define MESHQA_1D  1
#define MESHQA_2D  2
#define MESHQA_3D  3
#define MESHQA_4D  4
#define MESHQA_5D  5
#define MESHQA_6D  6
#define MESHQA_7D  7
#define MESHQA_8D  8
#define MESHQA_9D  9


#ifdef __cplusplus

/*-----------------------------------------------------------------------------
//
//  Templates
*/

/*!
  These are template definitions for min/max operations.
*/
template <typename T>
const T& MESHQAMin(const T& a, const T& b) {
  return (a < b ? a : b);
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b), c);
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b), MESHQAMin<T>(c,d));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d), e);
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d), MESHQAMin<T>(e,f));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d), MESHQAMin<T>(e,f,g));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d), MESHQAMin<T>(e,f,g,h));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),i);
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k,l));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k,l,m));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k,l,m,n));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n, const T& o) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k,l,m,n,o));
}

template <typename T>
const T& MESHQAMin(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n, const T& o, const T& p) {
  return MESHQAMin<T>(MESHQAMin<T>(a,b,c,d,e,f,g,h),MESHQAMin<T>(i,j,k,l,m,n,o,p));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b) {
  return (a > b ? a : b);
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b), c);
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b), MESHQAMax<T>(c,d));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d), e);
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d), MESHQAMax<T>(e,f));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d), MESHQAMax<T>(e,f,g));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d), MESHQAMax<T>(e,f,g,h));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),i);
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k,l));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k,l,m));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k,l,m,n));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n, const T& o) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k,l,m,n,o));
}

template <typename T>
const T& MESHQAMax(const T& a, const T& b, const T& c, const T& d,
               const T& e, const T& f, const T& g, const T& h,
               const T& i, const T& j, const T& k, const T& l,
               const T& m, const T& n, const T& o, const T& p) {
  return MESHQAMax<T>(MESHQAMax<T>(a,b,c,d,e,f,g,h),MESHQAMax<T>(i,j,k,l,m,n,o,p));
}


/*!
  This is a template definition for a swap operation.
*/
template <typename T>
void MESHQASwap(T& a, T& b)
{
  {
    T tmp = a;
    a     = b;
    b     = tmp;
  }
}

#endif   /* __cplusplus */

#endif	/* _MESHQACOMMON_H */

