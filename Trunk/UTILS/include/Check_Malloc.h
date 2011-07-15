/*******************************************************************************
* File:        Check_Malloc.h
* Author:      Ashish Gupta
* Revision:    4
******************************************************************************/

#include "License.h"

#ifndef _CHECK_MALLOC_H
#define _CHECK_MALLOC_H

#include <stdlib.h>

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| if DEBUG_CHECK_MALLOC is defined, all calls of check_*alloc() and check_free()
| are registered and freeing errors are detected.
------------------------------------------------------------------------------*/
#ifdef DEBUG_CHECK_MALLOC
#define check_malloc_print_memory_usage_file(filename) \
     private_print_memory_usage_file((filename), __FILE__, __LINE__)

extern void private_print_memory_usage_file(char *filename, const char *source,
        int line);
#endif

extern void ndm_print_memory_usage_file(char *filename); /* python wrapper */


/*******************************************************************************
 * general alloc, free functions including checks for pointers equal NULL
 *******************************************************************************/
#define check_calloc(no, sz)   private_check_calloc((no), (sz), __FILE__, \
                                                            __LINE__)

#define check_malloc(sz)       private_check_malloc((sz), __FILE__, __LINE__)

#define check_realloc(ptr, sz) private_check_realloc((ptr), (sz), __FILE__,\
                                                              __LINE__)

#define ALLOC(sz) (((sz) == 0) ? NULL : check_malloc(sz))

#define check_free(ptr) (private_free((void *)(ptr), __FILE__, __LINE__),\
                     (ptr) = NULL)

extern void *private_check_calloc(size_t number, size_t bytes,
        const char *sourcename, int line_no);

extern void *private_check_malloc(size_t bytes, const char *sourcename, int line_no);

extern void *private_check_realloc(void *old, size_t bytes,
        const char *sourcename, int line_no);

extern void private_free(void *ptr, const char *sourcename, int line_no);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _CHECK_MALLOC_H */
