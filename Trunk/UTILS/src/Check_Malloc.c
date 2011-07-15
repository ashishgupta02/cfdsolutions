/*******************************************************************************
 * File:        Check_Malloc.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Check_Malloc.h"
#include "Utils.h"
#include "Log.h"
#include "Warn_Error.h"

#ifdef DEBUG_CHECK_MALLOC
#include "Check_String.h"

#define MAX_SOURCE_NAME 50

typedef struct {
    void *ptr;
    char alloc_call[MAX_SOURCE_NAME];
    int alloc_line_no;
    char free_call[MAX_SOURCE_NAME];
    int free_line_no;
    size_t bytes;

} MemoryAllocationStack;

static MemoryAllocationStack *stack = NULL;
static int stack_size = 0;
static int call_no = -1;
static char print_format[] = "Pointer %p allocated in %s line %i,\n"
        "                is freed in %s line %i.\n\n";

static void register_alloc(void *ptr, size_t bytes, const char *sourcename,
        int line_no);
static void register_free(void *ptr, const char *sourcename, int line_no);

#endif /* DEBUG_CHECK_MALLOC */

/*******************************************************************************
 *
 ******************************************************************************/
void *private_check_realloc(void *old, size_t bytes,
        const char *sourcename, int line_no) {
    void *tmp;

    if (bytes <= (size_t) 0) {
        if (old != NULL) {
#ifdef DEBUG_CHECK_MALLOC
            register_free(old, sourcename, line_no);
#endif
            free(old);
        }
#ifdef DEBUG_CHECK_MALLOC
        ndm_errmsg("WARNING of check_realloc :\n"
                "        Allocating %i Bytes returns NULL to \n"
                "        File: %s  Line: %d\n", (int) bytes, sourcename,
                line_no);
#endif
        return NULL;
    }

    tmp = realloc(old, bytes);
#ifdef DEBUG_CHECK_MALLOC
    register_free(old, sourcename, line_no);
    register_alloc(tmp, bytes, sourcename, line_no);
#endif

    if (tmp == NULL) {
        ndm_errmsg("ERROR : Out of memory ?!\n"
                "        Trying to (re)allocate %i Bytes failed !\n"
                "        File: %s  Line: %d\n", (int) bytes, sourcename,
                line_no);
        ndm_exit("exit ndm");
    }
    return tmp;
} /** private_check_realloc() **/

/*******************************************************************************
 *
 ******************************************************************************/
void *private_check_calloc(size_t number, size_t bytes,
        const char *sourcename, int line_no) {
    void *tmp;

    if (number <= (size_t) 0 || bytes <= (size_t) 0) {
#ifdef DEBUG_CHECK_MALLOC
        ndm_errmsg("WARNING of check_calloc :\n"
                "        Allocating %i times %i Bytes returns NULL to \n"
                "        File: %s  Line: %d\n", (int) number, (int) bytes,
                sourcename, line_no);
#endif
        return NULL;
    }

    tmp = calloc(number, bytes);
#ifdef DEBUG_CHECK_MALLOC
    register_alloc(tmp, bytes * number, sourcename, line_no);
#endif

    if (tmp == NULL) {
        ndm_errmsg("ERROR : Out of memory ?!\n"
                "        Trying to (m)allocate %i Bytes failed !\n"
                "        File: %s  Line: %d\n", (int) (number * bytes),
                sourcename, line_no);
        ndm_exit("exit ndm");
    }
    return tmp;
} /** private_check_calloc() **/

/*******************************************************************************
 *
 ******************************************************************************/
void *private_check_malloc(size_t bytes, const char *sourcename, int line_no) {
    void *tmp;

    if (bytes <= (size_t) 0) {
#ifdef DEBUG_CHECK_MALLOC
        ndm_errmsg("WARNING of check_malloc :\n"
                "        Allocating %i Bytes returns NULL to \n"
                "        File: %s  Line: %d\n", (int) bytes, sourcename,
                line_no);
#endif
        return NULL;
    }

    tmp = malloc(bytes);
#ifdef DEBUG_CHECK_MALLOC
    register_alloc(tmp, bytes, sourcename, line_no);
#endif

#ifdef DEBUG_MEMPAT
    if (DEBUG_MEMPAT == 1)
        memset(tmp, 0, bytes);
    if (DEBUG_MEMPAT == 2)
        memset(tmp, 255, bytes);
#endif

    if (tmp == NULL) {
        ndm_errmsg("ERROR : Out of memory ?!\n"
                "        Trying to (m)allocate %i Bytes failed !\n"
                "        File: %s  Line: %d\n", (int) bytes, sourcename,
                line_no);
        ndm_exit("exit ndm");
    }
    return tmp;
} /** private_check_malloc() **/

/*******************************************************************************
 *
 ******************************************************************************/
void private_free(void *ptr, const char *sourcename, int line_no) {
    if (ptr != NULL) {
#ifdef DEBUG_CHECK_MALLOC
        register_free(ptr, sourcename, line_no);
#endif
        free(ptr);
    }

    return;

} /** private_free() **/


#ifdef DEBUG_CHECK_MALLOC

/*******************************************************************************
 *
 ******************************************************************************/
static void register_alloc(void *ptr, size_t bytes, const char *sourcename,
        int line_no) {
    if (ptr == NULL) /* out of memory */ {
        check_malloc_print_memory_usage_file("CHECK_MALLOC.debug");
        return;
    }

    call_no++;

    /*--------------------------------------------------------------------------
    | Prepare stack memory and init free_line_no = -1
    --------------------------------------------------------------------------*/
    if (call_no >= stack_size) {
        static const int blocksize = 2048;
        int i;

        stack_size += blocksize;
        stack = (MemoryAllocationStack *)
                realloc(stack, stack_size * sizeof (MemoryAllocationStack));

        for (i = stack_size - blocksize; i < stack_size; i++)
            stack[i].free_line_no = -1;
    }

    /*--------------------------------------------------------------------------
    | Store information
    --------------------------------------------------------------------------*/
    stack[call_no].ptr = ptr;
    stack[call_no].bytes = bytes;
    check_strcpy(stack[call_no].alloc_call, MAX_SOURCE_NAME, sourcename);
    stack[call_no].alloc_line_no = line_no;

} /** register_alloc() **/

/*******************************************************************************
 *
 ******************************************************************************/
static void register_free(void *ptr, const char *sourcename, int line_no) {
    int i;

    if (ptr == NULL)
        return;

    for (i = call_no; i >= 0; i--)
        if (stack[i].ptr == ptr)
            break;

    if (i >= 0 && stack[i].free_line_no != -1) {
        check_malloc_print_memory_usage_file("CHECK_MALLOC.debug");

        ndm_errmsg("ERROR in File: %s  Line: %d\n"
                "Try to free already freed pointer %p!\n\n"
                , sourcename, line_no, ptr);

        ndm_errmsg(print_format, stack[i].ptr,
                stack[i].alloc_call, stack[i].alloc_line_no,
                stack[i].free_call, stack[i].free_line_no);

        ndm_exit("exit ndm");
    }

    if (i < 0) {
        check_malloc_print_memory_usage_file("CHECK_MALLOC.debug");

        ndm_errmsg("ERROR : Try to free not allocated pointer %p!\n"
                "        File: %s  Line: %d\n", ptr, sourcename, line_no);
        ndm_exit("exit ndm");
    }

    check_strcpy(stack[i].free_call, MAX_SOURCE_NAME, sourcename);
    stack[i].free_line_no = line_no;

} /** register_free() **/

/*******************************************************************************
 *
 ******************************************************************************/
void private_print_memory_usage_file(char *filename, const char *source, int line) {
    char print_format1[] = "Pointer %p allocated in %s line %i (not freed,"
            " %i bytes)\n\n";
    int i, unfreed = 0;
    FILE *out;

    if ((out = fopen(filename, "w")) == NULL) {
        ndm_errmsg("ERROR : Can not open file %s!", filename);
        ndm_exit("exit ndm");
    }

    if (source == NULL)
        fprintf(out, "print_memory_usage_file()"
            " is called from python interface\n\n");

    else
        fprintf(out, "print_memory_usage_file() is called from %s in line %i\n\n",
            source, line);
    fprintf(out, "Statistics:\n");
    fprintf(out, "%i allocation calls registered!\n", call_no + 1);

    for (i = 0; i <= call_no; i++)
        if (stack[i].free_line_no == -1)
            unfreed++;
    if (unfreed > 0) {
        ndm_msg("\n%i allocation(s) not freed!\n", unfreed);
        fprintf(out, "%i allocation(s) not freed!\n\n", unfreed);
    } else {
        ndm_msg("\nAll memory is freed again!\n");
        fprintf(out, "All memory is freed again!\n\n");
    }

    for (i = 0; i <= call_no; i++)
        if (stack[i].free_line_no == -1)
            fprintf(out, print_format1, stack[i].ptr,
                stack[i].alloc_call, stack[i].alloc_line_no,
                stack[i].bytes);
        else
            fprintf(out, print_format, stack[i].ptr,
                stack[i].alloc_call, stack[i].alloc_line_no,
                stack[i].free_call, stack[i].free_line_no);

    fclose(out);

    if (unfreed == 0) {
        free(stack);
        stack = NULL;
        stack_size = 0;
        call_no = -1;
    }

    ndm_msg("Memory usage printed to %s!\n", filename);

} /** private_print_memory_usage_file() **/

#endif /* DEBUG_CHECK_MALLOC */


/*******************************************************************************
 *
 ******************************************************************************/
#ifdef DEBUG_CHECK_MALLOC

void ndm_print_memory_usage_file(char *filename) {
    private_print_memory_usage_file(filename, NULL, 0);
} /** ndm_print_memory_usage_file() **/

#else  /* DEBUG_CHECK_MALLOC */

void ndm_print_memory_usage_file(char *filename) {
}

#endif   /* DEBUG_CHECK_MALLOC */

