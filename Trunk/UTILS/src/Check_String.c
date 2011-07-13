/*******************************************************************************
 * File:        Check_String.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <string.h>
#include <stdarg.h>

#include "Check_String.h"
#include "Warn_Error.h"

/*******************************************************************************
 *
 ******************************************************************************/
char *private_check_strcpy(char *newstring, size_t max_len, const char *source,
        char *filename, int line, int warn_flag) {
    size_t string_end = 0;

    if (newstring == NULL || max_len < 1)
        fatal_error("Can not copy to NULL in %s line %i!", filename, line);

    if (source != NULL) {
        string_end = strlen(source);

        if (string_end >= max_len) {
            string_end = max_len - 1;

            if (warn_flag)
                warning("String: '%s'\n is reduced to '%i' characters\n"
                    " due to memory requirements in %s line %i\n",
                    source, max_len, filename, line);
        }
        strncpy(newstring, source, string_end);
    }

    strcpy(&newstring[string_end], "");

    return newstring;

} /** private_check_strcpy() **/


#ifndef __USE_ISOC99
/*******************************************************************************
 * define c99 construct in case of c89 compilation
 ******************************************************************************/
#define MAX_STRING_LENGTH 2048

extern int check_snprintf(char *string, size_t maxlen, const char *format, ...) {
    char ctmp[MAX_STRING_LENGTH] = "";
    size_t len = 0;
    va_list argument;

    if (string == NULL || maxlen == 0)
        return -1;

    va_start(argument, format);

    vsprintf(ctmp, format, argument);

    len = strlen(ctmp);

    if (len >= MAX_STRING_LENGTH)
        fatal_error("MAX_STRING_LENGTH too small!");

    va_end(argument);

    check_strcutcpy(string, ctmp, maxlen - 1);

    return len;
} /** check_snprintf() **/

#endif

/*******************************************************************************
 * strbftok() divides a string into tokens (like strtok()) but uses the external
 * buffer tok.
 * It copies the first token found in s to tok and returns a pointer to the
 * next token.
 * The standard function strtok() should NEVER be used in library
 * functions.
 ******************************************************************************/
char *strbftok(const char *s, char *tok, size_t toklen, char *delim) {
    char *lim;
    char *d;

    /*--------------------------------------------------------------------------
    | nothing to do if s starts with 0
    --------------------------------------------------------------------------*/
    if (!*s)
        return NULL;

    lim = tok + toklen - 1;

    /*--------------------------------------------------------------------------
    | search for the next character which is equal to a delimiter
    --------------------------------------------------------------------------*/
    while (*s && tok < lim) {
        for (d = delim; *d; d++) {
            if (*s == *d) {
                *tok = 0;
                /*--------------------------------------------------------------
                | skip all delimiters at the beginning of s and return the pointer
                --------------------------------------------------------------*/
                s += strspn(s, delim);
                return (char *) s;
            }
        }
        *tok++ = *s++;
    }
    *tok = 0;

    return (char *) s;

} /** strbftok() **/

