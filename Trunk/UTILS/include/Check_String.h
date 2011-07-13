/*******************************************************************************
 * File:        check_string.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _CHECK_STRING_H
#define _CHECK_STRING_H

#include <stdarg.h>

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| 'macro'-functions:
|
| check_strcpy() should be used instead of strcpy().
| To prevent memory mistakes the size of the destination 'dst' has to be
| specified in 'len'.
|
| check_strcutcpy() works similar to strncpy() but puts always a string end
| symbol at the end, even if 'len'-characters are copied.
| To prevent memory mistakes dst[] has to be larger or equal than dst[len + 1].
------------------------------------------------------------------------------*/
#define check_strcpy(dst, len, source) \
        private_check_strcpy((dst), (len), (source), __FILE__, __LINE__, TRUE)

#define check_strcutcpy(dst, source, len) \
        private_check_strcpy((dst), (len) + 1, (source), __FILE__, __LINE__, \
                             FALSE)

/*******************************************************************************
* The size of the destination is limited: dststring[max_len].
* Therefore, less than max_len characters are copied onto dststring.
* Depending on the warn_flag a warning is given,
*   if 'source' is not completely copied to 'dststring'.
* The string end symbol is appended to dststring to enable further string
* handling of dststring.
*******************************************************************************/
char *private_check_strcpy(char *dststring, size_t max_len, const char *source,
                           char *filename, int line, int warn_flag);


/*******************************************************************************
* declare c99 construct in case of c89 compilation
*******************************************************************************/
#ifdef __USE_ISOC99
# define check_snprintf snprintf
#else
extern int check_snprintf(char *string, size_t maxlen, const char *format, ...);
#endif

/*******************************************************************************
* strbftok() divides a string into tokens (similar to strtok()) but uses an
* external buffer tok.
* The function copies the first token found in s to tok and returns a pointer
* to the next token.
* The standard function strtok() should NEVER be used in library
* functions because it fails with nested calls.
*******************************************************************************/
char *strbftok(const char *s, char *tok, size_t toklen, char *delim);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _CHECK_STRING_H */
