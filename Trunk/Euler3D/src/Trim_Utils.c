/*******************************************************************************
 * File:        Trim_Utils.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#define PATH_DELIM ';'
#include <direct.h>
#include <io.h>
#else
#define PATH_DELIM ':'
#include <unistd.h>
#endif

/* Application Header Files */
#include "Trim_Utils.h"

/*---------- error -----------------------------------------------------
 * Print error message to stderr and exit
 *---------------------------------------------------------------------*/

void error(const char *fmt, ...) {
    va_list args;

    (void) fprintf(stderr, "ERROR: ");

    va_start(args, fmt);
    (void) vfprintf(stderr, fmt, args);
    va_end(args);

    (void) fprintf(stderr, "\n");

    /* To ensure log files are current */
    (void) fflush(stderr);
    exit(EXIT_FAILURE);
}

/*---------- warn ------------------------------------------------------
 * Print warning message to stderr and continue
 *---------------------------------------------------------------------*/

void warn(const char *fmt, ...) {
    va_list args;

    (void) fprintf(stderr, "WARNING: ");
    va_start(args, fmt);
    (void) vfprintf(stderr, fmt, args);
    va_end(args);

    (void) fprintf(stderr, "\n");
    /* To ensure log files are current */
    (void) fflush(stderr);
}

/*---------- info ------------------------------------------------------
 * Print Infomation message to stdout and continue
 *---------------------------------------------------------------------*/

void info(const char *fmt, ...) {
    va_list args;

    (void) fprintf(stdout, "INFO: ");
    va_start(args, fmt);
    (void) vfprintf(stdout, fmt, args);
    va_end(args);

    (void) fprintf(stdout, "\n");
    /* To ensure log files are current */
    (void) fflush(stdout);
}

/*---------- file_size ------------------------------------------------
 * Returns the size of file on disk
 *---------------------------------------------------------------------*/

double file_size(const char *fname) {
    struct stat st;

    if (stat(fname, &st)) return 0.0;
    return (double) st.st_size / 1048576.0;
}

/*---------- file_exists ----------------------------------------------
 * check if a file exists
 *---------------------------------------------------------------------*/

int file_exists(char *file) {
    struct stat st;

#ifdef _WIN32
    if (_access(file, 0) || stat(file, &st) ||
            S_IFREG != (st.st_mode & S_IFMT)) return 0;
#else
	if (access(file, 0) || stat(file, &st) ||
            S_IFREG != (st.st_mode & S_IFMT)) return 0;
#endif

    return 1;
}

/*---------- is_executable -----------------------------------------------
 * checks if pathname exists and is executable
 *------------------------------------------------------------------------*/

int is_executable(char *file) {
    struct stat st;

    /* needs to be executable and not a directory */
#ifdef _WIN32
    if (_access(file, 1) || stat(file, &st) ||
            S_IFDIR == (st.st_mode & S_IFMT))
        return (0);
#else
	if (access(file, 1) || stat(file, &st) ||
            S_IFDIR == (st.st_mode & S_IFMT))
        return (0);
#endif

    return (1);
}

#ifdef _WIN32

/*---------- check_extensions -------------------------------------------
 * check for DOS/Windows executable
 *-----------------------------------------------------------------------*/

static int check_extensions(char *pathname) {
    int n;
    char *p;
    static char *exts[] = {".com", ".exe", ".bat"};

    /* fix path */

    for (p = pathname; *p; p++) {
        if (*p == '/')
            *p = '\\';
    }
    if (is_executable(pathname))
        return (1);
    for (n = 0; n < 3; n++) {
        strcpy(p, exts[n]);
        if (is_executable(pathname))
            return (1);
    }
    *p = 0;
    return (0);
}

#endif

/*---------- find_executable --------------------------------------------
 * locate and build pathname to executable
 *-----------------------------------------------------------------------*/

char *find_executable(char *exename) {
    int n;
    char *p, *s;
    static char exepath[MAX_FILENAME_LEN + 1];

    if (exename == NULL || !*exename)
        return (NULL);

#ifdef _WIN32

    /* full path */

    if (*exename == '/' || *exename == '\\' || *(exename + 1) == ':') {
        strcpy(exepath, exename);
        return (check_extensions(exepath) ? exepath : NULL);
    }

    /* get current directory */

    if (NULL == _getcwd(exepath, MAX_FILENAME_LEN))
        error("find_executable", "couldn't get working directory");
    p = exepath + strlen(exepath);
    if (*(p - 1) == '\\')
        *--p = 0;

    /* relative path */

    if (0 == strncmp(exename, ".\\", 2) ||
            0 == strncmp(exename, "..\\", 3) ||
            0 == strncmp(exename, "./", 2) ||
            0 == strncmp(exename, "../", 3)) {
        if (exename[1] != '.')
            strcpy(p, &exename[1]);
        else {
            if (NULL == (p = strrchr(exepath, '\\')))
                p = exepath;
            strcpy(p, &exename[2]);
        }
        return (check_extensions(exepath) ? exepath : NULL);
    }

    /* current directory */

    *p++ = '\\';
    strcpy(p, exename);
    if (check_extensions(exepath))
        return (exepath);

#else

    /* full path */

    if (*exename == '/') {
        if (is_executable(exename))
            return (strcpy(exepath, exename));
        return (NULL);
    }

    /* relative path */

    if (0 == strncmp(exename, "./", 2) ||
            0 == strncmp(exename, "../", 3)) {
        if (NULL == getcwd(exepath, MAX_FILENAME_LEN))
            error("find_executable", "couldn't get working directory\n");
        p = exepath + strlen(exepath);
        if (*(p - 1) == '/')
            *--p = 0;
        if (exename[1] != '.')
            strcpy(p, &exename[1]);
        else {
            if (NULL == (p = strrchr(exepath, '/')))
                p = exepath;
            strcpy(p, &exename[2]);
        }
        return (is_executable(exepath) ? exepath : NULL);
    }

#endif

    /* scan $PATH environment variable */

    if (NULL == (p = getenv("PATH")))
        return (NULL);
    while (*p) {
        if (NULL == (s = strchr(p, PATH_DELIM))) {
            strcpy(exepath, p);
            n = strlen(exepath);
        } else {
            n = (int) (s++ - p);
            strncpy(exepath, p, n);
            exepath[n] = 0;
        }
        if (n) {
            p = exepath + n;
#ifdef _WIN32
            if (*(p - 1) != '\\')
                *p++ = '\\';
            strcpy(p, exename);
            if (check_extensions(exepath))
                return (exepath);
#else
            if (*(p - 1) != '/')
                *p++ = '/';
            strcpy(p, exename);
            if (is_executable(exepath))
                return (exepath);
#endif
        }
        if (NULL == (p = s)) break;
    }
    return (NULL);
}

/*---------- find_file ------------------------------------------------
 * get pathname to a file
 *---------------------------------------------------------------------*/

char *find_file(char *filename, char *exename) {
    char *p;
    static char pathname[MAX_FILENAME_LEN + 1];

    if (file_exists(filename))
        return strcpy(pathname, filename);
    if ((p = find_executable(exename)) == NULL)
        return NULL;
    strcpy(pathname, p);
    while ((p = strrchr(pathname, '/')) != NULL ||
            (p = strrchr(pathname, '\\')) != NULL) {
        strcpy(p + 1, filename);
        if (file_exists(pathname))
            return pathname;
        *p = 0;
    }
    return NULL;
}

/*---------- same_file ------------------------------------------------
 * check if 2 files are the same
 *---------------------------------------------------------------------*/

int same_file(char *file1, char *file2) {
    int n = (file_exists(file1) | (file_exists(file2) << 1));
#ifdef _WIN32
    char path1[257], path2[257];

    if (n == 1 || n == 2) return 0;
    if (_fullpath(path1, file1, 256) != NULL &&
            _fullpath(path2, file2, 256) != NULL)
        return (_stricmp(path1, path2) == 0);
    return (_stricmp(file1, file2) == 0);
#else
    if (n == 3) {
        struct stat st1, st2;
        stat(file1, &st1);
        stat(file2, &st2);
        return (st1.st_ino == st2.st_ino);
    }
    if (n == 0)
        return (strcmp(file1, file2) == 0);
    return 0;
#endif
}

/*---------- temporary_file -------------------------------------------
 * create a temporary file
 *---------------------------------------------------------------------*/

int temporary_file(char *name) {
    char *temp;
    size_t i, len = 0;

    /* Get the length of the name */
    len = strlen(name);
    if (len <= 6) return 1;

    /* Allocate the memory to hold name */
    temp = (char*) malloc((len + 1) * sizeof (char));
    strcpy(temp, name);
    for (i = (len - 7); i < len; i++) temp[i] = 'X';

    /*    strcpy (temp, "cgnsXXXXXX"); */
#ifdef _WIN32
    if (_mktemp(temp) == NULL)
        error("temporary_file", "failed to create temporary filename\n");
#else
    if (mkstemp(temp) == -1)
        error("temporary_file", "failed to create temporary filename\n");
#endif
    strcpy(name, temp);
    return 0;
}

/*---------- unlink_name -------------------------------------------
 * removes the link name
 *---------------------------------------------------------------------*/
int unlink_name(char *name) {
    int result;

    result = 0;
#ifdef _WIN32
    result = _unlink(name);
#else
	result = unlink(name);
#endif
    return result;
}

/*---------- copy_file ------------------------------------------------
 * make a copy of a file
 *---------------------------------------------------------------------*/

void copy_file(char *oldfile, char *newfile) {
    int c;
    FILE *oldfp, *newfp;

    if (NULL == (oldfp = fopen(oldfile, "rb")))
        error("copy_file", "error opening input file for reading\n");
    if (NULL == (newfp = fopen(newfile, "w+b"))) {
        fclose(oldfp);
        error("copy_file", "error opening output file for writing\n");
    }
    while (EOF != (c = getc(oldfp)))
        putc(c, newfp);
    fclose(oldfp);
    fclose(newfp);
}

/*---------- str_blank ------------------------------------------------
 * Blanks the strings with space
 *---------------------------------------------------------------------*/

void str_blank(char *str) {
    unsigned int i = 0;
    size_t length = 0;

    if (str != NULL) {
        length = strlen(str);
        for (i = 0; i < length; i++)
            str[i] = ' ';
        str[length] = '\0';
    }
}

/*---------- sort_low2high ----------------------------------------------
 * Sorts Data Unsorted Data in F Low to High Order and Store in List
 *---------------------------------------------------------------------*/

void sort_low2high(int n, int list[], double f[]) {
    int i, j, l, ir, ira;

    l = n / 2 + 1;
    ir = n;
    while (n > 1) {
        if (l > 1) {
            l = l - 1;
            ira = list[l - 1];
        } else {
            ira = list[ir - 1];
            list[ir - 1] = list[0];
            ir = ir - 1;
            if (ir == 1) {
                list[0] = ira;
                break;
            }
        }
        i = l;
        j = l * 2;
        while (j <= ir) {
            if (j < ir)
                if (f[list[j - 1]] < f[list[j]])
                    j = j + 1;
            if (f[ira] < f[list[j - 1]]) {
                list[i - 1] = list[j - 1];
                i = j;
                j = j + i;
            } else
                j = ir + 1;
        }
        list[i - 1] = ira;
    }
}

/*---------- sort_high2low --------------------------------------------
 * Sorts Data Unsorted Data in F High to Low Order and Store in List
 *---------------------------------------------------------------------*/

void sort_high2low(int n, int list[], double f[]) {
    int i, j, l, ir, ira;

    l = n / 2 + 1;
    ir = n;
    while (n > 1) {
        if (l > 1) {
            l = l - 1;
            ira = list[l - 1];
        } else {
            ira = list[ir - 1];
            list[ir - 1] = list[0];
            ir = ir - 1;
            if (ir == 1) {
                list[0] = ira;
                break;
            }
        }
        i = l;
        j = l * 2;
        while (j <= ir) {
            if (j < ir)
                if (f[list[j - 1]] > f[list[j]])
                    j = j + 1;
            if (f[ira] > f[list[j - 1]]) {
                list[i - 1] = list[j - 1];
                i = j;
                j = j + i;
            } else
                j = ir + 1;
        }
        list[i - 1] = ira;
    }
}

