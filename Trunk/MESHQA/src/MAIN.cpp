/*******************************************************************************
 * File:        MAIN.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef MESHQA_SHARED_LIB

/* Name of package */
#define PACKAGE "MESHQA"
/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "<ashish-gupta@utc.edu>"
/* Define to the full name of this package. */
#define PACKAGE_NAME ""
/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""
/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""
/* Define to the version of this package. */
#define PACKAGE_VERSION ""
/* Version number of package */
#define VERSION "1.0"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MESHQA.h"
#include "DBCGNS.h"

/* Command line options */
static char options[] = "ctsbhv";

static char *usgmsg[] = {
    "usage: meshqa [OPTIONS]... [FILE]",
    "options:",
    "	-c		= Mesh Quality Assissment for CGNS Data Format",
    "	-t		= Mesh Quality Assissment for NetCDF Data Format",
    "	-s		= Solution File (optional)",
    "	-b		= Boundary Map File (optional)",
    "	-h		= Print this message",
    "	-v		= Print Version",
    NULL
};

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static void usage(char **usgmsg) {
    int n;
    
    for (n = 0; NULL != usgmsg[n]; n++)
        fprintf(stderr, "%s\n", usgmsg[n]);
}

/*---------- arguments ---------------------------------------------------
 * Get option letter from argument vector or terminates on error
 *----------------------------------------------------------------------*/

static int arguments(int argc, char **argv) {
    unsigned int i, n, m;

    if (argc < 2) {
        fprintf(stderr, "ERROR: Too few arguments\n");
        usage(usgmsg);
        exit(1);
    }

    if (argc >= 2) {
        /* Check if help or version is enquired */
        if (argv[1][1] == options[4]) {
            usage(usgmsg);
            exit(0);
        } else if(argv[1][1] == options[5]) {
            printf("%s Utility, Version %s \n", PACKAGE, VERSION);
            printf("Copyright (C) 2010 Ashish Gupta. All rights reserved.\n");
            printf("Contact for Help or Bugs %s \n", PACKAGE_BUGREPORT);
            exit(0);
        } else if (argc < 3) {
            fprintf(stderr, "ERROR: Invalid Usage or No Input File\n");
            usage(usgmsg);
            exit(1);
        }
    }
    /* Check for options inputs */
    if (argc >= 3) {
        /* Check if output option is primary */
        if ((argc == 4) || (argc == 6) || (argv[1][1] == options[2]) || (argv[1][1] == options[3])) {
            fprintf(stderr, "ERROR: Invalid Usage\n");
            fprintf(stderr, "Example 1 = %s -[option] grid -s solution\n", PACKAGE);
            fprintf(stderr, "Example 2 = %s -[option] grid -s solution -b bmap\n", PACKAGE);
            fprintf(stderr, "Example 3 = %s -[option] grid\n", PACKAGE);
            exit(1);
        }
        /* Check if all arguments are proper */
        for (m = 1; m < (unsigned int)argc; m=m+2) {
            i = 0;
            for (n = 0; n < strlen(options); n++) {
                if (argv[m][1] == options[n]) {
                    i++;
                }
            }
            if (i == 0) {
                fprintf(stderr, "ERROR: Invalid option\n");
                usage(usgmsg);
                exit(1);
            }
        }
        /* Get the final option */
        for (n = 0; n < strlen(options); n++) {
            if (argv[1][1] == options[n]) {
                if ((n == 0) || (n ==1)) {
                    if (argc <= 4)
                        fprintf(stderr, "WARNING: Output File Name Missing, defaults: Input Name \n");
                }
                return n;
            }
        }
        fprintf(stderr, "ERROR: Invalid option\n");
        usage(usgmsg);
        exit(1);
    }
    return 0;
}


/*
 * 
 */
int main(int argc, char** argv) {
    int i, opt;
    char fnameg[257], fnames[257], fnameb[257];

    /* Get the argument options */
    opt = arguments(argc, argv);

    /* Initialize */
    for (i = 0; i < 257; i++) {
        fnameg[i] = ' ';
        fnames[i] = ' ';
        fnameb[i] = ' ';
    }
    fnameg[256] = '\0';
    fnames[256] = '\0';
    fnameb[256] = '\0';

    /* Get the input grid file name */
    strcpy(fnameg, argv[2]);
    /* Get the input solution/bmap file name */
    if (argc > 4) {
        for (i = 3; i < argc; i++) {
            if (argv[i][1] == options[2]) strcpy(fnames, argv[i+1]);
            if (argv[i][1] == options[3]) strcpy(fnameb, argv[i+1]);
            i++;
        }
    }

    DBCGNS object;
    MESHQARoot qualobject;
    /* Selecting the module */
    switch (opt) {
        case 0:
            /* CGNS File Format */
            printf("Reading Input CGNS File: %s \n", fnameg);
            object.Set_InputGrid_Filename(fnameg);
            object.Read_DB();
            if (qualobject.set_CFDDatabase(object.Get_DB())) {
            qualobject.update();
                printf("Start: Initialize MESH Quality Data Structure \n");
                qualobject.initialize();
                printf("Start: Analyze MESH Quality \n");
                qualobject.analyze();
            }
            break;
        case 1:
            /* TAU NetCDF File Format */
            printf("To Do\n");
            break;
    }
    
    return (EXIT_SUCCESS);
}

#else
#define MESHQA_VERSION 1.0
#endif
