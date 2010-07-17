/*******************************************************************************
 * File:        main.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
/* Name of package */
#define PACKAGE "EULER3D"
/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ashish-gupta@utc.edu"
/* Define to the full name of this package. */
#define PACKAGE_NAME ""
/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""
/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""
/* Define to the version of this package. */
#define PACKAGE_VERSION ""
/* Version number of package */
#define VERSION "0.1"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Custom header files */
#include "Commons.h"
#include "MeshIO.h"

/* Command line options */
static char options[] = "abcdefhv";

static char *usgmsg[] = {
    "usage: euler3d [OPTIONS]... [FILE]",
    "options:",
    " -a = Van Leer Flux Vector Splitting",
    " -b = Steger Warming Flux Vector Splitting",
    " -c = Advection Upstream Splitting Method",
    " -d = Low Diffusion Flux Splitting Scheme",
    " -e = Oshers Splitting",
    " -f = Roe Flux Splitting",
    " -h = Print this message",
    " -v = Print Version",
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
    unsigned int n;

    if (argc < 2) {
        fprintf(stderr, "ERROR: Too few arguments\n");
        usage(usgmsg);
        exit(1);
    }

    if (argc >= 2) {
        /* Check if help or version is enquired */
        if (argv[1][1] == options[6]) {
            usage(usgmsg);
            exit(0);
        } else if(argv[1][1] == options[7]) {
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
        for (n = 0; n < strlen(options); n++)
            if (argv[1][1] == options[n])
                return n;
        fprintf(stderr, "ERROR: Invalid option\n");
        usage(usgmsg);
        exit(1);
    }
    return 0;
}

// *****************************************************************************
// *****************************************************************************
int main(int argc, char *argv[]) {
    int opt;

    // Get the argument options
    opt = arguments(argc, argv);

    // Initialize the Common Data
    Commons_Init();

    // Read input Grid File
    UGrid_Reader(argv[2]);

    // Create All Connectivity Maps
    int reorder = 0;
    Create_Connectivity_Maps(reorder);
    
    /* Selecting the module */
    switch (opt) {
        case 0:
            // -a = Van Leer Flux Vector Splitting
            break;
        case 1:
            // -b = Steger Warming Flux Vector Splitting
            break;
        case 2:
            // -c = Advection Upstream Splitting Method
            break;
        case 3:
            // -d = Low Diffusion Flux Splitting Scheme
            break;
        case 4:
            // -e = Oshers Splitting
            break;
        case 5:
            // -f = Roe Flux Splitting
            break;
    }

    // Finalize the Common Data
    Commons_Finalize();
    
    return EXIT_SUCCESS;
}

