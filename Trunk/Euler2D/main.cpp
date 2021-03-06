/*******************************************************************************
 * File:        main.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
/* Name of package */
#define PACKAGE "EULER2D"
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
#define VERSION "2.0"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Custom header files */
#include "Utils.h"
#include "Euler2D_Solver_VanLeer.h"
#include "Euler2D_Solver_StegerWarming.h"
#include "Euler2D_Solver_Roe.h"
#include "Euler2D_Solver_Osher.h"
#include "Euler2D_Solver_LDFSS.h"
#include "Euler2D_Solver_AUSM.h"
#include "Euler2D_Design_VanLeer.h"

/* Command line options */
static char options[] = "abcdefghv";

static char *usgmsg[] = {
    "usage: euler2d [OPTIONS]... [FILE]",
    "options:",
    " -a = Van Leer Flux Vector Splitting",
    " -b = Steger Warming Flux Vector Splitting",
    " -c = Advection Upstream Splitting Method",
    " -d = Low Diffusion Flux Splitting Scheme",
    " -e = Oshers Splitting",
    " -f = Roe Flux Splitting",
    " -g = Design using Van Leer Flux Vector Splitting",
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
        if (argv[1][1] == options[7]) {
            usage(usgmsg);
            exit(0);
        } else if(argv[1][1] == options[8]) {
            printf("%s Utility, Version %s \n", PACKAGE, VERSION);
            printf("Copyright (C) 2010-11 Ashish Gupta. All rights reserved.\n");
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

    /* Get the argument options */
    opt = arguments(argc, argv);

    // Check if file exits
    if (!file_exists(argv[2]))
        error("EULER2D: %s\n", "Specified Mesh File Doesn't Exists");
    
    /* Selecting the module */
    switch (opt) {
        case 0:
            // -a = Van Leer Flux Vector Splitting
        {
            Euler2D_Solver_VanLeer VanLeer;
            VanLeer.Get_Solver_Inputs_VanLeer(argv[2]);
            VanLeer.Solver_Prepare();
            VanLeer.Solve();
            VanLeer.Solver_Finalize();
        }
            break;
        case 1:
            // -b = Steger Warming Flux Vector Splitting
        {
            Euler2D_Solver_StegerWarming StegerWarming;
            StegerWarming.Get_Solver_Inputs_StegerWarming(argv[2]);
            StegerWarming.Solver_Prepare();
            StegerWarming.Solve();
            StegerWarming.Solver_Finalize();
        }
            break;
        case 2:
            // -c = Advection Upstream Splitting Method
         {
            Euler2D_Solver_AUSM AUSM;
            AUSM.Get_Solver_Inputs_AUSM(argv[2]);
            AUSM.Solver_Prepare();
            AUSM.Solve();
            AUSM.Solver_Finalize();
        }
            break;
        case 3:
            // -d = Low Diffusion Flux Splitting Scheme
        {
            Euler2D_Solver_LDFSS LDFSS;
            LDFSS.Get_Solver_Inputs_LDFSS(argv[2]);
            LDFSS.Solver_Prepare();
            LDFSS.Solve();
            LDFSS.Solver_Finalize();
        }
            break;
        case 4:
            // -e = Oshers Splitting
        {
            Euler2D_Solver_Osher Osher;
            Osher.Get_Solver_Inputs_Osher(argv[2]);
            Osher.Solver_Prepare();
            Osher.Solve();
            Osher.Solver_Finalize();
        }
            break;
        case 5:
            // -f = Roe Flux Splitting
        {
            Euler2D_Solver_Roe Roe;
            Roe.Get_Solver_Inputs_Roe(argv[2]);
            Roe.Solver_Prepare();
            Roe.Solve();
            Roe.Solver_Finalize();
        }
            break;
        case 6:
            // -g = Design using Van Leer Flux Vector Splitting
        {
            Euler2D_Design_VanLeer Design;
            Design.Get_Design_Inputs_VanLeer(argv[2]);
            Design.Design_Prepare();
            Design.Design();
            Design.Design_Finalize();
        }
            break;
    }

    return EXIT_SUCCESS;
}

