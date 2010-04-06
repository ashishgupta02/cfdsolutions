/*
 * File:   main.cpp
 * Author: Ashish Gupta
 * 
 * Created on February 5, 2010, 7:22 PM
 * Modified on April 6, 2010
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
/* Name of package */
#define PACKAGE "MESH2D"
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
#include "Utils.h"
#include "TriangleMesher2D.h"
#include "TriangleMeshOptimizer.h"
#include "TriangleWinslowSmoother.h"
#include "TriangleLinearElasticSmoother.h"
#include "TriangleWinslowVirtualVolumeSmoother.h"
#include "TriangleAdaptiveRefinement.h"

/* Command line options */
static char options[] = "abcdefhv";

static char *usgmsg[] = {
    "usage: mesh2d [OPTIONS]... [FILE]",
    "options:",
    " -a = Generate Daulany Triangulation",
    " -b = Mesh Optimization",
    " -c = Mesh Winslow Smoothing",
    " -d = Mesh Linear Elasticity Smoothing",
    " -e = Mesh Winslow Virtual Control Volume Smoothing",
    " -f = Feature Based Adaptive Refinement",
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
    
    /* Get the argument options */
    opt = arguments(argc, argv);

    // Check if file exits
    if (!file_exists(argv[2]))
        error("MESH2D: %s\n", "Specified Mesh File Doesn't Exists");
    
    /* Selecting the module */
    switch (opt) {
        case 0:
            /* Generate Daulany Triangulation */
            TriangleMesher2D(argc, argv);
            break;
        case 1:
        {
            /* Mesh Optimization */
            TriangleMeshOptimizer MeshOpti;
            MeshOpti.SLK_MeshReader(argv[2]);
            MeshOpti.Rotate_Boundary(0);
            MeshOpti.Optimize();
            MeshOpti.SLK_MeshWriter("Triangle.mesh");
            MeshOpti.SLK_GnuplotWriter("Triangle.dat");
        }
            break;
        case 2:
        {
            /* Mesh Winslow Smoothing */
            TriangleWinslowSmoother MeshWSmooth;
            MeshWSmooth.SimCenterMeshReader(argv[2]);
            MeshWSmooth.Rotate_Boundary(0);
            MeshWSmooth.WinslowSmooth();
            MeshWSmooth.SimCenterMeshWriter("Triangle.mesh");
            MeshWSmooth.GnuplotWriter("Triangle.dat");
        }
            break;
        case 3:
        {
            /* Linear Elasticity Smoothing */
            TriangleLinearElasticSmoother MeshLESmooth;
            MeshLESmooth.SLK_MeshReader(argv[2]);
            MeshLESmooth.Rotate_Boundary(0);
            MeshLESmooth.LESmooth();
            MeshLESmooth.SLK_MeshWriter("Triangle.mesh");
            MeshLESmooth.SLK_GnuplotWriter("Triangle.dat");
        }
            break;
        case 4:
        {
            /* Mesh Winslow Virtual Control Volume Smoothing */
            TriangleWinslowVirtualVolumeSmoother MeshWVVSmooth;
            MeshWVVSmooth.SimCenterMeshReader(argv[2]);
            MeshWVVSmooth.Rotate_Boundary(0);
            MeshWVVSmooth.WinslowVirtualVolumeSmooth();
            MeshWVVSmooth.SimCenterMeshWriter("Triangle.mesh");
            MeshWVVSmooth.GnuplotWriter("Triangle.dat");
        }
            break;
        case 5:
        {
            /* Feature Based Adaptive Refinement */
            TriangleAdaptiveRefinement MeshRefine;
            MeshRefine.SLK_MeshReader(argv[2]);
            MeshRefine.Get_Input_Parameters();
            MeshRefine.AdaptiveRefinement();
            MeshRefine.SLK_MeshWriter("Triangle.mesh");
            MeshRefine.SLK_GnuplotWriter("Triangle.dat");
            break;
        }
    }
    
    return EXIT_SUCCESS;
}

