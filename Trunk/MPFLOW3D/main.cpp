/*******************************************************************************
 * File:        main.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
/* Name of package */
#define PACKAGE "MPFLOW3D"
/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ashishgupta02@gmail.com"
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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* Custom header files */
#include "Trim_Utils.h"
#include "Commons.h"
#include "MeshIO.h"
#include "Solver.h"

/* Command line options */
static char options[] = "-hv";

static char *usgmsg[] = {
    "usage: mpflow3d [PARAMFILE]",
    "options:",
    " -h = Print this message",
    " -v = Print Version",
    NULL
};

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static void usage(char **usgstr) {
    int n;

    for (n = 0; NULL != usgstr[n]; n++)
        fprintf(stderr, "%s\n", usgstr[n]);
}

/*---------- arguments ---------------------------------------------------
 * Get option letter from argument vector or terminates on error
 *----------------------------------------------------------------------*/

static int arguments(int argc, char **argv) {
    int rcode = EXIT_SUCCESS;
    
    if (argc < 2) {
        fprintf(stderr, "ERROR: Too few arguments\n");
        usage(usgmsg);
        rcode = EXIT_FAILURE;
    }

    if (argc >= 2) {
        /* Check if help or version is enquired */
        if (argv[1][0] == options[0] && argv[1][1] == options[1]) {
            usage(usgmsg);
            rcode = EXIT_FAILURE;
        } else if(argv[1][0] == options[0] && argv[1][1] == options[2]) {
            printf("%s Utility, Version %s \n", PACKAGE, VERSION);
            printf("Copyright (C) 2010-14 Ashish Gupta. All rights reserved.\n");
            printf("Contact for Help or Bugs %s \n", PACKAGE_BUGREPORT);
            rcode = EXIT_FAILURE;
        }
    }
    return rcode;
}

// *****************************************************************************
// *****************************************************************************
int main(int argc, char *argv[]) {
#ifdef HAVE_MPI
    int rank, nproc, ROOT;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    // Make process with Rank = 0 as Root Process
    ROOT = 0;
    
    // Run Solver on Root Process Only
    if (rank == ROOT) {
#endif
        // Check the arguments
        if (arguments(argc, argv) == EXIT_SUCCESS) {
            // Initialize the Common Data
            Commons_Init();

            // Initialize the Solver Parameters Data
            Solver_Parameters_Init();

            // Read the Solver Parameters
            Solver_Parameters_Read(argv[1]);

            // Read input Grid File
            UGrid_Reader(MeshInputFilename);

            // Create All Connectivity Maps
            Create_Connectivity_Maps(MeshReorder);

            // Calculate Areas and Control Volumes
            Calculate_Area_Volume();

            // Free Excess Memory Used for Connectivity Creation
            Trim_Connectivity_Memory();

            // Read Boundary Conditions
            Solver_BC_Parameters_Read(BCInputFilename);

            // Initialize the Solver Data
            Solver_Init();

            // Set Initial Condition for Solver
            Solver_Set_Initial_Conditions();

            // Solve and Write Output
            if (Solver() == EXIT_SUCCESS)
                VTK_Writer(SolutionOutputFilename, 1);

            // Finalize the Solver Data
            Solver_Finalize();

            // Finalize the Common Data
            Commons_Finalize();
        }
#ifdef HAVE_MPI
    }
    // Synchronize All Process
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}

