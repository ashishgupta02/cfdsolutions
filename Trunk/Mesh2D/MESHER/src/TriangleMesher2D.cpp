/*******************************************************************************
 * File:        TriangleMesher2D.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include "TriangleMesh.h"
#include "Utils.h"

// *****************************************************************************
// *****************************************************************************
int TriangleMesher2D(int argcs, char* pArgs[]) {
    int b, i, j, nn, nt, nb, tdim, n0, n1, n2;
    const int bdim = 132;
    char buff[bdim];
    char *cdum;
    FILE *fp;
    int *nbs, ***bs;
    int (*tri)[3];
    double *x, *y;
    time_t tm;
    clock_t time_0, time_1;
    char *t_char;

    printf("\n====================================================");
    printf("\n     Driver for 2d meshing module                   ");
    printf("\n====================================================\n");

    if (--argcs < 2) {
        printf("\nNo input points file specified!\n");
        exit(0);
    }

    if ((fp = fopen(pArgs[argcs], "r")) == 0) {
        printf("\nCouldn't open points file <%s>\n", pArgs[argcs]);
        exit(0);
    }

    // read number of nodes
    cdum = fgets(buff, bdim, fp);
    cdum = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &nn);
    printf("\nNumber of points = %d", nn);
    fflush(stdout);

    // allocate for coordinates
    x = new double[nn];
    y = new double[nn];

    // read in coordinates
    for (i = 0; i < nn; i++) {
        cdum = fgets(buff, bdim, fp);
        sscanf(buff, "%lf %lf", &(x[i]), &(y[i]));
    }

    // allocate for points per boundary
    cdum = fgets(buff, bdim, fp);
    cdum = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &nb);
    printf("\nNumber of boundaries = %d", nb);
    nbs = new int[nb];
    bs = (int***) malloc(nb * sizeof (int**));
    for (b = 0; b < nb; b++) {
        cdum = fgets(buff, bdim, fp);
        cdum = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &nbs[b]);

        printf("\nBoundary %d has %d segments", b, nbs[b]);

        bs[b] = (int**) malloc(nbs[b] * sizeof (int*));
        for (i = 0; i < nbs[b]; i++) {
            bs[b][i] = (int*) malloc(2 * sizeof (int));
            cdum = fgets(buff, bdim, fp);
            sscanf(buff, "%d %d", &bs[b][i][0], &bs[b][i][1]);
            bs[b][i][0]--;
            bs[b][i][1]--;
        }
    }
    
    // initial guess for number of triangles
    tdim = nn * 3;
    tri = new int[tdim][3];

    time(&tm);
    time_0 = clock();

    t_char = ctime(&tm);
    printf("\nTriangulation started at %s", t_char);
    
    // iterative loop to allocate space for triangle indices
    nt = -1;
    while (nt < 0) {
        nt = trimesh(nn, tdim, nb, nbs, bs, x, y, tri);

        if (nt < 0) {
            printf("\nExpanding triangle dimension from %d to %d", tdim, tdim + nn);
            fflush(stdout);
            // delete current memory for triangle
            delete[] tri;

            // increment triangle dimension by number of nodes
            tdim += nn;
            // allocate new space for triangles
            tri = new int[tdim][3];
        }
    }
    printf("\nNumber of triangles generated = %d\n", nt);
    fflush(stdout);

    time(&tm);
    t_char = ctime(&tm);
    printf("\nTriangulation completed at %s", t_char);
    time_1 = clock();
    printf("\nTotal time in seconds = %14.7e", (float) (time_1 - time_0) / CLOCKS_PER_SEC);

    // output resulting triangular mesh
    char filename[32];

    // output trimesh in generic format
    filename[0] = '\0';
    strcat(filename, "Triangle.mesh");
    printf("\nFilename = <%s>", filename);
    // Open file for write
    if ((fp = fopen(filename, "w")) == 0) {
        printf("\nError opening file <%s>.", filename);
        exit(0);
    }

    // Write out nodes
    fprintf(fp, "#Number of grid points\n");
    fprintf(fp, "%d\n", nn);
    for (i = 0; i < nn; i++)
        fprintf(fp, "%19.10e %19.10e\n", x[i], y[i]);

    fprintf(fp, "#Number of blocks\n");
    fprintf(fp, "1\n");

    fprintf(fp, "#Number of triangular elements\n");
    fprintf(fp, "%d\n", nt);
    for (i = 0; i < nt; i++)
        fprintf(fp, "%d %d %d\n", tri[i][0] + 1, tri[i][1] + 1, tri[i][2] + 1);

    fprintf(fp, "#Number of quadrilateral elements\n");
    fprintf(fp, "%d\n", 0);

    fprintf(fp, "#Number of boundaries\n");
    fprintf(fp, "%d\n", nb);

    for (i = 0; i < nb; i++) {
        fprintf(fp, "#Number of edges for boundary %d\n", i + 1);
        fprintf(fp, "%d\n", nbs[i]);

        for (j = 0; j < nbs[i]; j++)
            fprintf(fp, "%d %d\n", bs[i][j][0] + 1, bs[i][j][1] + 1);
    }

    fclose(fp);

    //
    // write GNUPLOT file
    //
    buff[0] = '\0';
    strcat(buff, "Triangle.dat");
    printf("\nFilename = <%s>\n", buff);
    // Open file for write
    if ((fp = fopen(buff, "w")) == 0) {
        printf("\nError opening file <%s>.", buff);
        exit(0);
    }
    for (i = 0; i < nt; i++) {
        n0 = tri[i][0];
        n1 = tri[i][1];
        n2 = tri[i][2];
        fprintf(fp, "%19.10e %19.10e 0.0\n", x[n0], y[n0]);
        fprintf(fp, "%19.10e %19.10e 0.0\n", x[n1], y[n1]);
        fprintf(fp, "%19.10e %19.10e 0.0\n", x[n2], y[n2]);
        fprintf(fp, "%19.10e %19.10e 0.0\n\n", x[n0], y[n0]);
    }
    fclose(fp);

    delete[] tri;
    delete[] bs;

    delete[] x;
    delete[] y;

    return 0;
}

