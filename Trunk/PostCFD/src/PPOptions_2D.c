/*[
 * Copyright 2006   Ashish Gupta
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
 * with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

/*
 * File		PPOptions_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <GL/glut.h>
#include "PPOptions_2D.h"
#include "ColorPlot_2D.h"
#include "ContourPlot_2D.h"
#include "VectorPlot_2D.h"
#include "Graphics_2D.h"
#include "PPInitializeSolution_2D.h"
#include "DrawMesh_2D.h"
#include "Stream_2D.h"
#include "Streamline_2D.h"
#include "VectorTopology_2D.h"
#include "ScalarTopology_2D.h"
#include "FluidProperties.h"
#include "ReferenceQuantities.h"
#include "PPScalarFlowVariables_2D.h"

int ColorPlotFlag = 0;
int ColorPlotArg;

int ContourPlotFlag = 0;
int ContourPlotArg;
int ContourLevel = 10;

int VectorPlotFlag = 0;
int VectorPlotArg1, VectorPlotArg2;

int StreamlineFlag = 0;

int MeshScalarArg;

int MeshActivateFlag = 0;

int LegendActivateFlag = 0;

int BoundaryActivateFlag = 0;

/*--------------------------------------------------------------*/

/* Mesh Menu Options */
static void MeshMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Activate */
            MeshActivateFlag = 1;
            printf("Mesh Activated\n");
            break;
        case 2:
            /* DeActivate */
            MeshActivateFlag = 0;
            printf("Mesh DeActivated\n");
            break;
    }
}

/*--------------------------------------------------------------*/

/* Boundary Menu Options */
static void BoundaryMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Activate */
            BoundaryActivateFlag = 1;
            printf("Boundary Activated\n");
            break;
        case 2:
            /* DeActivate */
            BoundaryActivateFlag = 0;
            printf("Boundary DeActivated\n");
            break;
    }
}

/*--------------------------------------------------------------*/

/* Legend Menu Options */
static void LegendMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Activate */
            LegendActivateFlag = 1;
            printf("Legend Activated\n");
            break;
        case 2:
            /* DeActivate */
            LegendActivateFlag = 0;
            printf("Legend DeActivated\n");
            break;
    }
}

/*--------------------------------------------------------------*/

/* Tensor Field Option Menu */
static void TensorFieldMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Velocity Gradient Tensor */
            break;
    }
}

/*--------------------------------------------------------------*/

/* Vector Field Option Menu */
static void VectorFieldMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Velocity */
            break;
        case 2:
            /* Momentum */
            break;
        case 3:
            /* Vorticity */
            break;
        case 4:
            /* Perturbation Velocity */
            break;
        case 5:
            /* Density Gradient */
            break;
        case 6:
            /* Pressure Gradient */
            break;
        case 7:
            /* Temperature Gradient */
            break;
    }
}

/*--------------------------------------------------------------*/

/* Scalar Field Option Menu */
static void ScalarFieldMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Density */
            CalculateDensity_2D();
            break;
        case 2:
            /* Pressure */
            CalculatePressure_2D();
            break;
        case 3:
            /* Temperature */
            CalculateTemperature_2D();
            break;
        case 9:
            /* Coefficient of Pressure */
            CalculatePressureCoefficient_2D();
            break;
        case 22:
            /* Velocity Magnitude */
            CalculateVelocityMagnitude_2D();
            break;
    }
}

/*--------------------------------------------------------------*/

/* Color Plot Option Menu */
static void ColorPlotMenuCallBacks_2D(int option) {
    option = 0;

    ColorPlotFlag = 1;
    ContourPlotFlag = 0;
    VectorPlotFlag = 0;
    StreamlineFlag = 0;

    GetSolutionList_2D();
    printf("Select Scalar = ");
    scanf("%d", &ColorPlotArg);

    glutDisplayFunc(ColorPlot_2D);
}

/*--------------------------------------------------------------*/

/* Contour Level Option Menu */
static void ContourLevelMenuCallBacks_2D(int option) {
    ContourLevel = (option * 10);

    if (option == 11) {
        printf("Contour Levels (>=10) = ");
        scanf("%d", &ContourLevel);
    }
}

/*--------------------------------------------------------------*/

/* Contour Plot Option Menu */
static void ContourPlotMenuCallBacks_2D(int option) {
    option = 0;
    ColorPlotFlag = 0;
    ContourPlotFlag = 1;
    VectorPlotFlag = 0;
    StreamlineFlag = 0;

    GetSolutionList_2D();
    printf("Contour Scalar = ");
    scanf("%d", &ContourPlotArg);

    glutDisplayFunc(ContourPlot_2D);
}

/*--------------------------------------------------------------*/

/* Vector Plot Option Menu */
static void VectorPlotMenuCallBacks_2D(int option) {
    option = 0;
    ColorPlotFlag = 0;
    ContourPlotFlag = 0;
    VectorPlotFlag = 1;
    StreamlineFlag = 0;

    /* Set Vector */
    GetSolutionList_2D();
    printf("VectorX = ");
    scanf("%d", &VectorPlotArg1);
    printf("VectorY = ");
    scanf("%d", &VectorPlotArg2);

    glutDisplayFunc(VectorPlot_2D);
}

/*--------------------------------------------------------------*/

/* Post Processing Option Sub Menu */
static void PPMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Fluid Properties */
            GetFluidProperties();
            break;
        case 2:
            /* Reference Quantities */
            GetReferenceQuantities();
            break;
    }
}

/*--------------------------------------------------------------*/

/* Post Processing Option Sub Menu */
static void StreamSeedMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Set Stream Seed Points */
            printf("Give No of Seed Points = ");
            scanf("%d", &NSeedPoint);
            StreamlineFlag = 1;
            glutDisplayFunc(DrawMesh_2D);
            GetSeedPoints_2D();
            break;
        case 2:
            /* Random Seeding */
            SeedRandom = 1;
            SeedLine = 0;
            break;
        case 3:
            /* Shade Line */
            SeedLine = 1;
            SeedRandom = 0;
            break;
    }
}

/*--------------------------------------------------------------*/

/* Post Processing Option Sub Menu */
static void StreamMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Display Seed Points */
            glutDisplayFunc(DrawSeedPoint_2D);
            break;
        case 2:
            /* Set Vectors */
            GetSolutionList_2D();
            printf("VectorX = ");
            scanf("%d", &StreamlineArg1);
            printf("VectorY = ");
            scanf("%d", &StreamlineArg2);
            break;
        case 3:
            /* Set Stream Time */
            printf("Give Stream Time Limit = ");
            scanf("%lf", &StreamTime);
            break;
        case 4:
            /* Generate Streamlines */
            /* Unset Stream Flag */
            StreamlineFlag = 0;
            StreamModule = 1;
            glutDisplayFunc(Streamline_2D);
            break;
    }
}

/*--------------------------------------------------------------*/

/* Scalar Topology Option Sub Menu */
static void ScalarTopologyMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Set Scalar */
            GetSolutionList_2D();
            printf("Scalar = ");
            scanf("%d", &ScalarTopologyArg);
            break;
        case 2:
            /* Set Stream Time */
            printf("Give Stream Time Limit = ");
            scanf("%lf", &StreamTime);
            break;
        case 3:
            /* Get Critical Points */
            StreamModule = 3;
            ScalarTopology_2D();
            break;
        case 4:
            /* Generate Streamlines */
            StreamModule = 2;
            glutDisplayFunc(Streamline_2D);
            break;
    }
}

/*--------------------------------------------------------------*/

/* Vector Topology Option Sub Menu */
static void VectorTopologyMenuCallBacks_2D(int option) {
    switch (option) {
        case 1:
            /* Set Vector */
            GetSolutionList_2D();
            printf("VectorX = ");
            scanf("%d", &VectorTopologyArg1);
            printf("VectorY = ");
            scanf("%d", &VectorTopologyArg2);
            break;
        case 2:
            /* Set Stream Time */
            printf("Give Stream Time Limit = ");
            scanf("%lf", &StreamTime);
            break;
        case 3:
            /* Get Critical Points */
            StreamModule = 2;
            VectorTopology_2D();
            break;
        case 4:
            /* Generate Streamlines */
            StreamModule = 2;
            glutDisplayFunc(Streamline_2D);
            break;
    }
}

/*--------------------------------------------------------------*/

/* Topology Option Sub Menu */
static void TopologyMenuCallBacks_2D(int option) {
    option++;
}

/*--------------------------------------------------------------*/

/* create a Main Menu */
static void MainMenuCallBacks_2D(int option) {

    switch (option) {
        case 1:
            /* Display FullScreen */
            glutFullScreen();
            glutPostRedisplay();
            break;
        case 2:
            /* Exit FullScreen */
            glutSetWindow(MainWindowID);
            glutReshapeWindow(600, 600);
            glutPositionWindow(winPositionX, winPositionY);
            glutPostRedisplay();
            break;
        case 3:
            /* Draw Grid */
            glutDisplayFunc(DrawMesh_2D);
            break;
        case 4:
            /* Draw Color Mesh */
            /* Set Scalar */
            GetSolutionList_2D();
            printf("Scalar = ");
            scanf("%d", &MeshScalarArg);
            glutDisplayFunc(MeshScalar_2D);
            break;
        case 5:
            /* Reset Display */
            glutDisplayFunc(DrawInitial_2D);
            break;
        case 6:
            /* About */
            break;
        case 7:
            /* Exit */
            glutDestroyWindow(MainWindowID);
            exit(0);
            break;
    }
}

/*--------------------------------------------------------------*/

/* Main Menu Drawing */
void MainMenu_2D(void) {
    int Menu;
    int PPMenu;
    int MeshMenu;
    int BoundaryMenu;
    int LegendMenu;
    /* Color Related */
    int ColorPlotMenu;
    /* Contour Related */
    int ContourLevelMenu, ContourPlotMenu;
    /* Vector Plot Related */
    int VectorPlotMenu;
    int ScalarFieldMenu;
    int VectorFieldMenu;
    int TensorFieldMenu;
    int StreamMenu;
    int StreamSeedMenu;
    /* Topology Related */
    int TopologyMenu, VectorTopologyMenu, ScalarTopologyMenu;

    /* Mesh Menu */
    MeshMenu = glutCreateMenu(MeshMenuCallBacks_2D);
    glutSetMenu(MeshMenu);
    glutAddMenuEntry("Activate Mesh", 1);
    glutAddMenuEntry("DeActivate Mesh", 2);

    /* Boundary Menu */
    BoundaryMenu = glutCreateMenu(BoundaryMenuCallBacks_2D);
    glutAddMenuEntry("Activate Boundary", 1);
    glutAddMenuEntry("DeActivate Boundary", 2);

    /* Legend Menu */
    LegendMenu = glutCreateMenu(LegendMenuCallBacks_2D);
    glutSetMenu(LegendMenu);
    glutAddMenuEntry("Activate Legend", 1);
    glutAddMenuEntry("DeActivate Legend", 2);

    /* Color Plot Menu */
    ColorPlotMenu = glutCreateMenu(ColorPlotMenuCallBacks_2D);
    glutSetMenu(ColorPlotMenu);
    glutAddMenuEntry("Select Scalar", 1);

    /* Contour Level Menu */
    ContourLevelMenu = glutCreateMenu(ContourLevelMenuCallBacks_2D);
    glutSetMenu(ContourLevelMenu);
    glutAddMenuEntry("10 Contours", 1);
    glutAddMenuEntry("20 Contours", 2);
    glutAddMenuEntry("30 Contours", 3);
    glutAddMenuEntry("40 Contours", 4);
    glutAddMenuEntry("50 Contours", 5);
    glutAddMenuEntry("60 Contours", 6);
    glutAddMenuEntry("70 Contours", 7);
    glutAddMenuEntry("80 Contours", 8);
    glutAddMenuEntry("90 Contours", 9);
    glutAddMenuEntry("100 Contours", 10);
    glutAddMenuEntry("User Defined", 11);

    /* Contour Plot Menu */
    ContourPlotMenu = glutCreateMenu(ContourPlotMenuCallBacks_2D);
    glutSetMenu(ContourPlotMenu);
    glutAddSubMenu("Contour Levels", ContourLevelMenu);
    glutAddMenuEntry("Select Scalar", 1);

    /* Vector Plot Menu */
    VectorPlotMenu = glutCreateMenu(VectorPlotMenuCallBacks_2D);
    glutSetMenu(VectorPlotMenu);
    glutAddMenuEntry("Set Vectors", 1);

    /* Scalar Field Menu */
    ScalarFieldMenu = glutCreateMenu(ScalarFieldMenuCallBacks_2D);
    glutSetMenu(ScalarFieldMenu);
    glutAddMenuEntry("Density", 1);
    glutAddMenuEntry("Pressure", 2);
    glutAddMenuEntry("Temperature", 3);
    glutAddMenuEntry("Speed of Sound", 4);
    glutAddMenuEntry("Mach Number", 5);
    glutAddMenuEntry("Stagnation Density", 6);
    glutAddMenuEntry("Stagnation Pressure", 7);
    glutAddMenuEntry("Stagnation Temperature", 8);
    glutAddMenuEntry("Pressure Coefficient", 9);
    glutAddMenuEntry("Stagnation Pressure Coefficient", 10);
    glutAddMenuEntry("Pitot Pressure", 11);
    glutAddMenuEntry("Pitot Pressure Ratio", 12);
    glutAddMenuEntry("Dynamic Pressure", 13);
    glutAddMenuEntry("Enthalpy", 14);
    glutAddMenuEntry("Stagnation Enthalpy", 15);
    glutAddMenuEntry("Internal Energy", 16);
    glutAddMenuEntry("Stagnation Energy", 17);
    glutAddMenuEntry("Stagnation Energy Density", 18);
    glutAddMenuEntry("Kinetic Energy", 19);
    glutAddMenuEntry("VelocityX", 20);
    glutAddMenuEntry("VelocityY", 21);
    glutAddMenuEntry("Velocity Magnitude", 22);
    glutAddMenuEntry("Equivalent Potential Velocity Ratio", 23);
    glutAddMenuEntry("MomentumX", 24);
    glutAddMenuEntry("MomentumY", 25);
    glutAddMenuEntry("Momentum Magnitude", 26);
    glutAddMenuEntry("Entropy", 27);
    glutAddMenuEntry("Entropy Measure S1", 28);
    glutAddMenuEntry("Swirl", 29);
    glutAddMenuEntry("Helicity", 30);
    glutAddMenuEntry("Relative Helicity", 31);
    glutAddMenuEntry("Filtered Relative Helicity", 32);
    glutAddMenuEntry("Shock", 33);
    glutAddMenuEntry("Filtered Shock", 34);
    glutAddMenuEntry("Pressure Gradient Magnitude", 35);
    glutAddMenuEntry("Density Gradient Magnitude", 36);
    glutAddMenuEntry("Temperature Gradient Magnitude", 37);
    glutAddMenuEntry("Divergence of Velocity", 38);
    glutAddMenuEntry("Sutherlands Law", 39);
    glutAddMenuEntry("Isentropic Density Ratio", 40);
    glutAddMenuEntry("Isentropic Pressure Ratio", 41);
    glutAddMenuEntry("Isentropic Temperature Ratio", 42);

    /* Vector Field Menu */
    VectorFieldMenu = glutCreateMenu(VectorFieldMenuCallBacks_2D);
    glutSetMenu(VectorFieldMenu);
    glutAddMenuEntry("Velocity", 1);
    glutAddMenuEntry("Momentum", 2);
    glutAddMenuEntry("Vorticity", 3);
    glutAddMenuEntry("Perturbation Velocity", 4);
    glutAddMenuEntry("Density Gradient", 5);
    glutAddMenuEntry("Pressure Gradient", 6);
    glutAddMenuEntry("Temperature Gradient", 7);

    /* Tensor Field Menu */
    TensorFieldMenu = glutCreateMenu(TensorFieldMenuCallBacks_2D);
    glutSetMenu(TensorFieldMenu);
    glutAddMenuEntry("Velocity Gradient Tensor", 1);

    /* Post Processing Menu Entries */
    PPMenu = glutCreateMenu(PPMenuCallBacks_2D);
    glutSetMenu(PPMenu);
    glutAddMenuEntry("Fluid Properties", 1);
    glutAddMenuEntry("Reference Quantities", 2);
    glutAddSubMenu("Scalar Fields", ScalarFieldMenu);
    glutAddSubMenu("Vector Fields", VectorFieldMenu);
    glutAddSubMenu("Tensor Fields", TensorFieldMenu);


    /* Stream Seed Strategy Menu Entries */
    StreamSeedMenu = glutCreateMenu(StreamSeedMenuCallBacks_2D);
    glutSetMenu(StreamSeedMenu);
    glutAddMenuEntry("No of Seed Points", 1);
    glutAddMenuEntry("Random Seed", 2);
    glutAddMenuEntry("Shade Line", 3);

    /* Stream Menu Entries */
    StreamMenu = glutCreateMenu(StreamMenuCallBacks_2D);
    glutSetMenu(StreamMenu);
    glutAddSubMenu("Seeding Strategy", StreamSeedMenu);
    glutAddMenuEntry("Display Seed Points", 1);
    glutAddMenuEntry("Set Vectors", 2);
    glutAddMenuEntry("Set Stream Time", 3);
    glutAddMenuEntry("Plot Streamlines", 4);

    /* Scalar Topology Menu Entries */
    ScalarTopologyMenu = glutCreateMenu(ScalarTopologyMenuCallBacks_2D);
    glutSetMenu(ScalarTopologyMenu);
    glutAddMenuEntry("Set Scalar", 1);
    glutAddMenuEntry("Set Stream Time", 2);
    glutAddMenuEntry("Get Critical Points", 3);
    glutAddMenuEntry("Generate Streamlines", 4);

    /* Vector Topology Menu Entries */
    VectorTopologyMenu = glutCreateMenu(VectorTopologyMenuCallBacks_2D);
    glutSetMenu(VectorTopologyMenu);
    glutAddMenuEntry("Set Vector", 1);
    glutAddMenuEntry("Set Stream Time", 2);
    glutAddMenuEntry("Get Critical Points", 3);
    glutAddMenuEntry("Generate Streamlines", 4);

    /* Topology Menu Entries */
    TopologyMenu = glutCreateMenu(TopologyMenuCallBacks_2D);
    glutSetMenu(TopologyMenu);
    glutAddSubMenu("Vector Topology", VectorTopologyMenu);
    glutAddSubMenu("Scalar Topology", ScalarTopologyMenu);

    /* Main Menu Entries */
    Menu = glutCreateMenu(MainMenuCallBacks_2D);
    glutSetMenu(Menu);
    glutAddMenuEntry("Full Screen", 1);
    glutAddMenuEntry("Full Screen Exit", 2);
    /* Add Post-Processing Menu */
    glutAddSubMenu("Post-Analyzer", PPMenu);
    glutAddMenuEntry("Draw Grid", 3);
    glutAddMenuEntry("Draw Color Mesh", 4);
    glutAddSubMenu("Mesh Option", MeshMenu);
    glutAddSubMenu("Boundary Option", BoundaryMenu);
    glutAddSubMenu("Legend Option", LegendMenu);
    /* Add Color-Plot Menu */
    glutAddSubMenu("Color Plot", ColorPlotMenu);
    glutAddSubMenu("Contour Plot", ContourPlotMenu);
    glutAddSubMenu("Vector Plot", VectorPlotMenu);
    glutAddSubMenu("Stream Operation", StreamMenu);
    glutAddSubMenu("Topology Operation", TopologyMenu);
    glutAddMenuEntry("Reset", 5);
    glutAddMenuEntry("About", 6);
    glutAddMenuEntry("Exit", 7);
    /* Attach Mouse Button to Main Menu */
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}
