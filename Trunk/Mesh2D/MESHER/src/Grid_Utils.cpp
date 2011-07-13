/*******************************************************************************
 * File:        Grid_Utils.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

// For Double Min-Max
#include <float.h>
// For Int Min-Max
#include <climits>
#include <math.h>
#include "Utils.h"
#include "Grid_Utils.h"


// *****************************************************************************
// Circumcircle (Alternate View) http://mathworld.wolfram.com/Circumcircle.html
// Return: 1 => Point inside the circle
//         0 => Point ouside the circle
// *****************************************************************************
int Point_In_Circle(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
    double a, bx, by, c;
    double s1, s2, s3;
    double x0, y0, Radius, Distance;

    s1 = x1 * x1 + y1*y1;
    s2 = x2 * x2 + y2*y2;
    s3 = x3 * x3 + y3*y3;

    a  = x1 * (y2 - y3) - y1 * (x2 - x3) + (x2 * y3 - y2 * x3);
    bx = s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2);
    by = s1 * (x2 - x3) + s2 * (x3 - x1) + s3 * (x1 - x2);
    c  = s1 * (x2 * y3 - y2 * x3) + s2 * (y1 * x3 - x1 * y3) + s3 * (x1 * y2 - x2 * y1);

    if (!(a == 0.0)) {
        x0 = bx / (2.0 * a);
        y0 = -by / (2.0 * a);
        Radius = (sqrt(bx * bx + by * by + 4.0 * a * c)) / (2.0 * fabs(a));
        Distance = sqrt(((x0 - x4)*(x0 - x4)) + ((y0 - y4)*(y0 - y4)));
    } else {
        // Co-linear Points
        Radius = DBL_MAX;
        Distance = 0.0;
    }

    if ((Distance - Radius) < DBL_MIN)
        return 1;
    else
        return 0;
}

// *****************************************************************************
// Search Algorithm to find Point in a Triangle from List of Triangles
// Search Starts from Seed Triangle and Ends at the Until the Condition of
// Search is satisfied.
// Condition: Dot Product of Normalized Normal Vector and Vector from Mid Point
//            of Edge to the Point of Interest is less then Zero for "All Edges"
// Returns: Found Triangle ID
//        : -1 => Triangle not Found
// *****************************************************************************
int Point_In_Triangle_Mesh(int seed, double xt, double yt, double x[], double y[], int tri[][3], int nbr[][3]) {
    int i, n0, n1, n2, seed_next, limit;
    double Nx, Ny, Vx, Vy, Nmag;
    double x0, x1, y0, y1;
    double DotProduct[3], MaxDotProduct;

    // Initialize
    x0 = 0.0;
    y0 = 0.0;
    x1 = 0.0;
    y1 = 0.0;
    
    // Start the Search
    limit = 0;
    while (limit < INT_MAX) {
        if (seed < 0)
            return -1;

        // Get the Nodes of Triangle
        n0 = tri[seed][0];
        n1 = tri[seed][1];
        n2 = tri[seed][2];
        
        // Compute Dot Product for each Edge Normal Vector and
        // Vector from Mid Point of Edge to Point of Interest
        for (i = 0; i < 3; i++) {
            switch (i) {
                // Edge 0
                case 0:
                    x0 = x[n0];
                    y0 = y[n0];
                    x1 = x[n1];
                    y1 = y[n1];
                    break;
                // Edge 1
                case 1:
                    x0 = x[n1];
                    y0 = y[n1];
                    x1 = x[n2];
                    y1 = y[n2];
                    break;
                // Edge 2
                case 2:
                    x0 = x[n2];
                    y0 = y[n2];
                    x1 = x[n0];
                    y1 = y[n0];
                    break;
            }

            // Normalized Normal Vector of Edge
            Nx = (y1 - y0);
            Ny = (x0 - x1);
            Nmag = sqrt(Nx*Nx + Ny*Ny);
            Nx /= Nmag;
            Ny /= Nmag;
            
            // Vector from Mid Point of Edge to Point of Interest
            Vx = xt - (x0 + 0.5*(x1 - x0));
            Vy = yt - (y0 + 0.5*(y1 - y0));

            // Compute Dot Project b/w Normal Vector and Vector to Point of Interest
            DotProduct[i] = (Nx * Vx) + (Ny * Vy);
        }

        // Condition to Check if Point is inside Triangle
        if ((DotProduct[0] <= DBL_ZERO) &&
                (DotProduct[1] <= DBL_ZERO) &&
                (DotProduct[2] <= DBL_ZERO)) {
            return seed;
        }

        // Get the Next Candidate Triangle to search
        seed_next = seed;
        MaxDotProduct = MAX3(DotProduct[0], DotProduct[1], DotProduct[2]);
        if (DotProduct[0] == MaxDotProduct)
            seed_next = nbr[seed][0];
        if (DotProduct[1] == MaxDotProduct)
            seed_next = nbr[seed][1];
        if (DotProduct[2] == MaxDotProduct)
            seed_next = nbr[seed][2];
        seed = seed_next;
        limit++;
    }
    return -1;
}

// *****************************************************************************
// Computes the Minimum and Maximum of Coordinates
// *****************************************************************************
void Get_Domain_Extent(int nn, double x[], double y[],
        double *Xmin, double *Ymin, double *Xmax, double *Ymax) {
    int i;

    if (nn <= 0)
        return;

    *Xmin = *Ymin = DBL_MAX;
    *Xmax = *Ymax = DBL_MIN;

    for (i = 0; i < nn; i++) {
        *Xmin = MIN(x[i], *Xmin);
        *Ymin = MIN(y[i], *Ymin);
        *Xmax = MAX(x[i], *Xmax);
        *Ymax = MAX(y[i], *Ymax);
    }
}

// *****************************************************************************
// Creates the Rectangular Domain from X and Y Min, Max
// *****************************************************************************
void Create_BoundingBox_Domain(double x[], double y[],
        double Xmin, double Ymin, double Xmax, double Ymax) {
    double DeltaX, DeltaY;

    DeltaX = Xmax - Xmin;
    DeltaY = Ymax - Ymin;
    Xmax = Xmax + DeltaX / 4.0;
    Ymax = Ymax + DeltaY / 4.0;
    Xmin = Xmin - DeltaX / 4.0;
    Ymin = Ymin - DeltaY / 4.0;

    x[0] = Xmin;
    y[0] = Ymin;
    x[1] = Xmax;
    y[1] = Ymin;
    x[2] = Xmax;
    y[2] = Ymax;
    x[3] = Xmin;
    y[3] = Ymax;
}

// *****************************************************************************
// Updates the Triangular Mesh Node ID -4 to remove Bounding Box Nodes
// *****************************************************************************
void Delete_BoundingBox_Domain(int ntri, int tri[][3]) {
    int i, j;

    // Update the Nodes ID so that Bounding Box Nodes are removed
    for (i = 0; i < ntri; i++) {
        for (j = 0; j < 3; j++)
            tri[i][j] -= 4;
    }
}

// *****************************************************************************
// Generate Node to Cell (Triangle) Connectivity
// Brute force method to create Node to Cell Connectivity
// *****************************************************************************
void Create_Node2Cell_Connectivity(int nn, int ntri, int tri[][3], List* node2cell[]) {
    int i, nt;

    // Reset the list => Delete all connectivity Information if "Exits"
    for (i = 0; i < nn; i++)
        node2cell[i]->Redimension(0);

    // Now Add Triangles connected to the Node
    for (nt = 0; nt < ntri; nt++) {
        for (i = 0; i < 3; i++) {
            if (tri[nt][i] >= 0)
                node2cell[tri[nt][i]]->Check_List(nt);
        }
    }
}

// *****************************************************************************
// Generate Cell to Cell (Triangle to Triangle) Connectivity
// Brute force method to create Cell to Cell Connectivity
// *****************************************************************************
void Create_Cell2Cell_Connectivity(int ntri, int tri[][3], int cell2cell[][3], List* node2cell[]) {
    int nt, nd, nc;
    int n0, n1, cellid;

    // For Each Triangle Edge Search if they have Common Cell
    // Except the cell which whose Edge is being searched
    for (nt = 0; nt < ntri; nt++) {
        for (nd = 0; nd < 3; nd++) {
            // Initialize the Neighbour Cells to -1
            cell2cell[nt][nd] = -1;
            // Get the Edge Nodes of the Triangle
            n0 = tri[nt][nd];
            n1 = tri[nt][(nd + 1) % 3];
            // For all Node0 connected Triangles
            for (nc = 0; nc < node2cell[n0]->max; nc++) {
                // Get the triangle id of Node0 connected cell from list
                cellid = node2cell[n0]->list[nc];
                // Check if Node0 and Node1 have common triangles
                // Max of two cells can be common for an Edge and
                // Min of One Cell if Edge is on Boundary
                if ((nt != cellid) && (node2cell[n1]->Is_In_List(cellid))) {
                    cell2cell[nt][nd] = cellid;
                }
            }
        }
    }
}


// *****************************************************************************
// Locally Updates the Cell to Cell (Triangle to Triangle Connecitivity)
// First generates the global list of triangles connected to the all nodes of
// Newly generated triangle and Recreate the Cell 2 Cell Connectivity of
// Effected triangles
// *****************************************************************************
void Update_Cell2Cell_Connectivity(int tri[][3], int cell2cell[][3],
        List* node2cell[], List* new_tri) {
    int nt, nd, nc;
    int n0, n1, cellid, nbr_tri;

    // Create the list to store effected triangles
    List* effected_tri = new List();
    // Search around the newly created triangle nodes
    for (nt = 0; nt < new_tri->max; nt++) {
        cellid = new_tri->list[nt];
        // Add newly created triangle to effected triangle list
        // Since newly created triangles have no Cell to Cell connectivity
        effected_tri->Add_To_List(cellid);
        // Search for triangles connected to nodes of newly created triangle
        for (nd = 0; nd < 3; nd++) {
            n0 = tri[cellid][nd];
            for (nc = 0; nc < node2cell[n0]->max; nc++) {
                nbr_tri = node2cell[n0]->list[nc];
                if (cellid != nbr_tri) {
                    effected_tri->Add_To_List(nbr_tri);
                }
            }
        }
    }

    // Now new_tri has new triangles and effected_tri has triangle
    // Sharing a node with a new triangles
    for (nt = 0; nt < effected_tri->max; nt++) {
        cellid = effected_tri->list[nt];
        for (nd = 0; nd < 3; nd++) {
            // Initialize the Neighbour Cells to -1
            cell2cell[cellid][nd] = -1;
             // Get the Edge Nodes of the Triangle
            n0 = tri[cellid][nd];
            n1 = tri[cellid][(nd + 1) % 3];
            // For all Node0 connected Triangles
            for (nc = 0; nc < node2cell[n0]->max; nc++) {
                // Get the triangle id of Node0 connected cell from list
                nbr_tri = node2cell[n0]->list[nc];
                // Check if Node0 and Node1 have common triangles
                // Max of two cells can be common for an Edge and
                // Min of One Cell if Edge is on Boundary
                if (cellid != nbr_tri && node2cell[n1]->Is_In_List(nbr_tri)) {
                    cell2cell[cellid][nd] = nbr_tri;
                }
            }
        }
    }
    delete effected_tri;
}


// *****************************************************************************
// Given a Test Point Recursive Function to Find Neighbour Triangles
// Voilating Circle Test and Appends to the List of Triangle Voilating
// Circle test
// Returns: List of Voilated Triangles
// *****************************************************************************
void Find_Neighbour_Triangles_Violating_CircleTest(int tid, int pid, double x[], double y[],
        int tri[][3], int (*cell2cell)[3], List *violated_tri) {

    int i, violated, nbr_tri;
    int n0, n1, n2;
    double x0, y0, x1, y1, x2, y2;

    // Check if triangle already exists in the Delete Triangle List
    if (violated_tri->Is_In_List(tid))
        return;

    // Get the Nodes of the Triangle
    n0 = tri[tid][0];
    n1 = tri[tid][1];
    n2 = tri[tid][2];

    // Get the Coordinates of Nodes
    x0 = x[n0];
    y0 = y[n0];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];

    // Perform Circle Test on the triangle
    violated = Point_In_Circle(x0, y0, x1, y1, x2, y2, x[pid], y[pid]);
    if (violated) {
        violated_tri->Add_To_List(tid);
        for (i = 0; i < 3; i++) {
            nbr_tri = cell2cell[tid][i];
            // Check if Neighbour Triangle is not a Boundary
            if (nbr_tri != -1) {
                Find_Neighbour_Triangles_Violating_CircleTest(nbr_tri, pid, x, y, tri, cell2cell, violated_tri);
            }
        }
    }
}

// *****************************************************************************
// Flood Fill Algorithm to separate Interior and Exterior triangles and Removing
// the outside triangles. Since the Boundary segment is stored such that the
// while walking along the boundary region on the left is interior and
// region on right is Outside the domain.
// To identify the triangle Interior, for boundary edge Node0-Node1 if common
// triangle has node oritentation along the boundary edge direction then it is
// interior else it is exterior cell
// *****************************************************************************
int Flood_Fill_Algorithm(int nb, int nbs[], int ***bs,
        int &ntri, int tri[][3], int cell2cell[][3], List* node2cell[]) {
    int i, j, nt, b, count;
    int cellid, commontri[2];
    int Node0, Node1, Node2;
    int BndNode0, BndNode1;

    // Check if Boundary Exits
    if (nb <= 0)
        return 0;

    // Allocate the memory to mark internal and external triangles
    int *mark = NULL;
    mark = (int *) malloc(ntri*sizeof(int));
    if (mark == NULL) {
        printf("ERROR: Flood_Fill_Algorithm - Memory Allocation Failure\n!");
        exit(EXIT_FAILURE);
    }

    // Initailize all triangles to be Neutral
    for (nt = 0; nt < ntri; nt++)
        mark[nt] = 0;

    // For all boundaries
    for (b = 0; b < nb; b++) {
        for (i = 0; i < nbs[b]; i++) {
            BndNode0 = bs[b][i][0];
            BndNode1 = bs[b][i][1];

            // Identify the common triangles shared by the boundary edge
            count = 0;
            commontri[0] = commontri[1] = -1;
            for (j = 0; ((j < node2cell[BndNode0]->max) && (count < 2)); j++) {
                cellid = node2cell[BndNode0]->list[j];
                    if (node2cell[BndNode1]->Is_In_List(cellid)) {
                        commontri[count] = cellid;
                        count++;
                    }
            }
            // Mark the Cells Interior and Exterior
            for (j = 0; j < 2; j++) {
                Node0 = tri[commontri[j]][0];
                Node1 = tri[commontri[j]][1];
                Node2 = tri[commontri[j]][2];
                // Triangle Edge Nodes Oriented along Boundary Direction
                if ((Node0 == BndNode0 && Node1 == BndNode1) ||
                        (Node1 == BndNode0 && Node2 == BndNode1) ||
                        (Node2 == BndNode0 && Node0 == BndNode1))
                    mark[commontri[j]] = 1;
                // Triangle Edge Nodes Oriented Opposite Boundary Direction
                else
                    mark[commontri[j]] = -1;
            }
        }
    }

    // At this point only Triangles at the Boundary are tagged interior
    // and exterior. Now tag the netural cells inside or outside
    // Tag all Interior Cells to 1 and all Exterior Cells to -1
    int tagged;
    for (;;) {
        for (nt = 0; nt < ntri; nt++) {
            if (mark[nt] == 0)
                continue;
            for (i = 0; i < 3; i++) {
                cellid = cell2cell[nt][i];
                // No Boundary Cells with -1
                if (cellid < 0)
                    continue;
                // Make the neutral neibouring cell interior or exterior if not
                // marked already previously
                if (mark[cellid] == 0)
                    mark[cellid] = mark[nt];
            }
        }
        // Check if all cells are marked either 1 or -1
        tagged = 0;
        for (nt = 0; nt < ntri; nt++)
            if (mark[nt] == 0)
                tagged++;
        if (tagged == 0)
            break;
    }

    // Get the Number of Interior Triangles
    tagged = 0;
    for (nt = 0; nt < ntri; nt++) {
        if (mark[nt] > 0) {
            mark[nt] = tagged;
            tagged++;
        }
    }

    // Update the Valid Triangulation with Mark 1
    for (nt = 0; nt < ntri; nt++) {
        if (mark[nt] >= 0) {
            tri[mark[nt]][0] = tri[nt][0];
            tri[mark[nt]][1] = tri[nt][1];
            tri[mark[nt]][2] = tri[nt][2];
        }
    }
    ntri = tagged;

    if (mark != NULL)
        free(mark);

    return 0;
}

// *****************************************************************************
// Boundry Edge Recovery Algorithm
// *****************************************************************************
int Boundary_Edge_Recovery(int nn, int nb, int nbs[], int ***bs,
        double x[], double y[], int ntri, int tri[][3],
        int cell2cell[][3], List* node2cell[]) {

    int i, j, b;
    int valid_edge, cellid;
    int BndNode0, BndNode1;
    double Magnitude;
    double BndEdgeVectorX, BndEdgeVectorY;
    double EdgeNormal1[2], EdgeNormal2[2];
    double DotProduct1, DotProduct2;

    // Check if Boundary Exits
    if (nb <= 0)
        return 0;

    // For all boundaries
    for (b = 0; b < nb; b++) {
        // For all segments of a boundary
        for (int ib = 0; ib < nbs[b]; ib++) {
            BndNode0 = bs[b][ib][0];
            BndNode1 = bs[b][ib][1];

            // Check if both nodes shared a common triangle
            valid_edge = 0;
            for (j = 0; j < node2cell[BndNode0]->max; j++) {
                cellid = node2cell[BndNode0]->list[j];
                if (node2cell[BndNode1]->Is_In_List(cellid)) {
                    valid_edge = 1;
                    break;
                }
            }

            // If Boundary Nodes are is not an Edges
            if (!valid_edge) {
                // Compute the Vector of Desired Boundary Edge
                BndEdgeVectorX = x[BndNode1] - x[BndNode0];
                BndEdgeVectorY = y[BndNode1] - y[BndNode0];
                Magnitude = sqrt(BndEdgeVectorX * BndEdgeVectorX + BndEdgeVectorY * BndEdgeVectorY);
                BndEdgeVectorX /= Magnitude;
                BndEdgeVectorY /= Magnitude;

                int bnode, Node0, Node1, Node2, OppTri;
                Node0  = -1;
                Node1  = -1;
                Node2  = -1;
                bnode  = -1;
                OppTri = -1;
                // Swap Edges Until Desired Boundary Edge is Achieved
                while (Node0 != BndNode1) {
                    Node0 = -1;
                    // Loop to Get the Opposite Triangle Node
                    for (j = 0; j < node2cell[BndNode0]->max; j++) {
                        if (Node0 >= 0)
                            break;

                        cellid = node2cell[BndNode0]->list[j];

                        // Get the other two nodes of the triangle and Opposite Triangle
                        for (i = 0; i < 3; i++) {
                            switch (i) {
                                case 0:
                                    bnode = tri[cellid][0];
                                    Node1 = tri[cellid][1];
                                    Node2 = tri[cellid][2];
                                    OppTri = cell2cell[cellid][1];
                                    break;
                                case 1:
                                    bnode = tri[cellid][1];
                                    Node1 = tri[cellid][2];
                                    Node2 = tri[cellid][0];
                                    OppTri = cell2cell[cellid][2];
                                    break;
                                case 2:
                                    bnode = tri[cellid][2];
                                    Node1 = tri[cellid][0];
                                    Node2 = tri[cellid][1];
                                    OppTri = cell2cell[cellid][0];
                                    break;
                            }
                            if (bnode == BndNode0)
                                break;
                        }
                        // Check if Opposite Triangle is not a Boundary
                        if (OppTri < 0)
                            continue;

                        // Compute the Normals of the two Edges
                        // Normalized BndNode0-Node1 Normal
                        EdgeNormal1[0] = y[BndNode0] - y[Node1];
                        EdgeNormal1[1] = x[Node1] - x[BndNode0];
                        Magnitude = sqrt(EdgeNormal1[0]*EdgeNormal1[0] + EdgeNormal1[1]*EdgeNormal1[1]);
                        EdgeNormal1[0] /= Magnitude;
                        EdgeNormal1[1] /= Magnitude;
                        // Normalized BndNode0-Node2 Normal
                        EdgeNormal2[0] = y[Node2] - y[BndNode0];
                        EdgeNormal2[1] = x[BndNode0] - x[Node2];
                        Magnitude = sqrt(EdgeNormal2[0]*EdgeNormal2[0] + EdgeNormal2[1]*EdgeNormal2[1]);
                        EdgeNormal2[0] /= Magnitude;
                        EdgeNormal2[1] /= Magnitude;

                        // Compute the Dot Product of Normals to the Desired Boundary Edge Vector
                        DotProduct1 = BndEdgeVectorX*EdgeNormal1[0] + BndEdgeVectorY*EdgeNormal1[1];
                        DotProduct2 = BndEdgeVectorX*EdgeNormal2[0] + BndEdgeVectorY*EdgeNormal2[1];

                        // Check if the Desired Boundary Edge Vector Lie inside the triangle
                        // and Pass/Cuts through Edge Node1-Node2
                        if ((DotProduct1 > DBL_ZERO) && (DotProduct2 > DBL_ZERO)) {
                            // Now Get the Node0 Opposite to Node1 and Node2 of OppTri
                            for (i = 0; i < 3; i++) {
                                switch (i) {
                                    case 0:
                                        if ((tri[OppTri][0] == Node2) && (tri[OppTri][1] == Node1))
                                            Node0 = tri[OppTri][2];
                                        break;
                                    case 1:
                                        if ((tri[OppTri][1] == Node2) && (tri[OppTri][2] == Node1))
                                            Node0 = tri[OppTri][0];
                                        break;
                                    case 2:
                                        if ((tri[OppTri][2] == Node2) && (tri[OppTri][0] == Node1))
                                            Node0 = tri[OppTri][1];
                                        break;
                                }
                            }

                            // Swap the Edge: Node1-Node2 with BndNode0-Node0 and Delete Create triangles
                            // Delete triangle current
                            node2cell[BndNode0]->Delete_From_List(cellid);
                            node2cell[Node1]->Delete_From_List(cellid);
                            node2cell[Node2]->Delete_From_List(cellid);
                            // Delete triangle OppTri
                            node2cell[Node1]->Delete_From_List(OppTri);
                            node2cell[Node2]->Delete_From_List(OppTri);
                            node2cell[Node0]->Delete_From_List(OppTri);

                            // Update the Nodes of New triangles and Add to Node2Cell Connectivity
                            tri[cellid][0] = BndNode0;
                            tri[cellid][1] = Node1;
                            tri[cellid][2] = Node0;
                            node2cell[BndNode0]->Add_To_List(cellid);
                            node2cell[Node1]->Add_To_List(cellid);
                            node2cell[Node0]->Add_To_List(cellid);
                            tri[OppTri][0] = BndNode0;
                            tri[OppTri][1] = Node0;
                            tri[OppTri][2] = Node2;
                            node2cell[BndNode0]->Add_To_List(OppTri);
                            node2cell[Node0]->Add_To_List(OppTri);
                            node2cell[Node2]->Add_To_List(OppTri);

                            // Create Cell to Cell Connectivity - Brute Force Method
                            Create_Node2Cell_Connectivity(nn, ntri, tri, node2cell);
                            // Create Cell to Cell Connectivity - Brute Force Method
                            Create_Cell2Cell_Connectivity(ntri, tri, cell2cell, node2cell);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

