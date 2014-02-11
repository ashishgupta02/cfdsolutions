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
 * File		PostProcessingInitialize_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "Algebra.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PrePostProcessing_2D.h"
#include "Boundary_2D.h"

/* Defination of externs of PostProcessing.h*/
int NoCells2D = 0;
int NoNodes2D = 0;
int NoEdges2D = 0;
int NoOfBoundaryEdges2D = 0;
Node_2D *Node2D;
Cell_2D *Cell2D;
Edge_2D *Edge2D;
Vector_2D *NodeAdjacentNodes2D;
Vector_2D *NodeAdjacentCells2D;
int MaxEdgeNum2D;
int *IA_2D;
int *JA_2D;

/*---------------------------------------------------------------*/
int compress_sort(int *A, int value, int size) {
    int i, j, count;
    int *data;

    count = 0;
    for (i = 0; i < size; i++) {
        if ((A[i] != A[i + 1])&&(A[i] != value))
            count += 1;
    }

    data = (int *) malloc(count * sizeof (int));
    if (data == NULL) {
        Warn("compress_sort: Memory Allocation failed for data");
        return 1;
    }

    j = 0;
    for (i = 0; i < size; i++) {
        if ((A[i] != A[i + 1])&&(A[i] != value)) {
            data[j] = A[i];
            j += 1;
        }
    }

    if (j != count) {
        Warn("compress_sort: Miss-Matched count and j");
        return 1;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int sort_bubble(int *A, int size) {
    int i, j, tmp;

    for (i = 0; i < size - 1; i++) {
        for (j = i + 1; j < size; j++) {
            if (A[i] > A[j]) {
                tmp = A[i];
                A[i] = A[j];
                A[j] = tmp;
            }
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetCellAdjacentCellsQuad_2D(int cellid) {
    int i, j, n1, n2, n3, n4, size, count;
    int *data;

    n1 = Cell2D[cellid].ConnectNode[0];
    n2 = Cell2D[cellid].ConnectNode[1];
    n3 = Cell2D[cellid].ConnectNode[2];
    n4 = Cell2D[cellid].ConnectNode[3];

    size = 0;
    size += NodeAdjacentCells2D[n1].Size;
    size += NodeAdjacentCells2D[n2].Size;
    size += NodeAdjacentCells2D[n3].Size;
    size += NodeAdjacentCells2D[n4].Size;

    data = (int *) malloc(size * sizeof (int));
    if (data == NULL) {
        Warn("GetCellAdjacentCellsQuad_2D: Memory Allocation Failed 1");
        return 1;
    }

    count = 0;

    for (i = 0; i < NodeAdjacentCells2D[n1].Size; i++) {
        data[count] = NodeAdjacentCells2D[n1].Data[i];
        count += 1;
    }

    for (i = 0; i < NodeAdjacentCells2D[n2].Size; i++) {
        data[count] = NodeAdjacentCells2D[n2].Data[i];
        count += 1;
    }

    for (i = 0; i < NodeAdjacentCells2D[n3].Size; i++) {
        data[count] = NodeAdjacentCells2D[n3].Data[i];
        count += 1;
    }

    for (i = 0; i < NodeAdjacentCells2D[n4].Size; i++) {
        data[count] = NodeAdjacentCells2D[n4].Data[i];
        count += 1;
    }

    sort_bubble(data, size);

    count = 0;
    for (i = 1; i < size; i++) {
        if ((data[i - 1] != data[i])&&(data[i - 1] != cellid))
            count += 1;
    }

    if (data[size - 1] != cellid)
        count += 1;

    Cell2D[cellid].NearCells = count;

    Cell2D[cellid].NearNeighbours = (int *) malloc(count * sizeof (int));
    if (Cell2D[cellid].NearNeighbours == NULL) {
        Warn("GetCellAdjacentCellsQuad_2D: Memory Allocation Failed 2");
        return 1;
    }

    j = 0;
    for (i = 1; i < size; i++) {
        if ((data[i - 1] != data[i])&&(data[i - 1] != cellid)) {
            Cell2D[cellid].NearNeighbours[j] = data[i - 1];
            j += 1;
        }
    }

    if (data[size - 1] != cellid)
        Cell2D[cellid].NearNeighbours[j] = data[size - 1];

    free(data);
    return 0;
}

/*---------------------------------------------------------------*/
int GetCellAdjacentCellsTri_2D(int cellid) {
    int i, j, n1, n2, n3, size, count;
    int *data;

    n1 = Cell2D[cellid].ConnectNode[0];
    n2 = Cell2D[cellid].ConnectNode[1];
    n3 = Cell2D[cellid].ConnectNode[2];

    size = 0;
    size += NodeAdjacentCells2D[n1].Size;
    size += NodeAdjacentCells2D[n2].Size;
    size += NodeAdjacentCells2D[n3].Size;

    data = (int *) malloc(size * sizeof (int));
    if (data == NULL) {
        Warn("GetCellAdjacentCellsTri_2D: Memory Allocation Failed 1");
        return 1;
    }

    count = 0;

    for (i = 0; i < NodeAdjacentCells2D[n1].Size; i++) {
        data[count] = NodeAdjacentCells2D[n1].Data[i];
        count += 1;
    }

    for (i = 0; i < NodeAdjacentCells2D[n2].Size; i++) {
        data[count] = NodeAdjacentCells2D[n2].Data[i];
        count += 1;
    }

    for (i = 0; i < NodeAdjacentCells2D[n3].Size; i++) {
        data[count] = NodeAdjacentCells2D[n3].Data[i];
        count += 1;
    }

    sort_bubble(data, size);

    count = 0;
    for (i = 1; i < size; i++) {
        if ((data[i - 1] != data[i])&&(data[i - 1] != cellid))
            count += 1;
    }

    if (data[size - 1] != cellid)
        count += 1;

    Cell2D[cellid].NearCells = count;

    Cell2D[cellid].NearNeighbours = (int *) malloc(count * sizeof (int));
    if (Cell2D[cellid].NearNeighbours == NULL) {
        Warn("GetCellAdjacentCellsTri_2D: Memory Allocation Failed 2");
        return 1;
    }

    j = 0;
    for (i = 1; i < size; i++) {
        if ((data[i - 1] != data[i])&&(data[i - 1] != cellid)) {
            Cell2D[cellid].NearNeighbours[j] = data[i - 1];
            j += 1;
        }
    }

    if (data[size - 1] != cellid)
        Cell2D[cellid].NearNeighbours[j] = data[size - 1];

    free(data);
    return 0;
}

/*---------------------------------------------------------------*/
int GetCellAdjacentCells_2D(void) {
    int icell, nnodes;

    for (icell = 0; icell < NoCells2D; icell++) {
        nnodes = Cell2D[icell].NodesPerCell;
        switch (nnodes) {
            case 3:
                GetCellAdjacentCellsTri_2D(icell);
                break;
            case 4:
                GetCellAdjacentCellsQuad_2D(icell);
                break;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetNodeAdjacentCells_2D(int *adjCells, int *noAdjCells, int nodeID) {
    int i;

    *noAdjCells = NodeAdjacentCells2D[nodeID].Size;

    for (i = 0; i < *noAdjCells; i++)
        adjCells[i] = NodeAdjacentCells2D[nodeID].Data[i];
    return 0;
}

/*---------------------------------------------------------------*/
int CountTotalNonZeroEntries_2D(void) {
    int nonzeros, inode;

    nonzeros = 0;
    for (inode = 0; inode < NoNodes2D; inode++)
        nonzeros += NodeAdjacentNodes2D[inode].Size;

    return (nonzeros);
}

/*---------------------------------------------------------------*/
int CreateSymmetricCSRFormat_2D(void) {
    int inode, icol, nodeID, halfBandWidth, nonZeroEntry, nonzeros;
    int degree;

    IA_2D = (int *) malloc((NoNodes2D + 1) * sizeof (int));
    if (IA_2D == NULL) {
        Warn("CreateSymmetricCSRFormat_2D: Memory Allocation Failed 1");
        return 1;
    }

    nonzeros = (CountTotalNonZeroEntries_2D() / 2) + NoNodes2D;

    JA_2D = (int *) malloc(nonzeros * sizeof (int));
    if (JA_2D == NULL) {
        Warn("CreateSymmetricCSRFormat_2D: Memory Allocation Failed 2");
        return 1;
    }

    IA_2D[0] = 1;
    nonZeroEntry = 0;
    for (inode = 0; inode < NoNodes2D; inode++) {
        degree = NodeAdjacentNodes2D[inode].Size;
        JA_2D[nonZeroEntry++] = inode;
        halfBandWidth = 1;
        for (icol = 0; icol < degree; icol++) {
            nodeID = NodeAdjacentNodes2D[inode].Data[icol];
            if (nodeID > inode) {
                JA_2D[nonZeroEntry++] = nodeID;
                halfBandWidth++;
            }
        }
        IA_2D[inode + 1] = IA_2D[inode] + halfBandWidth;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int FV_InsertSort_2D(int *A, int nsize) {
    int i, j, tmp;

    for (i = 1; i < nsize; i++) {
        tmp = A[i];
        for (j = i; j > 0 && A[j - 1] > tmp; j--)
            A[j] = A[j - 1];
        A[j] = tmp;
    }

    return 0;
}

/*---------------------------------------------------------------*/
int ResizeObject_2D(Vector_2D *vector) {
    int *data, nbytes, capacity;

    if (vector->Capacity == 0) {
        vector->Capacity = 1;
        vector->Data = (int *) malloc(sizeof (int));
    } else {
        /* First make a copy of original nodelist  */
        nbytes = (vector->Size) * sizeof (int);
        data = (int *) malloc(nbytes);
        memcpy(data, vector->Data, nbytes);

        /* Increase the capacity to double */
        capacity = vector->Capacity;
        capacity *= 2;

        /* Free Previous vector->Data */
        free(vector->Data);

        /* allocate new nodelist to double the capacity  */
        nbytes = capacity * sizeof (int);
        vector->Data = (int *) malloc(nbytes);
        vector->Capacity = capacity;

        /* Copy back the nodelist to new list */
        nbytes = (vector->Size) * sizeof (int);
        memcpy(vector->Data, data, nbytes);

        /* Free data */
        free(data);
    }

    return 0;
}

/*---------------------------------------------------------------*/
int Push_Back_2D(Vector_2D *vector, int item) {
    int size, capacity;

    capacity = vector->Capacity;
    size = vector->Size;

    if ((size + 1) > capacity)
        ResizeObject_2D(vector);

    vector->Data[size] = item;
    size++;
    FV_InsertSort_2D(vector->Data, size);

    (vector->Size)++;

    return 0;
}

/*---------------------------------------------------------------*/
int CreateNodeAdjacentCells_2D(void) {
    int icell, n1, n2, n3, n4;
    int inode;

    NodeAdjacentCells2D = NULL;
    NodeAdjacentCells2D = (Vector_2D *) malloc(NoNodes2D * sizeof (Vector_2D));

    if (NodeAdjacentCells2D == NULL) {
        Warn("CreateNodeAdjacentCells_2D: Memory Allocation Failed");
        return 1;
    }

    for (inode = 0; inode < NoNodes2D; inode++) {
        NodeAdjacentCells2D[inode].Size = 0;
        NodeAdjacentCells2D[inode].Capacity = 0;
        NodeAdjacentCells2D[inode].Data = NULL;
    }

    for (icell = 0; icell < NoCells2D; icell++) {
        switch (Cell2D[icell].EdgesPerCell) {
            case 3:
                n1 = Cell2D[icell].ConnectNode[0];
                n2 = Cell2D[icell].ConnectNode[1];
                n3 = Cell2D[icell].ConnectNode[2];
                Push_Back_2D(&NodeAdjacentCells2D[n1], icell);
                Push_Back_2D(&NodeAdjacentCells2D[n2], icell);
                Push_Back_2D(&NodeAdjacentCells2D[n3], icell);
                break;
            case 4:
                n1 = Cell2D[icell].ConnectNode[0];
                n2 = Cell2D[icell].ConnectNode[1];
                n3 = Cell2D[icell].ConnectNode[2];
                n4 = Cell2D[icell].ConnectNode[3];
                Push_Back_2D(&NodeAdjacentCells2D[n1], icell);
                Push_Back_2D(&NodeAdjacentCells2D[n2], icell);
                Push_Back_2D(&NodeAdjacentCells2D[n3], icell);
                Push_Back_2D(&NodeAdjacentCells2D[n4], icell);
                break;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetEdgeNumber_2D(int n1, int n2) {
    int startedge, start, degree, iedge, i;
    int edgenum, node1, node2;

    if (n1 < n2) {
        node1 = n1;
        node2 = n2;
    } else {
        node1 = n2;
        node2 = n1;
    }

    degree = IA_2D[node1 + 1] - IA_2D[node1];
    startedge = IA_2D[node1] - IA_2D[0] - node1;

    start = IA_2D[node1];
    iedge = 0;

    for (i = 0; i < degree; i++) {
        if (JA_2D[start + i] == node2) {
            edgenum = startedge + i;
            if (MaxEdgeNum2D < edgenum)
                MaxEdgeNum2D = edgenum;
            return (edgenum);
        }
    }

    printf("Warning ! There is no Edge between %d %d\n", node1, node2);
    exit(1);
}

/*---------------------------------------------------------------*/
int CreateQuadEdgeList_2D(int icell) {
    int n1, n2, n3, n4, edgeID;

    Cell2D[icell].ConnectEdge = (int *) malloc(4 * sizeof (int));

    n1 = Cell2D[icell].ConnectNode[0];
    n2 = Cell2D[icell].ConnectNode[1];
    n3 = Cell2D[icell].ConnectNode[2];
    n4 = Cell2D[icell].ConnectNode[3];

    edgeID = GetEdgeNumber_2D(n1, n2);
    Cell2D[icell].ConnectEdge[0] = edgeID;

    edgeID = GetEdgeNumber_2D(n2, n3);
    Cell2D[icell].ConnectEdge[1] = edgeID;

    edgeID = GetEdgeNumber_2D(n3, n4);
    Cell2D[icell].ConnectEdge[2] = edgeID;

    edgeID = GetEdgeNumber_2D(n4, n1);
    Cell2D[icell].ConnectEdge[3] = edgeID;

    return 0;
}

/*---------------------------------------------------------------*/
int CreateTriangleEdgeList_2D(int icell) {
    int n1, n2, n3, edgeID;

    Cell2D[icell].ConnectEdge = (int *) malloc(3 * sizeof (int));

    n1 = Cell2D[icell].ConnectNode[0];
    n2 = Cell2D[icell].ConnectNode[1];
    n3 = Cell2D[icell].ConnectNode[2];

    edgeID = GetEdgeNumber_2D(n1, n2);
    Cell2D[icell].ConnectEdge[0] = edgeID;

    edgeID = GetEdgeNumber_2D(n2, n3);
    Cell2D[icell].ConnectEdge[1] = edgeID;

    edgeID = GetEdgeNumber_2D(n3, n1);
    Cell2D[icell].ConnectEdge[2] = edgeID;

    return 0;
}

/*---------------------------------------------------------------*/
int CreateCellEdgeList_2D(void) {
    int icell, nedges;

    for (icell = 0; icell < NoCells2D; icell++) {
        nedges = Cell2D[icell].EdgesPerCell;
        switch (nedges) {
            case 3:
                CreateTriangleEdgeList_2D(icell);
                break;
            case 4:
                CreateQuadEdgeList_2D(icell);
                break;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int GetQuadrilateralEdges_2D(int cellID, int *EdgeNodeList) {
    int minnodeid, maxnodeid, n0, n1, n2, n3;

    n0 = Cell2D[cellID].ConnectNode[0];
    n1 = Cell2D[cellID].ConnectNode[1];
    n2 = Cell2D[cellID].ConnectNode[2];
    n3 = Cell2D[cellID].ConnectNode[3];

    minnodeid = MIN(n0, n1);
    maxnodeid = MAX(n0, n1);
    EdgeNodeList[0] = minnodeid;
    EdgeNodeList[1] = maxnodeid;

    minnodeid = MIN(n1, n2);
    maxnodeid = MAX(n1, n2);
    EdgeNodeList[2] = minnodeid;
    EdgeNodeList[3] = maxnodeid;

    minnodeid = MIN(n2, n3);
    maxnodeid = MAX(n2, n3);
    EdgeNodeList[4] = minnodeid;
    EdgeNodeList[5] = maxnodeid;

    minnodeid = MIN(n3, n0);
    maxnodeid = MAX(n3, n0);
    EdgeNodeList[6] = minnodeid;
    EdgeNodeList[7] = maxnodeid;

    return 0;
}

/*---------------------------------------------------------------*/
int GetTriangleEdges_2D(int cellID, int *EdgeNodeList) {
    int minnodeid, maxnodeid, n0, n1, n2;

    n0 = Cell2D[cellID].ConnectNode[0];
    n1 = Cell2D[cellID].ConnectNode[1];
    n2 = Cell2D[cellID].ConnectNode[2];

    minnodeid = MIN(n0, n1);
    maxnodeid = MAX(n0, n1);
    EdgeNodeList[0] = minnodeid;
    EdgeNodeList[1] = maxnodeid;

    minnodeid = MIN(n1, n2);
    maxnodeid = MAX(n1, n2);
    EdgeNodeList[2] = minnodeid;
    EdgeNodeList[3] = maxnodeid;

    minnodeid = MIN(n0, n2);
    maxnodeid = MAX(n0, n2);
    EdgeNodeList[4] = minnodeid;
    EdgeNodeList[5] = maxnodeid;

    return 0;
}

/*---------------------------------------------------------------*/
int GetCellEdgeNodes_2D(int CellID, int EdgesPerCell, int *EdgeNodeList) {
    switch (EdgesPerCell) {
        case 3:
            GetTriangleEdges_2D(CellID, EdgeNodeList);
            break;
        case 4:
            GetQuadrilateralEdges_2D(CellID, EdgeNodeList);
            break;
    }
    return 0;
}

/*---------------------------------------------------------------*/
int CreateEdgeCellList_2D(void) {
    int edgeID, icell, iedge, numNeighbours;
    int node1, node2, EdgeNodeList[8];
    int nedges;

    for (iedge = 0; iedge < NoEdges2D; iedge++)
        Edge2D[iedge].NoOfNeighbourCells = 1;

    for (icell = 0; icell < NoCells2D; icell++) {
        nedges = Cell2D[icell].EdgesPerCell;
        GetCellEdgeNodes_2D(icell, nedges, EdgeNodeList);

        for (iedge = 0; iedge < nedges; iedge++) {
            node1 = EdgeNodeList[2 * iedge];
            node2 = EdgeNodeList[2 * iedge + 1];
            edgeID = GetEdgeNumber_2D(node1, node2);
            numNeighbours = Edge2D[edgeID].NoOfNeighbourCells;

            if (numNeighbours == 1) {
                Edge2D[edgeID].ConnectedNodes[0] = node1;
                Edge2D[edgeID].ConnectedNodes[1] = node2;
                Edge2D[edgeID].Cell[1] = icell;
            } else {
                Edge2D[edgeID].Cell[0] = icell;
            }

            Edge2D[edgeID].NoOfNeighbourCells = numNeighbours + 1;
        }
    }

    for (iedge = 0; iedge < NoEdges2D; iedge++)
        Edge2D[iedge].NoOfNeighbourCells--;

    NoOfBoundaryEdges2D = 0;
    for (iedge = 0; iedge < NoEdges2D; iedge++) {
        if (Edge2D[iedge].NoOfNeighbourCells == 1) {
            NoOfBoundaryEdges2D++;
            node1 = Edge2D[iedge].ConnectedNodes[0];
            node2 = Edge2D[iedge].ConnectedNodes[1];
            Node2D[node1].Flag = EXTERNAL;
            Node2D[node2].Flag = EXTERNAL;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int SearchItem_2D(int objectID, int item) {
    int i, nitems;

    nitems = NodeAdjacentNodes2D[objectID].Size;
    for (i = 0; i < nitems; i++)
        if (NodeAdjacentNodes2D[objectID].Data[i] == item)
            return (1);
    return 0;
}

/*---------------------------------------------------------------*/
int SearchInsert_2D(int minnode, int maxnode) {
    int found;

    if (minnode > maxnode) {
        Warn("SearchInsert_2D: Node Order Is Wrong");
        return 1;
    }

    if (NodeAdjacentNodes2D[minnode].Size == 0)
        ResizeObject_2D(&NodeAdjacentNodes2D[minnode]);

    found = SearchItem_2D(minnode, maxnode);
    if (!found)
        Push_Back_2D(&NodeAdjacentNodes2D[minnode], maxnode);

    if (NodeAdjacentNodes2D[maxnode].Size == 0)
        ResizeObject_2D(&NodeAdjacentNodes2D[maxnode]);

    found = SearchItem_2D(maxnode, minnode);
    if (!found)
        Push_Back_2D(&NodeAdjacentNodes2D[maxnode], minnode);

    return 0;
}

/*---------------------------------------------------------------*/
int CreateNodeAdjacentNodes_2D(void) {
    int icell, node1, node2, EdgeNodeList[8];
    int iedge, inode, nedges;

    NodeAdjacentNodes2D = NULL;
    NodeAdjacentNodes2D = (Vector_2D *) malloc(NoNodes2D * sizeof (Vector_2D));
    if (NodeAdjacentNodes2D == NULL) {
        Warn("CreateNodeAdjacentNodes_2D: Memory Allocation Failed");
        return 1;
    }

    for (inode = 0; inode < NoNodes2D; inode++) {
        NodeAdjacentNodes2D[inode].Size = 0;
        NodeAdjacentNodes2D[inode].Capacity = 0;
        NodeAdjacentNodes2D[inode].Data = NULL;
    }

    for (icell = 0; icell < NoCells2D; icell++) {
        nedges = Cell2D[icell].EdgesPerCell;
        GetCellEdgeNodes_2D(icell, nedges, EdgeNodeList);
        for (iedge = 0; iedge < nedges; iedge++) {
            node1 = EdgeNodeList[2 * iedge];
            node2 = EdgeNodeList[2 * iedge + 1];
            SearchInsert_2D(node1, node2);
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CreatEdgeStructure_2D(void) {
    int iedge;

    CreateNodeAdjacentNodes_2D();
    CreateSymmetricCSRFormat_2D();
    NoEdges2D = CountTotalNonZeroEntries_2D() / 2;
    printf("No of Edges = %d\n", NoEdges2D);

    Edge2D = (Edge_2D *) malloc(NoEdges2D * sizeof (Edge_2D));
    if (Edge2D == NULL) {
        Warn("CreatEdgeStructure_2D: Memory Allocation Failed");
        return 1;
    }

    /*  Initialize all edges flags as internal Edges. */
    for (iedge = 0; iedge < NoEdges2D; iedge++) {
        Edge2D[iedge].Flag = INTERNAL;
        Edge2D[iedge].Cell[0] = 0;
        Edge2D[iedge].Cell[1] = 0;
        Edge2D[iedge].VisitFlag = 0;
    }
    CreateEdgeCellList_2D();
    CreateCellEdgeList_2D();
    CreateNodeAdjacentCells_2D();
    GetCellAdjacentCells_2D();

    return 0;
}

/*---------------------------------------------------------------*/
int InitializeUnstructuredGrid_2D(ZONE *P) {
    int i, k, e, nnodes, icell;
    ELEMSET *eset;

    NoNodes2D = P->dim[0];
    NoCells2D = P->dim[1];

    /* Allocating memory for nodes */
    Node2D = (Node_2D *) malloc(NoNodes2D * sizeof (Node_2D));
    if (Node2D == NULL) {
        Warn("InitializeUnstructuredGrid_2D: Memory Allocation Failed");
        return 1;
    }

    /* Copying Coordinates */
    for (i = 0; i < NoNodes2D; i++) {
        Node2D[i].Coordinate[0] = P->verts[i].x;
        Node2D[i].Coordinate[1] = P->verts[i].y;
    }

    /* Cell Connectivity */
    /* Below will work properly only when nsections = 1 */
    for (eset = P->esets, e = 1; e <= P->nesets; e++, eset++) {
        if (eset->start != 1 || eset->end != NoCells2D) {
            Warn("InitializeUnstructuredGrid_2D: Cell Connectivity Miss-Matched");
            return 1;
        }

        if (eset->type != TRI_3) {
            Warn("InitializeUnstructuredGrid_2D: Not A Triangluar Grid");
            return 1;
        }

        nnodes = 3;

        Cell2D = (Cell_2D *) malloc(NoCells2D * sizeof (Cell_2D));
        if (Cell2D == NULL) {
            Warn("InitializeUnstructuredGrid_2D: Memory Allocation Failed 2");
            return 1;
        }

        for (icell = 0; icell < NoCells2D; icell++) {
            Cell2D[icell].ConnectNode = (int *) malloc(nnodes * sizeof (int));
            Cell2D[icell].NodesPerCell = nnodes;
            Cell2D[icell].EdgesPerCell = nnodes;
        }

        for (k = 0, i = 0; k < NoCells2D; k++) {
            Cell2D[k].ConnectNode[0] = eset->conn[i];
            i++;
            Cell2D[k].ConnectNode[1] = eset->conn[i];
            i++;
            Cell2D[k].ConnectNode[2] = eset->conn[i];
            i++;
        }

        CreatEdgeStructure_2D();
        PrePostProcessing_2D();
    }

    /* Boundary Conditions */
    if (P->nbocos != 0) {
        InitializeBoundaryCondition_2D(P);
    } else {
        Warn("InitializeUnstructuredGrid_2D: Boundary Condition Missing");
    }

    return 0;
}

/*---------------------------------------------------------------*/
int InitializeStructuredGrid_2D(ZONE *P) {
    /* For Structured Grids */
    //printf ("Dimensions = %d x %d x %d\n",P->dim[0], P->dim[1], P->dim[2]);
    NoNodes2D = P->dim[0] * P->dim[1] * P->dim[2];
    Warn("No Support Now");

    return 0;

}

/*---------------------------------------------------------------*/
int InitializeGrid_2D(ZONE *P) {
    switch (P->type) {
        case 2:
            /* For Structured Grids */
            InitializeStructuredGrid_2D(P);
            break;
        case 3:
            /* For Unstructured Grids */
            InitializeUnstructuredGrid_2D(P);
            break;
        default:
            /* Unknown Type Grid */
            Error_Exit("Grid Type Not Supported");

    }
    return 0;
}
