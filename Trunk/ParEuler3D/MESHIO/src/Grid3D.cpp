/*******************************************************************************
 * File:        Grid3D.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

// Application Headers
#include "Utils.h"
#include "Node3D.h"
#include "Face3D.h"
#include "Ghost3D.h"
#include "Grid3D.h"
#include "CommMPI.h"
#include "Mesh3D_IO.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
Grid3D::Grid3D() {
    // Public
    myOffset        = 0;
    nodeCountOffset = 0;
    nodeCount       = 0;
    cellCount       = 0;
    faceCount       = 0;
    ghostCount      = 0;
    globalNodeCount = 0;
    globalCellCount = 0;
    globalFaceCount = 0;
    // Protected
    gridDB          = NULL;
    scaleFlag       = 0;
    rotateFlag      = 0;
    // Private
    adjIndex        = NULL;
    adjacency       = NULL;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
Grid3D::~Grid3D() {
    
}

//------------------------------------------------------------------------------
//! Reads Database and Creates the Basic Data Structure
//------------------------------------------------------------------------------
int Grid3D::Read_DB() {
    fstream file;

    // Check for Grid File
    file.open(Input_Grid_File.c_str());
    if (file.is_open())
        file.close();
    else
        CommMPI_Error(RankIO, "Unable to Read Grid File: %s", Input_Grid_File.c_str());

    // Check for Solution File
    file.open(Input_Solution_File.c_str());
    if (file.is_open())
        file.close();
    else
        CommMPI_Warn(RankIO, "Unable to Read Solution File: %s", Input_Solution_File.c_str());

    // Check for Parameter File
    file.open(Input_Parameter_File.c_str());
    if (file.is_open())
        file.close();
    else
        CommMPI_Warn(RankIO, "Unable to Read Parameter File: %s", Input_Parameter_File.c_str());
    
    // Now Read Database
    Mesh3D_IO objectDB;
    objectDB.Set_Input_DataBaseType(Input_MeshIO_Type);
    objectDB.Set_InputGrid_Filename(Input_Grid_File.c_str());
    objectDB.Set_InputSolution_Filename(Input_Solution_File.c_str());
    objectDB.Set_InputParam_Filename(Input_Parameter_File.c_str());
    objectDB.Read_DB();
    gridDB = objectDB.Get_DB();
    if (gridDB == NULL)
        CommMPI_Error(RankIO, "Grid3D::Read_DB: %s", "Database Read Failed!");

    // Create the Basic Connectivity to Start with
    Create_Global_DataStructure();

    // Reset the Database
    gridDB = NULL;
    objectDB.Reset_DB();
    
    return 1;
}

//------------------------------------------------------------------------------
//! Scans the CFDDB and Creates the Global Data Structure
//------------------------------------------------------------------------------
void Grid3D::Create_Global_DataStructure() {
    int i, j;
    int nNodes, nCells, nZones, nBocos, nSections;
    int POINT_LIST   = 1;
    int ELEMENT_LIST = 2;

    // Initialize the Global Counts
    globalNodeCount = 0;
    globalCellCount = 0;
    globalFaceCount = 0;

    // Check for Grid Database
    if (gridDB == NULL)
        CommMPI_Error(RankIO, "GridDB::Create_Connectivity_Base: %s", "No Data Found!");

    if (gridDB->nbases > 1)
        CommMPI_Error(RankIO, "GridDB::Create_Connectivity_Base: %s", "No of Bases > 1 Not Supported!");

    nZones = gridDB->bases[0].nzones;

    int zoneNodeCount[nZones];
    int zoneCellCount[nZones];
    // zoneCoordMap returns the global id of a given node in a given zone
    std::vector<int> zoneCoordMap[nZones];

    map<string, int>::iterator mit;
    
    // For each zone
    for (int zoneIndex = 1; zoneIndex <= nZones; zoneIndex++) {
        // These are the number of nodes and cells in the zone
        nNodes = gridDB->bases[0].zones[zoneIndex-1].dim[0];
        nCells = gridDB->bases[0].zones[zoneIndex-1].dim[1];
        
        zoneNodeCount[zoneIndex - 1] = nNodes;
        zoneCellCount[zoneIndex - 1] = nCells;

        // Initialize the Coordinates Mapper reserve memory to avoid rellocation
        zoneCoordMap[zoneIndex - 1].reserve(nNodes);
        for (i = 0; i < nNodes; i++) {
            // This map takes in zone index and node number within that zone
            // Returns the global index of that node
            // It is initialized to -1 for now
            zoneCoordMap[zoneIndex - 1].push_back(-1);
        }
        
        // In case there are multiple connected zones, collapse the repeated nodes and fix the node numbering
        // If the first zone
        if (zoneIndex == 1) {
            for (i = 0; i < nNodes; i++) {
                // Global node count is incremented as new nodes are found.
                // When in the first zone, every node is new.
                zoneCoordMap[0][i] = globalNodeCount;
                globalNodeCount++;
            }
        } else {
            // Scan the coordinates of all the other zones before this one for duplicates
            for (i = 0; i < nNodes; i++) {
                bool foundFlag = false;
                for (int z = 0; z < zoneIndex - 1; ++z) {
                    for (j = 0; j < zoneNodeCount[z]; j++) {
                        if ((fabs(gridDB->bases[0].zones[zoneIndex-1].verts[i].x - gridDB->bases[0].zones[z].verts[j].x) < DBL_TOLERANCE)
                                && (fabs(gridDB->bases[0].zones[zoneIndex-1].verts[i].y - gridDB->bases[0].zones[z].verts[j].y) < DBL_TOLERANCE)
                                && (fabs(gridDB->bases[0].zones[zoneIndex-1].verts[i].z - gridDB->bases[0].zones[z].verts[j].z) < DBL_TOLERANCE)) {
                            zoneCoordMap[zoneIndex - 1][i] = zoneCoordMap[z][j];
                            foundFlag = true;
                            break;
                        }
                    }
                    if (foundFlag) break;
                }
                if (!foundFlag) {
                    zoneCoordMap[zoneIndex - 1][i] = globalNodeCount;
                    globalNodeCount++;
                }
            }
        }

        // Boundary condition regions may be given as element or point ranges, or element or point lists
        // For each zone and each boundary condition region, store a point list, convert if another method is used
        vector<int> bc_method;
        vector<set<int> > bc_element_list;
        nBocos = gridDB->bases[0].zones[zoneIndex-1].nbocos;
        for (int bocoIndex = 1; bocoIndex <= nBocos; bocoIndex++) {
            // Find out the bc specification method
            int ptype;
            int npnts;
            
            ptype = gridDB->bases[0].zones[zoneIndex-1].bocos[bocoIndex-1].ptype;
            npnts = gridDB->bases[0].zones[zoneIndex-1].bocos[bocoIndex-1].npnts;

            // Check if this bc name matches to those found before
            bool new_bc = true;
            int bcIndex = -1;
            string bcName(gridDB->bases[0].zones[zoneIndex-1].bocos[bocoIndex-1].name);
            for (mit = global_bocoNameMap.begin(); mit != global_bocoNameMap.end(); mit++) {
                if ((*mit).first == bcName) {
                    new_bc = false;
                    bcIndex = (*mit).second;
                    break;
                }
            }

            // If no match, create new bc
            if (new_bc) {
                bcIndex = global_bocoNameMap.size();
                global_bocoNameMap.insert(pair<string, int>(bcName, bcIndex));
                global_bocoNodes.resize(bcIndex + 1);
            }


            if ((int)bc_method.size() < bcIndex + 1) {
                bc_method.resize(bcIndex + 1);
                bc_element_list.resize(bcIndex + 1);
            }

            int *list;
            // Check the bc specification method
            // PointList = 2; PointRange = 4; ElementRange = 6; ElementList = 7
            list = gridDB->bases[0].zones[zoneIndex-1].bocos[bocoIndex-1].pnts;
            switch (ptype) {
                case 2:
                    // PointList
                    bc_method[bcIndex] = POINT_LIST;
                    for (i = 0; i < npnts; i++)
                        global_bocoNodes[bcIndex].insert(zoneCoordMap[zoneIndex-1][list[i] - 1]);
                    break;
                case 4:
                    // PointRange
                    bc_method[bcIndex] = POINT_LIST;
                    for (i = list[0]; i <= list[1]; i++)
                        global_bocoNodes[bcIndex].insert(zoneCoordMap[zoneIndex-1][i - 1]);
                    break;
                case 6:
                    // ElementRange
                    // Convert element range to element list
                    bc_method[bcIndex] = ELEMENT_LIST;
                    for (i = list[0]; i <= list[1]; i++)
                        bc_element_list[bcIndex].insert(i);
                    break;
                case 7:
                    // ElementList
                    bc_method[bcIndex] = ELEMENT_LIST;
                    for (i = 0; i < npnts; i++)
                        bc_element_list[bcIndex].insert(list[i]);
                    break;
                default:
                    CommMPI_Error(RankIO, "GridDB::Create_Connectivity_Base: %s",
                            "Unsupported Boundary Conditions !");
            }  
        }
        nBocos = global_bocoNameMap.size();

        // Loop sections within the zone
        // These include connectivities of cells and boundary faces
        nSections = gridDB->bases[0].zones[zoneIndex-1].nesets;
        for (int sectionIndex = 1; sectionIndex <= nSections; sectionIndex++) {
            
            int elemType, elemNodeCount, elemStart, elemEnd;
            elemNodeCount = 0;
            elemType   = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].type;
            elemStart  = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].start;
            elemEnd    = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].end;

            // ElementType
            switch (elemType) {
                // Triangle
                case 5:
                    elemNodeCount = 3;
                    break;
                // Quadrilateral
                case 7:
                    elemNodeCount = 4;
                    break;
                // Tetrahydron
                case 10:
                    elemNodeCount = 4;
                    break;
                // Pyramid
                case 12:
                    elemNodeCount = 5;
                    break;
                // Prism
                case 14:
                    elemNodeCount = 6;
                    break;
                // Hexahydron
                case 17:
                    elemNodeCount = 8;
                    break;
            }

            // Mixed
            if (elemType == 20) {
                int count = 0;
                int nodeid;
                for (int elem = elemStart; elem <= elemEnd; elem++) {
                    // First Entry is Element Type
                    elemType = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].conn[count];
                    switch (elemType) {
                        // Triangle
                        case 5:
                            elemNodeCount = 3;
                            break;
                        // Quadrilateral
                        case 7:
                            elemNodeCount = 4;
                            break;
                        // Tetrahydron
                        case 10:
                            elemNodeCount = 4;
                            break;
                        // Pyramid
                        case 12:
                            elemNodeCount = 5;
                            break;
                        // Prism
                        case 14:
                            elemNodeCount = 6;
                            break;
                        // Hexahydron
                        case 17:
                            elemNodeCount = 8;
                            break;
                    }
                    count++;
                    // Get the Volume Elements
                    if (elemType == 10 || elemType == 12 || elemType == 14 || elemType == 17) {
                        global_cellConnIndex.push_back(global_cellConnectivity.size());
                        for (int n = 0; n < elemNodeCount; n++) {
                            nodeid = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].conn[count+n];
                            global_cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][nodeid - 1]);
                        }
                        globalCellCount++;
                    } else {
                        // Boundary Elements: Triangle and Quadrilaterals
                        for (int nbc = 0; nbc < nBocos; nbc++) {
                            if (bc_method[nbc] == ELEMENT_LIST) {
                                // Check if Element Belong to this Boundary
                                if (bc_element_list[nbc].find(elem) != bc_element_list[nbc].end()) {
                                    for (int n = 0; n < elemNodeCount; n++) {
                                        nodeid = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].conn[count+n];
                                        global_bocoNodes[nbc].insert(zoneCoordMap[zoneIndex - 1][nodeid - 1]);
                                    }
                                }
                            }
                        }
                    }
                    count += elemNodeCount;
                }
            } else {
                // Get the Volume Elements
                if (elemType == 10 || elemType == 12 || elemType == 14 || elemType == 17) {
                    // elements array serves as a start index for connectivity list elemConnectivity
                    int nodeid;
                    for (int elem = 0; elem <= (elemEnd - elemStart); elem++) {
                        global_cellConnIndex.push_back(global_cellConnectivity.size());
                        for (int n = 0; n < elemNodeCount; n++) {
                            nodeid = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].conn[elemNodeCount*elem+n];
                            global_cellConnectivity.push_back(zoneCoordMap[zoneIndex - 1][nodeid - 1]);
                        }
                    }
                    globalCellCount += (elemEnd - elemStart + 1);
                } else {
                    // Boundary Elements: Triangle and Quadrilaterals
                    for (int nbc = 0; nbc < nBocos; ++nbc) {
                        if (bc_method[nbc] == ELEMENT_LIST) {
                            int nodeid;
                            for (int elem = 0; elem <= (elemEnd - elemStart); elem++) {
                                if (bc_element_list[nbc].find(elemStart + elem) != bc_element_list[nbc].end()) {
                                    for (int n = 0; n < elemNodeCount; n++) {
                                        nodeid = gridDB->bases[0].zones[zoneIndex-1].esets[sectionIndex-1].conn[elemNodeCount*elem+n];
                                        global_bocoNodes[nbc].insert(zoneCoordMap[zoneIndex - 1][nodeid - 1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Now Merge the Coordinates from all zones
    global_x.reserve(globalNodeCount);
    global_y.reserve(globalNodeCount);
    global_z.reserve(globalNodeCount);
    // for zone 0
    for (i = 0; i < zoneNodeCount[0]; i++) {
        global_x.push_back(gridDB->bases[0].zones[0].verts[i].x);
        global_y.push_back(gridDB->bases[0].zones[0].verts[i].y);
        global_z.push_back(gridDB->bases[0].zones[0].verts[i].z);
    }
    for (int z = 1; z < nZones; z++) {
        for (i = 0; i < zoneNodeCount[z]; i++) {
            bool unique = true;
            for (j = 0; j < (int)zoneCoordMap[z - 1].size(); j++) {
                if (zoneCoordMap[z][i] <= zoneCoordMap[z - 1][j]) {
                    unique = false;
                    break;
                }
            }
            if (unique) {
                global_x.push_back(gridDB->bases[0].zones[z].verts[i].x);
                global_y.push_back(gridDB->bases[0].zones[z].verts[i].y);
                global_z.push_back(gridDB->bases[0].zones[z].verts[i].z);
            }
        }
    }

    if ((int)global_x.size() != globalNodeCount)
        CommMPI_Error(RankIO, "GridDB::Create_Connectivity_Base: %s",
                "Global Node Count Mis-Matched !");
}

//------------------------------------------------------------------------------
//! Deletes the Global Data Structure
//------------------------------------------------------------------------------
void Grid3D::Delete_Global_DataStructure() {
    global_x.clear();
    global_y.clear();
    global_z.clear();
    global_cellConnIndex.clear();
    global_cellConnectivity.clear();
    global_bocoNodes.clear();
    global_bocoNameMap.clear();
}

//------------------------------------------------------------------------------
//! Creates Local Cells and Nodes of the Partition Mesh
//! Inprocess creates Local Cell2Node Connectivity
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Nodes_And_Cells() {
    // Reserve the right memory of Cells after Partition
    cell.reserve(cellCount);

    // Reset the Local Partition Nodes
    nodeCount = 0;

    // Loop Over the Global Cell to Identify the Current Partition Cells
    for (int c = 0; c < globalCellCount; c++) {
        // Check if Cell Belongs to Current Process
        if (cellOwner[c] == Rank) {
            int cellNodeCount;
            // Find the number of nodes of the cell from global grid data
            if (c < globalCellCount - 1) {
                cellNodeCount = global_cellConnIndex[c + 1] - global_cellConnIndex[c];
            } else {
                cellNodeCount = global_cellConnectivity.size() - global_cellConnIndex[globalCellCount - 1];
            }
            int cellNodes[cellNodeCount];
            // Loop the cell  nodes
            for (int n = 0; n < cellNodeCount; n++) {
                // Get the Global Id of Node
                int gid = global_cellConnectivity[global_cellConnIndex[c] + n];
                // Check if Node Already Exist
                if (nodeGlobal2Local.find(gid) == nodeGlobal2Local.end()) {
                    // Create the node
                    Node3D newnode;
                    newnode.id       = nodeCount;
                    newnode.globalId = gid;
                    newnode.pos[0]   = global_x[gid];
                    newnode.pos[1]   = global_y[gid];
                    newnode.pos[2]   = global_z[gid];
                    // Insert the new element in map
                    nodeGlobal2Local[gid] = nodeCount;
                    node.push_back(newnode);
                    nodeCount++;
                }
                // Fill in cell nodes temp array with local node id's
                cellNodes[n] = nodeGlobal2Local[gid];
            }

            // Create the cell
            Cell3D newcell;
            newcell.nodeCount = cellNodeCount;
            switch (cellNodeCount) {
                // Tetrahyderon
                case 4:
                    newcell.faceCount = 4;
                    break;
                // Pyramid
                case 5:
                    newcell.faceCount = 5;
                    break;
                // Prism
                case 6: 
                    newcell.faceCount = 5;
                    break;
                // Hexahydron
                case 8: 
                    newcell.faceCount = 6;
                    break;
            }
            newcell.nodes.reserve(cellNodeCount);

            // Fill in the node list
            for (int n = 0; n < newcell.nodeCount; ++n) {
                newcell.nodes.push_back(cellNodes[n]);
            }

            newcell.globalId = c;
            cellGlobal2Local[newcell.globalId] = cell.size();

            cell.push_back(newcell);
        }
    }
}

//------------------------------------------------------------------------------
//! Create the Local Connectivity For Nodes Connected with Cells
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Connectivity_Node2Cell() {
    bool flag;
    
    for (int c = 0; c < cellCount; c++) {
        int n;
        // For each Node in the Cell
        for (int cn = 0; cn < cell[c].nodeCount; cn++) {
            // Get the Node ID
            n = cell[c].nodes[cn];
            flag = false;
            // Check if Cell is already added
            for (int i = 0; i < (int)node[n].cells.size(); i++) {
                if (node[n].cells[i] == c) {
                    flag = true;
                    break;
                }
            }
            if (!flag)
                node[n].cells.push_back(c);
        }
    }
}

//------------------------------------------------------------------------------
//! Create the Local Connectivity For Cells Connected to Cells
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Connectivity_Cell2Cell() {
    bool flag;
    int c2;
    for (int c = 0; c < cellCount; c++) {
        int n;
        // Loop nodes of the cell
        for (int cn = 0; cn < cell[c].nodeCount; cn++) {
            // Get the Node ID
            n = cell[c].nodes[cn];
            // Loop neighboring cells of the node
            for (int nc = 0; nc < (int)node[n].cells.size(); nc++) {
                // Get the Cell ID of the connected cell to the node
                c2 = node[n].cells[nc];
                flag = false;
                // Check if the cell was found before
                for (int cc = 0; cc < (int)cell[c].neighborCells.size(); cc++) {
                    if (cell[c].neighborCells[cc] == c2) {
                        flag = true;
                        break;
                    }
                }
                if (!flag)
                    cell[c].neighborCells.push_back(c2);
            }
        }
        cell[c].neighborCellCount = cell[c].neighborCells.size();
    }

    // Convert the Global Boundary Nodes to Local Boundary Nodes
    for (int nbc = 0; nbc < (int)global_bocoNodes.size(); nbc++) {
        set<int> temp;
        set<int>::iterator sit;
        for (sit = global_bocoNodes[nbc].begin(); sit != global_bocoNodes[nbc].end(); sit++) {
            if (nodeGlobal2Local.find(*sit) != nodeGlobal2Local.end()) {
                temp.insert(nodeGlobal2Local[*sit]);
            }
        }
        global_bocoNodes[nbc].swap(temp);
        temp.clear();
    }
}

//------------------------------------------------------------------------------
//! Create the Local Face Connectivity and
//! Categorize as Internal, Boundary and Unassigned
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Connectivity_Face() {
    // Set face connectivity lists using SIDS Convention
    // Tetrahedral
    int tetraFaces[4][3] = {
        {0, 2, 1},
        {0, 1, 3},
        {1, 2, 3},
        {2, 0, 3},
    };
    // Pyramid
    int pyraFaces[5][4] = {
        {0, 3, 2, 1},
        {0, 1, 4, 0},
        {1, 2, 4, 0},
        {2, 3, 4, 0},
        {3, 0, 4, 0}
    };
    // Prism
    int prismFaces[5][4] = {
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {2, 0, 3, 5},
        {0, 2, 1, 0},
        {3, 4, 5, 0}
    };
    // Hexahedral
    int hexaFaces[6][4] = {
        {0, 3, 2, 1},
        {0, 1, 5, 4},
        {1, 2, 6, 5},
        {2, 3, 7, 6},
        {0, 4, 7, 3},
        {4, 5, 6, 7}
    };

    // Initialize the Face Count
    faceCount = 0;

    // Start Searching the Faces
    // Loop through all the cells
    for (int c = 0; c < cellCount; c++) {
        // Loop through the faces of the current cell
        for (int cf = 0; cf < cell[c].faceCount; cf++) {
            Face3D newface;
            int *tempNodes = NULL;
            switch (cell[c].nodeCount) {
                // Tetrahedral
                case 4:
                    newface.nodeCount = 3;
                    tempNodes = new int[3];
                    break;
                // Pyramid
                case 5:
                    // Quad Face
                    if (cf < 1) {
                        newface.nodeCount = 4;
                        tempNodes = new int[4];
                    } else {
                        // Tri Face
                        newface.nodeCount = 3;
                        tempNodes = new int[3];
                    }
                    break;
                // Prism
                case 6:
                    // Quad Face
                    if (cf < 3) {
                        newface.nodeCount = 4;
                        tempNodes = new int[4];
                    } else {
                        // Tri Face
                        newface.nodeCount = 3;
                        tempNodes = new int[3];
                    }
                    break;
                // Hexahedral
                case 8:
                    newface.nodeCount = 4;
                    tempNodes = new int[4];
                    break;
            }

            if (tempNodes == NULL)
                CommMPI_Error(Rank, "Grid3D::Create_Local_Connectivity_Face: %s", "Face Creation Failed !");
            
            // Face count is incremented everytime a new face is found
            newface.id = faceCount;
            // Assign current cell as the parent cell
            newface.parent = c;
            // Assign boundary type as internal by default, will be overwritten later
            newface.bc = INTERNAL;
            // Store the node local ids of the current face
            for (int fn = 0; fn < newface.nodeCount; fn++) {
                switch (cell[c].nodeCount) {
                    // Tetrahedral
                    case 4:
                        tempNodes[fn] = cell[c].nodes[tetraFaces[cf][fn]];
                        break;
                    // Pyramid
                    case 5:
                        tempNodes[fn] = cell[c].nodes[pyraFaces[cf][fn]];
                        break;
                    // Prism
                    case 6:
                        tempNodes[fn] = cell[c].nodes[prismFaces[cf][fn]];
                        break;
                    // Hexahedral
                    case 8:
                        tempNodes[fn] = cell[c].nodes[hexaFaces[cf][fn]];
                        break;
                }
            }

            // Find the neighbor cell
            bool internal = false;
            bool unique = true;
            // Loop cells neighboring the first node of the current face
            for (int nc = 0; nc < (int)node[tempNodes[0]].cells.size(); nc++) {
                // i is the neighbor cell index
                int i = node[tempNodes[0]].cells[nc];
                // If neighbor cell is not the current cell itself, and it has the same nodes as the face
                if (i != c && cell[i].HaveNodes(newface.nodeCount, tempNodes)) {
                    // If the neighbor cell index is smaller then the current cell index,
                    // it has already been processed so skip it
                    if (i > c) {
                        newface.neighbor = i;
                        internal = true;
                    } else {
                        unique = false;
                    }
                    break;
                }
            }

            // Add only New Faces
            if (unique) {
                // Insert the node list
                for (int fn = 0; fn < newface.nodeCount; fn++)
                    newface.nodes.push_back(tempNodes[fn]);
                // Check if Face is Boundary or Inter-Partition Face
                if (!internal) {
                    newface.bc = UNASSIGNED;
                    vector<int> face_matched_bcs;
                    int cell_matched_bc = -1;
                    bool match;
                    // For each boundary condition region
                    for (int nbc = 0; nbc < (int)global_bocoNameMap.size(); nbc++) {
                        match = true;
                        // For each node of the current face
                        for (int i = 0; i < newface.nodeCount; i++) {
                            if (global_bocoNodes[nbc].find(tempNodes[i]) == global_bocoNodes[nbc].end()) {
                                match = false;
                                break;
                            }
                        }
                        // This means that all the face nodes are on the current bc node list
                        if (match) {
                            face_matched_bcs.push_back(nbc);
                        }
                        // There can be situations like back and front symmetry BC's in which
                        // face nodes will match more than one boundary condition
                        // Check if the owner cell has all its nodes on one of those bc's
                        // and eliminate those
                        if (cell_matched_bc == -1) {
                            match = true;
                            for (int i = 0; i < cell[c].nodeCount; i++) {
                                if (global_bocoNodes[nbc].find(cell[c].nodes[i]) == global_bocoNodes[nbc].end()) {
                                    match = false;
                                    break;
                                }
                            }
                            // This means that all the cell nodes are on the current bc node list
                            if (match) {
                                cell_matched_bc = nbc;
                            }
                        }
                    }
                    if (face_matched_bcs.size() > 1) {
                        for (int fbc = 0; fbc < (int)face_matched_bcs.size(); fbc++) {
                            if (face_matched_bcs[fbc] != cell_matched_bc) {
                                newface.bc = face_matched_bcs[fbc];
                                break;
                            }
                        }
                    } else if (face_matched_bcs.size() == 1) {
                        newface.bc = face_matched_bcs[0];
                    }
                }
                newface.parentIndex = cf;
                face.push_back(newface);
                cell[c].faces.push_back(newface.id);
                if (internal)
                    cell[newface.neighbor].faces.push_back(newface.id);
                // Increment the Face Count
                faceCount++;
            }
            delete [] tempNodes;
        }
    }
}

//------------------------------------------------------------------------------
//! Create the Local Ghost Cell Connectitity with Other Partion
//! Mark Inter Patition Faces
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Connectivity_Ghost() {
    ghostCount = 0;

    if (NoProc == 1) {
        for (int c = 0; c < cellCount; c++)
            cell[c].ghostCount = 0;
    } else {
        int counter = 0;
        int cellCountOffset[NoProc];

        // Find out other partition's cell counts
        int otherCellCounts[NoProc];
        for (int i = 0; i < NoProc; i++)
            otherCellCounts[i] = 0;
        for (int c = 0; c < globalCellCount; c++)
            otherCellCounts[cellOwner[c]] += 1;

        for (int p = 0; p < NoProc; p++) {
            cellCountOffset[p] = counter;
            counter += otherCellCounts[p];
        }

        // Now find metis2global index mapping
        int metis2global[globalCellCount];
        int counter2[NoProc];
        for (int p = 0; p < NoProc; p++)
            counter2[p] = 0;
        for (int c = 0; c < globalCellCount; c++) {
            metis2global[cellCountOffset[cellOwner[c]] + counter2[cellOwner[c]]] = c;
            counter2[cellOwner[c]]++;
        }

        int foundFlag[globalCellCount];
        for (int c = 0; c < globalCellCount; c++)
            foundFlag[c] = 0;

        int parent, metisIndex, gg, matchCount;
        // Stores global id's of ghost cells near each node
        map<int, set<int> > nodeGhostSet;

        // Loop faces
        for (int f = 0; f < faceCount; f++) {
            // Check if Face is Unassigned or Boundary: Not defined Internal
            if (face[f].bc == UNASSIGNED || face[f].bc >= 0) {
                parent = face[f].parent;
                // Loop through the cells that are adjacent to the current face's parent
                for (int adjCount = 0; adjCount < (adjIndex[parent + 1] - adjIndex[parent]); adjCount++) {
                    metisIndex = adjacency[adjIndex[parent] + adjCount];
                    // Get global id of the adjacent cell
                    gg = metis2global[metisIndex];
                    // Find the number of nodes of the cell from raw grid data
                    int cellNodeCount;
                    if (gg < globalCellCount - 1) {
                        cellNodeCount = global_cellConnIndex[gg + 1] - global_cellConnIndex[gg];
                    } else {
                        cellNodeCount = global_cellConnectivity.size() - global_cellConnIndex[globalCellCount - 1];
                    }
                    // If that cell is not on the current partition
                    if (metisIndex < cellCountOffset[Rank] || metisIndex >= (cellCount + cellCountOffset[Rank])) {
                        // Count number of matches in node lists of the current face and the adjacent cell
                        matchCount = 0;
                        for (int fn = 0; fn < face[f].nodeCount; fn++) {
                            set<int> tempSet;
                            nodeGhostSet.insert(pair<int, set<int> >(face[f].nodes[fn], tempSet));
                            for (int gn = 0; gn < cellNodeCount; gn++) {
                                if (global_cellConnectivity[global_cellConnIndex[gg] + gn] == node[face[f].nodes[fn]].globalId) {
                                    nodeGhostSet[face[f].nodes[fn]].insert(gg);
                                    matchCount++;
                                }
                            }
                        }
                        // foundFlag is 0 by default
                        // 0 means that particular adjacent cell wasn't discovered as a ghost before
                        // If so, create a new ghost
                        if (matchCount > 0 && foundFlag[gg] == 0) {
                            foundFlag[gg] = matchCount;
                            Ghost3D newghost;
                            newghost.globalId = gg;
                            newghost.partition = cellOwner[gg];
                            ghostGlobal2Local[newghost.globalId] = ghost.size();
                            ghost.push_back(newghost);
                            ghostCount++;
                        }
                        // If that ghost was found before, now we discovered another face also neighbors the same ghost
                        if (matchCount == face[f].nodeCount) {
                            foundFlag[gg] = matchCount;
                            face[f].bc = GHOST;
                            face[f].neighbor = -1 * ghostGlobal2Local[gg] - 1;
                        }
                    }
                }
            }
        }

        // Store the local id's of ghosts touching each node
        map<int, set<int> >::iterator mit;
        set<int>::iterator sit;
        for (mit = nodeGhostSet.begin(); mit != nodeGhostSet.end(); mit++) {
            for (sit = (*mit).second.begin(); sit != (*mit).second.end(); sit++) {
                node[(*mit).first].ghosts.push_back(ghostGlobal2Local[*sit]);
            }
        }

        // Construct the list of neighboring ghosts for each cell
        int g;
        bool flag;
        for (int c = 0; c < cellCount; c++) {
            int n;
            for (int cn = 0; cn < cell[c].nodeCount; cn++) {
                n = cell[c].nodes[cn];
                for (int ng = 0; ng < (int)node[n].ghosts.size(); ng++) {
                    g = node[n].ghosts[ng];
                    // Check if ghost is already added
                    flag = false;
                    for (int cg = 0; cg < (int)cell[c].ghosts.size(); cg++) {
                        if (cell[c].ghosts[cg] == g) {
                            flag = true;
                            break;
                        }
                    }
                    // Add the ghost Cell
                    if (flag == false) {
                        cell[c].ghosts.push_back(g);
                        ghost[g].cells.push_back(c);
                    }
                }
            }
            cell[c].ghostCount = cell[c].ghosts.size();
        }
    }

    nodeGlobal2Output.resize(globalNodeCount);
    for (int ng = 0; ng < globalNodeCount; ng++)
        nodeGlobal2Output[ng] = -1;
    nodeCountOffset = 0;
    int nodeCountOffsetPlus = 0;

    // Receive nodeCountOffset from processor Rank-1
    // Set tag to destination
    if (Rank != 0)
        MPI_Recv(&nodeCountOffset, 1, MPI_INT, Rank - 1, Rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (Rank != 0)
        MPI_Recv(&nodeGlobal2Output[0], globalNodeCount, MPI_INT, Rank - 1, Rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int n = 0; n < nodeCount; n++) {
        if (nodeGlobal2Output[node[n].globalId] == -1) {
            nodeGlobal2Output[node[n].globalId] = nodeCountOffset + nodeCountOffsetPlus;
            nodeCountOffsetPlus++;
        }
    }
    nodeCountOffsetPlus += nodeCountOffset;

    if (Rank != NoProc - 1)
        MPI_Send(&nodeCountOffsetPlus, 1, MPI_INT, Rank + 1, Rank + 1, MPI_COMM_WORLD);
    if (Rank != NoProc - 1)
        MPI_Send(&nodeGlobal2Output[0], globalNodeCount, MPI_INT, Rank + 1, Rank + 1, MPI_COMM_WORLD);

}

//------------------------------------------------------------------------------
//! Creates All Connectivity Graphs
//------------------------------------------------------------------------------
void Grid3D::Create_Local_Connectivity_Maps() {
    Create_Local_Nodes_And_Cells();
    Create_Local_Connectivity_Node2Cell();
    Create_Local_Connectivity_Cell2Cell();
    Parmetis_Create_Mesh2Dual();
    Create_Local_Connectivity_Face();
    Create_Local_Connectivity_Ghost();

    // Optimize the Local Connectivity Memory Usage
    for (int n = 0; n < nodeCount; n++) {
        vector<int> (node[n].cells).swap(node[n].cells);
        vector<int> (node[n].ghosts).swap(node[n].ghosts);
    }

    for (int f = 0; f < faceCount; f++)
        vector<int> (face[f].nodes).swap(face[f].nodes);

    for (int c = 0; c < cellCount; c++) {
        vector<int> (cell[c].nodes).swap(cell[c].nodes);
        vector<int> (cell[c].faces).swap(cell[c].faces);
        vector<int> (cell[c].neighborCells).swap(cell[c].neighborCells);
        vector<int> (cell[c].ghosts).swap(cell[c].ghosts);
    }

    for (int g = 0; g < ghostCount; g++)
        vector<int> (ghost[g].cells).swap(ghost[g].cells);

    vector<Node3D >  (node).swap(node);
    vector<Face3D >  (face).swap(face);
    vector<Cell3D >  (cell).swap(cell);
    vector<Ghost3D > (ghost).swap(ghost);
    vector<int>      (partitionOffset).swap(partitionOffset);
}

//------------------------------------------------------------------------------
//! Sets the Scaling Flag and Value of Scale
//------------------------------------------------------------------------------
void Grid3D::Set_Scale(Vector3D value) {
    scaleFlag = 1;
    scaleFactor[0] = value[0];
    scaleFactor[1] = value[1];
    scaleFactor[2] = value[2];
}

//------------------------------------------------------------------------------
//! Scales the Grid Coordinates in X, Y and Z directions
//------------------------------------------------------------------------------
void Grid3D::Scale_DB() {
    if (scaleFlag) {
        for (int n = 0; n < globalNodeCount; n++) {
            global_x[n] *= scaleFactor[0];
            global_y[n] *= scaleFactor[1];
            global_z[n] *= scaleFactor[2];
        }
        CommMPI_Info(RankIO, "Grid Data Base Scaled: [%lf, %lf, %lf]",
                scaleFactor[0], scaleFactor[1], scaleFactor[2]);
    }
}

//------------------------------------------------------------------------------
//! Sets the Rotation Center and Angles
//------------------------------------------------------------------------------
void Grid3D::Set_Rotation(Vector3D center, Vector3D angle) {
    rotateFlag = 1;
    rotationCenter[0] = center[0];
    rotationCenter[1] = center[1];
    rotationCenter[2] = center[2];
    rotationAngle[0]  = angle[0];
    rotationAngle[1]  = angle[1];
    rotationAngle[2]  = angle[2];
}

//------------------------------------------------------------------------------
//! Rotates the Grid Coordinates using Given Center and Angles
//------------------------------------------------------------------------------
void Grid3D::Rotate_DB() {
    if (rotateFlag) {
        // Convert angles to radian
        rotationAngle *= 4. * atan(1.) / 180.;
        double x, y, z;
        for (int n = 0; n < globalNodeCount; ++n) {
            // Translate to the new center
            global_x[n] -= rotationCenter[0];
            global_y[n] -= rotationCenter[1];
            global_z[n] -= rotationCenter[2];
            x = global_x[n];
            y = global_y[n];
            z = global_z[n];
            // Rotate around x axis
            global_y[n] = y * cos(rotationAngle[0]) - z * sin(rotationAngle[0]);
            global_z[n] = y * sin(rotationAngle[0]) + z * cos(rotationAngle[0]);
            // Rotate around y axis
            global_x[n] = x * cos(rotationAngle[1]) + z * sin(rotationAngle[1]);
            global_z[n] = -x * sin(rotationAngle[1]) + z * cos(rotationAngle[1]);
            // Rotate around z axis
            global_x[n] = x * cos(rotationAngle[2]) - y * sin(rotationAngle[2]);
            global_y[n] = x * sin(rotationAngle[2]) + y * cos(rotationAngle[2]);
            // Translate back to the original center
            global_x[n] += rotationCenter[0];
            global_y[n] += rotationCenter[1];
            global_z[n] += rotationCenter[2];
        }
        CommMPI_Info(RankIO, "Grid Data Base Rotated: Center[%lf, %lf, %lf], Angle[%lf, %lf, %lf]",
                rotationCenter[0], rotationCenter[1], rotationCenter[2],
                rotationAngle[0], rotationAngle[1], rotationAngle[2]);
    }
}

//------------------------------------------------------------------------------
//! Grid Partition Using Parmetis
//------------------------------------------------------------------------------
void Grid3D::Parmetis_Partition_DB() {
    // Initialize the partition sizes
    // This is just a simple manual partitioning to be able to use parmetis afterwards
    cellCount = floor(double(globalCellCount) / double(NoProc));
    
    int baseCellCount = cellCount;
    int offset = Rank*cellCount;
    // Adjust the cellCount for last partition to contain remaining cells
    if (Rank == NoProc - 1)
        cellCount = cellCount + globalCellCount - NoProc * cellCount;
    
    // This Array Describes how elements of mesh are distributed among processors
    idxtype elmdist[NoProc + 1];
    // This Array Specify the elements that are stored locally at each processor
    // like xadj
    idxtype *eptr;
    eptr = new idxtype[cellCount + 1];
    // This Array Specify the elements that are stored locally at each processor
    // like adjncy
    idxtype *eind;
    int eindSize = 0;
    if ((offset + cellCount) == globalCellCount) {
        eindSize = global_cellConnectivity.size() - global_cellConnIndex[offset];
    } else {
        eindSize = global_cellConnIndex[offset + cellCount] - global_cellConnIndex[offset] + 1;
    }
    
    eind = new idxtype[eindSize];
    // This Array stores the weight of the elements
    idxtype* elmwgt = NULL;
    // This is used to indicate if the graph is weighted
    // no weights associated with elem or edges
    int wgtflag = 0;
    // This used to indicate the numbering scheme
    // 0 C-style numbers, 1 Fortran-style numbers
    int numflag = 0;
    // This is used to specify the number of weights that each vertex has
    int ncon = 1;
    // This parameter determines the degree of connectivty among the vertices in the dual graph
    // set to 3 for tetrahedra or mixed type
    int ncommonnodes = 3;
    // This is used to specify the number of sub-domain that are desired.
    // (Note: BE CAREFUL if != to # of proc)
    int nparts = NoProc;
    // This array is used to specify the fraction of vertex weight that should be distributed
    // to each sub-domain for each balance constraint. Setting same weight for each vertex = 1/nparts
    float tpwgts[NoProc];
    for (int p = 0; p < NoProc; p++)
        tpwgts[p] = 1. / float(nparts);
    // This Array specify the imbalance tolerance for each vertex weight
    // 1.0 = perfect balance; nparts = perfect imbalance, recommended = 1.05
    float ubvec = 1.02;
    // This Array is used to pass additional paramters for the routine
    // Using the default settings: [0 1 15] for default
    int options[3];
    options[0] = 0;
    options[1] = 1;
    options[2] = 15;
    // Output: Stores the number of edges that are cut by the partitioning
    int edgecut;
    // Output: This Array of size number of locally-stored verticies
    idxtype* part = new idxtype[cellCount];
    // This a pointer to the MPI Communicator of the process
    MPI_Comm commWorld = MPI_COMM_WORLD;
    
    for (int p = 0; p < NoProc; p++)
        elmdist[p] = p * floor(globalCellCount / NoProc);
    // Note this is because #elements mod(NoProc) are all on last proc
    elmdist[NoProc] = globalCellCount;

    for (int c = 0; c < cellCount; c++) {
        eptr[c] = global_cellConnIndex[offset + c] - global_cellConnIndex[offset];
    }
    if ((offset + cellCount) == globalCellCount) {
        eptr[cellCount] = global_cellConnectivity.size() - global_cellConnIndex[offset];
    } else {
        eptr[cellCount] = global_cellConnIndex[offset + cellCount] - global_cellConnIndex[offset];
    }

    for (int i = 0; i < eindSize; ++i) {
        eind[i] = global_cellConnectivity[global_cellConnIndex[offset] + i];
    }
    
    ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt,
            &wgtflag, &numflag, &ncon, &ncommonnodes,
            &nparts, tpwgts, &ubvec, options, &edgecut,
            part, &commWorld);
    
    delete[] eptr;
    delete[] eind;

    // Distribute the part list to each proc
    // Each proc has an array of length globalCellCount which says the processor number that cell belongs to [cellMap]
    int recvCounts[NoProc];
    int displs[NoProc];
    for (int p = 0; p < NoProc; ++p) {
        recvCounts[p] = baseCellCount;
        displs[p] = p*baseCellCount;
    }
    recvCounts[NoProc - 1] = baseCellCount + globalCellCount - NoProc*baseCellCount;

    cellOwner.resize(globalCellCount);
    //cellMap of a cell returns which processor it is assigned to
    MPI_Allgatherv(part, cellCount, MPI_INT, &cellOwner[0], recvCounts, displs, MPI_INT, MPI_COMM_WORLD);

    // Find new local cellCount after ParMetis distribution
    cellCount = 0.;
    int otherCellCounts[NoProc];
    for (int p = 0; p < NoProc; p++) {
        otherCellCounts[p] = 0;
        partitionOffset.push_back(0);
    }

    for (int c = 0; c < globalCellCount; ++c) {
        otherCellCounts[cellOwner[c]] += 1;
        if (cellOwner[c] == Rank)
            cellCount++;
    }
    CommMPI_Info(Rank, "Number of Cells = %d", cellCount);

    myOffset = 0;
    partitionOffset[0] = 0;
    for (int p = 1; p < NoProc; ++p)
        partitionOffset[p] = partitionOffset[p - 1] + otherCellCounts[p - 1];
    myOffset = partitionOffset[Rank];

    delete[] part;
}

//------------------------------------------------------------------------------
//! Create Mesh Duals Using Parmetis
//------------------------------------------------------------------------------
void Grid3D::Parmetis_Create_Mesh2Dual() {
    // Find out other partition's cell counts
    int otherCellCounts[NoProc];
    for (int i = 0; i < NoProc; i++)
        otherCellCounts[i] = 0;
    for (int c = 0; c < globalCellCount; c++)
        otherCellCounts[cellOwner[c]] += 1;

    // Create the Mesh2Dual inputs
    idxtype elmdist[NoProc + 1];
    idxtype *eptr;
    eptr = new idxtype[cellCount + 1];
    idxtype *eind;
    int eindSize = 0;
    int ncommonnodes = 1;
    // C-style numbering
    int numflag = 0;
    MPI_Comm commWorld = MPI_COMM_WORLD;

    for (int c = 0; c < cellCount; c++)
        eindSize += cell[c].nodeCount;
    eind = new idxtype[eindSize];

    elmdist[0] = 0;
    for (int p = 1; p <= NoProc; p++)
        elmdist[p] = otherCellCounts[p - 1] + elmdist[p - 1];

    eptr[0] = 0;
    for (int c = 1; c <= cellCount; c++)
        eptr[c] = eptr[c - 1] + cell[c - 1].nodeCount;
    
    int eindIndex = 0;
    for (int c = 0; c < cellCount; c++) {
        for (int cn = 0; cn < cell[c].nodeCount; cn++) {
            eind[eindIndex] = node[cell[c].nodes[cn]].globalId;
            eindIndex++;
        }
    }
    
    ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommonnodes,
            &adjIndex, &adjacency, &commWorld);

    delete[] eptr;
    delete[] eind;
}

//------------------------------------------------------------------------------
//! Partition the Mesh and Create Local Connectivity Graphs
//------------------------------------------------------------------------------
void Grid3D::Partition_And_Create_Connectivity() {
    // Partition the Mesh using Parmetis
    Parmetis_Partition_DB();
    // Create the Local Connectitity Graphs
    Create_Local_Connectivity_Maps();
    // Free the Memory
    Delete_Global_DataStructure();
}

//------------------------------------------------------------------------------
//! Computes various required Metrices: Area, Volume,
//------------------------------------------------------------------------------
void Grid3D::Compute_Grid_Metrics() {
    int nodeid[8];
    Vector3D vtemp1, vtemp2;
    Point3D  centroid;
    Point3D  patchCentroid;
    
    // Now loop through faces and calculate centroids and areas
    for (int f = 0; f < faceCount; f++) {
        Vector3D areaVec  = 0.0;
        Vector3D patchArea;

        // Find an approxiamate centroid (true centroid for triangle)
        centroid = 0.0;
        for (int n = 0; n < face[f].nodeCount; n++) {
            nodeid[n] = face[f].nodes[n];
            centroid += node[nodeid[n]];
        }
        centroid /= double(face[f].nodeCount);

        // Sum the area as a patch of triangles formed by connecting two nodes
        // and an interior point
        face[f].centroid = 0.0;
        areaVec = 0.0;

        // Calculate face normal
        vtemp1.vec[0] = (node[nodeid[2]] - node[nodeid[1]]).pos[0];
        vtemp1.vec[1] = (node[nodeid[2]] - node[nodeid[1]]).pos[1];
        vtemp1.vec[2] = (node[nodeid[2]] - node[nodeid[1]]).pos[2];
        vtemp2.vec[0] = (node[nodeid[0]] - node[nodeid[1]]).pos[0];
        vtemp2.vec[1] = (node[nodeid[0]] - node[nodeid[1]]).pos[1];
        vtemp2.vec[2] = (node[nodeid[0]] - node[nodeid[1]]).pos[2];
        face[f].normal = vtemp1.cross(vtemp2);
        face[f].normal.normalize();

        int next;
        for (int n = 0; n < face[f].nodeCount; n++) {
            next = n + 1;
            if (next == face[f].nodeCount)
                next = 0;
            vtemp1.vec[0] = node[nodeid[n]].pos[0] - centroid.pos[0];
            vtemp1.vec[1] = node[nodeid[n]].pos[1] - centroid.pos[1];
            vtemp1.vec[2] = node[nodeid[n]].pos[2] - centroid.pos[2];
            vtemp2.vec[0] = node[nodeid[next]].pos[0] - centroid.pos[0];
            vtemp2.vec[1] = node[nodeid[next]].pos[1] - centroid.pos[1];
            vtemp2.vec[2] = node[nodeid[next]].pos[2] - centroid.pos[2];
            patchArea = 0.5 * vtemp1.cross(vtemp2);
            patchCentroid = (1.0/3.0) * (node[nodeid[n]] + node[nodeid[next]] + centroid);
            face[f].centroid += patchCentroid * patchArea.dot(face[f].normal);
            areaVec += patchArea;
        }
        face[f].area = areaVec.magnitude();
        face[f].centroid /= face[f].area;
    }

    // Loop through the cells and calculate the volumes
    double totalVolume = 0.0;
    for (int c = 0; c < cellCount; c++) {
        double volume = 0.0;
        double patchVolume = 0.0;
        Vector3D height, basePatchArea;
        int f, next;
        // Calculate cell centroid
        // First calculate an approximate one
        centroid = 0.0;
        for (int cn = 0; cn < cell[c].nodeCount; cn++) {
            nodeid[cn] = cell[c].nodes[cn];
            centroid += node[nodeid[cn]];
        }
        centroid /= double(cell[c].nodeCount);

        // Break-up the volume into tetrahedras and add the volumes
        // Calculate the centroid of the cell by taking moments of each tetra
        cell[c].centroid = 0.0;
        double sign;
        for (int cf = 0; cf < cell[c].faceCount; cf++) {
            f = cell[c].faces[cf];
            // Every cell face is broken to triangles
            for (int n = 0; n < face[f].nodeCount; n++) {
                next = n + 1;
                if (next == face[f].nodeCount)
                    next = 0;
                // Triangle area
                vtemp1.vec[0] = node[face[f].nodes[n]].pos[0] - face[f].centroid.pos[0];
                vtemp1.vec[1] = node[face[f].nodes[n]].pos[1] - face[f].centroid.pos[1];
                vtemp1.vec[2] = node[face[f].nodes[n]].pos[2] - face[f].centroid.pos[2];
                vtemp2.vec[0] = node[face[f].nodes[next]].pos[0] - face[f].centroid.pos[0];
                vtemp2.vec[1] = node[face[f].nodes[next]].pos[1] - face[f].centroid.pos[1];
                vtemp2.vec[2] = node[face[f].nodes[next]].pos[2] - face[f].centroid.pos[2];
                basePatchArea = 0.5 * vtemp1.cross(vtemp2);
                // Height of the tetrahedra
                vtemp1.vec[0] = face[f].centroid.pos[0] - centroid.pos[0];
                vtemp1.vec[1] = face[f].centroid.pos[1] - centroid.pos[1];
                vtemp1.vec[2] = face[f].centroid.pos[2] - centroid.pos[2];
                height = vtemp1.dot(face[f].normal) * face[f].normal;
                // Fix face orientation issue
                sign = -1.;
                if (face[f].parent == c)
                    sign = 1.0;
                patchVolume = sign * basePatchArea.dot(height) / 3.0;
                patchCentroid = 0.25 * (face[f].centroid + node[face[f].nodes[n]] + node[face[f].nodes[next]] + centroid);
                cell[c].centroid += patchVolume*patchCentroid;
                if (patchVolume < 0.0)
                    CommMPI_Error(Rank, "Grid3D::Compute_Grid_Metrics: cell[%d] %s", c, "Negative Volume Encountered !");
                volume += patchVolume;
            }
        }
        cell[c].volume = volume;
        cell[c].centroid /= volume;
        totalVolume += volume;
    }

    double globalTotalVolume = 0.;
    MPI_Allreduce(&totalVolume, &globalTotalVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (NoProc > 1)
        CommMPI_Info(RankIO, "Global Total Volume = %lf", globalTotalVolume);

    for (int f = 0; f < faceCount; f++) {
        vtemp1.vec[0] = face[f].centroid.pos[0] - cell[face[f].parent].centroid.pos[0];
        vtemp1.vec[1] = face[f].centroid.pos[1] - cell[face[f].parent].centroid.pos[1];
        vtemp1.vec[2] = face[f].centroid.pos[2] - cell[face[f].parent].centroid.pos[2];
        if (face[f].normal.dot(vtemp1) <= 0.0) {
            CommMPI_Warn(Rank, "Face[%d]: Negative Normal Detected", f);
            face[f].normal *= -1.0;
            // Reverse the Node Order
            face[f].nodes.assign(face[f].nodes.rbegin(), face[f].nodes.rend());
        }
    }
}

