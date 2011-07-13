/*******************************************************************************
 * File:        Grid3D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _GRID3D_H
#define	_GRID3D_H

#include <vector>
#include <map>
#include <set>
#include <parmetis.h>

#include "Point3D.h"
#include "Vector3D.h"
#include "corestruct.h"

using namespace std;

class Node3D;
class Face3D;
class Cell3D;
class Ghost3D;

class Grid3D {
public:
    int             myOffset;        // TODO
    int             nodeCountOffset;
    int             nodeCount;
    int             cellCount;
    int             faceCount;
    int             ghostCount;
    int             globalNodeCount;
    int             globalCellCount;
    int             globalFaceCount;
    vector<int>     partitionOffset; // TODO
    vector<Node3D>  node;
    vector<Face3D>  face;
    vector<Cell3D>  cell;
    vector<Ghost3D> ghost;
    vector<int>         cellOwner;
    map<int, int>       nodeGlobal2Local;
    map<int, int>       cellGlobal2Local;
    map<int, int>       ghostGlobal2Local;
    vector<int>         nodeGlobal2Output;
    idxtype             *adjIndex;
    idxtype             *adjacency;
protected:
    ROOT        *gridDB;
    int         scaleFlag;
    int         rotateFlag;
    Vector3D    scaleFactor;
    Vector3D    rotationCenter;
    Vector3D    rotationAngle;
private:
    vector<double>      global_x;
    vector<double>      global_y;
    vector<double>      global_z;
    vector<int>         global_cellConnIndex;
    vector<int>         global_cellConnectivity;
    vector< set<int> >  global_bocoNodes;
    map<string, int>    global_bocoNameMap;
public:
    Grid3D();       // TODO
    ~Grid3D();      // TODO
    
    int Read_DB();   // Done
    void Set_Scale(Vector3D value);
    void Set_Rotation(Vector3D center, Vector3D angle);
    void Partition_And_Create_Connectivity();   // Done
    void Compute_Grid_Metrics();                // Done
    
//    void gradMaps();            // Check
//    void gradients();           // Check
//    void limit_gradients(void); // Check
//    void lengthScales(void);    // Check
//    void nodeAverages();                // TODO
//    void sortStencil(Node3D& n);        // TODO
//    void interpolate_tetra(Node3D& n);  // TODO
//    void interpolate_tri(Node3D& n);    // TODO
//    void interpolate_line(Node3D& n);   // TODO
//    void nodeAverages_idw();            // TODO
//    void faceAverages();                // TODO
protected:
    void Create_Global_DataStructure();         // Done
    void Delete_Global_DataStructure();         // Done
    void Create_Local_Nodes_And_Cells();        // Done
    void Create_Local_Connectivity_Node2Cell(); // Done
    void Create_Local_Connectivity_Cell2Cell(); // Done
    void Create_Local_Connectivity_Face();      // Done
    void Create_Local_Connectivity_Ghost();     // Done
    void Create_Local_Connectivity_Maps();      // Done

    void Scale_DB();                            // Done
    void Rotate_DB();                           // Done
    void Parmetis_Partition_DB();               // Done
    void Parmetis_Create_Mesh2Dual();           // Done
};

#endif	/* _GRID3D_H */

