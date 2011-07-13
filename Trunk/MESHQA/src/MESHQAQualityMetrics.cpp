/*******************************************************************************
 * File:        MESHQAQualityMetrics.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQAQualityMetrics.h"

//! Egde Quality
double MESHQAQualityMetrics::analyze_Edge2(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 2;
    
    switch (qt) {
        case MESHQAEnums::QT_Length:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_edge_length(num_nodes, coordinates);
            #else
            qvalue = v_edge_length(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    return qvalue;
}

//! Triangle Quality
double MESHQAQualityMetrics::analyze_Tri3(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 3;
    
    switch (qt) {
        case MESHQAEnums::QT_Aspect:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_aspect_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tri_aspect_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Area:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_area(num_nodes, coordinates);
            #else
            qvalue = v_tri_area(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_SmallestAngle:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_minimum_angle(num_nodes, coordinates);
            #else
            qvalue = v_tri_minimum_angle(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_LargestAngle:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_maximum_angle(num_nodes, coordinates);
            #else
            qvalue = v_tri_maximum_angle(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Condition:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_condition(num_nodes, coordinates);
            #else
            qvalue = v_tri_condition(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_NormalizedJacobian:
            #ifdef HAVE_VTK_VERDICT
             qvalue = vtk_v_tri_scaled_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_tri_scaled_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shape:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_shape(num_nodes, coordinates);
            #else
            qvalue = v_tri_shape(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RelativeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_relative_size_squared(num_nodes, coordinates);
            #else
            qvalue = v_tri_relative_size_squared(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShapeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_shape_and_size(num_nodes, coordinates);
            #else
            qvalue = v_tri_shape_and_size(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_EdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tri_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_AspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_tri_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Distortion:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_distortion(num_nodes, coordinates);
            #else
            qvalue = v_tri_distortion(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RadiusRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tri_radius_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tri_radius_ratio(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    return qvalue;
}

//! Qualdrilateral Quality
double MESHQAQualityMetrics::analyze_Quad4(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 4;
    
    switch (qt) {
        case MESHQAEnums::QT_Aspect:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_aspect_ratio(num_nodes, coordinates);
            #else
            qvalue = v_quad_aspect_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Area:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_area(num_nodes, coordinates);
            #else
            qvalue = v_quad_area(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_SmallestAngle:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_minimum_angle(num_nodes, coordinates);
            #else
            qvalue = v_quad_minimum_angle(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_LargestAngle:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_maximum_angle(num_nodes, coordinates);
            #else
            qvalue = v_quad_maximum_angle(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Condition:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_condition(num_nodes, coordinates);
            #else
            qvalue = v_quad_condition(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Jacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_quad_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_NormalizedJacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_scaled_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_quad_scaled_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shear:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_shear(num_nodes, coordinates);
            #else
            qvalue = v_quad_shear(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shape:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_shape(num_nodes, coordinates);
            #else
            qvalue = v_quad_shape(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RelativeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_relative_size_squared(num_nodes, coordinates);
            #else
            qvalue = v_quad_relative_size_squared(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShapeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_shape_and_size(num_nodes, coordinates);
            #else
            qvalue = v_quad_shape_and_size(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Skew:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_skew(num_nodes, coordinates);
            #else
            qvalue = v_quad_skew(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Taper:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_taper(num_nodes, coordinates);
            #else
            qvalue = v_quad_taper(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Warpage:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_warpage(num_nodes, coordinates);
            #else
            qvalue = v_quad_warpage(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Stretch:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_stretch(num_nodes, coordinates);
            #else
            qvalue = v_quad_stretch(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Oddy:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_oddy(num_nodes, coordinates);
            #else
            qvalue = v_quad_oddy(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_EdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_quad_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MaxEdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_max_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_quad_max_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MedAspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_med_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_quad_med_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MaxAspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_max_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_quad_max_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Distortion:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_distortion(num_nodes, coordinates);
            #else
            qvalue = v_quad_distortion(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShearSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_shear_and_size(num_nodes, coordinates);
            #else
            qvalue = v_quad_shear_and_size(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RadiusRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_quad_radius_ratio(num_nodes, coordinates);
            #else
            qvalue = v_quad_radius_ratio(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    
    return qvalue;
}

//! Tetra Quality
double MESHQAQualityMetrics::analyze_Tetra4(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 4;

    switch (qt) {
        case MESHQAEnums::QT_Aspect:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_aspect_beta(num_nodes, coordinates);
            #else
            qvalue = v_tet_aspect_beta(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_AspectGamma:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_aspect_gamma(num_nodes, coordinates);
            #else
            qvalue = v_tet_aspect_gamma(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_SmallestAngle:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_minimum_angle(num_nodes, coordinates);
            #else
            qvalue = v_tet_minimum_angle(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Condition:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_condition(num_nodes, coordinates);
            #else
            qvalue = v_tet_condition(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Jacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_tet_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_NormalizedJacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_scaled_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_tet_scaled_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shape:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_shape(num_nodes, coordinates);
            #else
            qvalue = v_tet_shape(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RelativeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_relative_size_squared(num_nodes, coordinates);
            #else
            qvalue = v_tet_relative_size_squared(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShapeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_shape_and_size(num_nodes, coordinates);
            #else
            qvalue = v_tet_shape_and_size(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Volume:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_volume(num_nodes, coordinates);
            #else
            qvalue = v_tet_volume(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_EdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tet_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_AspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_tet_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Distortion:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_distortion(num_nodes, coordinates);
            #else
            qvalue = v_tet_distortion(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RadiusRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_radius_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tet_radius_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_AspectBeta:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_aspect_beta(num_nodes, coordinates);
            #else
            qvalue = v_tet_aspect_beta(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_CollapseRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_tet_collapse_ratio(num_nodes, coordinates);
            #else
            qvalue = v_tet_collapse_ratio(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    
    return qvalue;
}

//! Pyramid Quality
double MESHQAQualityMetrics::analyze_Pyra5(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 5;
    
    switch (qt) {
        case MESHQAEnums::QT_Volume:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_pyramid_volume(num_nodes, coordinates);
            #else
            qvalue = v_pyramid_volume(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    return qvalue;
}
//! Prism Quality
double MESHQAQualityMetrics::analyze_Prism6(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 6;

    switch (qt) {
        case MESHQAEnums::QT_Volume:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_wedge_volume(num_nodes, coordinates);
            #else
            qvalue = v_wedge_volume(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    return qvalue;
}

//! Hexa Quality
double MESHQAQualityMetrics::analyze_Hexa8(MESHQAEnums::QualityType qt, double coordinates[][3])
{
    double qvalue = 0.0;
    int num_nodes = 6;

    switch (qt) {
        case MESHQAEnums::QT_Condition:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_condition(num_nodes, coordinates);
            #else
            qvalue = v_hex_condition(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Jacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_hex_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_NormalizedJacobian:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_scaled_jacobian(num_nodes, coordinates);
            #else
            qvalue = v_hex_scaled_jacobian(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shear:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_shear(num_nodes, coordinates);
            #else
            qvalue = v_hex_shear(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Shape:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_shape(num_nodes, coordinates);
            #else
            qvalue = v_hex_shape(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_RelativeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_relative_size_squared(num_nodes, coordinates);
            #else
            qvalue = v_hex_relative_size_squared(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShapeSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_shape_and_size(num_nodes, coordinates);
            #else
            qvalue = v_hex_shape_and_size(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Skew:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_skew(num_nodes, coordinates);
            #else
            qvalue = v_hex_skew(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Taper:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_taper(num_nodes, coordinates);
            #else
            qvalue = v_hex_taper(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Stretch:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_stretch(num_nodes, coordinates);
            #else
            qvalue = v_hex_stretch(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Oddy:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_oddy(num_nodes, coordinates);
            #else
            qvalue = v_hex_oddy(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Volume:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_volume(num_nodes, coordinates);
            #else
            qvalue = v_hex_volume(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Diagonal:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_diagonal(num_nodes, coordinates);
            #else
            qvalue = v_hex_diagonal(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Dimension:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_dimension(num_nodes, coordinates);
            #else
            qvalue = v_hex_dimension(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_EdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_hex_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MaxEdgeRatio:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_max_edge_ratio(num_nodes, coordinates);
            #else
            qvalue = v_hex_max_edge_ratio(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MedAspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_med_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_hex_med_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_MaxAspectFrobenius:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_max_aspect_frobenius(num_nodes, coordinates);
            #else
            qvalue = v_hex_max_aspect_frobenius(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_Distortion:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_distortion(num_nodes, coordinates);
            #else
            qvalue = v_hex_distortion(num_nodes, coordinates);
            #endif
            break;
        case MESHQAEnums::QT_ShearSize:
            #ifdef HAVE_VTK_VERDICT
            qvalue = vtk_v_hex_shear_and_size(num_nodes, coordinates);
            #else
            qvalue = v_hex_shear_and_size(num_nodes, coordinates);
            #endif
            break;
        default:
            qvalue = 0.0;
    }
    return qvalue;
}

