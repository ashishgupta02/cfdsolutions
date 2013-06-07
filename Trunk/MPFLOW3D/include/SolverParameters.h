/*******************************************************************************
 * File:        SolverParameters.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"
#include "SolverUtils.h"

#ifndef _SOLVERPARAMETERS_H
#define	_SOLVERPARAMETERS_H

/*!
 * \brief Solver Boolean Options
 */
enum SOLVER_BOOLEAN {
    NO  = 0,
    YES = 1
};
static const map<string, SOLVER_BOOLEAN> SolverBooleanMap = CCreateMap<string, SOLVER_BOOLEAN>
("NO",  NO)
("YES", YES);

/*!
 * \brief Boundary Condition Method
 */
enum BC_METHOD {
    BC_METHOD_NONE               = -1,  /*!< \brief Boundary Condition Method None. */
    BC_METHOD_CHARCTERISTIC      =  0,  /*!< \brief Boundary Condition Method Charactersistic. */
    BC_METHOD_CHARCTERISTIC_WEAK =  1,  /*!< \brief Boundary Condition Method Weak Characteristic. */
    BC_METHOD_PRESSURE_WALL      =  2   /*!< \brief Boundary Condition With Pressure Wall */
};
static const map<string, BC_METHOD> BCMethodMap = CCreateMap<string, BC_METHOD>
("BC_METHOD_NONE",               BC_METHOD_NONE)
("BC_METHOD_CHARCTERISTIC",      BC_METHOD_CHARCTERISTIC)
("BC_METHOD_CHARCTERISTIC_WEAK", BC_METHOD_CHARCTERISTIC_WEAK)
("BC_METHOD_PRESSURE_WALL",      BC_METHOD_PRESSURE_WALL);

/*!
 * \brief Boundary Conditions Type
 */
enum BC_TYPE {
    BC_TYPE_NONE                         = -1,  /*!< \brief No Boundary Condition definition. */
    BC_TYPE_EULER_WALL                   =  0,  /*!< \brief Boundary Euler Solid Wall definition. */
    BC_TYPE_FAR_FIELD                    =  1,  /*!< \brief Boundary Far-Field / Free-Stream definition. */
    BC_TYPE_INFLOW                       =  2,  /*!< \brief Boundary Inlet Flow definition. */
    BC_TYPE_OUTFLOW                      =  3,  /*!< \brief Boundary Outlet Flow definition. */
    BC_TYPE_SUBSONIC_INFLOW              =  4,  /*!< \brief Boundary Subsonic Inlet Flow definition. */
    BC_TYPE_SUBSONIC_OUTFLOW             =  5,  /*!< \brief Boundary Subsonic Outlet Flow definition. */
    BC_TYPE_SUPERSONIC_INFLOW            =  6,  /*!< \brief Boundary Supersonic Inlet Flow definition. */
    BC_TYPE_SUPERSONIC_OUTFLOW           =  7,  /*!< \brief Boundary Supersonic Outlet Flow definition. */
    BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW    =  8,  /*!< \brief Boundary Subsonic Inlet Far-Field / Free-Stream Flow definition. */
    BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW   =  9,  /*!< \brief Boundary Subsonic Outlet Far-Field / Free-Stream Flow definition. */
    BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW  = 10,  /*!< \brief Boundary Supersonic Inlet Far-Field / Free-Stream Flow definition. */
    BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW = 11,  /*!< \brief Boundary Supersonic Outlet Far-Field / Free-Stream Flow definition. */
    BC_TYPE_NO_SLIP_WALL                 = 12,  /*!< \brief Boundary No Slip Wall definition. */
    BC_TYPE_SYMMETRY_PLANE               = 13,  /*!< \brief Boundary Symmetry Plane definition. */
    BC_TYPE_PERIODIC                     = 14,  /*!< \brief Boundary Periodic definition. */
    BC_TYPE_DIRICHLET                    = 15,  /*!< \brief Boundary Dirichlet definition. */
    BC_TYPE_NEUMANN                      = 16   /*!< \brief Boundary Neumann definition. */
};
static const map<string, BC_TYPE> BCTypeMap = CCreateMap<string, BC_TYPE>
("BC_TYPE_NONE",                         BC_TYPE_NONE)
("BC_TYPE_EULER_WALL",                   BC_TYPE_EULER_WALL)
("BC_TYPE_FAR_FIELD",                    BC_TYPE_FAR_FIELD)
("BC_TYPE_INFLOW",                       BC_TYPE_INFLOW)
("BC_TYPE_OUTFLOW",                      BC_TYPE_OUTFLOW)
("BC_TYPE_SUBSONIC_INFLOW",              BC_TYPE_SUBSONIC_INFLOW)
("BC_TYPE_SUBSONIC_OUTFLOW",             BC_TYPE_SUBSONIC_OUTFLOW)
("BC_TYPE_SUPERSONIC_INFLOW",            BC_TYPE_SUPERSONIC_INFLOW)
("BC_TYPE_SUPERSONIC_OUTFLOW",           BC_TYPE_SUPERSONIC_OUTFLOW)
("BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW",    BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW)
("BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW",   BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW)
("BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW",  BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW)
("BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW", BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW)
("BC_TYPE_NO_SLIP_WALL",                 BC_TYPE_NO_SLIP_WALL)
("BC_TYPE_SYMMETRY_PLANE",               BC_TYPE_SYMMETRY_PLANE)
("BC_TYPE_PERIODIC",                     BC_TYPE_PERIODIC)
("BC_TYPE_DIRICHLET",                    BC_TYPE_DIRICHLET)
("BC_TYPE_NEUMANN",                      BC_TYPE_NEUMANN);

/*!
 * \brief Solver Method
 */
enum SOLVER_METHOD {
    SOLVER_METHOD_NONE     = -1,    /*!< \brief Solver Method None. */
    SOLVER_METHOD_STEADY   =  0,    /*!< \brief Solver Method Steady. */
    SOLVER_METHOD_UNSTEADY =  1     /*!< \brief Solver Method Steady. */
};
static const map<string, SOLVER_METHOD> SolverMethodMap = CCreateMap<string, SOLVER_METHOD>
("SOLVER_METHOD_NONE",     SOLVER_METHOD_NONE)
("SOLVER_METHOD_STEADY",   SOLVER_METHOD_STEADY)
("SOLVER_METHOD_UNSTEADY", SOLVER_METHOD_UNSTEADY);

/*!
 * \brief Solver Scheme
 */
enum SOLVER_SCHEME {
    SOLVER_SCHEME_NONE     = -1,    /*!< \brief Solver Scheme None. */
    SOLVER_SCHEME_EXPLICIT =  0,    /*!< \brief Solver Scheme Explicit. */
    SOLVER_SCHEME_IMPLICIT =  1     /*!< \brief Solver Scheme Implicit. */
};
static const map<string, SOLVER_SCHEME> SolverSchemeMap = CCreateMap<string, SOLVER_SCHEME>
("SOLVER_SCHEME_NONE",     SOLVER_SCHEME_NONE)
("SOLVER_SCHEME_EXPLICIT", SOLVER_SCHEME_EXPLICIT)
("SOLVER_SCHEME_IMPLICIT", SOLVER_SCHEME_IMPLICIT);

/*!
 * \brief Flux Scheme Type
 */
enum FLUX_SCHEME {
    FLUX_SCHEME_NONE          = -1, /*!< \brief Flux Scheme None. */
    FLUX_SCHEME_ROE           =  0, /*!< \brief Flux Roe Scheme. */
    FLUX_SCHEME_HLLC          =  1, /*!< \brief Flux HLLC Scheme. */
    FLUX_SCHEME_AUSM          =  2, /*!< \brief Flux AUSM Scheme. */
    FLUX_SCHEME_VANLEER       =  3, /*!< \brief Flux VanLeer Scheme. */
    FLUX_SCHEME_LDFSS         =  4, /*!< \brief Flux LDFSS Scheme. */
    FLUX_SCHEME_OSHER         =  5, /*!< \brief Flux Osher Scheme. */
    FLUX_SCHEME_STEGERWARMING =  6, /*!< \brief Flux Steger-Warming Scheme. */
    FLUX_SCHEME_JST           =  7  /*!< \brief Flux JST Scheme. */
};
static const map<string, FLUX_SCHEME> FluxSchemeMap = CCreateMap<string, FLUX_SCHEME>
("FLUX_SCHEME_NONE",          FLUX_SCHEME_NONE)
("FLUX_SCHEME_ROE",           FLUX_SCHEME_ROE)
("FLUX_SCHEME_HLLC",          FLUX_SCHEME_HLLC)
("FLUX_SCHEME_AUSM",          FLUX_SCHEME_AUSM)
("FLUX_SCHEME_VANLEER",       FLUX_SCHEME_VANLEER)
("FLUX_SCHEME_LDFSS",         FLUX_SCHEME_LDFSS)
("FLUX_SCHEME_OSHER",         FLUX_SCHEME_OSHER)
("FLUX_SCHEME_STEGERWARMING", FLUX_SCHEME_STEGERWARMING)
("FLUX_SCHEME_JST",           FLUX_SCHEME_JST);

/*!
 * \brief Precondition Method
 */
enum PRECOND_METHOD {
    PRECOND_METHOD_NONE     = -1,   /*!< \brief Solver Preconditioning None. */
    PRECOND_METHOD_LMFIX    =  0,   /*!< \brief Solver Preconditioning Roe LMFix */
    PRECOND_METHOD_THORNBER =  1,   /*!< \brief Solver Preconditioning Roe THORNBER */ 
    PRECOND_METHOD_BTW      =  2,   /*!< \brief Solver Preconditioning Roe Briley Taylor Whitfield */
    PRECOND_METHOD_ERIKSSON =  3,   /*!< \brief Solver Preconditioning Roe Eriksson */
    PRECOND_METHOD_MERKEL   =  4,   /*!< \brief Solver Preconditioning Roe Merkel */
    PRECOND_METHOD_TURKEL   =  5    /*!< \brief Solver Preconditioning Roe Turkel*/
};
static const map<string, PRECOND_METHOD> PrecondMethodMap = CCreateMap<string, PRECOND_METHOD>
("PRECOND_METHOD_NONE",     PRECOND_METHOD_NONE)
("PRECOND_METHOD_LMFIX",    PRECOND_METHOD_LMFIX)
("PRECOND_METHOD_THORNBER", PRECOND_METHOD_THORNBER)
("PRECOND_METHOD_BTW",      PRECOND_METHOD_BTW)
("PRECOND_METHOD_ERIKSSON", PRECOND_METHOD_ERIKSSON)
("PRECOND_METHOD_MERKEL",   PRECOND_METHOD_MERKEL)
("PRECOND_METHOD_TURKEL",   PRECOND_METHOD_TURKEL);

/*!
 * \brief Precondition Type
 */
enum PRECOND_TYPE {
    PRECOND_TYPE_NONE     = -1,   /*!< \brief Solver Preconditioning Type None. */
    PRECOND_TYPE_LOCAL    =  0,   /*!< \brief Solver Preconditioning Type Local. */
    PRECOND_TYPE_GLOBAL   =  1,    /*!< \brief Solver Preconditioning Type Global. */
    PRECOND_TYPE_URLOCAL  =  2    /*!< \brief Solver Preconditioning Type Un-Restricted Local. */     
};
static const map<string, PRECOND_TYPE> PrecondTypeMap = CCreateMap<string, PRECOND_TYPE>
("PRECOND_TYPE_NONE",   PRECOND_TYPE_NONE)
("PRECOND_TYPE_LOCAL",  PRECOND_TYPE_LOCAL)
("PRECOND_TYPE_GLOBAL", PRECOND_TYPE_GLOBAL)
("PRECOND_TYPE_URLOCAL", PRECOND_TYPE_URLOCAL);

/*!
 * \brief Residual Smoother Method
 */
enum RESIDUAL_SMOOTH_METHOD {
    RESIDUAL_SMOOTH_METHOD_NONE              = -1,  /*!< \brief Residual Smoother Method None. */
    RESIDUAL_SMOOTH_METHOD_EXPLICIT          =  0,  /*!< \brief Residual Smoother Method Explicit. */
    RESIDUAL_SMOOTH_METHOD_EXPLICIT_WEIGHTED =  1,  /*!< \brief Residual Smoother Method Explicit Weighted. */
    RESIDUAL_SMOOTH_METHOD_IMPLICIT          =  2,  /*!< \brief Residual Smoother Method Implicit. */
    RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT =  3   /*!< \brief Residual Smoother Method Implicit Constant . */
};
static const map<string, RESIDUAL_SMOOTH_METHOD> ResidualSmoothMethodMap = CCreateMap<string, RESIDUAL_SMOOTH_METHOD>
("RESIDUAL_SMOOTH_METHOD_NONE",              RESIDUAL_SMOOTH_METHOD_NONE)
("RESIDUAL_SMOOTH_METHOD_EXPLICIT",          RESIDUAL_SMOOTH_METHOD_EXPLICIT)
("RESIDUAL_SMOOTH_METHOD_EXPLICIT_WEIGHTED", RESIDUAL_SMOOTH_METHOD_EXPLICIT_WEIGHTED)
("RESIDUAL_SMOOTH_METHOD_IMPLICIT",          RESIDUAL_SMOOTH_METHOD_IMPLICIT)
("RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT", RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT);

/*!
 * \brief Residual Smoother Type
 */
enum RESIDUAL_SMOOTH_TYPE {
    RESIDUAL_SMOOTH_TYPE_NONE   = -1,   /*!< \brief Residual Smoother Type None. */
    RESIDUAL_SMOOTH_TYPE_LOCAL  =  0,   /*!< \brief Residual Smoother Type Local. */
    RESIDUAL_SMOOTH_TYPE_GLOBAL =  1    /*!< \brief Residual Smoother Type Global. */
};
static const map<string, RESIDUAL_SMOOTH_TYPE> ResidualSmoothTypeMap = CCreateMap<string, RESIDUAL_SMOOTH_TYPE>
("RESIDUAL_SMOOTH_TYPE_NONE",   RESIDUAL_SMOOTH_TYPE_NONE)
("RESIDUAL_SMOOTH_TYPE_LOCAL",  RESIDUAL_SMOOTH_TYPE_LOCAL)
("RESIDUAL_SMOOTH_TYPE_GLOBAL", RESIDUAL_SMOOTH_TYPE_GLOBAL);

/*!
 * \brief Residual Smoother Scheme
 */
enum RESIDUAL_SMOOTH_SCHEME {
    RESIDUAL_SMOOTH_SCHEME_NONE                   = -1, /*!< \brief Residual Smoother Scheme None. */
    RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION =  0, /*!< \brief Residual Smoother Scheme both Convective and Dissipation. */
    RESIDUAL_SMOOTH_SCHEME_ONLY_DISSIPATION       =  1  /*!< \brief Residual Smoother Scheme only Dissipation. */
};
static const map<string, RESIDUAL_SMOOTH_SCHEME> ResidualSmoothSchemeMap = CCreateMap<string, RESIDUAL_SMOOTH_SCHEME>
("RESIDUAL_SMOOTH_SCHEME_NONE",                   RESIDUAL_SMOOTH_SCHEME_NONE)
("RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION", RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION)
("RESIDUAL_SMOOTH_SCHEME_ONLY_DISSIPATION",       RESIDUAL_SMOOTH_SCHEME_ONLY_DISSIPATION);

/*!
 * \brief Solver Time Integration Method
 */
enum TIME_INTEGRATION_METHOD {
    TIME_INTEGRATION_METHOD_NONE                  = -1, /*!< \brief Time Integration Method None. */
    TIME_INTEGRATION_METHOD_EXPLICIT_EULER        =  0, /*!< \brief Time Integration Method Explicit Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK3          =  1, /*!< \brief Time Integration Method Explicit RK4. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK4          =  2, /*!< \brief Time Integration Method Explicit RK4. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5          =  3, /*!< \brief Time Integration Method Explicit RK5. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ       =  4, /*!< \brief Time Integration Method Explicit RK5 Martinelli and Jameson */
    TIME_INTEGRATION_METHOD_EXPLICIT_EULER_EULER  =  5, /*!< \brief Time Integration Method Explicit Dual Time Euler Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK3_EULER    =  6, /*!< \brief Time Integration Method Explicit Dual Time RK3 Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK4_EULER    =  7, /*!< \brief Time Integration Method Explicit Dual Time RK4 Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_EULER    =  8, /*!< \brief Time Integration Method Explicit Dual Time RK5 Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_EULER =  9, /*!< \brief Time Integration Method Explicit Dual Time RK5 Martinelli and Jameson Euler. */
    TIME_INTEGRATION_METHOD_EXPLICIT_EULER_BDF    = 10, /*!< \brief Time Integration Method Explicit Dual Time Euler BDF. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK3_BDF      = 11, /*!< \brief Time Integration Method Explicit Dual Time RK3 BDF. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK4_BDF      = 12, /*!< \brief Time Integration Method Explicit Dual Time RK4 BDF. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_BDF      = 13, /*!< \brief Time Integration Method Explicit Dual Time RK5 BDF. */
    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_BDF   = 14, /*!< \brief Time Integration Method Explicit Dual Time RK5 Martinelli and Jameson BDF. */
    TIME_INTEGRATION_METHOD_IMPLICIT_EULER        = 15, /*!< \brief Time Integration Method Implicit Euler Newton. */
    TIME_INTEGRATION_METHOD_IMPLICIT_EULER_EULER  = 16, /*!< \brief Time Integration Method Implicit Dual Time Euler Euler Newton. */
    TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF    = 17  /*!< \brief Time Integration Method Implicit Dual Time Euler BDF Newton. */
};
static const map<string, TIME_INTEGRATION_METHOD> TimeIntegrationMethodMap = CCreateMap<string, TIME_INTEGRATION_METHOD>
("TIME_INTEGRATION_METHOD_NONE",                  TIME_INTEGRATION_METHOD_NONE)
("TIME_INTEGRATION_METHOD_EXPLICIT_EULER",        TIME_INTEGRATION_METHOD_EXPLICIT_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK3",          TIME_INTEGRATION_METHOD_EXPLICIT_RK3)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK4",          TIME_INTEGRATION_METHOD_EXPLICIT_RK4)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5",          TIME_INTEGRATION_METHOD_EXPLICIT_RK5)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ",       TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ)
("TIME_INTEGRATION_METHOD_EXPLICIT_EULER_EULER",  TIME_INTEGRATION_METHOD_EXPLICIT_EULER_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK3_EULER",    TIME_INTEGRATION_METHOD_EXPLICIT_RK3_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK4_EULER",    TIME_INTEGRATION_METHOD_EXPLICIT_RK4_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5_EULER",    TIME_INTEGRATION_METHOD_EXPLICIT_RK5_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_EULER", TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_EULER)
("TIME_INTEGRATION_METHOD_EXPLICIT_EULER_BDF",    TIME_INTEGRATION_METHOD_EXPLICIT_EULER_BDF)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK3_BDF",      TIME_INTEGRATION_METHOD_EXPLICIT_RK3_BDF)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK4_BDF",      TIME_INTEGRATION_METHOD_EXPLICIT_RK4_BDF)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5_BDF",      TIME_INTEGRATION_METHOD_EXPLICIT_RK5_BDF)
("TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_BDF",   TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ_BDF)
("TIME_INTEGRATION_METHOD_IMPLICIT_EULER",        TIME_INTEGRATION_METHOD_IMPLICIT_EULER)
("TIME_INTEGRATION_METHOD_IMPLICIT_EULER_EULER",  TIME_INTEGRATION_METHOD_IMPLICIT_EULER_EULER)
("TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF",    TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF);


/*!
 * \brief Time Integration Type
 */
enum TIME_INTEGRATION_TYPE {
    TIME_INTEGRATION_TYPE_NONE   = -1,  /*!< \brief Time Integration Type None. */
    TIME_INTEGRATION_TYPE_LOCAL  =  0,  /*!< \brief Time Integration Type Local. */
    TIME_INTEGRATION_TYPE_GLOBAL =  1   /*!< \brief Time Integration Type Global. */
};
static const map<string, TIME_INTEGRATION_TYPE> TimeIntegrationTypeMap = CCreateMap<string, TIME_INTEGRATION_TYPE>
("TIME_INTEGRATION_TYPE_NONE",   TIME_INTEGRATION_TYPE_NONE)
("TIME_INTEGRATION_TYPE_LOCAL",  TIME_INTEGRATION_TYPE_LOCAL)
("TIME_INTEGRATION_TYPE_GLOBAL", TIME_INTEGRATION_TYPE_GLOBAL);

/*!
 * \brief Solver Jacobian Computation Method
 */
enum JACOBIAN_METHOD {
    JACOBIAN_METHOD_NONE      = -1, /*!< \brief Jacobian Compuatation Method None. */
    JACOBIAN_METHOD_CENTRAL   =  0, /*!< \brief Jacobian Compuatation Method Central Difference. */
    JACOBIAN_METHOD_FORWARD   =  1, /*!< \brief Jacobian Compuatation Method Forward Difference. */
    JACOBIAN_METHOD_BACKWARD  =  2, /*!< \brief Jacobian Compuatation Method Backward Difference. */
    JACOBIAN_METHOD_ALTERNATE =  3, /*!< \brief Jacobian Compuatation Method Forward then Backward Difference. */
    JACOBIAN_METHOD_APPROX    =  4, /*!< \brief Jacobian Compuatation Method Approximate. */
    JACOBIAN_METHOD_EXACT     =  5  /*!< \brief Jacobian Compuatation Method Central Exact. */
};
static const map<string, JACOBIAN_METHOD> JacobianMethodMap = CCreateMap<string, JACOBIAN_METHOD>
("JACOBIAN_METHOD_NONE",      JACOBIAN_METHOD_NONE)
("JACOBIAN_METHOD_CENTRAL",   JACOBIAN_METHOD_CENTRAL)
("JACOBIAN_METHOD_FORWARD",   JACOBIAN_METHOD_FORWARD)
("JACOBIAN_METHOD_BACKWARD",  JACOBIAN_METHOD_BACKWARD)
("JACOBIAN_METHOD_ALTERNATE", JACOBIAN_METHOD_ALTERNATE)
("JACOBIAN_METHOD_APPROX",    JACOBIAN_METHOD_APPROX)
("JACOBIAN_METHOD_EXACT",     JACOBIAN_METHOD_EXACT);

/*!
 * \brief Variable Type
 */
enum VARIABLE_TYPE {
    VARIABLE_NONE = -1, /*!< \brief Variable Type None. */
    VARIABLE_CON  =  0, /*!< \brief Variable Type Conservative. */
    VARIABLE_RUP  =  1, /*!< \brief Variable Type Density Vecocity Pressure. */
    VARIABLE_PUT  =  2, /*!< \brief Variable Type Presure Velocity Temperature. */
    VARIABLE_PUS  =  3, /*!< \brief Variable Type Pressure Velocity Entropy. */
    VARIABLE_RUT  =  4  /*!< \brief Variable Type Density Velocity Temperature. */
};
static const map<string, VARIABLE_TYPE> VariableTypeMap = CCreateMap<string, VARIABLE_TYPE>
("VARIABLE_NONE", VARIABLE_NONE)
("VARIABLE_CON",  VARIABLE_CON)
("VARIABLE_RUP",  VARIABLE_RUP)
("VARIABLE_PUT",  VARIABLE_PUT)
("VARIABLE_PUS",  VARIABLE_PUS)
("VARIABLE_RUT",  VARIABLE_RUT);

/*!
 * \brief Solver Order
 */
enum SOLVER_ORDER {
    SOLVER_ORDER_NONE   = -1,   /*!< \brief Solver Order None. */
    SOLVER_ORDER_FIRST  =  0,   /*!< \brief Solver First Order. */
    SOLVER_ORDER_SECOND =  1    /*!< \brief Solver Second Order. */
};
static const map<string, SOLVER_ORDER> SolverOrderMap = CCreateMap<string, SOLVER_ORDER>
("SOLVER_ORDER_NONE",   SOLVER_ORDER_NONE)
("SOLVER_ORDER_FIRST",  SOLVER_ORDER_FIRST)
("SOLVER_ORDER_SECOND", SOLVER_ORDER_SECOND);

/*!
 * \brief Gradient Computation Method
 */
enum GRADIENT_METHOD {
    GRADIENT_METHOD_NONE                   = -1,    /*!< \brief Gradient Method None. */
    GRADINET_METHOD_GREEN_GAUSS            =  0,    /*!< \brief Gradient Method Green Gauss. */
    GRADIENT_METHOD_LEAST_SQUARES          =  1,    /*!< \brief Gradient Method Least Squares. */
    GRADIENT_METHOD_LEAST_SQUARES_WEIGHTED =  2     /*!< \brief Gradient Method Weighted Least Squares. */
};
static const map<string, GRADIENT_METHOD> GradientMethodMap = CCreateMap<string, GRADIENT_METHOD>
("GRADIENT_METHOD_NONE",                   GRADIENT_METHOD_NONE)
("GRADINET_METHOD_GREEN_GAUSS",            GRADINET_METHOD_GREEN_GAUSS)
("GRADIENT_METHOD_LEAST_SQUARES",          GRADIENT_METHOD_LEAST_SQUARES)
("GRADIENT_METHOD_LEAST_SQUARES_WEIGHTED", GRADIENT_METHOD_LEAST_SQUARES_WEIGHTED);

/*!
 * \brief Flux Limiter Method
 */
enum LIMITER_METHOD {
    LIMITER_METHOD_NONE            = -1,    /*!< \brief Limiter Method None. */
    LIMITER_METHOD_BERTH_JESPERSEN =  0,    /*!< \brief Limiter Method Berth and Jespersen. */
    LIMITER_METHOD_VENKATAKRISHNAN =  1     /*!< \brief Limiter Method Venkatakrishnan. */
};
static const map<string, LIMITER_METHOD> LimiterMethodMap = CCreateMap<string, LIMITER_METHOD>
("LIMITER_METHOD_NONE",            LIMITER_METHOD_NONE)
("LIMITER_METHOD_BERTH_JESPERSEN", LIMITER_METHOD_BERTH_JESPERSEN)
("LIMITER_METHOD_VENKATAKRISHNAN", LIMITER_METHOD_VENKATAKRISHNAN);

/*!
 * \brief Limiter Order
 */
enum LIMITER_ORDER {
    LIMITER_ORDER_NONE   = -1,  /*!< \brief Limiter Order None. */
    LIMITER_ORDER_FIRST  =  0,  /*!< \brief Limiter First Order. */
    LIMITER_ORDER_SECOND =  1   /*!< \brief Limiter Second Order. */
};
static const map<string, LIMITER_ORDER> LimiterOrderMap = CCreateMap<string, LIMITER_ORDER>
("LIMITER_ORDER_NONE",   LIMITER_ORDER_NONE)
("LIMITER_ORDER_FIRST",  LIMITER_ORDER_FIRST)
("LIMITER_ORDER_SECOND", LIMITER_ORDER_SECOND);

/*!
 * \brief Material Type
 */
enum MATERIAL_TYPE {
    MATERIAL_TYPE_NONE     = -1,    /*!< \brief Material Type None. */
    MATERIAL_TYPE_IDEALGAS =  0,    /*!< \brief Material Type Ideal Gas. */
    MATERIAL_TYPE_NIST     =  1     /*!< \brief Material Type NIST. */
};
static const map<string, MATERIAL_TYPE> MaterialTypeMap = CCreateMap<string, MATERIAL_TYPE>
("MATERIAL_TYPE_NONE",     MATERIAL_TYPE_NONE)
("MATERIAL_TYPE_IDEALGAS", MATERIAL_TYPE_IDEALGAS)
("MATERIAL_TYPE_NIST",     MATERIAL_TYPE_NIST);

/*!
 * \brief Material Computation Type
 */
enum MATERIAL_COMP_TYPE {
    MATERIAL_COMP_TYPE_NONE = -1,   /*!< \brief Material Dimensional Computation Type None. */
    MATERIAL_COMP_TYPE_NDIM =  0,   /*!< \brief Material Non-Dimensional Computation. */
    MATERIAL_COMP_TYPE_DIM  =  1    /*!< \brief Material Dimensional Computation */
};
static const map<string, MATERIAL_COMP_TYPE> MaterialCompTypeMap = CCreateMap<string, MATERIAL_COMP_TYPE>
("MATERIAL_COMP_TYPE_NONE", MATERIAL_COMP_TYPE_NONE)
("MATERIAL_COMP_TYPE_NDIM", MATERIAL_COMP_TYPE_NDIM)
("MATERIAL_COMP_TYPE_DIM",  MATERIAL_COMP_TYPE_DIM);

/*!
 * \brief Averaging Type
 */
enum AVERAGE_TYPE {
    AVERAGE_TYPE_NONE          = -1,    /*!< \brief Average Computation Type None. */
    AVERAGE_TYPE_SIMPLE        =  0,    /*!< \brief Average Computation Simple */
    AVERAGE_TYPE_ROE           =  1,    /*!< \brief Average Computation Roe */
    AVERAGE_TYPE_SIMPLE_APPROX =  2,    /*!< \brief Average Computation Simple Approximate */
    AVERAGE_TYPE_ROE_APPROX    =  3     /*!< \brief Average Computation Roe Approximate */
};
static const map<string, AVERAGE_TYPE> AverageTypeMap = CCreateMap<string, AVERAGE_TYPE>
("AVERAGE_TYPE_NONE",          AVERAGE_TYPE_NONE)
("AVERAGE_TYPE_SIMPLE",        AVERAGE_TYPE_SIMPLE)
("AVERAGE_TYPE_ROE",           AVERAGE_TYPE_ROE)
("AVERAGE_TYPE_SIMPLE_APPROX", AVERAGE_TYPE_SIMPLE_APPROX)
("AVERAGE_TYPE_ROE_APPROX",    AVERAGE_TYPE_ROE_APPROX);

// No of Equations
#define NEQUATIONS                              5

// Mesh Parameters
extern char   MeshInputFilename[256];
extern int    MeshReorder;      /* Cuthill Mckee */

// Boundary Condition Parameters
extern char   BCInputFilename[256];

// Solution Parameters
extern char   SolutionOutputFilename[256];
extern char   RMSOutputFilename[256];

// Variable Type
extern int    VariableType;

// Solver Parameters   
extern int    SolverMethod;             /* Steady, Unsteady */
extern int    SolverScheme;             /* Explicit, Implicit */
extern int    FluxScheme;               /* Roe, LMRoe */
extern int    TimeIntegrationMethod;    /* Euler, Runge-Kutta 4 */
extern int    TimeIntegrationType;      /* Global Time, Local Time */
extern int    BCMethod;                 /* Characteristic, Weak Characteristic */
extern int    JacobianMethod;           /* Central, Forward, Backward, Approx, Exact */
extern int    JacobianUpdate;           
extern int    SolverOrder;              /* First, Second */
extern int    SolverNIteration;         /* Number of Main Iteration */ 
extern int    InnerNIteration;          /* Inner Iterations: Newton or Dual Time */
extern int    LinearSolverNIteration;   /* Linear Solver Iterations */
extern int    FirstOrderNIteration;
extern double PhysicalDeltaTime;        /* Physical Time Step For Unsteady Calculations */
extern double Relaxation;               
extern int    AverageType;              /* Interface Average Computation Type */

// Solution Smoother
extern int    SolutionSmooth;
extern int    SolutionSmoothNIteration;
extern double SolutionSmoothRelaxation;

// Residual Smoother Variables
extern int    ResidualSmoothMethod;
extern int    ResidualSmoothType;
extern int    ResidualSmoothScheme;
extern int    ResidualSmoothNIteration;
extern double ResidualSmoothRelaxation;

// Flux Limiters
extern int    LimiterMethod;
extern int    LimiterSmooth;
extern int    LimiterOrder;
extern int    LimiterStartSolverIteration;
extern int    LimiterEndSolverIteration;
extern double Venkat_KThreshold;

// Entropy Fix
extern int EntropyFix;

// Precondition Variables
extern int    PrecondMethod;
extern int    PrecondType;
extern int    PrecondSmooth;
extern int    PrecondVariableType;
extern double PrecondGlobalMach;

// Restart Parameters
extern int  RestartInput;
extern int  RestartOutput;
extern int  RestartIteration;
extern int  RestartCycle;
extern int  RestartInputIOType;
extern int  RestartOutputIOType;
extern char RestartInputFilename[256];
extern char RestartOutputFilename[256];

// Solver Tuning Parameters
// CFL Conditions
extern int    CFL_Ramp;
extern double CFL_MAX;
extern double CFL_MIN;
extern double CFL;

// Zero Pressure Gradient No of Iterations
extern int    ZPGIteration;

// Function Definitions
void Solver_Parameters_Init(void);
void Solver_Parameters_Read(const char *filename);
void Solver_BC_Parameters_Read(const char *filename);

void   RMS_Writer_Init(void);
void   RMS_Writer_Finalize(void);
double RMS_Writer(int Iteration);

#endif	/* _SOLVERPARAMETERS_H */

