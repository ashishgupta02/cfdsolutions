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
 * File		Algebra.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

#ifndef __Algebra_H__
#define __Algebra_H__

#include "Vector.h"

/* 2-by-2 matrix */
typedef struct Matrix2Struct {
    double element[2][2];
} Matrix2;

/* 3-by-3 matrix */
typedef struct Matrix3Struct {
    double element[3][3];
} Matrix3;

/* 4-by-4 matrix */
typedef struct Matrix4Struct {
    double element[4][4];
} Matrix4;

/********************/
/* useful constants */
/********************/
#define PI			3.141592	/* the venerable pi */
#define PITIMES2	6.283185	/* 2 * pi */
#define PIOVER2	1.570796	/* pi / 2 */
#define E			2.718282	/* the venerable e */
#define SQRT2		1.414214	/* sqrt(2) */
#define SQRT3		1.732051	/* sqrt(3) */
#define GOLDEN	1.618034	/* the golden ratio */
#define DTOR		0.017453	/* convert degrees to radians */
#define RTOD		57.29578	/* convert radians to degrees */
#define EQN_EPS	1e-9
#define DEG_TO_RAD   (PI / 180.0)
#define RAD_TO_DEG   (180.0 / PI)

/***********************/
/* one-argument macros */
/***********************/
/* radian to degrees */
#define DEGREES(x)	((x)/(PI)*180)

/* degrees to radian */
#define RAD(x)	((x)/180.0*(PI))

/* absolute value of a */
#define ABS(a)		(((a)<0) ? -(a) : (a))

/* round a to nearest int */
#define ROUND(a)	floor((a)+0.5)

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a)	(((a)<0) ? -1 : (a)>0 ? 1 : 0)

/* take binary sign of a, either -1, or 1 if >= 0 */
#define SGN(a)		(((a)<0) ? -1 : 1)

/* shout if something that should be true isn't */
#define ASSERT(x) \
if (!(x)) fprintf(stderr," Assert failed: x\n");

/* square a */
#define SQR(a)		((a)*(a))

/* Magnitude of a times sign of b */
#define SIGN(a, b)	((b) > 0.0 ? fabs(a) : -fabs(a))

/***********************/
/* two-argument macros */
/***********************/
/* find minimum of a and b */
#define MIN(a,b)	(((a)<(b))?(a):(b))

/* find maximum of a and b */
#define MAX(a,b)	(((a)>(b))?(a):(b))

/* swap int a and b */
#define SWAP(a,b)	{ a^=b; b^=a; a^=b; }

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) */
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

/* clamp the input to the specified range */
#define CLAMP(v,l,h)		((v)<(l) ? (l) : (v) > (h) ? (h) : v)

/****************************/
/* memory allocation macros */
/****************************/
/* create a new instance of a structure (see Gem by Hultquist) */
#define NEWSTRUCT(x)	(struct x *)(malloc((unsigned)sizeof(struct x)))

/* create a new instance of a type */
#define NEWTYPE(x)		(x *)(malloc((unsigned)sizeof(x)))

/************/
/* booleans */
/************/
#define TRUE		1
#define FALSE		0
#define ON		1
#define OFF		0
#define ERROR		0
#define SUCCESS		1
#define FAIL		0
typedef int boolean; /* boolean data type */
typedef boolean flag; /* flag data type */

/* Functions Declearations */
/******************/
/*   2d Library	  */
/******************/
Vector2 *Triangle2Centroid(Vector2 *, Vector2 *, Vector2 *, Vector2 *);
Point2 *V2MulPointByProjMatrix(Point2 *, Matrix3 *, Point2 *);
Matrix3 *V2MatMul(Matrix3 *, Matrix3 *, Matrix3 *);
Matrix3 *TransposeMatrix3(Matrix3 *, Matrix3 *);

/******************/
/*   3d Library	  */
/******************/
Vector3 *Triangle3Centroid(Vector3 *, Vector3 *, Vector3 *, Vector3 *);
Point3 *V3MulPointByMatrix(Point3 *, Matrix3 *, Point3 *);
Point3 *V3MulPointByProjMatrix(Point3 *, Matrix4 *, Point3 *);
Matrix4 *V3MatMul(Matrix4 *, Matrix4 *, Matrix4 *);

/***********************/
/*   Useful Routines   */
/***********************/
double acosh2(double);
double asinh2(double);
double atanh2(double);
double atan4(double, double);
double cot(double);
double agud(double);
int i_gcf(int, int);
int i_lcm(int, int);
double Matrix2_Det(Matrix2 *);
double Matrix2_Inverse(Matrix2 *, Matrix2 *);
double Matrix3_Det(Matrix3 *);
double Matrix3_Inverse(Matrix3 *, Matrix3 *);
double Matrix4_Det(Matrix4 *);
double Matrix4_Inverse(Matrix4 *, Matrix4 *);

#endif
