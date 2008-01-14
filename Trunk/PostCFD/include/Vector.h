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
 * File		Vector.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

#ifndef __Vector_H__
#define __Vector_H__

/*********************/
/* 2d geometry types */
/*********************/
/* 2D point */
typedef struct Point2Struct {
	double x, y;
} Point2;

typedef Point2 Vector2;

/* 2D integer point */
typedef struct IntPoint2Struct {
	int x, y;
} IntPoint2;

/* 2D box */
typedef struct Box2dStruct {
	Point2 min, max;
} Box2;


/* Functions Declearations */
/******************/
/*   2d Library	  */
/******************/
double V2SquaredLength(Vector2 *);
double V2Length(Vector2 *);
double V2Dot(Vector2 *, Vector2 *);
double V2Cross(Vector2 *, Vector2 *);
double DotO_2D(Point2 *, Point2 *, Point2 *);
double CrossO_2D(Point2 *, Point2 *, Point2 *); 
double V2DistanceBetween2Points(Point2 *, Point2 *);
double V2ENorm(Vector2 *);
double ENormO_2D(Point2 *, Point2 *);
Vector2 *V2Negate(Vector2 *);
Vector2 *V2Normalize(Vector2 *);
Vector2 *V2Scale(Vector2 *, double);
Vector2 *V2Add(Vector2 *, Vector2 *, Vector2 *);
Vector2 *V2Sub(Vector2 *, Vector2 *, Vector2 *);
Vector2 *V2Lerp(Vector2 *, Vector2 *, double, Vector2 *);
Vector2 *V2Combine(Vector2 *, Vector2 *, Vector2 *, double, double);
Vector2 *V2Mul(Vector2 *, Vector2 *, Vector2 *);
Vector2 *V2MakePerpendicular(Vector2 *, Vector2 *);
Vector2 *V2New(double, double);
Vector2 *V2Duplicate(Vector2 *);

/*********************/
/* 3d geometry types */
/*********************/
/* 3D point */
typedef struct Point3Struct {
	double x, y, z;
} Point3;

typedef Point3 Vector3;

/* 3D integer point */
typedef struct IntPoint3Struct {
	int x, y, z;
} IntPoint3;

/* 3D box */
typedef struct Box3dStruct {
	Point3 min, max;
} Box3;


/* Functions Declearations */
/******************/
/*   3d Library	  */
/******************/
double V3SquaredLength(Vector3 *);
double V3Length(Vector3 *);
double V3Dot(Vector3 *, Vector3 *);
double DotO_3D(Point3 *, Point3 *, Point3 *);
double V3DistanceBetween2Points(Point3 *, Point3 *);
double V3ENorm(Vector3 *);
double ENormO_3D(Point3 *, Point3 *);
Vector3 *V3Negate(Vector3 *);
Vector3 *V3Normalize(Vector3 *);
Vector3 *V3Scale(Vector3 *, double);
Vector3 *V3Add(Vector3 *, Vector3 *, Vector3 *);
Vector3 *V3Sub(Vector3 *, Vector3 *, Vector3 *);
Vector3 *V3Lerp(Vector3 *, Vector3 *, double, Vector3 *);
Vector3 *V3Combine(Vector3 *, Vector3 *, Vector3 *, double, double);
Vector3 *V3Mul(Vector3 *, Vector3 *, Vector3 *);
Vector3 *V3Cross(Vector3 *, Vector3 *, Vector3 *);
Vector3 *CrossO_3D(Point3 *, Point3 *, Point3 *, Vector3 *);
Vector3 *V3New(double, double, double);
Vector3 *V3Duplicate(Vector3 *);

#endif

