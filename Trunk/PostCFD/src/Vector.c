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
 * File		Vector.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Vector.h"
#include "Algebra.h"

/******************/
/*   2d Library	  */
/******************/
/* returns squared length of input vector */
double V2SquaredLength(Vector2 *a) {
	return((a->x * a->x)+(a->y * a->y));
}

/* returns length of input vector */
double V2Length(Vector2 *a) {
	return(sqrt(V2SquaredLength(a)));
}

/* negates the input vector and returns it */
Vector2 *V2Negate(Vector2 *v) {
	v->x = -v->x;
	v->y = -v->y;
	return(v);
}

/* normalizes the input vector and returns it */
Vector2 *V2Normalize(Vector2 *v) {
	double len = V2Length(v);
	if (len != 0.0) {
		v->x /= len;
		v->y /= len;
	}
	return(v);
}

/* scales the input vector to the new length and returns it */
Vector2 *V2Scale(Vector2 *v, double newlen) {
	double len = V2Length(v);
	if (len != 0.0) {
		v->x *= newlen/len;
		v->y *= newlen/len;
	}
	return(v);
}

/* return vector sum c = a+b */
Vector2 *V2Add(Vector2 *a, Vector2 *b, Vector2 *c) {
	c->x = a->x+b->x;
	c->y = a->y+b->y;
	return(c);
}

/* return vector difference c = a-b */
Vector2 *V2Sub(Vector2 *a, Vector2 *b, Vector2 *c) {
	c->x = a->x-b->x;
	c->y = a->y-b->y;
	return(c);
}

/* return the dot product of vectors a and b */
double V2Dot(Vector2 *a, Vector2 *b) {
	return((a->x*b->x)+(a->y*b->y));
}

/* return the dot product of vectors oa and ob */
double DotO_2D(Point2 *o, Point2 *a, Point2 *b) {
	double value;
	
	value = (a->x - o->x) * (b->x - o->x) + (a->y - o->y) * (b->y - o->y);

	return(value);
}

/* return the cross product of a and b */
double V2Cross(Vector2 *a, Vector2 *b) {
	return((a->x*b->y)-(a->y*b->x));
}

/* return the cross product of oa and ob */
double CrossO_2D(Point2 *o, Point2 *a, Point2 *b) {
	double value;
	
	value = (a->x - o->x) * (b->y - o->y) - (a->y - o->y) * (b->x - o->x);
	
	return(value);
}

/* linearly interpolate between vectors by an amount alpha */
/* and return the resulting vector. */
/* When alpha=0, result=lo.  When alpha=1, result=hi. */
Vector2 *V2Lerp(Vector2 *lo, Vector2* hi, double alpha, Vector2 *result) {
	result->x = LERP(alpha, lo->x, hi->x);
	result->y = LERP(alpha, lo->y, hi->y);
	return(result);
}

/* make a linear combination of two vectors and return the result. */
/* result = (a * ascl) + (b * bscl) */
Vector2 *V2Combine (Vector2 *a, Vector2 *b, Vector2 *result, double ascl, double bscl) {
	result->x = (ascl * a->x) + (bscl * b->x);
	result->y = (ascl * a->y) + (bscl * b->y);
	return(result);
}

/* multiply two vectors together component-wise */
Vector2 *V2Mul (Vector2 *a, Vector2 *b, Vector2 *result) {
	result->x = a->x * b->x;
	result->y = a->y * b->y;
	return(result);
}

/* return the distance between two points */
double V2DistanceBetween2Points(Point2 *a, Point2 *b) {
	double dx = a->x - b->x;
	double dy = a->y - b->y;
	return(sqrt((dx*dx)+(dy*dy)));
}

/* return the vector perpendicular to the input vector a */
Vector2 *V2MakePerpendicular(Vector2 *a, Vector2 *ap) {
	ap->x = -a->y;
	ap->y = a->x;
	return(ap);
}

/* create, initialize, and return a new vector */
Vector2 *V2New(double x, double y) {
	Vector2 *v = NEWTYPE(Vector2);
	v->x = x;
	v->y = y;
	return(v);
}

/* create, initialize, and return a duplicate vector */
Vector2 *V2Duplicate(Vector2 *a) {
	Vector2 *v = NEWTYPE(Vector2);
	v->x = a->x;
	v->y = a->y;
	return(v);
}

/* Euclidean norm of a vector a */
double V2ENorm(Vector2 *a) {
	double value;
	
	value = sqrt (a->x * a->x + a->y * a->y);
	
	return (value);
}

/* Euclidean norm of a vector oa */
double ENormO_2D(Point2 *o, Point2 *a) {
	double value;
	
	value = sqrt((a->x - o->x) * (a->x - o->x) + (a->y - o->y) * (a->y - o->y));

	return (value);
}

/******************/
/*   3d Library	  */
/******************/
/* returns squared length of input vector */
double V3SquaredLength(Vector3 *a) {
	return((a->x * a->x)+(a->y * a->y)+(a->z * a->z));
}

/* returns length of input vector */
double V3Length(Vector3 *a) {
	return(sqrt(V3SquaredLength(a)));
}

/* negates the input vector and returns it */
Vector3 *V3Negate(Vector3 *v) {
	v->x = -v->x;
	v->y = -v->y;
	v->z = -v->z;
	return(v);
}

/* normalizes the input vector and returns it */
Vector3 *V3Normalize(Vector3 *v) {
	double len = V3Length(v);
	if (len != 0.0) {
		v->x /= len;
		v->y /= len;
		v->z /= len;
	}
	return(v);
}

/* scales the input vector to the new length and returns it */
Vector3 *V3Scale(Vector3 *v, double newlen) {
	double len = V3Length(v);
	if (len != 0.0) {
		v->x *= newlen/len;
		v->y *= newlen/len;
		v->z *= newlen/len;
	}
	return(v);
}

/* return vector sum c = a+b */
Vector3 *V3Add(Vector3 *a, Vector3 *b, Vector3 *c) {
	c->x = a->x+b->x;
	c->y = a->y+b->y;
	c->z = a->z+b->z;
	return(c);
}

/* return vector difference c = a-b */
Vector3 *V3Sub(Vector3 *a, Vector3 *b, Vector3 *c) {
	c->x = a->x-b->x;
	c->y = a->y-b->y;
	c->z = a->z-b->z;
	return(c);
}

/* return the dot product of vectors a and b */
double V3Dot(Vector3 *a, Vector3 *b) {
	return((a->x*b->x)+(a->y*b->y)+(a->z*b->z));
}

/* return the dot product of vectors oa and ob */
double DotO_3D(Point3 *o, Point3 *a, Point3 *b) {
	double value;
	
	value = (a->x - o->x) * (b->x - o->x)
		+ (a->y - o->y) * (b->y - o->y) 
		+ (a->z - o->z) * (b->z - o->z);
	
	return(value);
}

/* linearly interpolate between vectors by an amount alpha */
/* and return the resulting vector. */
/* When alpha=0, result=lo.  When alpha=1, result=hi. */
Vector3 *V3Lerp(Vector3 *lo, Vector3 *hi, double alpha, Vector3 *result) {
	result->x = LERP(alpha, lo->x, hi->x);
	result->y = LERP(alpha, lo->y, hi->y);
	result->z = LERP(alpha, lo->z, hi->z);
	return(result);
}

/* make a linear combination of two vectors and return the result. */
/* result = (a * ascl) + (b * bscl) */
Vector3 *V3Combine (Vector3 *a, Vector3 *b, Vector3 *result, double ascl, double bscl) {
	result->x = (ascl * a->x) + (bscl * b->x);
	result->y = (ascl * a->y) + (bscl * b->y);
	result->z = (ascl * a->z) + (bscl * b->z);
	return(result);
}

/* multiply two vectors together component-wise and return the result */
Vector3 *V3Mul (Vector3 *a, Vector3 *b, Vector3 *result) {
	result->x = a->x * b->x;
	result->y = a->y * b->y;
	result->z = a->z * b->z;
	return(result);
}

/* return the distance between two points */
double V3DistanceBetween2Points(Point3 *a, Point3 *b) {
	double dx = a->x - b->x;
	double dy = a->y - b->y;
	double dz = a->z - b->z;
	return(sqrt((dx*dx)+(dy*dy)+(dz*dz)));
}

/* return the cross product c = a cross b */
Vector3 *V3Cross(Vector3 *a, Vector3 *b, Vector3 *c) {
	c->x = (a->y*b->z) - (a->z*b->y);
	c->y = (a->z*b->x) - (a->x*b->z);
	c->z = (a->x*b->y) - (a->y*b->x);
	return(c);
}

/* return the cross product c = oa cross ob */
Vector3 *CrossO_3D(Point3 *o, Point3 *a, Point3 *b, Vector3 *c) {
	c->x = (a->y - o->y) * (b->z - o->z) - (a->z - o->z) * (b->y - o->y);
	c->y = (a->z - o->z) * (b->x - o->x) - (a->x - o->x) * (b->z - o->z);
	c->z = (a->x - o->x) * (b->y - o->y) - (a->y - o->y) * (b->x - o->x);
	return(c);
}

/* create, initialize, and return a new vector */
Vector3 *V3New(double x, double y, double z) {
	Vector3 *v = NEWTYPE(Vector3);
	v->x = x;
	v->y = y;
	v->z = z;
	return(v);
}

/* create, initialize, and return a duplicate vector */
Vector3 *V3Duplicate(Vector3 *a) {
	Vector3 *v = NEWTYPE(Vector3);
	v->x = a->x;
	v->y = a->y;
	v->z = a->z;
	return(v);
}

/* Euclidean norm of a vector a */
double V3ENorm(Vector3 *a) {
	double value;
	
	value = sqrt (a->x * a->x + a->y * a->y + a->z * a->z);
	
	return (value);
}

/* Euclidean norm of a vector oa */
double ENormO_3D(Point3 *o, Point3 *a) {
	double value;
	
	value = sqrt((a->x - o->x) * (a->x - o->x) 
		+ (a->y - o->y) * (a->y - o->y) 
		+ (a->z - o->z) * (a->z - o->z));
 
	return (value);
}

