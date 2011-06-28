/*******************************************************************************
 * File:        Vector3D.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef _VECTOR3D_H
#define	_VECTOR3D_H

#include "Point3D.h"

class Vector3D {
public:
    double vec[3];
public:
    Vector3D(double dx = 0.0, double dy = 0.0, double dz = 0.0);
    Vector3D(const Point3D &from, const Point3D &to);
    ~Vector3D();

    Vector3D operator +(const Vector3D &v);
    Vector3D operator -(const Vector3D &v);
    Vector3D & operator +=(const Vector3D &v);
    Vector3D & operator -=(const Vector3D &v);
    Vector3D & operator  =(const Vector3D &v);
    Vector3D & operator +=(const double &s);
    Vector3D & operator -=(const double &s);
    Vector3D & operator *=(const double &s);
    Vector3D & operator /=(const double &s);
    Vector3D & operator  =(const double &s);
    Vector3D operator +(const double &s);
    Vector3D operator -(const double &s);
    Vector3D operator *(const double &s);
    Vector3D operator /(const double &s);

    double operator *(const Vector3D &v);
    Vector3D operator %(const Vector3D &v);
    bool operator ==(const Vector3D &v);
    bool operator !=(const Vector3D &v);
    double operator () (int i) const;
    double & operator () (int i);
    double & operator[] (int i);

    double dot(const Vector3D &v);
    Vector3D cross(const Vector3D &v);
    Vector3D normalvec(void);
    void normalize();
    double magnitude();
    void print();
};

Vector3D operator+(const double &s, const Vector3D &v);
Vector3D operator-(const double &s, const Vector3D &v);
Vector3D operator*(const double &s, const Vector3D &v);

#endif	/* _VECTOR3D_H */
