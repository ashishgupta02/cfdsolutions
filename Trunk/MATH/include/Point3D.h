/*******************************************************************************
 * File:        Point3D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _POINT3D_H
#define	_POINT3D_H

class Point3D {
public:
    double pos[3];
public:
    Point3D(double x = 0.0, double y = 0.0, double z = 0.0);
    ~Point3D();

    Point3D operator +(const Point3D &p);
    Point3D operator -(const Point3D &p);
    Point3D & operator +=(const Point3D &p);
    Point3D & operator -=(const Point3D &p);
    Point3D & operator  =(const Point3D &p);
    Point3D & operator +=(const double &s);
    Point3D & operator -=(const double &s);
    Point3D & operator *=(const double &s);
    Point3D & operator /=(const double &s);
    Point3D & operator  =(const double &s);
    Point3D operator +(const double &s);
    Point3D operator -(const double &s);
    Point3D operator *(const double &s);
    Point3D operator /(const double &s);

    bool operator ==(const Point3D &p);
    bool operator !=(const Point3D &p);
    double operator () (int) const;
    double & operator () (int);
    double & operator [] (int);
    void print();
};

Point3D operator+(const double &s, const Point3D &p);
Point3D operator-(const double &s, const Point3D &p);
Point3D operator*(const double &s, const Point3D &p);

#endif	/* _POINT3D_H */

