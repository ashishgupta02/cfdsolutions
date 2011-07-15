/*******************************************************************************
 * File:        Point2D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _POINT2D_H
#define	_POINT2D_H

class Point2D {
public:
    double pos[2];
public:
    Point2D(double x = 0.0, double y = 0.0);
    ~Point2D();

    Point2D operator +(const Point2D &p);
    Point2D operator -(const Point2D &p);
    Point2D & operator +=(const Point2D &p);
    Point2D & operator -=(const Point2D &p);
    Point2D & operator  =(const Point2D &p);
    Point2D & operator +=(const double &s);
    Point2D & operator -=(const double &s);
    Point2D & operator *=(const double &s);
    Point2D & operator /=(const double &s);
    Point2D & operator  =(const double &s);
    Point2D operator +(const double &s);
    Point2D operator -(const double &s);
    Point2D operator *(const double &s);
    Point2D operator /(const double &s);

    bool operator ==(const Point2D &p);
    bool operator !=(const Point2D &p);
    double operator () (int) const;
    double & operator () (int);
    double & operator [] (int);
    void print();
};

Point2D operator+(const double &s, const Point2D &p);
Point2D operator-(const double &s, const Point2D &p);
Point2D operator*(const double &s, const Point2D &p);

#endif	/* _POINT2D_H */

