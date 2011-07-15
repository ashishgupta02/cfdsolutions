/*******************************************************************************
 * File:        Vector2D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _VECTOR2D_H
#define	_VECTOR2D_H

#include "Point2D.h"

class Vector2D {
public:
    double vec[2];
public:
    Vector2D(double dx = 0.0, double dy = 0.0);
    Vector2D(const Point2D &from, const Point2D &to);
    ~Vector2D();

    Vector2D operator +(const Vector2D &v);
    Vector2D operator -(const Vector2D &v);
    Vector2D & operator +=(const Vector2D &v);
    Vector2D & operator -=(const Vector2D &v);
    Vector2D & operator  =(const Vector2D &v);
    Vector2D & operator +=(const double &s);
    Vector2D & operator -=(const double &s);
    Vector2D & operator *=(const double &s);
    Vector2D & operator /=(const double &s);
    Vector2D & operator  =(const double &s);
    Vector2D operator +(const double &s);
    Vector2D operator -(const double &s);
    Vector2D operator *(const double &s);
    Vector2D operator /(const double &s);

    double operator *(const Vector2D &v);
    double operator %(const Vector2D &v);
    bool operator ==(const Vector2D &v);
    bool operator !=(const Vector2D &v);
    double operator () (int i) const;
    double &operator () (int i);
    double &operator [] (int i);

    double dot(const Vector2D &v);
    double cross(const Vector2D &v);
    Vector2D normalvec(void);
    void normalize();
    double magnitude();
    void print();
};

Vector2D operator+(const double &s, const Vector2D &v);
Vector2D operator-(const double &s, const Vector2D &v);
Vector2D operator*(const double &s, const Vector2D &v);

#endif	/* _VECTOR2D_H */
