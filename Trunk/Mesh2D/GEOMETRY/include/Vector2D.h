/* 
 * File:   Vector2D.h
 * Author: Ashish Gupta
 *
 * Created on April 20, 2010, 5:52 PM
 */

#include "Point2D.h"

#ifndef _VECTOR2D_H
#define	_VECTOR2D_H

class Vector2D {
private:
    double vec[2];
public:
    Vector2D(double dx = 0.0, double dy = 0.0);
    Vector2D(const Point2D &from, const Point2D &to);
    ~Vector2D();

    Vector2D operator +(const Vector2D &v);
    Vector2D operator -(const Vector2D &v);
    Vector2D & operator +=(const Vector2D &v);
    Vector2D & operator -=(const Vector2D &v);
    Vector2D & operator +=(const double &s);
    Vector2D & operator -=(const double &s);
    Vector2D & operator *=(const double &s);
    Vector2D & operator /=(const double &s);

    double operator *(const Vector2D &v);
    double operator %(const Vector2D &v);
    Vector2D operator *(const double &s);
    Vector2D operator /(const double &s);
    double operator () (int i) const;
    double &operator () (int i);
    double &operator [] (int i);
    double magnitude();
    void normalize();
    void print();
};

#endif	/* _VECTOR2D_H */
