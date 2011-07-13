/*******************************************************************************
 * File:        Point2D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _POINT2D_H
#define	_POINT2D_H

class Point2D {
public:
    double pos[2];
public:
    Point2D(double x = 0.0, double y = 0.0);
    ~Point2D();

    Point2D operator +(const Point2D &) const;
    Point2D operator -(const Point2D &) const;
    Point2D & operator +=(const Point2D &);
    Point2D & operator -=(const Point2D &);
    Point2D & operator +=(double);
    Point2D & operator -=(double);
    Point2D & operator *=(double);
    Point2D & operator /=(double);
    Point2D operator *(double) const;
    Point2D operator /(double) const;
    double & operator () (int);
    double operator () (int) const;
    double & operator [] (int);
    void print();
};

#endif	/* _POINT2D_H */

