/*******************************************************************************
 * File:        Vector3D.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <math.h>

#include "Vector3D.h"

// *****************************************************************************
// *****************************************************************************
Vector3D::Vector3D(double x, double y, double z) {
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
}

// *****************************************************************************
// *****************************************************************************
Vector3D::Vector3D(const Point3D &from, const Point3D &to) {
    vec[0] = to.pos[0] - from.pos[0];
    vec[1] = to.pos[1] - from.pos[1];
    vec[2] = to.pos[2] - from.pos[2];
}

// *****************************************************************************
// *****************************************************************************
Vector3D::~Vector3D() {

}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator +(const Vector3D &v) {
    return Vector3D(vec[0] + v.vec[0], vec[1] + v.vec[1], vec[2] + v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator -(const Vector3D &v) {
    return Vector3D(vec[0] - v.vec[0], vec[1] - v.vec[1], vec[2] - v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator +=(const Vector3D &v) {
    vec[0] += v.vec[0];
    vec[1] += v.vec[1];
    vec[2] += v.vec[2];
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator -=(const Vector3D &v) {
    vec[0] -= v.vec[0];
    vec[1] -= v.vec[1];
    vec[2] -= v.vec[2];
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator =(const Vector3D &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
    vec[2] = v.vec[2];
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator +=(const double &s) {
    vec[0] += s;
    vec[1] += s;
    vec[2] += s;
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator -=(const double &s) {
    vec[0] -= s;
    vec[1] -= s;
    vec[2] -= s;
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator *=(const double &s) {
    vec[0] *= s;
    vec[1] *= s;
    vec[2] *= s;
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator /=(const double &s) {
    vec[0] /= s;
    vec[1] /= s;
    vec[2] /= s;
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D &Vector3D::operator =(const double &s) {
    vec[0] = s;
    vec[1] = s;
    vec[2] = s;
    return *this;
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator +(const double &s) {
    return Vector3D(vec[0] + s, vec[1] + s, vec[2] + s);
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator -(const double &s) {
    return Vector3D(vec[0] - s, vec[1] - s, vec[2] - s);
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator *(const double &s) {
    return Vector3D(vec[0] * s, vec[1] * s, vec[2] * s);
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::operator /(const double &s) {
    return Vector3D(vec[0] / s, vec[1] / s, vec[2] / s);
}

// *****************************************************************************
// Dot Product
// *****************************************************************************
double Vector3D::operator *(const Vector3D &v) {
    return double(vec[0] * v.vec[0] + vec[1] * v.vec[1] + vec[2] * v.vec[2]);
}

// *****************************************************************************
// Cross Product
// *****************************************************************************
Vector3D Vector3D::operator %(const Vector3D &v) {
    Vector3D temp;
    temp.vec[0] = +vec[1] * v.vec[2] - vec[2] * v.vec[1];
    temp.vec[1] = -vec[0] * v.vec[2] + vec[2] * v.vec[0];
    temp.vec[2] = +vec[0] * v.vec[1] - vec[1] * v.vec[0];
    return temp;
}

// *****************************************************************************
// *****************************************************************************
bool Vector3D::operator ==(const Vector3D &v) {
    return (vec[0] == v.vec[0] && vec[1] == v.vec[1] && vec[2] == v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
bool Vector3D::operator !=(const Vector3D &v) {
    return (vec[0] != v.vec[0] || vec[1] != v.vec[1] || vec[2] != v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
double Vector3D::operator () (int i) const {
    return vec[i];
}

// *****************************************************************************
// *****************************************************************************
double &Vector3D::operator () (int i) {
    return vec[i];
}

// *****************************************************************************
// *****************************************************************************
double &Vector3D::operator[] (int i) {
    return vec[i];
}

// *****************************************************************************
// *****************************************************************************
double Vector3D::dot(const Vector3D &v) {
    return (vec[0] * v.vec[0] + vec[1] * v.vec[1] + vec[2] * v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::cross(const Vector3D &v) {
    Vector3D temp;
    temp.vec[0] = +vec[1] * v.vec[2] - vec[2] * v.vec[1];
    temp.vec[1] = -vec[0] * v.vec[2] + vec[2] * v.vec[0];
    temp.vec[2] = +vec[0] * v.vec[1] - vec[1] * v.vec[0];
    return temp;
}

// *****************************************************************************
// *****************************************************************************
Vector3D Vector3D::normalvec(void) {
    return (*this) /= magnitude();
}

// *****************************************************************************
// *****************************************************************************
void Vector3D::normalize() {
    double mag = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (mag > 1.0e-20) {
        vec[0] = vec[0] / mag;
        vec[1] = vec[1] / mag;
        vec[2] = vec[2] / mag;
    }
}

// *****************************************************************************
// *****************************************************************************
double Vector3D::magnitude() {
    return double(sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

// *****************************************************************************
// *****************************************************************************
void Vector3D::print() {
    printf("Vector3D (x,y,z)= (%g, %g, %g)\n", vec[0], vec[1], vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D operator +(const double &s, const Vector3D &v) {
    return Vector3D(s + v.vec[0], s + v.vec[1], s + v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D operator -(const double &s, const Vector3D &v) {
    return Vector3D(s - v.vec[0], s - v.vec[1], s - v.vec[2]);
}

// *****************************************************************************
// *****************************************************************************
Vector3D operator *(const double &s, const Vector3D &v) {
    return Vector3D(s * v.vec[0], s * v.vec[1], s * v.vec[2]);
}

