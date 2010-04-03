/*
 * File:   MUtils.c
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#include "MUtils.h"

/* Matrix Operations */

/* Matrix Determinant 2x2 */
double DetMat2x2(double a, double b, double c, double d)
{
  return (a*d - b*c);
}

/* Matrix Inverse 2x2 */
void InvMat2x2(double* M, double* N)
{
  double det;
  det = DetMat2x2(M[0], M[1], M[2], M[3]);
  N[0] = M[3]/det;
  N[1] = -M[1]/det;
  N[2] = -M[2]/det;
  N[3] = M[0]/det;
}

/* Matrix Multiply 2x2 */
void MulMat2x2(double* M, double* N, double* O)
{
  O[0] = M[0]*N[0] + M[1]*N[2];
  O[1] = M[0]*N[1] + M[1]*N[3];
  O[2] = M[2]*N[0] + M[3]*N[2];
  O[3] = M[2]*N[1] + M[3]*N[3];
}

