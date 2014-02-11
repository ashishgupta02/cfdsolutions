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
 * File		Algebra.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Algebra.h"
#include "Vector.h"

/******************/
/*   2d Library	  */
/******************/

/* returns the centriod of the triangle */
Vector2 *Triangle2Centroid(Vector2 *a, Vector2 *b, Vector2 *c, Vector2 *centroid) {
    centroid->x = (a->x + b->x + c->x) / 3;
    centroid->y = (a->y + b->y + c->y) / 3;
    return (centroid);
}

/* multiply a point by a projective matrix and return the transformed point */
Point2 *V2MulPointByProjMatrix(Point2 *pin, Matrix3 *m, Point2 *pout) {
    double w;
    pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) + m->element[2][0];
    pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) + m->element[2][1];
    w = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) + m->element[2][2];
    if (w != 0.0) {
        pout->x /= w;
        pout->y /= w;
    }
    return (pout);
}

/* multiply together matrices c = ab */

/* note that c must not point to either of the input matrices */
Matrix3 *V2MatMul(Matrix3 *a, Matrix3 *b, Matrix3 *c) {
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            c->element[i][j] = 0;
            for (k = 0; k < 3; k++)
                c->element[i][j] += a->element[i][k] * b->element[k][j];
        }
    }
    return (c);
}

/* transpose matrix a, return b */
Matrix3 *TransposeMatrix3(Matrix3 *a, Matrix3 *b) {
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            b->element[i][j] = a->element[j][i];
    }
    return (b);
}

/******************/
/*   3d Library	  */
/******************/

/* returns the centroid of triangle 3d */
Vector3 *Triangle3Centroid(Vector3 *a, Vector3 *b, Vector3 *c, Vector3 *centroid) {
    centroid->x = (a->x + b->x + c->x) / 3;
    centroid->y = (a->y + b->y + c->y) / 3;
    centroid->z = (a->z + b->z + c->z) / 3;
    return (centroid);
}

/* multiply a point by a matrix and return the transformed point */
Point3 *V3MulPointByMatrix(Point3 *pin, Matrix3 *m, Point3 *pout) {
    pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) + (pin->z * m->element[2][0]);
    pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) + (pin->z * m->element[2][1]);
    pout->z = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) + (pin->z * m->element[2][2]);
    return (pout);
}

/* multiply a point by a projective matrix and return the transformed point */
Point3 *V3MulPointByProjMatrix(Point3 *pin, Matrix4 *m, Point3 *pout) {
    double w;
    pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) + (pin->z * m->element[2][0]) + m->element[3][0];
    pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) + (pin->z * m->element[2][1]) + m->element[3][1];
    pout->z = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) + (pin->z * m->element[2][2]) + m->element[3][2];
    w = (pin->x * m->element[0][3]) + (pin->y * m->element[1][3]) + (pin->z * m->element[2][3]) + m->element[3][3];
    if (w != 0.0) {
        pout->x /= w;
        pout->y /= w;
        pout->z /= w;
    }
    return (pout);
}

/* multiply together matrices c = ab */

/* note that c must not point to either of the input matrices */
Matrix4 *V3MatMul(Matrix4 *a, Matrix4 *b, Matrix4 *c) {
    int i, j, k;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            c->element[i][j] = 0;
            for (k = 0; k < 4; k++)
                c->element[i][j] += a->element[i][k] * b->element[k][j];
        }
    }
    return (c);
}

/***********************/
/*   Useful Routines   */
/***********************/
/*--------------------------------------------------------------*/

/*
Purpose:

  ACOSH2 returns the inverse hyperbolic cosine of a number.

  Definition:

    Applying the inverse function

      Y = ACOSH2(X)

    implies that

      X = COSH(Y)
        = 0.5 * ( EXP(Y) + EXP(-Y) ).

    For every X greater than or equal to 1, there are two possible
    choices Y such that X = COSH(Y), differing only in sign.  It
    is usual to resolve this choice by taking the value of ACOSH(X)
    to be nonnegative.

  Method:

    One formula is:

      ACOSH2 = LOG ( X + SQRT ( X**2 - 1.0 ) )

    but this formula suffers from roundoff and overflow problems.
  Parameters:

    Input, double X, the number whose inverse hyperbolic cosine is desired.
    X should be greater than or equal to 1.

    Output, double ACOSH2, the inverse hyperbolic cosine of X.  The
    principal value (that is, the positive value of the two ) is returned.
 */
double acosh2(double x) {
    double a;

    if (x < 1.0) {
        printf("ACOSH2 - Fatal error!\n");
        printf("  Argument X must be >= 1.\n");
        printf("  The input X = %f\n", x);
        exit(1);
    }

    a = 2.0 * log(sqrt(0.5 * (x + 1.0)) + sqrt(0.5 * (x - 1.0)));

    return a;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    ASINH2 returns the inverse hyperbolic sine of a number.

  Definition:

    Y = ASINH2(X) implies that 
    X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).

  Parameters:

    Input, double X, the number whose inverse hyperbolic sine is desired.

    Output, double ASINH2, the inverse hyperbolic sine of X.
 */
double asinh2(double x) {
    double value;

    value = log(x + sqrt(x * x + 1.0));

    return value;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    ATANH2 returns the inverse hyperbolic tangent of a number.

  Definition:

    Y = ATANH2(X) implies that
    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )

  Discussion:

    Since a library function ATANH may be available on some systems,
    this routine is named ATANH2 to avoid name conflicts.

  Parameters:

    Input, double X, the number whose inverse hyperbolic tangent is desired.
    The absolute value of X should be less than or equal to 1.

    Output, double ATANH2, the inverse hyperbolic tangent of X.
 */
double atanh2(double x) {
    double value;

    if (fabs(x) >= 1.0) {
        printf("ATANH2 - Fatal error\n");
        printf("  ABS(X) must be < 1.\n");
        printf("  Your input is X = %f\n", x);
        exit(1);
    }

    value = 0.5 * log((1.0 + x) / (1.0 - x));

    return value;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    ATAN4 computes the inverse tangent of the ratio Y / X.

  Discussion:

    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
    the built in functions ATAN and ATAN2 already do.  

    However:

 * ATAN4 always returns a positive angle, between 0 and 2 PI, 
      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
      and [-PI,+PI] respectively;

 * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN 
     function by contrast always returns an angle in the first or fourth
     quadrants.

  Parameters:

    Input, double Y, X, two quantities which represent the tangent of
    an angle.  If Y is not zero, then the tangent is (Y/X).

    Output, double ATAN4, a positive angle whose tangent is (Y/X), and
    which lies in the appropriate quadrant so that the signs of its
    cosine and sine match those of X and Y.
 */
double atan4(double y, double x) {
    double abs_x, abs_y;
    double theta = 0, theta_0;

    /* Special cases: */
    if (x == 0.0) {
        if (y > 0.0) {
            theta = PI / 2.0;
        } else if (y < 0.0) {
            theta = 3.0 * PI / 2.0;
        } else {
            theta = 0.0;
        }
    } else if (y == 0.0) {
        if (x > 0.0) {
            theta = 0.0;
        } else if (x < 0.0) {
            theta = PI;
        }
    }
        /*
        We assume that ATAN2 is correct when both arguments are positive.
         */
    else {
        abs_y = fabs(y);
        abs_x = fabs(x);

        theta_0 = atan2(abs_y, abs_x);

        if (x > 0.0 && y > 0.0) {
            theta = theta_0;
        } else if (x < 0.0 && y > 0.0) {
            theta = PI - theta_0;
        } else if (x < 0.0 && y < 0.0) {
            theta = PI + theta_0;
        } else if (x > 0.0 && y < 0.0) {
            theta = 2.0 * PI - theta_0;
        }
    }

    return theta;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    COT returns the cotangent of an angle.

  Parameters:

    Input, double ANGLE, the angle, in radians.

    Output, double COT, the cotangent of the angle.
 */
double cot(double angle) {
    double value;

    value = cos(angle) / sin(angle);

    return value;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    AGUD evaluates the inverse Gudermannian function.

  Definition:

    The Gudermannian function relates the hyperbolic and trigonomentric
    functions.  For any argument X, there is a corresponding value
    GAMMA so that

      SINH(X) = TAN(GAMMA).

    This value GAMMA(X) is called the Gudermannian of X.  The inverse
    Gudermannian function is given as input a value GAMMA and computes
    the corresponding value X.

  Parameters:

    Input, double GAMMA, the value of the Gudermannian.

    Output, double AGUD, the argument of the Gudermannian.
 */
double agud(double gamma) {
    double value;

    value = log(tan(0.25 * PI + 0.5 * gamma));

    return value;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    I_GCF finds the greatest common factor of I and J.

  Parameters:

    Input, int I, J, two numbers whose greatest common factor
    is desired.

    Output, int I_GCF, the greatest common factor of I and J.

    Note that only the absolute values of I and J are
    considered, so that I_GCF is always nonnegative.

    If I or J is 0, I_GCF is returned as max ( 1, abs ( I ), abs ( J ) ).

    If I and J have no common factor, I_GCF is returned as 1.

    Otherwise, using the Euclidean algorithm, I_GCF is the
    largest common factor of I and J.
 */
int i_gcf(int i, int j) {
    int ip, iq, ir;

    i = abs(i);
    j = abs(j);

    /* Return immediately if either I or J is zero. */
    if (i == 0) {
        return j;
    } else if (j == 0) {
        return i;
    }

    /* 
    Set IP to the larger of I and J, IQ to the smaller.
    This way, we can alter IP and IQ as we go.
     */
    if (i > j) {
        ip = i;
        iq = j;
    } else {
        ip = j;
        iq = i;
    }

    /* Carry out the Euclidean algorithm. */
    for (;;) {
        ir = ip % iq;

        if (ir == 0) {
            return iq;
        }

        ip = iq;
        iq = ir;
    }
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    I_LCM computes the least common multiple of two ints.

  Definition:

    The least common multiple may be defined as

      LCM(I,J) = ABS( I * J ) / GCF(I,J)

    where GCF(I,J) is the greatest common factor of I and J.

  Parameters:

    Input, int I, J, the integers whose I_LCM is desired.

    Output, int I_LCM, the least common multiple of I and J.
    I_LCM is never negative.  I_LCM is 0 if either I or J is zero.
 */
int i_lcm(int i, int j) {
    return abs(i * (j / i_gcf(i, j)));
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix2_Det computes the determinant of a 2 by 2 matrix.

  Formula:

    The determinant of a 2 by 2 matrix is

      m11 * m22 - m12 * m21.

  Parameters:

    Input, Matrix2 *M, the matrix whose determinant is desired.

    Output, double Matrix2_Det, the determinant of the matrix.
 */
double Matrix2_Det(Matrix2 *m) {
    double det;

    det = (m->element[0][0] * m->element[1][1])
            - (m->element[0][1] * m->element[1][0]);

    return det;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix2_Inverse inverts a 2 by 2 float matrix using Cramer's rule.

  Parameters:

    Input, Matrix2 *M, the matrix to be inverted.

    Output, Matrix2 *N, the inverse of the matrix M.

    Output, double Matrix2_Inverse, the determinant of the matrix M.

    If the determinant is zero, then M is singular, and does not have an
    inverse.  In that case, N is simply set to zero, and a
    message is printed.

    If the determinant is nonzero, then its value is roughly an estimate
    of how nonsingular the matrix M is.
 */
double Matrix2_Inverse(Matrix2 *m, Matrix2 *n) {
    double det;
    int i, j;

    /* Compute the determinant of M.*/
    det = (m->element[0][0] * m->element[1][1])
            - (m->element[0][1] * m->element[1][0]);

    /* If the determinant is zero, bail out.*/
    if (det == 0.0) {
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 2; j++) {
                n->element[i][j] = 0.0;
            }
        }

        return det;
    }

    /*
    Compute the entries of the inverse matrix using an explicit formula.
     */
    n->element[0][0] = +m->element[1][1] / det;
    n->element[0][1] = -m->element[0][1] / det;
    n->element[1][0] = -m->element[1][0] / det;
    n->element[1][1] = +m->element[0][0] / det;

    return det;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix3_Det computes the determinant of a 3 by 3 matrix.

  Formula:

    The determinant of a 3 by 3 matrix is

    a11 * a22 * a33 - a11 * a23 * a32
  + a12 * a23 * a31 - a12 * a21 * a33
  + a13 * a21 * a32 - a13 * a22 * a31

  Parameters:

    Input, Matrix3 *A, the matrix whose determinant is desired.

    Output, double Matrix3_Det, the determinant of the matrix.
 */
double Matrix3_Det(Matrix3 *a) {
    double det;

    det = a->element[0][0] * (a->element[1][1] * a->element[2][2]
            - a->element[1][2] * a->element[2][1])
            + a->element[0][1] * (a->element[1][2] * a->element[2][0]
            - a->element[1][0] * a->element[2][2])
            + a->element[0][2] * (a->element[1][0] * a->element[2][1]
            - a->element[1][1] * a->element[2][0]);

    return det;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix3_Inverse inverts a 3 by 3 matrix using Cramer's rule.

  Parameters:

    Input, Matrix3 *A, the matrix to be inverted.

    Output, Matrix3 *B, the inverse of the matrix A.

    Output, double Matrix3_Inverse, the determinant of the matrix A.

    If the determinant is zero, A is singular, and does not have an
    inverse.  In that case, B is simply set to zero.

    If the determinant is nonzero, its value is an estimate
    of how nonsingular the matrix A is.
 */
double Matrix3_Inverse(Matrix3 *a, Matrix3 *b) {
    double det;
    int i, j;

    /* Compute the determinant of A.*/
    det =
            a->element[0][0] * (a->element[1][1] * a->element[2][2]
            - a->element[1][2] * a->element[2][1])
            + a->element[0][1] * (a->element[1][2] * a->element[2][0]
            - a->element[1][0] * a->element[2][2])
            + a->element[0][2] * (a->element[1][0] * a->element[2][1]
            - a->element[1][1] * a->element[2][0]);

    if (det == 0.0) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                b->element[i][j] = 0.0;
            }
        }
    } else {
        b->element[0][0] = (a->element[1][1] * a->element[2][2]
                - a->element[1][2] * a->element[2][1]) / det;
        b->element[0][1] = -(a->element[0][1] * a->element[2][2]
                - a->element[0][2] * a->element[2][1]) / det;
        b->element[0][2] = (a->element[0][1] * a->element[1][2]
                - a->element[0][2] * a->element[1][1]) / det;

        b->element[1][0] = -(a->element[1][0] * a->element[2][2]
                - a->element[1][2] * a->element[2][0]) / det;
        b->element[1][1] = (a->element[0][0] * a->element[2][2]
                - a->element[0][2] * a->element[2][0]) / det;
        b->element[1][2] = -(a->element[0][0] * a->element[1][2]
                - a->element[0][2] * a->element[1][0]) / det;

        b->element[2][0] = (a->element[1][0] * a->element[2][1]
                - a->element[1][1] * a->element[2][0]) / det;
        b->element[2][1] = -(a->element[0][0] * a->element[2][1]
                - a->element[0][1] * a->element[2][0]) / det;
        b->element[2][2] = (a->element[0][0] * a->element[1][1]
                - a->element[0][1] * a->element[1][0]) / det;
    }

    return det;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix4_Det computes the determinant of a 4 by 4 matrix.

  Parameters:

    Input, Matrix4 *A, the matrix whose determinant is desired.

    Output, double Matrix4_Det, the determinant of the matrix.
 */
double Matrix4_Det(Matrix4 *a) {
    double det;

    det =
            a->element[0][0] * (
            a->element[1][1] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            - a->element[1][2] * (a->element[2][1] * a->element[3][3] - a->element[2][3] * a->element[3][1])
            + a->element[1][3] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1]))
            - a->element[0][1] * (
            a->element[1][0] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            - a->element[1][2] * (a->element[2][0] * a->element[3][3] - a->element[2][3] * a->element[3][0])
            + a->element[1][3] * (a->element[2][0] * a->element[3][2] - a->element[2][2] * a->element[3][0]))
            + a->element[0][2] * (
            a->element[1][0] * (a->element[2][1] * a->element[3][3] - a->element[2][3] * a->element[3][1])
            - a->element[1][1] * (a->element[2][0] * a->element[3][3] - a->element[2][3] * a->element[3][0])
            + a->element[1][3] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0]))
            - a->element[0][3] * (
            a->element[1][0] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1])
            - a->element[1][1] * (a->element[2][0] * a->element[3][2] - a->element[2][2] * a->element[3][0])
            + a->element[1][2] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0]));

    return det;
}

/*--------------------------------------------------------------*/

/*
  Purpose:

    Matrix4_Inverse inverts a 4 by 4 real matrix using Cramer's rule.

  Parameters:

    Input, Matrix4 *A, the matrix to be inverted.

    Output, Matrix4 *B, the inverse of the matrix A.

    Output, double Matrix4_Inverse, the determinant of the matrix A.
 */
double Matrix4_Inverse(Matrix4 *a, Matrix4 *b) {
    int i, j;
    double det;

    /* Compute the determinant of A.*/
    det = Matrix4_Det(a);

    /* If the determinant is zero, bail out.*/
    if (det == 0.0) {
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                b->element[i][j] = 0.0;
            }
        }

        return det;
    }

    /*
    Compute the entries of the inverse matrix using an explicit formula.
     */
    b->element[0][0] =
            +(
            +a->element[1][1] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            + a->element[1][2] * (a->element[2][3] * a->element[3][1] - a->element[2][1] * a->element[3][3])
            + a->element[1][3] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1])
            ) / det;

    b->element[1][0] =
            -(
            +a->element[1][0] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            + a->element[1][2] * (a->element[2][3] * a->element[3][0] - a->element[2][0] * a->element[3][3])
            + a->element[1][3] * (a->element[2][0] * a->element[3][2] - a->element[2][2] * a->element[3][0])
            ) / det;

    b->element[2][0] =
            +(
            +a->element[1][0] * (a->element[2][1] * a->element[3][3] - a->element[2][3] * a->element[3][1])
            + a->element[1][1] * (a->element[2][3] * a->element[3][0] - a->element[2][0] * a->element[3][3])
            + a->element[1][3] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0])
            ) / det;

    b->element[3][0] =
            -(
            +a->element[1][0] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1])
            + a->element[1][1] * (a->element[2][2] * a->element[3][0] - a->element[2][0] * a->element[3][2])
            + a->element[1][2] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0])
            ) / det;

    b->element[0][1] =
            -(
            +a->element[0][1] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            + a->element[0][2] * (a->element[2][3] * a->element[3][1] - a->element[2][1] * a->element[3][3])
            + a->element[0][3] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1])
            ) / det;

    b->element[1][1] =
            +(
            +a->element[0][0] * (a->element[2][2] * a->element[3][3] - a->element[2][3] * a->element[3][2])
            + a->element[0][2] * (a->element[2][3] * a->element[3][0] - a->element[2][0] * a->element[3][3])
            + a->element[0][3] * (a->element[2][0] * a->element[3][2] - a->element[2][2] * a->element[3][0])
            ) / det;

    b->element[2][1] =
            -(
            +a->element[0][0] * (a->element[2][1] * a->element[3][3] - a->element[2][3] * a->element[3][1])
            + a->element[0][1] * (a->element[2][3] * a->element[3][0] - a->element[2][0] * a->element[3][3])
            + a->element[0][3] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0])
            ) / det;

    b->element[3][1] =
            +(
            +a->element[0][0] * (a->element[2][1] * a->element[3][2] - a->element[2][2] * a->element[3][1])
            + a->element[0][1] * (a->element[2][2] * a->element[3][0] - a->element[2][0] * a->element[3][2])
            + a->element[0][2] * (a->element[2][0] * a->element[3][1] - a->element[2][1] * a->element[3][0])
            ) / det;

    b->element[0][2] =
            +(
            +a->element[0][1] * (a->element[1][2] * a->element[3][3] - a->element[1][3] * a->element[3][2])
            + a->element[0][2] * (a->element[1][3] * a->element[3][1] - a->element[1][1] * a->element[3][3])
            + a->element[0][3] * (a->element[1][1] * a->element[3][2] - a->element[1][2] * a->element[3][1])
            ) / det;

    b->element[1][2] =
            -(
            +a->element[0][0] * (a->element[1][2] * a->element[3][3] - a->element[1][3] * a->element[3][2])
            + a->element[0][2] * (a->element[1][3] * a->element[3][0] - a->element[1][0] * a->element[3][3])
            + a->element[0][3] * (a->element[1][0] * a->element[3][2] - a->element[1][2] * a->element[3][0])
            ) / det;

    b->element[2][2] =
            +(
            +a->element[0][0] * (a->element[1][1] * a->element[3][3] - a->element[1][3] * a->element[3][1])
            + a->element[0][1] * (a->element[1][3] * a->element[3][0] - a->element[1][0] * a->element[3][3])
            + a->element[0][3] * (a->element[1][0] * a->element[3][1] - a->element[1][1] * a->element[3][0])
            ) / det;

    b->element[3][2] =
            -(
            +a->element[0][0] * (a->element[1][1] * a->element[3][2] - a->element[1][2] * a->element[3][1])
            + a->element[0][1] * (a->element[1][2] * a->element[3][0] - a->element[1][0] * a->element[3][2])
            + a->element[0][2] * (a->element[1][0] * a->element[3][1] - a->element[1][1] * a->element[3][0])
            ) / det;

    b->element[0][3] =
            -(
            +a->element[0][1] * (a->element[1][2] * a->element[2][3] - a->element[1][3] * a->element[2][2])
            + a->element[0][2] * (a->element[1][3] * a->element[2][1] - a->element[1][1] * a->element[2][3])
            + a->element[0][3] * (a->element[1][1] * a->element[2][2] - a->element[1][2] * a->element[2][1])
            ) / det;

    b->element[1][3] =
            +(
            +a->element[0][0] * (a->element[1][2] * a->element[2][3] - a->element[1][3] * a->element[2][2])
            + a->element[0][2] * (a->element[1][3] * a->element[2][0] - a->element[1][0] * a->element[2][3])
            + a->element[0][3] * (a->element[1][0] * a->element[2][2] - a->element[1][2] * a->element[2][0])
            ) / det;

    b->element[2][3] =
            -(
            +a->element[0][0] * (a->element[1][1] * a->element[2][3] - a->element[1][3] * a->element[2][1])
            + a->element[0][1] * (a->element[1][3] * a->element[2][0] - a->element[1][0] * a->element[2][3])
            + a->element[0][3] * (a->element[1][0] * a->element[2][1] - a->element[1][1] * a->element[2][0])
            ) / det;

    b->element[3][3] =
            +(
            +a->element[0][0] * (a->element[1][1] * a->element[2][2] - a->element[1][2] * a->element[2][1])
            + a->element[0][1] * (a->element[1][2] * a->element[2][0] - a->element[1][0] * a->element[2][2])
            + a->element[0][2] * (a->element[1][0] * a->element[2][1] - a->element[1][1] * a->element[2][0])
            ) / det;

    return det;
}

