/*******************************************************************************
 * File:        Complex.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef COMPLEX_H
#define	COMPLEX_H

#include <stdio.h>
#include <iostream>
#include <math.h>

template <class BaseType>
class Complex {
public:
    // Data
    BaseType a[2];

    // Functions
    //--------------------------------------------------------------------------
    void assign(BaseType r, BaseType i) {
        a[0] = r;
        a[1] = i;
    };

    // Constructors
    //--------------------------------------------------------------------------
    Complex() {
    };
    //--------------------------------------------------------------------------
    Complex(BaseType r) {
        assign(r, 0.0);
    };
    //--------------------------------------------------------------------------
    Complex(BaseType r, BaseType i) {
        assign(r, i);
    };
    
    //--------------------------------------------------------------------------
    const BaseType& real(void) const {
        return a[0];
    };
    //--------------------------------------------------------------------------
    const BaseType& imag(void) const {
        return a[1];
    };

    // Unary Operations
    //--------------------------------------------------------------------------
    Complex operator-() const {
        Complex lval(-a[0], -a[1]);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex operator+() const {
        return *this;
    }

    // Assignment Operations
    //--------------------------------------------------------------------------
    Complex& operator=(const BaseType d) {
        a[0] = d;
        a[1] = 0.0;
        return *this;
    }

    // Addition Operations
    //--------------------------------------------------------------------------
    Complex operator+(const Complex &rval) const {
        return Complex(a[0] + rval.a[0], a[1] + rval.a[1]);
    }
    //--------------------------------------------------------------------------
    Complex operator+(const BaseType num) const {
        return Complex(a[0] + num, a[1]);
    }
    //--------------------------------------------------------------------------
    friend Complex operator+(const BaseType d, const Complex<BaseType>&rval) {
        const BaseType r = rval.a[0] + d;
        const BaseType i = rval.a[1];
        return Complex<BaseType > (r, i);
    };
    //--------------------------------------------------------------------------
    Complex& operator+=(const Complex &rval) {
        a[0] += rval.a[0];
        a[1] += rval.a[1];
        return *this;
    }
    //--------------------------------------------------------------------------
    Complex& operator+=(const BaseType num) {
        a[0] += num;
        return *this;
    }

    // Subtraction Operations
    //--------------------------------------------------------------------------
    Complex operator-(const Complex &rval) const {
        return Complex<BaseType > (a[0] - rval.a[0], a[1] - rval.a[1]);
    }
    //--------------------------------------------------------------------------
    Complex operator-(const BaseType num) const {
        return Complex<BaseType > (a[0] - num, a[1]);
    }
    //--------------------------------------------------------------------------
    friend Complex operator-(const BaseType d, const Complex<BaseType>&rval) {
        return Complex<BaseType > (d - rval.a[0], -rval.a[1]);
    };
    //--------------------------------------------------------------------------
    Complex& operator-=(const Complex &rval) {
        a[0] -= rval.a[0];
        a[1] -= rval.a[1];
        return *this;
    }
    //--------------------------------------------------------------------------
    Complex& operator-=(const BaseType r) {
        a[0] -= r;
        return *this;
    }

    // Multiplication Operations
    //--------------------------------------------------------------------------
    Complex operator*(const Complex &rval) const {
        Complex lval(
                a[0] * rval.a[0] - a[1] * rval.a[1],
                a[0] * rval.a[1] + a[1] * rval.a[0]
                );
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex operator*(const BaseType r) const {
        Complex lval(
                a[0] * r,
                a[1] * r
                );
        return lval;
    }
    //--------------------------------------------------------------------------
    friend Complex operator*(const BaseType d, const Complex<BaseType>&rval) {
        const BaseType r = rval.a[0] * d;
        const BaseType i = rval.a[1] * d;
        Complex<BaseType> lval(r, i);
        return lval;
    };
    //--------------------------------------------------------------------------
    Complex& operator*=(const Complex &rval) {
        const BaseType r = a[0] * rval.a[0] - a[1] * rval.a[1];
        const BaseType i = a[0] * rval.a[1] + a[1] * rval.a[0];
        a[0] = r;
        a[1] = i;
        return *this;
    }
    //--------------------------------------------------------------------------
    Complex& operator*=(const BaseType r) {
        a[0] *= r;
        a[1] *= r;
        return *this;
    }
    
    // Division Operations
    //--------------------------------------------------------------------------
    Complex operator/(const Complex &rval) const {
        const BaseType denom = 1.0 / (rval.a[0] * rval.a[0] + rval.a[1] * rval.a[1]);
        Complex lval(
                (a[0] * rval.a[0] + a[1] * rval.a[1]) * denom,
                (a[1] * rval.a[0] - a[0] * rval.a[1]) * denom
                );
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex operator/(const BaseType r) const {
        const BaseType denom = 1.0 / r;
        Complex lval(a[0] * denom, a[1] * denom);
        return lval;
    }
    //--------------------------------------------------------------------------
    friend Complex operator/(const BaseType d, const Complex<BaseType>&rval) {
        const BaseType denom = 1.0 / (rval.a[0] * rval.a[0] + rval.a[1] * rval.a[1]);
        const BaseType dd = d*denom;
        Complex<BaseType> lval(dd * rval.a[0], -dd * rval.a[1]);
        return lval;
    };
    //--------------------------------------------------------------------------
    Complex& operator/=(const Complex &rval) {
        const BaseType denom = 1.0 / (rval.a[0] * rval.a[0] + rval.a[1] * rval.a[1]);
        const BaseType r = (a[0] * rval.a[0] + a[1] * rval.a[1]) * denom;
        const BaseType i = (a[1] * rval.a[0] - a[0] * rval.a[1]) * denom;
        a[0] = r;
        a[1] = i;
        return *this;
    }
    //--------------------------------------------------------------------------
    Complex& operator/=(const BaseType r) {
        const BaseType denom = 1.0 / r;
        a[0] *= denom;
        a[1] *= denom;
        return *this;
    }

    // Logical Operations
    // In these logical ops, only the real part of the complex number is used
    // in the comparison
    //--------------------------------------------------------------------------
    int operator>(const Complex &rval) const {
        return (a[0] > rval.a[0]);
    }
    int operator>(const BaseType d) const {
        return (a[0] > d);
    }
    friend int operator>(const BaseType d, const Complex<BaseType>&rval) {
        return (d > rval.a[0]);
    }
    //--------------------------------------------------------------------------
    int operator>=(const Complex &rval) const {
        return (a[0] >= rval.a[0]);
    }
    int operator>=(const BaseType d) const {
        return (a[0] >= d);
    }
    friend int operator>=(const BaseType d, const Complex<BaseType>&rval) {
        return (d >= rval.a[0]);
    }
    //--------------------------------------------------------------------------
    int operator<(const Complex &rval) const {
        return (a[0] < rval.a[0]);
    }
    int operator<(const BaseType d) const {
        return (a[0] < d);
    }
    friend int operator<(const BaseType d, const Complex<BaseType>&rval) {
        return (d < rval.a[0]);
    }
    //--------------------------------------------------------------------------
    int operator<=(const Complex &rval) const {
        return (a[0] <= rval.a[0]);
    }
    int operator<=(const BaseType d) const {
        return (a[0] <= d);
    }
    friend int operator<=(const BaseType d, const Complex<BaseType>&rval) {
        return (d <= rval.a[0]);
    }
    //--------------------------------------------------------------------------
    int operator==(const Complex &rval) const {
        return (a[0] == rval.a[0] && a[1] == rval.a[1]);
    }
    int operator==(const BaseType d) const {
        return (a[0] == d && a[1] == 0.0);
    }
    friend int operator==(const BaseType d, const Complex<BaseType>&rval) {
        return (d == rval.a[0] && 0.0 == rval.a[1]);
    }
    //--------------------------------------------------------------------------
    int operator!=(const Complex &rval) const {
        return (a[0] != rval.a[0] || a[1] != rval.a[1]);
    }
    int operator!=(const BaseType d) const {
        return (a[0] != d || a[1] != 0.0);
    }
    friend int operator!=(const BaseType d, const Complex<BaseType>&rval) {
        return (d != rval.a[0] || 0.0 != rval.a[1]);
    }
    
    // Math Library Functions
    //--------------------------------------------------------------------------
    BaseType fabs(void) const {
        return ::sqrt(a[0] * a[0] + a[1] * a[1]);
    }
    //--------------------------------------------------------------------------
    BaseType norm(void) const {
        // is this right? It looks like in the code we want 
        // the length squared for the norm function.
        return a[0] * a[0] + a[1] * a[1];
    }
    //--------------------------------------------------------------------------
    Complex sqrt(void) const {
        const BaseType r = ::sqrt(a[0] * a[0] + a[1] * a[1]);
        const BaseType term1 = ::sqrt(r + a[0]);
        const BaseType term2 = ::sqrt(r - a[0]);
        Complex lval(M_SQRT1_2*term1,
                a[1] > 0.0 ? M_SQRT1_2 * term2 : -M_SQRT1_2 * term2
                );
        return lval;
    }
    //--------------------------------------------------------------------------
    // Exponentiation of a complex number is reasonably well
    // defined, however, the periodicity of the trigonometric
    // functions make values repeat for different imaginary
    // parts, i.e.,
    // cos(x) = cos(x + 2*pi) = cos(x - 2*pi)....
    // sin(x) = sin(x + 2*pi) = sin(x - 2*pi)....
    // What this means is the if the real part is unchanged
    // and the imaginary part is offset by an integer multiple
    // of 2*pi, we get the same return value
    Complex exp(void) const {
        const BaseType r = ::exp(a[0])*::cos(a[1]);
        const BaseType i = ::exp(a[0])*::sin(a[1]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    // log has issues because of the inherent periodicity of
    // trigonometric functions, i.e.,
    // atan2(x) = atan2(x + 2*pi) = atan2(x + 4*pi) ....
    // therefore the log is not unique
    Complex log(void) const {
        const BaseType r = 0.5 * ::log(a[0] * a[0] + a[1] * a[1]);
        const BaseType i = ::atan2(a[1], a[0]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex log10(void) const {
        return this->log() / ::log(10.0);
    }
    //--------------------------------------------------------------------------
    Complex pow(BaseType n) const {
        // z^n = (r*exp(i*theta))^n = r^n*exp(i*theta*n)
        const BaseType r = ::pow(a[0] * a[0] + a[1] * a[1], n / 2.0); // r^n
        const BaseType theta = ::atan2(a[1], a[0]) * n; // theta*n
        Complex lval(r*::cos(theta), r*::sin(theta));
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex pow(const Complex& rval) const {
        // x = a + i*b
        // y = c + i*d
        // x^y = ρ^c*exp(-d*θ)*[cos(c*θ + d*ln ρ) + i*sin(c*θ + d*ln ρ)]
        // ρ = sqrt(a^2 + b^2)
        // θ = atan2(b,a)
        const BaseType rho = ::sqrt(a[0] * a[0] + a[1] * a[1]); // ρ
        const BaseType theta = atan2(a[1], a[0]); // θ
        const BaseType c = rval.a[0];
        const BaseType d = rval.a[1];
        const BaseType temp1 = c * theta + d*::log(rho); // c*θ + d*ln ρ
        const BaseType temp2 = ::pow(rho, c)*::exp(-d * theta); // ρ^c*exp(-d*θ)
        Complex lval(temp2*::cos(temp1), temp2*::sin(temp1));
        return lval;
    }
    //--------------------------------------------------------------------------
    // Trigonometric functions
    Complex sin(void) const {
        const BaseType r = ::sin(a[0])*::cosh(a[1]);
        const BaseType i = ::cos(a[0])*::sinh(a[1]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex cos(void) const {
        const BaseType r = ::cos(a[0])*::cosh(a[1]);
        const BaseType i = -::sin(a[0])*::sinh(a[1]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex tan(void) const {
        const BaseType sinhy = ::sinh(a[1]);
        const BaseType cosx = ::cos(a[0]);
        const BaseType denom = 1.0 / (cosx * cosx + sinhy * sinhy);
        const BaseType r = ::sin(a[0]) * cosx*denom;
        const BaseType i = sinhy*::cosh(a[1]) * denom;
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex asin(void) const {
        // sin^-1(z) = -i*log(i*z + sqrt(1 - z^2))
        Complex temp1(-a[1], a[0]); // i*z
        Complex temp2 = 1.0 - (*this)*(*this); // 1 - z^2
        Complex temp3 = (temp1 + temp2.sqrt()).log(); // log(i*z + sqrt(1 - z^2))
        Complex lval(temp3.a[1], -temp3.a[0]); // -i*log(i*z + sqrt(1 - z^2))
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex acos(void) const {
        // cos^-1(z) = pi/2 + i*log(i*z + sqrt(1 - z^2))
        Complex temp1(-a[1], a[0]); // i*z
        Complex temp2 = 1.0 - (*this)*(*this); // 1 - z^2
        Complex temp3 = (temp1 + temp2.sqrt()).log(); // log(i*z + sqrt(1 - z^2))
        Complex lval(M_PI_2 - temp3.a[1], temp3.a[0]); // pi/2 + i*log(i*z + sqrt(1 - z^2))
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex atan(void) const {
        // tan^-1(z) = i*(log(1 - i*z) - log(1 + i*z))/2 = -i*tanh^-1(i*z)
        Complex temp1(-a[1], a[0]); // i*z
        Complex temp2 = temp1.atanh();
        Complex lval(temp2.a[1], -temp2.a[0]);
        return lval;
    }
    //--------------------------------------------------------------------------
    // Hyperbolic trigonometric functions
    Complex sinh(void) const {
        const BaseType r = ::sinh(a[0])*::cos(a[1]);
        const BaseType i = ::cosh(a[0])*::sin(a[1]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex cosh(void) const {
        const BaseType r = ::cosh(a[0])*::cos(a[1]);
        const BaseType i = ::sinh(a[0])*::sin(a[1]);
        Complex lval(r, i);
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex tanh(void) const {
        Complex lval(this->sinh() / this->cosh());
        return lval;
    }
    //--------------------------------------------------------------------------
    Complex asinh(void) const {
        // sinh^-1(z) = log[z + sqrt(z^2 + 1)] = i*sin^-1(-i*z)
        const Complex temp(a[1], -a[0]); // -i*z
        const Complex temp2 = temp.asin(); // sin^-1(-i*z)
        Complex lval(-temp2.a[1], temp2.a[0]); // i*sin^-1(-i*z)
        return lval;
    }
    //--------------------------------------------------------------------------
    // cosh^-1(z) = log(z + sqrt(z - 1)*sqrt(z + 1)) = +- i*cos^-1(z)
    // This has a +- issue, so implementation currently deferred
    //--------------------------------------------------------------------------
    Complex atanh(void) const {
        // tanh^-1(z) = 0.5*[log(1 + z) - log(1 - z)]
        const Complex temp1(a[0] + 1.0, a[1]);
        const Complex temp2(1.0 - a[0], -a[1]);
        Complex lval = 0.5 * (temp1.log() - temp2.log());
        return lval;
    }
    //--------------------------------------------------------------------------
    void print() {
        printf("%.15E %.15E\n", a[0], a[1]);
    }
    //--------------------------------------------------------------------------
    void printnl() {
        printf("%.2f + %.2f * I", a[0], a[1]);
    }
};

//------------------------------------------------------------------------------
// helper functions
template <class BaseType> inline
Complex<BaseType> sqrt(const Complex<BaseType> &a) {
    return a.sqrt();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> exp(const Complex<BaseType> &a) {
    return a.exp();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> log(const Complex<BaseType> &a) {
    return a.log();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> log10(const Complex<BaseType> &a) {
    return a.log10();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> sin(const Complex<BaseType> &a) {
    return a.sin();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> cos(const Complex<BaseType> &a) {
    return a.cos();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> tan(const Complex<BaseType> &a) {
    return a.tan();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> asin(const Complex<BaseType> &a) {
    return a.asin();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> acos(const Complex<BaseType> &a) {
    return a.acos();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> atan(const Complex<BaseType> &a) {
    return a.atan();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> sinh(const Complex<BaseType> &a) {
    return a.sinh();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> cosh(const Complex<BaseType> &a) {
    return a.cosh();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> tanh(const Complex<BaseType> &a) {
    return a.tanh();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> asinh(const Complex<BaseType> &a) {
    return a.asinh();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> atanh(const Complex<BaseType> &a) {
    return a.atanh();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> pow(const Complex<BaseType> &a, const double p) {
    return a.pow(p);
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> pow(const double a, const Complex<BaseType> &p) {
    const BaseType c = p.a[0];
    const BaseType d = p.a[1];
    const BaseType temp1 = d*::log(a);
    const BaseType temp2 = ::pow(a, c);
    Complex<BaseType> lval(temp2*::cos(temp1), temp2*::sin(temp1));
    return lval;
}
//------------------------------------------------------------------------------
template <class BaseType> inline
Complex<BaseType> pow(const Complex<BaseType> &a, const Complex<BaseType> &p) {
    return a.pow(p);
}
//------------------------------------------------------------------------------
template <class BaseType> inline
BaseType real(const Complex<BaseType> &a) {
    return a.real();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
BaseType imag(const Complex<BaseType> &a) {
    return a.imag();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
BaseType fabs(const Complex<BaseType> &a) {
    return a.fabs();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
BaseType abs(const Complex<BaseType> &a) {
    return a.fabs();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
BaseType norm(const Complex<BaseType> &a) {
    return a.norm();
}
//------------------------------------------------------------------------------
template <class BaseType> inline
std::ostream &operator <<(std::ostream &out_file, const Complex<BaseType> &a) {
    out_file << '(' << a.real() << ',' << a.imag() << ')';
    return out_file;
}
//------------------------------------------------------------------------------
#endif	/* COMPLEX_H */
