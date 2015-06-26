/******************************************
 * C++ Quaternions
 * Version: 1.0.9
 * Author:  Douglas Wilhelm Harder
 * Date:    2008/03/03
 *
 * Copyright (c) 2006-2008 by Douglas Wilhelm Harder.
 * All rights reserved.
 ******************************************/

#include "Complex.h"
#include "Quaternion.h"
#include "Support.h"
#include <iostream>
#include <cmath>
#include <string>

/******************************************
 * Constructors
 ******************************************/

template <typename T> Quaternion<T>::Quaternion( T re, T im, T jm, T km ):r(re), i(im), j(jm), k(km) {
}

template <typename T> Quaternion<T>::Quaternion( T real ):r(real), i(0.0), j(0.0), k(0.0) {
}


/******************************************
 * Assignment Operator
 ******************************************/

template <typename T> const Quaternion<T> & Quaternion<T>::operator = ( const T & x ) {
    r = x;
    i = 0.0;
    j = 0.0;
    k = 0.0;

    return *this;
}

/******************************************
 * Mutating Arithmetic Operators
 ******************************************/

template <typename T> Quaternion<T> & Quaternion<T>::operator += ( const Quaternion<T> & q ) {
    r += q.r;
    i += q.i;
    j += q.j;
    k += q.k;

    return *this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator -= ( const Quaternion<T> & q ) {
    r -= q.r;
    i -= q.i;
    j -= q.j;
    k -= q.k;

    return *this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator *= ( const Quaternion<T> & q ) {
    T RE = r, I = i, J = j;

    r = RE*q.r - I*q.i - J*q.j - k*q.k;
    i = RE*q.i + I*q.r + J*q.k - k*q.j;
    j = RE*q.j - I*q.k + J*q.r + k*q.i;
    k = RE*q.k + I*q.j - J*q.i + k*q.r;

    return * this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator /= ( const Quaternion<T> & q ) {
    T denom = q.norm();
    T RE = r, I = i, J = j;

    r = ( RE*q.r + I*q.i + J*q.j + k*q.k)/denom;
    i = (-RE*q.i + I*q.r - J*q.k + k*q.j)/denom;
    j = (-RE*q.j + I*q.k + J*q.r - k*q.i)/denom;
    k = (-RE*q.k - I*q.j + J*q.i + k*q.r)/denom;

    return * this;
}


template <typename T> Quaternion<T> & Quaternion<T>::operator += ( T x ) {
    r += x;

    return *this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator -= ( T x ) {
    r -= x;

    return *this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator *= ( T x ) {
    if ( Support<T>::is_inf( x ) && norm() > 0 ) {
        r *= ( r == 0 )?Support<T>::sign( x ):x;
        i *= ( i == 0 )?Support<T>::sign( x ):x;
        j *= ( j == 0 )?Support<T>::sign( x ):x;
        k *= ( k == 0 )?Support<T>::sign( x ):x;
    } else {
        r *= x;
        i *= x;
        j *= x;
        k *= x;
    }

    return * this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator /= ( T x ) {
    if ( x == 0.0 && norm() > 0 ) {
        r /= ( r == 0 ) ? Support<T>::sign( x ):x;
        i /= ( i == 0 ) ? Support<T>::sign( x ):x;
        j /= ( j == 0 ) ? Support<T>::sign( x ):x;
        k /= ( k == 0 ) ? Support<T>::sign( x ):x;
    } else {
        r /= x;
        i /= x;
        j /= x;
        k /= x;
    }

    return * this;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator ++() {
    ++r;

    return *this;
}

template <typename T> Quaternion<T> Quaternion<T>::operator ++( int ) {
    Quaternion copy = *this;
    ++r;

    return copy;
}

template <typename T> Quaternion<T> & Quaternion<T>::operator --() {
    --r;

    return *this;
}

template <typename T> Quaternion<T> Quaternion<T>::operator --( int ) {
    Quaternion copy = *this;
    --r;

    return copy;
}

/******************************************
 * Real-valued Functions
 ******************************************/

template <typename T> T Quaternion<T>::real() const {
    return r;
}

template <typename T> T Quaternion<T>::operator []( int n ) const {
    return reinterpret_cast<const T *>( this )[n];
}

template <typename T> T &Quaternion<T>::operator []( int n ) {
    return reinterpret_cast<T *>( this )[n];
}

template <typename T> T Quaternion<T>::imag_i() const {
    return i;
}

template <typename T> T Quaternion<T>::imag_j() const {
    return j;
}

template <typename T> T Quaternion<T>::imag_k() const {
    return k;
}

template <typename T> T Quaternion<T>::csgn() const {
    return is_zero() ? 0.0 : Support<T>::sign( r );
}

template <typename T> T Quaternion<T>::abs() const {
    return is_inf()? Support<T>::POS_INF : std::sqrt( r*r + i*i + j*j + k*k);
}

template <typename T> T Quaternion<T>::norm() const {
    return is_inf() ? Support<T>::POS_INF : r*r + i*i + j*j + k*k;
}

template <typename T> T Quaternion<T>::abs_imag() const {
    return is_inf() ? Support<T>::POS_INF : std::sqrt( i*i + j*j + k*k );
}

template <typename T> T Quaternion<T>::norm_imag() const {
    return is_inf() ? Support<T>::POS_INF : i*i + j*j + k*k;
}

template <typename T> T Quaternion<T>::arg() const {
    return std::atan2( abs_imag(), r );
}

/******************************************
 * Quaternion<T>-valued Functions
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::imag() const {
    return Quaternion<T>( 0.0, i, j, k );
}

template <typename T> Quaternion<T> Quaternion<T>::conj() const {
    return Quaternion<T>( r, -i, -j, -k );
}

template <typename T> Quaternion<T> Quaternion<T>::operator * () const {
    return Quaternion<T>( r, -i, -j, -k );
}

template <typename T> Quaternion<T> Quaternion<T>::signum() const {
    T absq = abs();

    if ( absq == 0.0 || Support<T>::is_nan( absq ) || Support<T>::is_inf( absq ) ) {
        return *this;
    } else {
        return Quaternion<T>( r/absq, i/absq, j/absq, k/absq );
    }
}

template <typename T> Quaternion<T> Quaternion<T>::sqr() const {
    return Quaternion<T>( r*r - i*i - j*j - k*k, 2*r*i, 2*r*j, 2*r*k );
}

template <typename T> Quaternion<T> Quaternion<T>::sqrt() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).sqrt();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::rotate( const Quaternion<T> & q ) const {
    // the assumption is that |q| = 1
    // in this case, q.inverse() == q.conj()

    T rr = q.r*q.r;
    T ii = q.i*q.i;
    T jj = q.j*q.j;
    T kk = q.k*q.k;

    T ri = q.r*q.i;
    T ij = q.i*q.j;
    T ik = q.i*q.k;

    T rj = q.r*q.j;
    T jk = q.j*q.k;

    T rk = q.r*q.k;

    return Quaternion<T>(
        ( ( r == 0 ) ? r : r*(rr + ii + jj + kk) ),
        i*(rr + ii - jj - kk) + 2*(               j*(-rk + ij) + k*( rj + ik)),
        j*(rr - ii + jj - kk) + 2*(i*( rk + ij)                + k*(-ri + jk)),
        k*(rr - ii - jj + kk) + 2*(i*(-rj + ik) + j*( ri + jk)               )
    );
}

/******************************************
 * Boolean-valued Functions
 ******************************************/


template <typename T> bool Quaternion<T>::is_imaginary() const {
    return ( r == 0.0 );
}

template <typename T> bool Quaternion<T>::is_inf() const {
    return
        ( r == Support<T>::POS_INF ) ||
        ( r == Support<T>::NEG_INF ) ||
        ( i == Support<T>::POS_INF ) ||
        ( i == Support<T>::NEG_INF ) ||
        ( j == Support<T>::POS_INF ) ||
        ( j == Support<T>::NEG_INF ) ||
        ( k == Support<T>::POS_INF ) ||
        ( k == Support<T>::NEG_INF );
}

template <typename T> bool Quaternion<T>::is_nan() const {
    return ( r != r ) || ( i != i ) || ( j != j ) || ( k != k );
}

template <typename T> bool Quaternion<T>::is_neg_inf() const {
    return ( r == Support<T>::NEG_INF ) && ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

template <typename T> bool Quaternion<T>::is_pos_inf() const {
    return ( r == Support<T>::POS_INF ) && ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

template <typename T> bool Quaternion<T>::is_real() const {
    return ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

template <typename T> bool Quaternion<T>::is_real_inf() const {
    return ( r == Support<T>::POS_INF || r == Support<T>::NEG_INF ) && ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

template <typename T> bool Quaternion<T>::is_zero() const {
    return ( r == 0.0 ) && ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

/******************************************
 * Multiplier Function
 ******************************************/

template <typename T> inline Quaternion<T> Quaternion<T>::multiplier( T r, T mltplr, const Quaternion<T> & q ) {
    if ( Support<T>::is_nan( mltplr ) || Support<T>::is_inf( mltplr ) ) {
        if ( q.i == 0 && q.j == 0 && q.k == 0 ) {
            return Quaternion<T>( r, mltplr*q.i, mltplr*q.j, mltplr*q.k );
        } else {
            return Quaternion<T>(
                r,
                (q.i == 0) ? Support<T>::sign( mltplr )*q.i : mltplr*q.i,
                (q.j == 0) ? Support<T>::sign( mltplr )*q.j : mltplr*q.j,
                (q.k == 0) ? Support<T>::sign( mltplr )*q.k : mltplr*q.k
            );
        }
    } else {
        return Quaternion<T>( r, mltplr*q.i, mltplr*q.j, mltplr*q.k );
    }
}

template <typename T> inline Quaternion<T> Quaternion<T>::make_inf( T r, T i ) {
    return Quaternion<T>( r, i, i, i );
}

template <typename T> inline Quaternion<T> Quaternion<T>::make_i( T r, T i ) {
    return Quaternion<T>( r, i, 0.0, 0.0 );
}

/******************************************
 * Exponential and Logarithmic Functions
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::exp() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).exp();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::log() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).log();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::log10() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).log10();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::pow( const Quaternion<T> & q ) const {
    return ( log() * q ).exp();
}

template <typename T> Quaternion<T> Quaternion<T>::pow(T x) const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).pow( x );

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::inverse() const {
    if ( is_zero() ) {
        return Quaternion(
             Support<T>::sign(  r )*Support<T>::POS_INF,
            -Support<T>::sign(  i )*Support<T>::POS_INF,
            -Support<T>::sign(  j )*Support<T>::POS_INF,
            -Support<T>::sign(  k )*Support<T>::POS_INF
        );
    } else if ( is_inf() ) {
        return Quaternion(
             Support<T>::sign(  r )*0.0,
            -Support<T>::sign(  i )*0.0,
            -Support<T>::sign(  j )*0.0,
            -Support<T>::sign(  k )*0.0
        );
    } else if ( is_nan() ) {
        return Quaternion( Support<T>::NaN, Support<T>::NaN, Support<T>::NaN, Support<T>::NaN );
    } else {
        T denom = norm();

        return Quaternion( r/denom, -i/denom, -j/denom, -k/denom );
    }
}

/********************************************************************
 * Trigonometric and Hyperbolic Functions
 *
 * For each function f:H -> H, we define:
 *
 *            ~                       Im(q)     ~
 *  f(q) = Re(f(Re(q) + i|Im(q)|)) + ------- Im(f(Re(q) + i|Im(q)|))
 *                                   |Im(q)|
 *       ~
 * where f:C -> C is the complex equivalent of the
 * function f.
 ********************************************************************/

/**********************************************************
 * Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::sin() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).sin();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Complementary Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::cos() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).cos();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Tangent Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::tan() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).tan();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::sec() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).sec();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Complementary Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::csc() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).csc();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Complementary Tangent Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::cot() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).cot();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::sinh() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).sinh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Complementary Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::cosh() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).cosh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Tangent Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::tanh() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).tanh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::sech() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).sech();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Complementary Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::csch() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).csch();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Hyperbolic Complementary Tangent Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::coth() const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).coth();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

// Real Branch Cut:    (-oo, -1) U (1, oo)

/**********************************************************
 * Inverse Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::asin( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    // Branch Cuts:   (-oo, -1) U (1, oo)

    if ( absIm == 0 ) {
        if ( r > 1 ) {
            // Branch cut (1, oo)

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( Support<T>::PI2, std::log( r + std::sqrt( r*r - 1 ) ) );
            } else {
                return multiplier( Support<T>::PI2, std::log( r + std::sqrt( r*r - 1 ) )/absq, q );
            }
        } else if ( r < -1 ) {
            // Branch cut (-oo, -1)

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( -Support<T>::PI2, std::log( -r + std::sqrt( r*r - 1 ) ) );
            } else {
                return multiplier( -Support<T>::PI2, std::log( -r + std::sqrt( r*r - 1 ) )/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).asin();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );

}

// Real Branch Cut:    (-oo, -1) U (1, oo)

/**********************************************************
 * Inverse Complementary Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acos( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    // Branch Cuts:   (-oo, -1) U (1, oo)

    if ( absIm == 0 ) {
        if ( r > 1 ) {
            // Branch cut (1, oo)

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, -std::log( r + std::sqrt( r*r - 1 ) ) );
            } else {
                return multiplier( 0.0, -std::log( r + std::sqrt( r*r - 1 ) )/absq, q );
            }
        } else if ( r < -1 ) {
            // Branch cut (-oo, -1)

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( Support<T>::PI, -std::log( -r + std::sqrt( r*r - 1 ) ) );
            } else {
                return multiplier( Support<T>::PI, -std::log( -r + std::sqrt( r*r - 1 ) )/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).acos();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

// Complex Branch Cut:    (-ooi, -i] U [i, ooi)

/**********************************************************
 * Inverse Tangent Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::atan() const {
    T absIm = abs_imag();

    if ( r == 0 ) {
        if ( absIm == 1 ) {
                return multiplier( Support<T>::NaN, Support<T>::POS_INF, *this );
        } else if ( absIm > 1 ) {
            // Branch cut [ui, Inf)
            //  - ui is a unit purely-imaginary quaternion

            T p = absIm + 1;
            T m = absIm - 1;

            T mltplr = 0.25*std::log( (p*p)/(m*m) )/absIm;

            return multiplier( Support<T>::sign( r )*Support<T>::PI2, mltplr, *this );
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).atan();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Inverse Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::asec( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    // Branch Cut: (-1, 1)

    if ( absIm == 0 ) {
        if ( r == 0.0 ) {
            return multiplier( Support<T>::POS_INF, Support<T>::POS_INF, q );
        } else if ( r > 0.0 && r < 1.0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, std::log( 1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) ) );
            } else {
                return multiplier( 0.0, std::log( 1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) )/absq, q );
            }
        } else if ( r > -1.0 && r < 0.0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( Support<T>::PI, std::log( -1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) ) );
            } else {
                return multiplier( Support<T>::PI, std::log( -1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) )/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).asec();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

// Real Branch Cut:    (-1, 1)
/**********************************************************
 * Inverse Complementary Secant Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acsc( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    if ( absIm == 0.0 ) {
        if ( r == 0.0 ) {
            return multiplier( Support<T>::POS_INF, Support<T>::POS_INF, q );
        } else if ( r > 0.0 && r <  1.0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( Support<T>::PI2, -std::log( 1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) ) );
            } else {
                return multiplier( Support<T>::PI2, -std::log( 1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) )/absq, q );
            }
        } else if ( r > -1.0 && r < 0.0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( -Support<T>::PI2, -std::log( -1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) ) );
            } else {
                return multiplier( -Support<T>::PI2, -std::log( -1.0/r + std::sqrt( 1.0/(r*r) - 1.0 ) )/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).acsc();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Complementary Tangent Function
 * Complex Branch Cut:    (-i, i)
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acot() const {
    T absIm = abs_imag();

    if ( r == 0 ) {
        if ( absIm == 0 ) {
            return Quaternion<T>( Support<T>::PI2, -i, -j, -k );
        } else if ( absIm < 1 ) {
            // Branch cut [ui, Inf)
            //  - ui is a unit purely-imaginary quaternion

            T p = absIm + 1;
            T m = absIm - 1;

            T mltplr = -0.25*std::log( (p*p)/(m*m) )/absIm;

            return multiplier( Support<T>::PI2, mltplr, *this );
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).acot();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

// Complex Branch Cut:    (-ooi, -i) U (i, ooi)

/**********************************************************
 * Inverse Hyperbolic Sine Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::asinh() const {
    T absIm = abs_imag();

    if ( r == 0 ) {
        if ( absIm > 1 ) {
            return multiplier( Support<T>::sign( r )*std::log( absIm + std::sqrt(absIm*absIm - 1) ), Support<T>::PI2/absIm, *this );
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).asinh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Hyperbolic Complementary Sine Function
 * Real Branch Cut:    (-oo, 1)
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acosh( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    if ( absIm == 0 ) {
        if ( r < -1 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( std::log(-r + std::sqrt(r*r - 1)), Support<T>::PI*Support<T>::sign(i) );
            } else {
                return multiplier( std::log(-r + std::sqrt(r*r - 1)), Support<T>::PI/absq, q );
            }
        } else if ( r == -1 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, Support<T>::PI*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, Support<T>::PI/absq, q );
            }
        } else if ( r < 0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, std::acos(r)*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, std::acos(r)/absq, q );
            }
        } else if ( r == 0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, Support<T>::PI2*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, Support<T>::PI2/absq, q );
            }
        } else if ( r < 1 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, std::acos( r )*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, std::acos(r)/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).acosh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Hyperbolic Tangent Function
 * Real Branch Cut:    (-oo, -1] U [1, oo)
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::atanh( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    if ( absIm == 0 ) {
        if ( r == -1 ) {
            return make_inf( Support<T>::NEG_INF, Support<T>::NaN );
        } else if ( r == 1 ) {
            return make_inf( Support<T>::POS_INF, Support<T>::NaN );
        } else if ( r < -1 || r > 1 ) {
            T p = r + 1;
            T m = r - 1;

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.25*std::log( (p*p)/(m*m) ), -Support<T>::sign(r)*Support<T>::PI2 );
            } else {
                return multiplier( 0.25*std::log( (p*p)/(m*m) ), -Support<T>::sign(r)*Support<T>::PI2/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).atanh();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Hyperbolic Secant Function
 * Real Branch Cut:    (-oo, 0] U (1, oo)
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::asech( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    if ( absIm == 0 ) {
        if ( r < -1 || r > 1 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, -std::acos( 1/r )*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, -std::acos( 1/r )/absq, q );
            }
        } else if ( r == -1 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, Support<T>::PI*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, Support<T>::PI/absq, q );
            }
        } else if ( r < 0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( std::log( -1/r + std::sqrt( 1/(r*r) - 1 ) ), -Support<T>::PI );
            } else {
                return multiplier( std::log( -1/r + std::sqrt( 1/(r*r) - 1 ) ), -Support<T>::PI/absq, q );
            }
        } else if ( r == 0 ) {
            return make_inf( Support<T>::POS_INF, Support<T>::NaN );
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).asech();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Hyperbolic Complementary Secant Function
 * Complex Branch Cut:    (-i, i)
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acsch() const {
    T absIm = abs_imag();

    if ( r == 0 ) {
            if ( absIm == 0 ) {
                return make_inf( Support<T>::NEG_INF, Support<T>::NaN );
            } else if ( absIm < 1 ) {
                return multiplier( Support<T>::sign( r )*std::log( 1/absIm + std::sqrt(1/(absIm*absIm) - 1) ), -Support<T>::PI2/absIm, *this );
            }
    }

    Complex<T> z = Complex<T>( r, absIm ).acsch();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}


/**********************************************************
 * Inverse Hyperbolic Complementary Tangent Function
 * Real Branch Cut:    [-1, 1]
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::acoth( const Quaternion<T> & q ) const {
    T absIm = abs_imag();

    if ( absIm == 0 ) {
        if ( r == -1 ) {
            return make_inf( Support<T>::NEG_INF, Support<T>::NaN );
        } else if ( r == 1 ) {
            return make_inf( Support<T>::POS_INF, Support<T>::NaN );
        } else if ( r == 0 ) {
            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.0, -Support<T>::PI2*Support<T>::sign(i) );
            } else {
                return multiplier( 0.0, -Support<T>::PI2/absq, q );
            }
        } else if ( r > -1 && r < 0 ) {
            T p = r + 1;
            T m = r - 1;

            T absq = q.abs_imag();

            if ( q == I || absq == 0 ) {
                return make_i( 0.25*std::log( (p*p)/(m*m) ), -Support<T>::sign(r)*Support<T>::PI2 );
            } else {
                return multiplier( 0.25*std::log( (p*p)/(m*m) ), -Support<T>::sign(r)*Support<T>::PI2/absq, q );
            }
        }
    }

    Complex<T> z = Complex<T>( r, absIm ).acoth();

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/**********************************************************
 * Bessel J Function
 **********************************************************/

template <typename T> Quaternion<T> Quaternion<T>::bessel_J( int n ) const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).bessel_J( n );

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/******************************************
 * Integer Functions
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::ceil() const {
    return Quaternion<T>(
        std::ceil(r), std::ceil(i), std::ceil(j), std::ceil(k)
    );
}

template <typename T> Quaternion<T> Quaternion<T>::floor() const {
    return Quaternion<T>(
        std::floor(r), std::floor(i), std::floor(j), std::floor(k)
    );
}


/******************************************
 * Horner's Rule
 *
 *   The polynomial is defined by giving the highest
 *   coefficient first:
 *
 *            n - 1         n - 2
 *      v[0]*q      + v[1]*q      + ... + v[n-2]*q + v[n-1]
 *
 *   This is the same as with Matlab.  Because quaternions are
 *   not commutative, this only makes sense if the coefficients
 *   and offsets are real.
 *
 *        Re(q) + i |Imag(q)|
 *
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::horner( T * v, unsigned int n ) const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).horner( v, n );

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

template <typename T> Quaternion<T> Quaternion<T>::horner( T * v, T * c, unsigned int n ) const {
    T absIm = abs_imag();
    Complex<T> z = Complex<T>( r, absIm ).horner( v, c, n );

    T mltplr;

    if ( absIm == 0.0 ) {
        mltplr = z.imag_i();
    } else {
        mltplr = z.imag_i() / absIm;
    }

    return multiplier( z.real(), mltplr, *this );
}

/******************************************
 * Random Factories
 ******************************************/


template <typename T> Quaternion<T> Quaternion<T>::random() {
    return Quaternion<T>(
        (static_cast<T>( rand() ))/RAND_MAX,
        (static_cast<T>( rand() ))/RAND_MAX,
        (static_cast<T>( rand() ))/RAND_MAX,
        (static_cast<T>( rand() ))/RAND_MAX
    );
}

template <typename T> Quaternion<T> Quaternion<T>::random_imag() {
    return Quaternion<T>( 0.0, (static_cast<T>( rand() ))/RAND_MAX, (static_cast<T>( rand() ))/RAND_MAX, (static_cast<T>( rand() ))/RAND_MAX );
}

template <typename T> Quaternion<T> Quaternion<T>::random_real() {
    return Quaternion<T>( (static_cast<T>( rand() ))/RAND_MAX, 0.0, 0.0, 0.0 );
}

/******************************************
 * Binary Arithmetic Operators
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::operator + ( const Quaternion<T> & z ) const {
    return Quaternion<T>( r + z.r, i + z.i, j + z.j, k + z.k );
}

template <typename T> Quaternion<T> Quaternion<T>::operator + ( T x ) const {
    return Quaternion<T>( r + x, i, j, k );
}

template <typename T> Quaternion<T> operator + ( T x, const Quaternion<T> & z ) {
    return Quaternion<T>( x + z.real(), z.imag_i(), z.imag_j(), z.imag_k() );
}

template <typename T> Quaternion<T> operator + ( long x, const Quaternion<T> & z ) {
    return Quaternion<T>( static_cast<T>( x ) + z.real(), z.imag_i(), z.imag_j(), z.imag_k() );
}

template <typename T> Quaternion<T> Quaternion<T>::operator - ( const Quaternion<T> & z ) const {
    return Quaternion<T>( r - z.r, i - z.i, j - z.j, k - z.k );
}

template <typename T> Quaternion<T> Quaternion<T>::operator - ( T x ) const {
    return Quaternion<T>( r - x, i, j, k );
}

template <typename T> Quaternion<T> operator - ( T x, const Quaternion<T> & z ) {
    return Quaternion<T>( x - z.real(), -z.imag_i(), -z.imag_j(), -z.imag_k() );
}

template <typename T> Quaternion<T> operator - ( long x, const Quaternion<T> & z ) {
    return Quaternion<T>( static_cast<T>( x ) - z.real(), -z.imag_i(), -z.imag_j(), -z.imag_k() );
}

template <typename T> Quaternion<T> Quaternion<T>::operator * ( const Quaternion<T> & q ) const {
    return Quaternion<T>(
        r*q.r - i*q.i - j*q.j - k*q.k,
        r*q.i + i*q.r + j*q.k - k*q.j,
        r*q.j - i*q.k + j*q.r + k*q.i,
        r*q.k + i*q.j - j*q.i + k*q.r
    );
}

template <typename T> Quaternion<T> Quaternion<T>::operator * ( T x ) const {
    if ( Support<T>::is_inf( x ) && norm() > 0 ) {
        return Quaternion<T>(
            ( (  r == 0 )?Support<T>::sign(x):x )*r,
            ( (  i == 0 )?Support<T>::sign(x):x )*i,
            ( (  j == 0 )?Support<T>::sign(x):x )*j,
            ( (  k == 0 )?Support<T>::sign(x):x )*k
        );
    } else {
        return Quaternion<T>( x*r, x*i, x*j, x*k );
    }
}

template <typename T> Quaternion<T> operator * ( T x, const Quaternion<T> & q ) {
    return q.operator * ( x );
}

template <typename T> Quaternion<T> operator * ( long x, const Quaternion<T> & q ) {
    return q.operator * ( static_cast<T>( x ) );
}

template <typename T> Quaternion<T> Quaternion<T>::operator / ( const Quaternion<T> & q ) const {
    T denom = q.norm();

    return Quaternion<T>(
        ( r*q.r + i*q.i + j*q.j + k*q.k)/denom,
        (-r*q.i + i*q.r - j*q.k + k*q.j)/denom,
        (-r*q.j + i*q.k + j*q.r - k*q.i)/denom,
        (-r*q.k - i*q.j + j*q.i + k*q.r)/denom
    );
}

template <typename T> Quaternion<T> Quaternion<T>::operator / ( T x ) const {
    if ( x == 0.0 && norm() > 0 ) {
        return Quaternion<T>(
            r / ( ( r == 0 )?Support<T>::sign( x ):x ),
            i / ( ( i == 0 )?Support<T>::sign( x ):x ),
            j / ( ( j == 0 )?Support<T>::sign( x ):x ),
            k / ( ( k == 0 )?Support<T>::sign( x ):x )
        );
    } else {
        return Quaternion<T>( r/x, i/x, j/x, k/x );
    }
}

template <typename T> Quaternion<T> operator / ( T x, const Quaternion<T> & q ) {
    T mltplr = x/q.norm();

    return Quaternion<T>( mltplr*q.real(), -mltplr*q.imag_i(), -mltplr*q.imag_j(), -mltplr*q.imag_k() );
}

template <typename T> Quaternion<T> operator / ( long x, const Quaternion<T> & q ) {
    T mltplr = static_cast<T>( x )/q.norm();

    return Quaternion<T>( mltplr*q.real(), -mltplr*q.imag_i(), -mltplr*q.imag_j(), -mltplr*q.imag_k() );
}

/******************************************
 * Unary Arithmetic Operators
 ******************************************/

template <typename T> Quaternion<T> Quaternion<T>::operator - () const {
    return Quaternion<T>( -r, -i, -j, -k );
}

/******************************************
 * Binary Boolean Operators
 ******************************************/

template <typename T> bool Quaternion<T>::operator == ( const Quaternion<T> & q ) const {
    return ( r == q.r ) && ( i == q.i ) && ( j == q.j ) && ( k == q.k );
}

template <typename T> bool Quaternion<T>::operator == ( T x ) const {
    return ( r == x ) && ( i == 0.0 ) && ( j == 0.0 ) && ( k == 0.0 );
}

template <typename T> bool operator == ( T x, const Quaternion<T> & q ) {
    return q.operator == ( x );
}

template <typename T> bool operator == ( long x, const Quaternion<T> & q ) {
    return q.operator == ( static_cast<T>( x ) );
}

template <typename T> bool Quaternion<T>::operator != ( const Quaternion<T> & q ) const {
    return ( r != q.r ) || ( i != q.i ) || ( j != q.j ) || ( k != q.k );
}

template <typename T> bool Quaternion<T>::operator != ( T x ) const {
    return ( r != x ) || ( i != 0.0 ) || ( j != 0.0 ) || ( k != 0.0 );
}

template <typename T> bool operator != ( T x, const Quaternion<T> & q ) {
    return q.operator != ( x );
}

template <typename T> bool operator != ( long x, const Quaternion<T> & q ) {
    return q.operator != ( static_cast<T>( x ) );
}

/******************************************
 * IO Stream Operators
 ******************************************/

template <typename T> std::ostream & operator << ( std::ostream & out, const Quaternion<T> & z ) {
    Support<T>::print_real( z.real(), out );
    Support<T>::print_imaginary( z.imag_i(), 'i', out );
    Support<T>::print_imaginary( z.imag_j(), 'j', out );
    Support<T>::print_imaginary( z.imag_k(), 'k', out );

    return out;
}

/******************************************
 * ************************************** *
 * *                                    * *
 * *      Procedural Functions          * *
 * *                                    * *
 * ************************************** *
 ******************************************/

/******************************************
 * Real-valued Functions
 ******************************************/

template <typename T> T real( const Quaternion<T> & q ) {
    return q.real();
}

template <typename T> T imag_i( const Quaternion<T> & q ) {
    return q.imag_i();
}

template <typename T> T imag_j( const Quaternion<T> & q ) {
    return q.imag_j();
}

template <typename T> T imag_k( const Quaternion<T> & q ) {
    return q.imag_k();
}

template <typename T> T csgn( const Quaternion<T> & q ) {
    return q.csgn();
}

template <typename T> T abs( const Quaternion<T> & q ) {
    return q.abs();
}

template <typename T> T norm( const Quaternion<T> & q ) {
    return q.norm();
}

template <typename T> T abs_imag( const Quaternion<T> & z ) {
    return z.abs_imag();
}

template <typename T> T norm_imag( const Quaternion<T> & z ) {
    return z.norm_imag();
}

/******************************************
 * Quaternion-valued Functions
 ******************************************/

template <typename T> Quaternion<T> imag( const Quaternion<T> & q ) {
    return q.imag();
}

template <typename T> Quaternion<T> conj( const Quaternion<T> & q ) {
    return q.conj();
}

template <typename T> Quaternion<T> signum( const Quaternion<T> & q ) {
    return q.signum();
}

template <typename T> Quaternion<T> sqr( const Quaternion<T> & q ) {
    return q.sqr();
}

template <typename T> Quaternion<T> sqrt( const Quaternion<T> & q ) {
    return q.sqrt();
}

template <typename T> Quaternion<T> rotate( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.rotate( p );
}

/******************************************
 * Exponential and Logarithmic Functions
 ******************************************/

template <typename T> Quaternion<T> exp( const Quaternion<T> & q ) {
    return q.exp();
}

template <typename T> Quaternion<T> log( const Quaternion<T> & q ) {
    return q.log();
}

template <typename T> Quaternion<T> log10( const Quaternion<T> & q ) {
    return q.log10();
}

template <typename T> Quaternion<T> pow( const Quaternion<T> & q, const Quaternion<T> & w ) {
    return q.pow( w );
}

template <typename T> Quaternion<T> pow( const Quaternion<T> & q, T x ) {
    return q.pow( x );
}

template <typename T> Quaternion<T> inverse( const Quaternion<T> & q ) {
    return q.inverse();
}

/******************************************
 * Trigonometric and Hyperbolic Functions
 ******************************************/

template <typename T> Quaternion<T> sin( const Quaternion<T> & q ) {
    return q.sin();
}

template <typename T> Quaternion<T> cos( const Quaternion<T> & q ) {
    return q.cos();
}

template <typename T> Quaternion<T> tan( const Quaternion<T> & q ) {
    return q.tan();
}

template <typename T> Quaternion<T> sec( const Quaternion<T> & q ) {
    return q.sec();
}

template <typename T> Quaternion<T> csc( const Quaternion<T> & q ) {
    return q.csc();
}

template <typename T> Quaternion<T> cot( const Quaternion<T> & q ) {
    return q.cot();
}

template <typename T> Quaternion<T> sinh( const Quaternion<T> & q ) {
    return q.sinh();
}

template <typename T> Quaternion<T> cosh( const Quaternion<T> & q ) {
    return q.cosh();
}

template <typename T> Quaternion<T> tanh( const Quaternion<T> & q ) {
    return q.tanh();
}

template <typename T> Quaternion<T> sech( const Quaternion<T> & q ) {
    return q.sech();
}

template <typename T> Quaternion<T> csch( const Quaternion<T> & q ) {
    return q.csch();
}

template <typename T> Quaternion<T> coth( const Quaternion<T> & q ) {
    return q.coth();
}

template <typename T> Quaternion<T> asin( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.asin( p );
}

template <typename T> Quaternion<T> acos( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.acos( p );
}

template <typename T> Quaternion<T> atan( const Quaternion<T> & q ) {
    return q.atan();
}

template <typename T> Quaternion<T> asec( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.asec( p );
}

template <typename T> Quaternion<T> acsc( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return p.acsc( p );
}

template <typename T> Quaternion<T> acot( const Quaternion<T> & q ) {
    return q.acot();
}

template <typename T> Quaternion<T> asinh( const Quaternion<T> & q ) {
    return q.asinh();
}

template <typename T> Quaternion<T> acosh( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.acosh( p );
}

template <typename T> Quaternion<T> atanh( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.atanh( p );
}

template <typename T> Quaternion<T> asech( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.asech( p );
}

template <typename T> Quaternion<T> acsch( const Quaternion<T> & q ) {
    return q.acsch();
}

template <typename T> Quaternion<T> acoth( const Quaternion<T> & q, const Quaternion<T> & p ) {
    return q.acoth( p );
}

template <typename T> Quaternion<T> bessel_J( int n, const Quaternion<T> & z ) {
    return z.bessel_J( n );
}

/******************************************
 * Integer Functions
 ******************************************/

template <typename T> Quaternion<T> floor( const Quaternion<T> & q ) {
    return q.floor();
}

template <typename T> Quaternion<T> ceil( const Quaternion<T> & q ) {
    return q.ceil();
}

/******************************************
 * Horner's Rule
 ******************************************/

template <typename T> Quaternion<T> horner( const Quaternion<T> & q, T * v, unsigned int n ) {
    return q.horner( v, n );
}

template <typename T> Quaternion<T> horner( const Quaternion<T> & q, T * v, T * c, unsigned int n ) {
    return q.horner( v, c, n );
}

/**************************************************
 * ********************************************** *
 * *                                            * *
 * *    Double-precision Floating-point         * *
 * *    Instance of Template                    * *
 * *                                            * *
 * ********************************************** *
 **************************************************/

template class Quaternion<double>;
template std::ostream & operator << ( std::ostream &, const Quaternion<double> & );

template Quaternion<double> operator + ( double, const Quaternion<double> & );
template Quaternion<double> operator + ( long, const Quaternion<double> & );
template Quaternion<double> operator - ( double, const Quaternion<double> & );
template Quaternion<double> operator - ( long, const Quaternion<double> & );
template Quaternion<double> operator * ( double, const Quaternion<double> & );
template Quaternion<double> operator * ( long, const Quaternion<double> & );
template Quaternion<double> operator / ( double, const Quaternion<double> & );
template Quaternion<double> operator / ( long, const Quaternion<double> & );

template bool operator == ( double, const Quaternion<double> & );
template bool operator == ( long, const Quaternion<double> & );
template bool operator != ( double, const Quaternion<double> & );
template bool operator != ( long, const Quaternion<double> & );

template <> const Quaternion<double> Quaternion<double>::ZERO = Quaternion<double>( 0, 0, 0, 0 );
template <> const Quaternion<double>  Quaternion<double>::ONE = Quaternion<double>( 1, 0, 0, 0 );
template <> const Quaternion<double>    Quaternion<double>::I = Quaternion<double>( 0, 1, 0, 0 );
template <> const Quaternion<double>    Quaternion<double>::J = Quaternion<double>( 0, 0, 1, 0 );
template <> const Quaternion<double>    Quaternion<double>::K = Quaternion<double>( 0, 0, 0, 1 );

template <> const Quaternion<double> Quaternion<double>::UNITS[4] = {
    Quaternion<double>::ONE,
    Quaternion<double>::I,
    Quaternion<double>::J,
    Quaternion<double>::K
};

template double real( const Quaternion<double> & );
template double imag_i( const Quaternion<double> & );
template double imag_j( const Quaternion<double> & );
template double imag_k( const Quaternion<double> & );
template double csgn( const Quaternion<double> & );
template double abs( const Quaternion<double> & );
template double norm( const Quaternion<double> & );
template double abs_imag( const Quaternion<double> & );
template double norm_imag( const Quaternion<double> & );
template Quaternion<double> imag( const Quaternion<double> & );
template Quaternion<double> conj( const Quaternion<double> & );
template Quaternion<double> signum( const Quaternion<double> & );
template Quaternion<double> sqr( const Quaternion<double> & );
template Quaternion<double> sqrt( const Quaternion<double> & );
template Quaternion<double> rotate( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> exp( const Quaternion<double> & );
template Quaternion<double> log( const Quaternion<double> & );
template Quaternion<double> log10( const Quaternion<double> & );
template Quaternion<double> pow( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> pow( const Quaternion<double> &, double );
template Quaternion<double> inverse( const Quaternion<double> & );
template Quaternion<double> sin( const Quaternion<double> & );
template Quaternion<double> cos( const Quaternion<double> & );
template Quaternion<double> tan( const Quaternion<double> & );
template Quaternion<double> sec( const Quaternion<double> & );
template Quaternion<double> csc( const Quaternion<double> & );
template Quaternion<double> cot( const Quaternion<double> & );
template Quaternion<double> sinh( const Quaternion<double> & );
template Quaternion<double> cosh( const Quaternion<double> & );
template Quaternion<double> tanh( const Quaternion<double> & );
template Quaternion<double> sech( const Quaternion<double> & );
template Quaternion<double> csch( const Quaternion<double> & );
template Quaternion<double> coth( const Quaternion<double> & );
template Quaternion<double> asin( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> acos( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> atan( const Quaternion<double> & );
template Quaternion<double> asec( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> acsc( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> acot( const Quaternion<double> & );
template Quaternion<double> asinh( const Quaternion<double> & );
template Quaternion<double> acosh( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> atanh( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> asech( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> acsch( const Quaternion<double> & );
template Quaternion<double> acoth( const Quaternion<double> &, const Quaternion<double> & );
template Quaternion<double> bessel_J( int, const Quaternion<double> & );
template Quaternion<double> floor( const Quaternion<double> & );
template Quaternion<double> ceil( const Quaternion<double> & );
template Quaternion<double> horner( const Quaternion<double> &, double *, unsigned int );
template Quaternion<double> horner( const Quaternion<double> &, double *, double *, unsigned int );

/**************************************************
 * ********************************************** *
 * *                                            * *
 * *    Floating-point Instance of Template     * *
 * *                                            * *
 * ********************************************** *
 **************************************************/

template class Quaternion<float>;
template std::ostream & operator << ( std::ostream &, const Quaternion<float> & );

template Quaternion<float> operator + ( float, const Quaternion<float> & );
template Quaternion<float> operator + ( long, const Quaternion<float> & );
template Quaternion<float> operator - ( float, const Quaternion<float> & );
template Quaternion<float> operator - ( long, const Quaternion<float> & );
template Quaternion<float> operator * ( float, const Quaternion<float> & );
template Quaternion<float> operator * ( long, const Quaternion<float> & );
template Quaternion<float> operator / ( float, const Quaternion<float> & );
template Quaternion<float> operator / ( long, const Quaternion<float> & );

template bool operator == ( float, const Quaternion<float> & );
template bool operator == ( long, const Quaternion<float> & );
template bool operator != ( float, const Quaternion<float> & );
template bool operator != ( long, const Quaternion<float> & );

template <> const Quaternion<float> Quaternion<float>::ZERO = Quaternion<float>( 0, 0, 0, 0 );
template <> const Quaternion<float>  Quaternion<float>::ONE = Quaternion<float>( 1, 0, 0, 0 );
template <> const Quaternion<float>    Quaternion<float>::I = Quaternion<float>( 0, 1, 0, 0 );
template <> const Quaternion<float>    Quaternion<float>::J = Quaternion<float>( 0, 0, 1, 0 );
template <> const Quaternion<float>    Quaternion<float>::K = Quaternion<float>( 0, 0, 0, 1 );

template <> const Quaternion<float> Quaternion<float>::UNITS[4] = {
    Quaternion<float>::ONE,
    Quaternion<float>::I,
    Quaternion<float>::J,
    Quaternion<float>::K
};

template float real( const Quaternion<float> & );
template float imag_i( const Quaternion<float> & );
template float imag_j( const Quaternion<float> & );
template float imag_k( const Quaternion<float> & );
template float csgn( const Quaternion<float> & );
template float abs( const Quaternion<float> & );
template float norm( const Quaternion<float> & );
template float abs_imag( const Quaternion<float> & );
template float norm_imag( const Quaternion<float> & );
template Quaternion<float> imag( const Quaternion<float> & );
template Quaternion<float> conj( const Quaternion<float> & );
template Quaternion<float> signum( const Quaternion<float> & );
template Quaternion<float> sqr( const Quaternion<float> & );
template Quaternion<float> sqrt( const Quaternion<float> & );
template Quaternion<float> rotate( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> exp( const Quaternion<float> & );
template Quaternion<float> log( const Quaternion<float> & );
template Quaternion<float> log10( const Quaternion<float> & );
template Quaternion<float> pow( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> pow( const Quaternion<float> &, float );
template Quaternion<float> inverse( const Quaternion<float> & );
template Quaternion<float> sin( const Quaternion<float> & );
template Quaternion<float> cos( const Quaternion<float> & );
template Quaternion<float> tan( const Quaternion<float> & );
template Quaternion<float> sec( const Quaternion<float> & );
template Quaternion<float> csc( const Quaternion<float> & );
template Quaternion<float> cot( const Quaternion<float> & );
template Quaternion<float> sinh( const Quaternion<float> & );
template Quaternion<float> cosh( const Quaternion<float> & );
template Quaternion<float> tanh( const Quaternion<float> & );
template Quaternion<float> sech( const Quaternion<float> & );
template Quaternion<float> csch( const Quaternion<float> & );
template Quaternion<float> coth( const Quaternion<float> & );
template Quaternion<float> asin( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> acos( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> atan( const Quaternion<float> & );
template Quaternion<float> asec( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> acsc( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> acot( const Quaternion<float> & );
template Quaternion<float> asinh( const Quaternion<float> & );
template Quaternion<float> acosh( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> atanh( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> asech( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> acsch( const Quaternion<float> & );
template Quaternion<float> acoth( const Quaternion<float> &, const Quaternion<float> & );
template Quaternion<float> bessel_J( int, const Quaternion<float> & );
template Quaternion<float> floor( const Quaternion<float> & );
template Quaternion<float> ceil( const Quaternion<float> & );
template Quaternion<float> horner( const Quaternion<float> &, float *, unsigned int );
template Quaternion<float> horner( const Quaternion<float> &, float *, float *, unsigned int );

/************************************************************************
 * ******************************************************************** *
 * *                                                                  * *
 * *    Long Double-precision Floating-point Instance of Template     * *
 * *                                                                  * *
 * ******************************************************************** *
 ************************************************************************/

template class Quaternion<long double>;
template std::ostream & operator << ( std::ostream &, const Quaternion<long double> & );

template Quaternion<long double> operator + ( long double, const Quaternion<long double> & );
template Quaternion<long double> operator + ( long, const Quaternion<long double> & );
template Quaternion<long double> operator - ( long double, const Quaternion<long double> & );
template Quaternion<long double> operator - ( long, const Quaternion<long double> & );
template Quaternion<long double> operator * ( long double, const Quaternion<long double> & );
template Quaternion<long double> operator * ( long, const Quaternion<long double> & );
template Quaternion<long double> operator / ( long double, const Quaternion<long double> & );
template Quaternion<long double> operator / ( long, const Quaternion<long double> & );

template bool operator == ( long double, const Quaternion<long double> & );
template bool operator == ( long, const Quaternion<long double> & );
template bool operator != ( long double, const Quaternion<long double> & );
template bool operator != ( long, const Quaternion<long double> & );

template <> const Quaternion<long double> Quaternion<long double>::ZERO = Quaternion<long double>( 0, 0, 0, 0 );
template <> const Quaternion<long double>  Quaternion<long double>::ONE = Quaternion<long double>( 1, 0, 0, 0 );
template <> const Quaternion<long double>    Quaternion<long double>::I = Quaternion<long double>( 0, 1, 0, 0 );
template <> const Quaternion<long double>    Quaternion<long double>::J = Quaternion<long double>( 0, 0, 1, 0 );
template <> const Quaternion<long double>    Quaternion<long double>::K = Quaternion<long double>( 0, 0, 0, 1 );

template <> const Quaternion<long double> Quaternion<long double>::UNITS[4] = {
    Quaternion<long double>::ONE,
    Quaternion<long double>::I,
    Quaternion<long double>::J,
    Quaternion<long double>::K
};

template long double real( const Quaternion<long double> & );
template long double imag_i( const Quaternion<long double> & );
template long double imag_j( const Quaternion<long double> & );
template long double imag_k( const Quaternion<long double> & );
template long double csgn( const Quaternion<long double> & );
template long double abs( const Quaternion<long double> & );
template long double norm( const Quaternion<long double> & );
template long double abs_imag( const Quaternion<long double> & );
template long double norm_imag( const Quaternion<long double> & );
template Quaternion<long double> imag( const Quaternion<long double> & );
template Quaternion<long double> conj( const Quaternion<long double> & );
template Quaternion<long double> signum( const Quaternion<long double> & );
template Quaternion<long double> sqr( const Quaternion<long double> & );
template Quaternion<long double> sqrt( const Quaternion<long double> & );
template Quaternion<long double> rotate( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> exp( const Quaternion<long double> & );
template Quaternion<long double> log( const Quaternion<long double> & );
template Quaternion<long double> log10( const Quaternion<long double> & );
template Quaternion<long double> pow( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> pow( const Quaternion<long double> &, long double );
template Quaternion<long double> inverse( const Quaternion<long double> & );
template Quaternion<long double> sin( const Quaternion<long double> & );
template Quaternion<long double> cos( const Quaternion<long double> & );
template Quaternion<long double> tan( const Quaternion<long double> & );
template Quaternion<long double> sec( const Quaternion<long double> & );
template Quaternion<long double> csc( const Quaternion<long double> & );
template Quaternion<long double> cot( const Quaternion<long double> & );
template Quaternion<long double> sinh( const Quaternion<long double> & );
template Quaternion<long double> cosh( const Quaternion<long double> & );
template Quaternion<long double> tanh( const Quaternion<long double> & );
template Quaternion<long double> sech( const Quaternion<long double> & );
template Quaternion<long double> csch( const Quaternion<long double> & );
template Quaternion<long double> coth( const Quaternion<long double> & );
template Quaternion<long double> asin( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> acos( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> atan( const Quaternion<long double> & );
template Quaternion<long double> asec( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> acsc( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> acot( const Quaternion<long double> & );
template Quaternion<long double> asinh( const Quaternion<long double> & );
template Quaternion<long double> acosh( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> atanh( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> asech( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> acsch( const Quaternion<long double> & );
template Quaternion<long double> acoth( const Quaternion<long double> &, const Quaternion<long double> & );
template Quaternion<long double> bessel_J( int, const Quaternion<long double> & );
template Quaternion<long double> floor( const Quaternion<long double> & );
template Quaternion<long double> ceil( const Quaternion<long double> & );
template Quaternion<long double> horner( const Quaternion<long double> &, long double *, unsigned int );
template Quaternion<long double> horner( const Quaternion<long double> &, long double *, long double *, unsigned int );
