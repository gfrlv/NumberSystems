/******************************************
 * C++ Complex Numbers
 * Version: 1.0.9
 * Author:  Douglas Wilhelm Harder
 * Date:    2008/03/03
 *
 * Copyright (c) 2006-2008 by Douglas Wilhelm Harder.
 * All rights reserved.
 *
 * Acknowledgments:
 *      Andre Kostro
 ******************************************/

#include "Complex.h"
#include "Support.h"
#include <iostream>
#include <cmath>
#include <string>

template<typename T> inline char Complex<T>::use_symbol( char sym ) {
    char tmp = imaginary_symbol;

    imaginary_symbol = sym;

    return tmp;
}

template<typename T> inline char Complex<T>::get_symbol() {
    return imaginary_symbol;
}

/******************************************
 * Constructors
 ******************************************/

// (a, b) -> a + bi

template <typename T> Complex<T>::Complex( T re, T im ):r(re), i(im) {
    // empty
}

template <typename T> Complex<T>::Complex( T re ):r(re), i(0.0) {
    // empty
}

/******************************************
 * Assignment Operator
 ******************************************/

// z = 3 -> z == 3 + 0i

template <typename T> inline const Complex<T> & Complex<T>::operator = ( T x ) {
    r = x;
    i = 0.0;

    return *this;
}

/******************************************
 * Mutating Arithmetic Operators
 ******************************************/

// (a + bi) + (c + di) = (a + c) + (b + d)i

template <typename T> inline Complex<T> & Complex<T>::operator += ( const Complex<T> & z ) {
    r += z.r;
    i += z.i;

    return *this;
}

// (a + bi) + x = (a + c) + bi

template <typename T> inline Complex<T> & Complex<T>::operator += ( T x ) {
    r += x;

    return *this;
}

// (a + bi) - (c + di) = (a - c) + (b - d)i

template <typename T> inline Complex<T> & Complex<T>::operator -= ( const Complex<T> & z ) {
    r -= z.r;
    i -= z.i;

    return *this;
}

// (a + bi) - x = (a - x) + bi

template <typename T> inline Complex<T> & Complex<T>::operator -= ( T x ) {
    r -= x;

    return *this;
}

// (a + bi) * (c + di) = (a*c - b*d) + (a*c + b*d)i

template <typename T> inline Complex<T> & Complex<T>::operator *= ( const Complex<T> & z ) {
    T RE = r;

    r =  r * z.r - i * z.i;
    i = RE * z.i + z.r * i;

    return *this;
}

// if ( x != oo )
//     (a + bi) * x = ax + bxi
// else
//     (a + 0i) * oo =  oo +   0i
//     (0 + bi) * oo =   0 +  ooi
//     (0 + 0i) * oo = NaN + NaNi

template <typename T> inline Complex<T> & Complex<T>::operator *= ( T x ) {
    if ( Support<T>::is_inf( x ) && !is_zero() ) {
        r *= ( r == 0.0 ) ? Support<T>::sign( x ) : x;
        i *= ( i == 0.0 ) ? Support<T>::sign( x ) : x;
    } else {
        r *= x;
        i *= x;
    }

    return * this;
}

// n = c^2 + d^2
// (a + bi) / (c + di) = (a*c + b*d)/n + ((a*c - b*d)/n)i

template <typename T> inline Complex<T> & Complex<T>::operator /= ( const Complex<T> & z ) {
    T denom = z.norm();
    T RE = r;

    r = ( RE*z.r + i * z.i)/denom;
    i = (-RE*z.i + z.r * i)/denom;

    return * this;
}

// if ( x != 0 )
//     (a + bi) / x = (a/x) + (b/x)i
// else
//     (a + 0i) / 0 =  oo +   0i
//     (0 + bi) / 0 =   0 +  ooi
//     (0 + 0i) / 0 = NaN + NaNi

template <typename T> inline Complex<T> & Complex<T>::operator /= ( T x ) {
    if ( ( x == 0.0 ) && !is_zero() ) {
        r /= ( r == 0.0 ) ? Support<T>::sign( x ) : x;
        i /= ( i == 0.0 ) ? Support<T>::sign( x ) : x;
    } else {
        r /= x;
        i /= x;
    }

    return * this;
}

// ++(a + bi) = (a + 1) + bi
// Returns (a + 1) + bi

template <typename T> inline Complex<T> & Complex<T>::operator ++() {
    ++r;

    return *this;
}

// (a + bi)++ = (a + 1) + bi
// Returns a + bi

template <typename T> inline Complex<T> Complex<T>::operator ++( int ) {
    Complex<T> copy = *this;
    ++r;

    return copy;
}

// --(a + bi) = (a - 1) + bi
// Returns (a - 1) + bi

template <typename T> inline Complex<T> & Complex<T>::operator --() {
    --r;

    return *this;
}

// (a + bi)-- = (a - 1) + bi
// Returns a + bi

template <typename T> inline Complex<T> Complex<T>::operator --( int ) {
    Complex<T> copy = *this;
    --r;

    return copy;
}

/******************************************
 * Real-valued Functions
 ******************************************/

template <typename T> inline T Complex<T>::real() const {
    return r;
}

template <typename T> inline T Complex<T>::operator []( int n ) const {
    return reinterpret_cast<const T *>( this )[n];
}

template <typename T> inline T & Complex<T>::operator []( int n ) {
    return reinterpret_cast<T *>( this )[n];
}

template <typename T> inline T Complex<T>::imag_i() const {
    return i;
}

template <typename T> inline T Complex<T>::csgn() const {
    return is_zero() ? 0 : Support<T>::sign( r );
}

// if a != oo && b != oo
//     (a + bi) -> sqrt( a^2 + b^2 )
// else
//     ( a + ooi) -> oo
//     (oo +  bi) -> oo
//     (oo + ooi) -> oo

template <typename T> inline T Complex<T>::abs() const {
    return is_inf()? Support<T>::POS_INF : std::sqrt( r*r + i*i );
}

// if a != oo && b != oo
//     (a + bi) -> a^2 + b^2
// else
//     ( a + ooi) -> oo
//     (oo +  bi) -> oo
//     (oo + ooi) -> oo

template <typename T> inline T Complex<T>::norm() const {
    return is_inf() ? Support<T>::POS_INF : r*r + i*i;
}

// a + bi -> |b|

template <typename T> inline T Complex<T>::abs_imag() const {
    return ( i >= 0.0 ) ? i : -i;
}

// a + bi -> b^2

template <typename T> inline T Complex<T>::norm_imag() const {
    return i*i;
}

// (a + bi) -> arctan( b, a )

template <typename T> inline T Complex<T>::arg() const {
    return std::atan2( i, r );
}

/******************************************
 * Complex<T>-valued Functions
 ******************************************/

// (a + bi) -> (0 + bi)

template <typename T> inline Complex<T> Complex<T>::imag() const {
    return Complex<T>( 0.0, i );
}

// (a + bi) -> (a - bi)

template <typename T> inline Complex<T> Complex<T>::conj() const {
    return Complex<T>( r, -i );
}

template <typename T> inline Complex<T> Complex<T>::operator * () const {
    return Complex<T>( r, -i );
}

template <typename T> inline Complex<T> Complex<T>::signum() const {
    T absq = abs();

    if ( absq == 0.0 || Support<T>::is_nan( absq ) || Support<T>::is_inf( absq ) ) {
        return *this;
    } else {
        return Complex<T>( r/absq, i/absq );
    }
}

// (a + bi)^2 = (a^2 - b^2) + 2abi

template <typename T> inline Complex<T> Complex<T>::sqr() const {
    return Complex<T>( r*r - i*i, 2*r*i );
}

// Branch cut:    (-oo, 0)
// Thanks:   Chris Saunders

template <typename T> inline Complex<T> Complex<T>::sqrt() const {
    T srss = std::sqrt( r*r + i*i );

    return Complex<T>(
        0.5*std::sqrt( 2*(srss + r) ),
        0.5*Complex<T>( i, -r ).csgn()*std::sqrt( 2*( srss - r ) )
    );
}

/******************************************
 * Boolean-valued Functions
 ******************************************/

// (0 + bi) -> true
// false otherwise

template <typename T> inline bool Complex<T>::is_imaginary() const {
    return ( r == 0.0 );
}

// if a == +/-oo || b == +/- oo
//     (a + bi) -> true
// false otherwise

template <typename T> inline bool Complex<T>::is_inf() const {
    return
        ( r == Support<T>::POS_INF ) ||
        ( r == Support<T>::NEG_INF ) ||
        ( i == Support<T>::POS_INF ) ||
        ( i == Support<T>::NEG_INF );
}

// if a == NaN or b == NaN
//     (a + bi) -> true
// false otherwise

template <typename T> inline bool Complex<T>::is_nan() const {
    return ( r != r ) || ( i != i );
}

// -oo + 0i -> true

template <typename T> inline bool Complex<T>::is_neg_inf() const {
    return ( r == Support<T>::NEG_INF ) && i == 0.0;
}

// oo + 0i -> true

template <typename T> inline bool Complex<T>::is_pos_inf() const {
    return ( r == Support<T>::POS_INF ) && ( i == 0.0 );
}

// +/-oo + 0i -> true

template <typename T> inline bool Complex<T>::is_real_inf() const {
    return ( r == Support<T>::POS_INF || r == Support<T>::NEG_INF ) && i == 0.0;
}

// a + 0i -> true

template <typename T> inline bool Complex<T>::is_real() const {
    return ( i == 0.0 );
}

// 0 + 0i -> true

template <typename T> inline bool Complex<T>::is_zero() const {
    return ( r == 0.0 ) && ( i == 0.0 );
}

/******************************************
 * Exponential and Logarithmic Functions
 ******************************************/

/******************************************
 * Exponential Function
 *
 * Checked for errors.
 ******************************************/

template <typename T> inline Complex<T> Complex<T>::exp() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::exp( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( std::cos( i ), std::sin( i ) );
    } else {
        T expr = std::exp( r );

        return Complex<T>(
            expr * std::cos( i ),
            expr * std::sin( i )
        );
    }
}

/******************************************
 * Logarithmic Function
 *
 * Checked for errors.
 ******************************************/

template <typename T> inline Complex<T> Complex<T>::log() const {
    if ( i == 0 ) {
        if ( r > 0 ) {
            return Complex<T>( std::log( r ), i );
        } else if ( r == 0 ) {
            return Complex<T>( Support<T>::NEG_INF, Support<T>::NaN );
        } else {
            return Complex<T>( std::log( -r ), Support<T>::sign( i )*Support<T>::PI );
        }
    } else if ( r == 0 ) {
        if ( i > 0 ) {
            return Complex<T>( std::log( i ), Support<T>::PI2 );
        } else {
            return Complex<T>( std::log( -i ), -Support<T>::PI2 );
        }
    } else {
        return Complex<T>(
            0.5 * std::log( r*r + i*i ),
            std::atan2( i, r )
        );
    }
}

// Branch cut:    (-oo, 0]

template <typename T> inline Complex<T> Complex<T>::log10() const {
    {
        T ln10 = std::log( 10.0 );

        return Complex<T>(
            0.5 * std::log( r*r + i*i ) / ln10,
            std::atan2( i, r ) / ln10
        );
    }
}

// (a + bi)^(c + di) -> exp( log(a + bi)*(c + di) )

template <typename T> inline Complex<T> Complex<T>::pow( const Complex<T> & z ) const {
    if ( is_zero() && z.is_zero() ) {
        return Complex<T>( Support<T>::NaN, Support<T>::NaN );
    } else if ( is_zero() ) {
        if ( Support<T>::sign( r ) == 1 ) {
            return Complex<T>( 0.0, 0.0 );
        } else {
            return Complex<T>( Support<T>::POS_INF, Support<T>::POS_INF );
        }
    } else if ( z.is_zero() ) {
            return Complex<T>( 1.0, 0.0 );
    }

    if ( i == 0 ) {
        if ( z.i == 0 ) {
            if ( r > 0 ) {
                return Complex<T>( std::pow( r, z.r ), 0.0 );
            } else {
                std::cout << "B" << std::endl;

                return Complex<T>(
                    std::exp(z.r*std::log(-r))*std::cos(z.r*Support<T>::PI),
                    -std::exp(z.r*std::log(-r))*std::sin(z.r*Support<T>::PI)
                );
            }
        } else if ( z.r == 0 ) {
            if ( r > 0 ) {
                return Complex<T>(
                    std::cos( z.i*std::log( r ) ),
                    std::sin( z.i*std::log( r ) )
                );
            } else {
                std::cout << "D" << std::endl;
                if ( Support<T>::sign( i ) == 1 ) {
                    return Complex<T>(
                        std::exp( -z.i*Support<T>::PI )*std::cos( z.i*std::log( -r ) ),
                        std::exp( -z.i*Support<T>::PI )*std::sin( z.i*std::log( -r ) )
                    );
                } else {
                    return Complex<T>(
                        std::exp( z.i*Support<T>::PI )*std::cos( z.i*std::log( -r ) ),
                        std::exp( z.i*Support<T>::PI )*std::sin( z.i*std::log( -r ) )
                    );
                }
            }
        } else {
            if ( r > 0 ) {
                return Complex<T>(
                    std::exp( z.r*std::log(r) )*std::cos( z.i*std::log(r)),
                    std::exp( z.r*std::log(r) )*std::sin( z.i*std::log(r) )
                );
            } else {
                if ( Support<T>::sign( i ) == 1 ) {
                    return Complex<T>(
                        std::exp( z.r*std::log(-r) - z.i*Support<T>::PI )*
                        std::cos( z.i*std::log(-r) + z.r*Support<T>::PI ),
                        std::exp( z.r*std::log(-r) - z.i*Support<T>::PI )*
                        std::sin( z.i*std::log(-r) + z.r*Support<T>::PI )
                    );
                } else {
                    return Complex<T>(
                        std::pow( -r, z.r )*
                        std::exp(z.i*Support<T>::PI)*std::cos(-z.i*std::log(-r) + z.r*Support<T>::PI),
                        -std::pow( -r, z.r )*
                        std::exp(z.i*Support<T>::PI)*std::sin(-z.i*std::log(-r) + z.r*Support<T>::PI)
                    );
                }
            }
        }
    } else if ( r == 0 ) {
        if ( z.i == 0 ) {
            return Complex<T>(
                std::exp(0.5*z.r*std::log(i*i))*std::cos(0.5*z.r*Support<T>::sign(i)*Support<T>::PI),
                std::exp(0.5*z.r*std::log(i*i))*std::sin(0.5*z.r*Support<T>::sign(i)*Support<T>::PI)
            );
        } else if ( z.r == 0 ) {
            return Complex<T>(
                std::exp(-0.5*z.i*Support<T>::sign(i)*Support<T>::PI)*std::cos(0.5*z.i*std::log(i*i)),
                std::exp(-0.5*z.i*Support<T>::sign(i)*Support<T>::PI)*std::sin(0.5*z.i*std::log(i*i))
            );
        } else {
            return Complex<T>(
                std::exp(0.5*z.r*std::log(i*i) - 0.5*z.i*Support<T>::sign(i)*Support<T>::PI)*std::cos(0.5*z.i*std::log(i*i) + 0.5*z.r*Support<T>::sign(i)*Support<T>::PI),
                std::exp(0.5*z.r*std::log(i*i) - 0.5*z.i*Support<T>::sign(i)*Support<T>::PI)*std::sin(0.5*z.i*std::log(i*i) + 0.5*z.r*Support<T>::sign(i)*Support<T>::PI)
            );
        }
    } else {
        if ( z.i == 0 ) {
            return Complex<T>(
                std::exp(0.5*z.r*std::log(r*r + i*i))*std::cos(z.r*std::atan2(i,r)),
                std::exp(0.5*z.r*std::log(r*r + i*i))*std::sin(z.r*std::atan2(i,r))
            );
        } else if ( z.r == 0 ) {
            return Complex<T>(
                std::exp(-z.i*std::atan2(i,r))*std::cos(0.5*z.i*std::log(r*r + i*i)),
                std::exp(-z.i*std::atan2(i,r))*std::sin(0.5*z.i*std::log(r*r + i*i))
            );
        } else {
            return Complex<T>(
                std::exp(0.5*z.r*std::log(r*r + i*i)-z.i*std::atan2(i,r))*std::cos(0.5*z.i*std::log(r*r + i*i) + z.r*std::atan2(i,r)),
                std::exp(0.5*z.r*std::log(r*r + i*i)-z.i*std::atan2(i,r))*std::sin(0.5*z.i*std::log(r*r + i*i) + z.r*std::atan2(i,r))
            );
        }
    }
}

template <typename T> inline Complex<T> Complex<T>::pow(T x) const {
    if ( i == 0 ) {
        return Complex<T>(
            std::exp( x*std::log( std::fabs(r) ) )*std::cos( x*(0.5 - 0.5*Support<T>::sign(r) )*Support<T>::PI ),
            std::exp( x*std::log( std::fabs(r) ) )*std::sin( x*(0.5 - 0.5*Support<T>::sign(r) )*Support<T>::PI )
        );
    } else if ( r == 0 ) {
        return Complex<T>(
            std::exp( 0.5*x*std::log( i*i ) )*std::cos( 0.5*x*Support<T>::sign(i)*Support<T>::PI ),
            std::exp( 0.5*x*std::log( i*i ) )*std::sin( 0.5*x*Support<T>::sign(i)*Support<T>::PI )
        );
    } else {
        return Complex<T>(
            std::exp( 0.5*x*std::log( r*r + i*i ) )*std::cos( x*std::atan2( i, r ) ),
            std::exp( 0.5*x*std::log( r*r + i*i ) )*std::sin( x*std::atan2( i, r ) )
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::inverse() const {
    if ( is_zero() ) {
        return Complex( Support<T>::sign( r )*Support<T>::POS_INF, -Support<T>::sign( i )*Support<T>::POS_INF );
    } else if ( is_inf() ) {
        return Complex( Support<T>::sign( r )*0.0, -Support<T>::sign( i )*0.0 );
    } else if ( is_nan() ) {
        return Complex( Support<T>::NaN, Support<T>::NaN );
    } else {
        T denom = norm();

        return Complex( r/denom, -i/denom );
    }
}

/******************************************
 * Trigonometric and Hyperbolic Functions
 ******************************************/

// a + 0i -> sin(a) + 0i
// 0 + bi ->     0  + sinh(b)i
// a + bi -> [sin(a) cosh(b)] + [cos(a) sinh(b)]i


template <typename T> inline Complex<T> Complex<T>::sin() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::sin( r ), std::cos( r )*i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, std::sinh( i ) );
    } else {
        return Complex<T>(
            std::sin( r ) * std::cosh( i ),
            std::cos( r ) * std::sinh( i )
        );
    }
}

/******************************************
 * Cosine Function
 *
 * a + 0i -> cos(a) + 0i
 * 0 + bi ->     0  + cosh(b)i
 * a + bi -> [cos(a) cosh(b)] - [sin(a) sinh(b)]i
 *
 * Checked for errors.
 ******************************************/

template <typename T> inline Complex<T> Complex<T>::cos() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::cos( r ), -i*std::sin( r ) );
    } else if ( r == 0.0 ) {
        return Complex<T>( std::cosh( i ), -i*r );
    } else {
        return Complex<T>(
             std::cos( r ) * std::cosh( i ),
            -std::sin( r ) * std::sinh( i )
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::tan() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::tan( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, std::tanh( i ) );
    } else {
        T cosr =  std::cos(r);
        T sinhi = std::sinh(i);

        T denom = cosr*cosr + sinhi*sinhi;

        return Complex<T>(
            std::sin( r )  * cosr  / denom,
            sinhi * std::cosh( i ) / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::sec() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::sec( r ), -i );
    } else if ( r == 0.0 ) {
        return Complex<T>( Support<T>::sech( i ), r*Support<T>::sign( i ) );
    } else {
        T cosr =  std::cos(r);
        T sinhi = std::sinh(i);

        T denom = cosr*cosr + sinhi*sinhi;

        return Complex<T>(
            cosr  * std::cosh( i )  / denom,
            std::sin( r ) * sinhi / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::csc() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::csc( r ), -i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, -Support<T>::csch( i ) );
    } else {
        T sinr =  std::sin(r);
        T sinhi = std::sinh(i);

        T denom = sinr*sinr + sinhi*sinhi;

        return Complex<T>(
             sinr  * std::cosh( i ) / denom,
            -std::cos( r ) * sinhi  / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::cot() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::cot( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, -Support<T>::coth( i ) );
    } else {
        T sinr =  std::sin(r);
        T sinhi = std::sinh(i);

        T denom = sinr*sinr + sinhi*sinhi;

        return Complex<T>(
             sinr  * std::cos( r )  / denom,
            -sinhi * std::cosh( i ) / denom
        );
    }
}

// a + 0i ->  sinh(a) + 0i
// 0 + bi ->      0  + sin(b)i
// a + bi -> [sinh(a)  cos(b)] + [cosh(a) sin(b)]i

template <typename T> inline Complex<T> Complex<T>::sinh() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::sinh( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, std::sin( i ) );
    } else {
        return Complex<T>(
            std::sinh( r ) * std::cos( i ),
            std::cosh( r ) * std::sin( i )
        );
    }
}

// a + 0i ->  cosh(a) + 0i
// 0 + bi ->   cos(b) + 0i
// a + bi -> [cosh(a) cos(b)] + [sinh(a) sin(b)]i

template <typename T> inline Complex<T> Complex<T>::cosh() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::cosh( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( std::cos( i ), -r );
    } else {
        return Complex<T>(
            std::cosh( r ) * std::cos( i ),
            std::sinh( r ) * std::sin( i )
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::tanh() const {
    if ( i == 0.0 ) {
        return Complex<T>( std::tanh( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, std::tan( i ) );
    } else {
        T sinhr = std::sinh(r);
        T cosi =  std::cos(i);

        T denom = sinhr*sinhr + cosi*cosi;

        return Complex<T>(
            sinhr  * std::cosh( r )  / denom,
            std::sin( i ) * cosi / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::sech() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::sech( r ), -i );
    } else if ( r == 0.0 ) {
        return Complex<T>( Support<T>::sec( i ), r );
    } else {
        T sinhr = std::sinh(r);
        T cosi =  std::cos(i);

        T denom = sinhr*sinhr + cosi*cosi;

        return Complex<T>(
             std::cosh( r )  * cosi  / denom,
            -sinhr * std::sin( i ) / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::csch() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::csch( r ), -i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, -Support<T>::csc( i ) );
    } else {
        T sinhr =  std::sinh(r);
        T sini = std::sin(i);

        T denom = sinhr*sinhr + sini*sini;

        return Complex<T>(
             sinhr  * std::cos( i )  / denom,
            -std::cosh( r ) * sini / denom
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::coth() const {
    if ( i == 0.0 ) {
        return Complex<T>( Support<T>::coth( r ), i );
    } else if ( r == 0.0 ) {
        return Complex<T>( r, -Support<T>::cot( i ) );
    } else {
        T sinhr =  std::sinh(r);
        T sini = std::sin(i);

        T denom = sinhr*sinhr + sini*sini;

        return Complex<T>(
             sinhr  * std::cosh( r )  / denom,
            -std::cos( i ) * sini / denom
        );
    }
}

// Branch cut:    (-oo, -1) U (1, oo)

template <typename T> inline Complex<T> Complex<T>::asin() const {
    if ( i == 0 ) {
        if ( r >= -1 && r <= 1 ) {
            return Complex<T>( std::asin( r ), i );
        } else {
            if ( r > 1 ) {
                return Complex<T>(
                    Support<T>::PI2,
                    Support<T>::sign( i )*std::log(
                        r + std::sqrt(r*r  - 1)
                    )
                );
            } else {
                return Complex<T>(
                    -Support<T>::PI2,
                    Support<T>::sign( i )*std::log(
                        -r + std::sqrt(r*r  - 1)
                    )
                );
            }
        }
    } else if ( r == 0 ) {
        return Complex<T>( r, std::log( i + std::sqrt( i*i + 1 ) ) );
    } else {
        T ss = r*r + i*i + 1;
        T ssp2r = std::sqrt( ss + 2*r );
        T ssm2r = std::sqrt( ss - 2*r );
        T sum = 0.5*( ssp2r + ssm2r );

        return Complex<T>(
            std::asin( ssp2r/2 - ssm2r/2 ),
            Complex<T>( i, -r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            )
        );
    }
}

// Branch cut:    (-oo, -1) U (1, oo)

template <typename T> inline Complex<T> Complex<T>::acos() const {
    if ( i == 0 ) {
        if ( r >= -1 && r <= 1 ) {
            return Complex<T>( std::acos( r ), -i );
        } else if ( r > 1 ) {
            return Complex<T>(
                0.0,
                -Support<T>::sign( i )*std::log(
                    r + std::sqrt(r*r  - 1)
                )
            );
        } else {
            return Complex<T>(
                Support<T>::PI,
                -Support<T>::sign( i )*std::log(
                    -r + std::sqrt(r*r  - 1)
                )
            );
        }
    } else {
        T ss = r*r + i*i + 1;
        T ssp2r = std::sqrt( ss + 2*r );
        T ssm2r = std::sqrt( ss - 2*r );
        T sum = 0.5*( ssp2r + ssm2r );

        return Complex<T>(
            std::acos( 0.5*(ssp2r - ssm2r) ),
            -Complex<T>( i, -r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            )
        );
    }
}

/**********************************************************
 * Arc Tangent Function
 *
 * Branch cut:    (-ooi, -i] U [i, ooi)
 *
 * There do not appear to be any errors in this function
 * for finite arguments.
 * Should simplify for r == 0, however.
 **********************************************************/

template <typename T> inline Complex<T> Complex<T>::atan() const {
    if ( r == 0 ) {
        if ( i < -1 ) {
            T ip = i + 1;
            T im = i - 1;

            return Complex<T>(
                Support<T>::sign(r)*Support<T>::PI2,
                0.25*std::log( ip*ip/(im*im) )
            );
        } else if ( i == -1 ) {
            return Complex<T>( Support<T>::NaN, Support<T>::NEG_INF );
        } else if ( i == 0 ) {
            return *this;
        } else if ( i < 1 ) {
            T ip = i + 1;
            T im = i - 1;

            return Complex<T>(
                r,
                0.25*std::log( ip*Support<T>::PI/(im*im) )
            );
        } else if ( i == 1 ) {
            return Complex<T>( Support<T>::NaN, Support<T>::POS_INF );
        } else {
            T ip = i + 1;
            T im = i - 1;

            return Complex<T>(
                Support<T>::sign(r)*Support<T>::PI2,
                0.25*std::log( ip*ip/(im*im) )
            );
        }
    } else if ( i == 0 ) {
        return Complex<T>( std::atan( r ), i );
    } else {
      T opi = 1.0 + i;
      T omi = 1.0 - i;
        T rr = r*r;

        return Complex<T>(
            0.5*( std::atan2( r, omi ) - std::atan2( -r, opi ) ),
            0.25*std::log( (rr + opi*opi)/(rr + omi*omi) )
        );
    }
}

/**********************************************************
 * Arc Secant Function
 *
 * Branch cut:    (-1, 1)
 *
 * There do not appear to be any errors in this function
 * for finite arguments.
 **********************************************************/

template <typename T> inline Complex<T> Complex<T>::asec() const {
    if ( i == 0 ) {
        if ( ( r <= -1 ) || ( r >= 1 ) ) {
            return Complex<T>( std::acos( 1/r ), i );
        } else if ( r == 0 ) {
            return Complex<T>( Support<T>::POS_INF, Support<T>::POS_INF );
        } else {
            if ( r < 0 ) {
                return Complex<T>(
                    Support<T>::PI,
                    Support<T>::sign( i )*std::log(
                        -1/r + std::sqrt(1/(r*r)  - 1)
                    )
                );
            } else {
                return Complex<T>(
                    0.0,
                    Support<T>::sign( i )*std::log(
                        1/r + std::sqrt(1/(r*r)  - 1)
                    )
                );
            }
        }
    } else {
        T ss = r*r + i*i;
        T p1 = r/ss;

        T p1p1 = p1 + 1;
        p1p1 = p1p1 * p1p1;

        T p1m1 = p1 - 1;
        p1m1 = p1m1 * p1m1;

        T p2 = i/ss;
        p2 = p2*p2;

        T p1p1p2 = std::sqrt( p1p1 + p2 );
        T p1m1p2 = std::sqrt( p1m1 + p2 );

        T sum = 0.5*( p1p1p2 + p1m1p2 );

        return Complex<T>(
            std::acos( 0.5*(p1p1p2 - p1m1p2 ) ),
            Complex<T>( i, r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            )
        );
    }
}

/**********************************************************
 * Arc Cosecant Function
 *
 * Branch cut:    (-1, 1)
 *
 * There do not appear to be any errors in this function
 * for finite arguments.
 **********************************************************/

template <typename T> inline Complex<T> Complex<T>::acsc() const {
    if ( i == 0 ) {
        if ( ( r <= -1 ) || ( r >= 1 ) ) {
            // (-inf, 1] U [1, inf)

            return Complex<T>( std::asin( 1/r ), -i );
        } else if ( r == 0 ) {
            // 0

            return Complex<T>( Support<T>::POS_INF, Support<T>::POS_INF );
        } else if ( r < 0 ) {
            // (-1, 0)

            return Complex<T>(
                -Support<T>::PI2,
                -Support<T>::sign( i )*std::log(
                    -1/r + std::sqrt(1/(r*r)  - 1)
                )
            );
        } else {
            // (0, 1)

            return Complex<T>(
                Support<T>::PI2,
                -Support<T>::sign( i )*std::log(
                    1/r + std::sqrt(1/(r*r)  - 1)
                )
            );
        }
    } else if ( r == 0 ) {
        return Complex<T>(
            r,
            -std::log( 1/i + std::sqrt( 1 + 1/(i*i) ) )
        );
    } else {
        T ss = r*r + i*i;
        T p1 = r/ss;

        T p1p1 = p1 + 1;
        p1p1 = p1p1 * p1p1;

        T p1m1 = p1 - 1;
        p1m1 = p1m1 * p1m1;

        T p2 = i/ss;
        p2 = p2*p2;

        T p1p1p2 = std::sqrt( p1p1 + p2 );
        T p1m1p2 = std::sqrt( p1m1 + p2 );

        T sum = 0.5*( p1p1p2 + p1m1p2 );

        return Complex<T>(
            std::asin( 0.5*(p1p1p2 - p1m1p2) ),
            -Complex<T>( i, r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            )
        );
    }
}

// Branch cut:    (-i, i)

template <typename T> inline Complex<T> Complex<T>::acot() const {
    if ( i == 0 ) {
        return Complex<T>( Support<T>::PI2 - std::atan( r ), -i );
    } else if ( r == 0 ) {
        if ( i < -1 ) {
            // (-inf i, -i)

            return Complex<T>(
                Support<T>::sign( r ) == 1 ? 0.0 : Support<T>::PI,
                0.5*std::log( (-i + 1)/(-i - 1) )
            );
        } else if ( i == -1 ) {
            return Complex<T>( Support<T>::NaN, Support<T>::POS_INF );
        } else if ( i < 1 ) {
            // (-i, i)

            return Complex<T>(
                Support<T>::PI2,
                -0.25*std::log( (i + 1)*(i + 1)/((i - 1)*(i - 1)) )
            );
        } else if ( i == 1 ) {
            return Complex<T>( Support<T>::NaN, Support<T>::NEG_INF );
        } else {
            // (i, inf i)

            return Complex<T>(
                Support<T>::sign( r ) == 1 ? 0.0 : Support<T>::PI,
                0.5*std::log( (-i + 1)/(-i - 1) )
            );
        }
    } else {
      T opi = 1.0 + i;
      T omi = 1.0 - i;
        T rr = r*r;

        return Complex<T>(
            Support<T>::PI2 + 0.5*( std::atan2( -r, opi ) - std::atan2( r, omi ) ),
            -0.25*std::log( (rr + opi*opi)/(rr + omi*omi) )
        );
    }
}

/**********************************************************
 * Arc Hyperbolic Sine Function
 *
 * Branch cut:    (-ooi, -i) U (i, ooi)
 *
 * There do not appear to be any errors in this function
 * for finite arguments.
 **********************************************************/

template <typename T> inline Complex<T> Complex<T>::asinh() const {
    if ( r == 0 ) {
        if ( ( i >= -1 ) && ( i <= 1 ) ) {
            // [-i, i]
            return Complex<T>( r, std::asin( i ) );
        } else if ( i > 1 ) {
            // (i, ooi)

            return Complex<T>(
                Support<T>::sign( r )*std::log( i + std::sqrt( i*i - 1 ) ),
                Support<T>::PI2
            );
        } else {
            // (-ooi, -i)

            return Complex<T>(
                Support<T>::sign( r )*std::log( -i + std::sqrt( i*i - 1 ) ),
                -Support<T>::PI2
            );
        }
    } else if ( i == 0 ) {
        return Complex<T>( std::log( r + std::sqrt( r*r + 1 ) ), i );
    } else {
        T ss = r*r + i*i + 1;
        T ssp2i = std::sqrt( ss + 2*i );
        T ssm2i = std::sqrt( ss - 2*i );
        T sum = 0.5*( ssp2i + ssm2i );

        return Complex<T>(
            Complex<T>( r, i ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            ),
            std::asin( ssp2i/2 - ssm2i/2 )
        );
    }
}

// Branch cut:    (-oo, 1)

template <typename T> inline Complex<T> Complex<T>::acosh() const {
    if ( i == 0 ) {
        if ( r > 1 ) {
            // (1, oo)

            return Complex<T>( std::log( r + std::sqrt( (r - 1)*(r + 1) ) ), i );
        } else if ( r == 1 ) {
            // 1

            return Complex<T>( 0.0, i );
        } else if ( r == 0 ) {
            // 0

            return Complex<T>( 0.0, Support<T>::sign( i )*Support<T>::PI2 );
        } else if ( r > -1 ) {
            // (-1, 0) U (0, 1)

            return Complex<T>( 0.0, Support<T>::sign( i )*std::acos( r ) );
        } else if ( r == -1 ) {
            // -1

            return Complex<T>( 0.0, Support<T>::sign( i )*Support<T>::PI );
        } else {
            // (-oo, -1)
            return Complex<T>( std::log( -r + std::sqrt( (-1 - r)*(1 - r) ) ), Support<T>::sign( i )*Support<T>::PI );
        }
    } else if ( r == 0 ) {
        return Complex<T>( std::log( std::sqrt( i*i + 1 ) + Support<T>::sign( i )*i ), Support<T>::sign( i )*Support<T>::PI2 );
    } else {
        T ss = r*r + i*i + 1;
        T ssp2r = std::sqrt( ss + 2*r );
        T ssm2r = std::sqrt( ss - 2*r );
        T sum = 0.5*( ssp2r + ssm2r );

        return Complex<T>(
            Complex<T>( -i, r - 1.0 ).csgn()*
            Complex<T>( -i, r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            ),
            -Complex<T>( -i, r - 1.0 ).csgn()*std::acos( ssp2r/2 - ssm2r/2 )
        );
    }
}

// Branch cut:    (-oo, -1] U [1, oo)

template <typename T> inline Complex<T> Complex<T>::atanh() const {
    if ( i == 0 ) {
        if ( r > 1 ) {
            // [1, oo)

            return Complex<T>(
                0.25*std::log( (r + 1)*(r + 1)/((r - 1)*(r - 1)) ),
                Support<T>::sign( i )*Support<T>::PI2
            );
        } else if ( r == 1 ) {
            return Complex<T>( Support<T>::POS_INF, Support<T>::NaN );
        } else if ( r > 0 ) {
            return Complex<T>(
                0.5*std::log( (1 + r)/(1 - r) ),
                i
            );
        } else if ( r == 0 ) {
            return *this;
        } else if ( r > -1 ) {
            return Complex<T>(
                -0.5*std::log( (1 - r)/(1 + r) ),
                i
            );
        } else if ( r == -1 ) {
            return Complex<T>( Support<T>::NEG_INF, Support<T>::NaN );
        } else {
            T rp = r + 1;
            T rm = r - 1;
            rp *= rp;
            rm *= rm;
            T logrprm = 0.25*std::log( rp/rm );

            return Complex<T>( logrprm, Support<T>::sign( i )*Support<T>::PI2 );
        }
    } else if ( r == 0 ) {
        return Complex<T>( r, std::atan( i ) );
    } else {
      T opr = 1.0 + r;
      T omr = 1.0 - r;
        T ii = i*i;

        return Complex<T>(
            0.25*std::log( (ii + opr*opr)/(ii + omr*omr) ),
            0.5*( std::atan2( i, opr ) - std::atan2( -i, omr ) )
        );
    }
}

// Branch cut:    (-oo, 0] U (1, oo)

template <typename T> inline Complex<T> Complex<T>::asech() const {
    if ( i == 0 ) {
        if ( r < -1 ) {
            return Complex<T>( 0.0, -Support<T>::sign( i )*std::acos( 1/r ) );
        } else if ( r == -1 ) {
            return Complex<T>( 0.0, -Support<T>::sign( i )*Support<T>::PI );
        } else if ( r < 0 ) {
            return Complex<T>(
                std::log( -1/r + std::sqrt( -(1 - 1/r)*(1 + 1/r) ) ),
                -Support<T>::sign( i )*Support<T>::PI
            );
        } else if ( r == 0 ) {
            return Complex<T>( Support<T>::POS_INF, Support<T>::NaN );
        } else if ( r <= 1 ) {
            return Complex<T>(
                std::log( 1/r + std::sqrt( (1/r + 1)*(1/r - 1) ) ),
                -Support<T>::sign( i )*0.0
            );
        } else {
            return Complex<T>( 0.0, -Support<T>::sign( i )*std::acos( 1/r ) );
        }
    } else if ( r == 0 ) {
        if ( i > 0 ) {
            return Complex<T>(
                std::log( std::sqrt( 1 + 1/(i*i) ) + 1/i ), -Support<T>::PI2
            );
        } else {
            return Complex<T>(
                std::log( std::sqrt( 1 + 1/(i*i) ) - 1/i ), Support<T>::PI2
            );
        }
    } else {
        T ss = r*r + i*i;
        T Rr = r/ss;

        T Rrp1 = Rr + 1;
        Rrp1 = Rrp1 * Rrp1;

        T Rrm1 = Rr - 1;
        Rrm1 = Rrm1 * Rrm1;

        T Ri = i/ss;
        Ri = Ri*Ri;

        T Rrp1pRi = std::sqrt( Rrp1 + Ri );
        T Rrm1pRi = std::sqrt( Rrm1 + Ri );

        T sum = 0.5*( Rrp1pRi + Rrm1pRi );

        return Complex<T>(
            -Complex<T>( -i, ss - r ).csgn()*Complex<T>( i, r ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            ),
            Complex<T>( -i, ss - r ).csgn()*(
                std::acos( 0.5*(Rrp1pRi - Rrm1pRi) )
            )
        );
    }
}

// Branch cut:    (-i, i)
// NOT YET FINISHED

template <typename T> inline Complex<T> Complex<T>::acsch() const {
    if ( r == 0 ) {
        if ( i <= -1 ) {
            return Complex<T>( r, -std::asin( 1/i ) );
        } else if ( i == -1 ) {
            return Complex<T>( r, -Support<T>::PI2 );
        } else if ( i < 0 ) {
            return Complex<T>(
                Support<T>::sign( r )*std::log( -1/i + std::sqrt( 1/(i*i) - 1 ) ),
                Support<T>::PI2
            );
        } else if ( i == 0 ) {
            return Complex<T>( Support<T>::POS_INF, Support<T>::NEG_INF );
        } else if ( i < 1 ) {
            return Complex<T>(
                Support<T>::sign( r )*std::log( 1/i + std::sqrt( 1/(i*i) - 1 ) ),
                -Support<T>::PI2
            );
        } else if ( i == 1 ) {
            return Complex<T>( r, -Support<T>::PI2 );
        } else {
            return Complex<T>( r, -std::asin( 1/i ) );
        }
    } else if ( i == 0 ) {
        return Complex<T>(
            std::log( 1/r + std::sqrt( 1/(r*r) + 1 ) ),
            -i
        );
    } else {
        T ss = r*r + i*i;
        T p1 = i/ss;

        T p1p1 = p1 + 1;
        p1p1 = p1p1 * p1p1;

        T p1m1 = p1 - 1;
        p1m1 = p1m1 * p1m1;

        T p2 = r/ss;
        p2 = p2*p2;

        T p1p1p2 = std::sqrt( p1p1 + p2 );
        T p1m1p2 = std::sqrt( p1m1 + p2 );

        T sum = 0.5*( p1p1p2 + p1m1p2 );

        return Complex<T>(
            Complex<T>( r, -i ).csgn()*std::log(
                sum + std::sqrt(sum*sum  - 1)
            ),
            std::asin( 0.5*(p1m1p2 - p1p1p2) )
        );
    }
}

// Branch cut:    [-1, 1]
// NOT YET FINISHED

template <typename T> inline Complex<T> Complex<T>::acoth() const {
    if ( i == 0 ) {
        if ( r < -1 ) {
            return Complex<T>(
                0.5*std::log( (r + 1)/(r - 1) ),
                -i
            );
        } else if ( r == -1 ) {
            return Complex<T>( Support<T>::NEG_INF, Support<T>::NaN );
        } else if ( r < 0 ) {
            return Complex<T>(
                0.5*std::log( (r + 1)/(1 - r) ),
                -Support<T>::sign( i )*Support<T>::PI2
            );
        } else if ( r == 0 ) {
            return Complex( r, -Support<T>::sign( i )*Support<T>::PI2 );
        } else if ( r < 1 ) {
            return Complex<T>(
                0.5*std::log( (r + 1)/(1 - r) ),
                -Support<T>::sign( i )*Support<T>::PI2
            );
        } else if ( r == 1 ) {
            return Complex<T>( Support<T>::POS_INF, Support<T>::NaN );
        } else {
            return Complex<T>(
                0.5*std::log( (r + 1)/(r - 1) ),
                -i
            );
        }
    } else if ( r == 0 ) {
        if ( i > 0 ) {
            return Complex<T>( r, -Support<T>::PI2 + std::atan( i ) );
        } else {
            return Complex<T>( r, Support<T>::PI2 + std::atan( i ) );
        }
    } else {
      T rp1 = r + 1.0;
      T rm1 = r - 1.0;
        T ii = i*i;

        return Complex<T>(
            0.25*std::log( (ii + rp1*rp1)/(ii + rm1*rm1) ),
            0.5*( std::atan2( i, rp1 ) - std::atan2( i, rm1 ) )
        );
    }
}

template <typename T> inline Complex<T> Complex<T>::bessel_J( int n ) const {
    if ( n < 0 ) {
        return (n & 1) ? -bessel_J( -n ) : bessel_J( -n );
    }

    Complex<T> z = *this / 2.0;

    Complex<T> term = 1.0;

    for ( int i = 1; i <= n; ++i ) {
        term *= z/i;
    }

    Complex<T> z2 = z.sqr();

    Complex<T> result = term;
    Complex<T> result_old = term;

    for ( T i = 1; ; i += 2.0 ) {
        term *= z2 / (i*(n + i));
        result -= term;

        term *= z2 / ((i + 1.0)*(n + i + 1.0));
        result += term;

        if ( result_old == result ) {
            return result;
        }

        result_old = result;
    }
}

/******************************************
 * Integer Functions
 ******************************************/

// (a + bi) -> ceil(a) + ceil(b)i

template <typename T> inline Complex<T> Complex<T>::ceil() const {
    return Complex<T>(
        std::ceil(r), std::ceil(i)
    );
}

// (a + bi) -> floor(a) + floor(b)i

template <typename T> inline Complex<T> Complex<T>::floor() const {
    return Complex<T>(
        std::floor(r), std::floor(i)
    );
}


/******************************************
 * Horner's Rule
 *
 *   The polynomial is defined by giving the highest
 *   coefficient first:
 *
 *            n - 1         n - 2
 *      v[0]*z      + v[1]*z      + ... + v[n-2]*z + v[n-1]
 *
 *   This is the same as with Matlab.
 ******************************************/

template <typename T> inline Complex<T> Complex<T>::horner( T * v, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = *v;

        for ( unsigned int j = 1; j < n; ++j ) {
            ar = ar*r + *(++v);
        }

        return Complex<T>( ar, 0 );
    } else if ( r == 0 ) {              // imaginary
        T ar = *v;
        T ai = 0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T aR = ar;

            ar = -ai*i + *(++v);
            ai = aR*i;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = *v;
        T ai = 0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T aR = ar;

            ar = ar*r - ai*i + *(++v);
            ai = aR*i + ai*r;
        }

        return Complex<T>( ar, ai );
    }
}

template <typename T> inline Complex<T> Complex<T>::horner( T * v, T * c, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = *v;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - *(++c);

            ar = ar*dr + *(++v);
        }

        return Complex<T>( ar, 0 );
    } else if ( r == 0 ) {              // imaginary
        T ar = *v;
        T ai = 0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = -*(++c);
            T aR = ar;

            ar = ar*dr - ai*i + *(++v);
            ai = aR*i + ai*dr;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = *v;
        T ai = 0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - *(++c);
            T aR = ar;

            ar = ar*dr - ai*i + *(++v);
            ai = aR*i + ai*dr;
        }

        return Complex<T>( ar, ai );
    }

}

template <typename T> inline Complex<T> Complex<T>::horner( Complex<T> * v, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            ar = ar*r + (*(++v)).r;
            ai = ai*r + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else if ( r == 0 ) {              // imaginary
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T aR = ar;

            ar = -ai*i + (*(++v)).r;
            ai = aR*i + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T aR = ar;

            ar = ar*r - ai*i + (*(++v)).r;
            ai = aR*i + ai*r + (*v).i;
        }

        return Complex<T>( ar, ai );
    }
}

template <typename T> inline Complex<T> Complex<T>::horner( Complex<T> * v, T * c, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - *(++c);

            ar = ar*dr + (*(++v)).r;
            ai = ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else if ( r == 0 ) {              // imaginary
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = -*(++c);
            T aR = ar;

            ar = ar*dr - ai*i + (*(++v)).r;
            ai = aR*i + ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - *(++c);
            T aR = ar;

            ar = ar*dr - ai*i + (*(++v)).r;
            ai = aR*i + ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    }
}

template <typename T> inline Complex<T> Complex<T>::horner( T * v, Complex<T> * c, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = *v;
        T ai = 0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - (*(++c)).r;
            T di = -(*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + *(++v);
            ai = aR*di + ai*dr;
        }

        return Complex<T>( ar, ai );
    } else if ( r == 0 ) {              // imaginary
        T ar = *v;
        T ai = 0.0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = -(*(++c)).r;
            T di = i - (*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + *(++v);
            ai = aR*di + ai*dr;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = *v;
        T ai = 0.0;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - (*(++c)).r;
            T di = i - (*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + *(++v);
            ai = aR*di + ai*dr;
        }

        return Complex<T>( ar, ai );
    }
}

template <typename T> inline Complex<T> Complex<T>::horner( Complex<T> * v, Complex<T> * c, unsigned int n ) const {
    if ( n == 0 ) {
        return ZERO;
    }

    if ( i == 0 ) {                     // real
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - (*(++c)).r;
            T di = -(*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + (*(++v)).r;
            ai = aR*di + ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else if ( r == 0 ) {              // imaginary
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = -(*(++c)).r;
            T di = i - (*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + (*(++v)).r;
            ai = aR*di + ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    } else {                            // complex
        T ar = (*v).r;
        T ai = (*v).i;

        for ( unsigned int j = 1; j < n; ++j ) {
            T dr = r - (*(++c)).r;
            T di = i - (*c).i;
            T aR = ar;

            ar = ar*dr - ai*di + (*(++v)).r;
            ai = aR*di + ai*dr + (*v).i;
        }

        return Complex<T>( ar, ai );
    }
}

/******************************************
 * Random Factories
 ******************************************/

// -> a + bi

template <typename T> inline Complex<T> Complex<T>::random() {
    return Complex<T>(
        (static_cast<T>( rand() ))/RAND_MAX,
        (static_cast<T>( rand() ))/RAND_MAX
    );
}

// -> 0 + bi

template <typename T> inline Complex<T> Complex<T>::random_imag() {
    return Complex<T>( 0.0, (static_cast<T>( rand() ))/RAND_MAX );
}

// -> a + 0i

template <typename T> inline Complex<T> Complex<T>::random_real() {
    return Complex<T>( (static_cast<T>( rand() ))/RAND_MAX, 0.0 );
}

/******************************************
 * Binary Arithmetic Operators
 ******************************************/

// (a + bi) + (c + di) -> (a + c) + (b + d)i

template <typename T> inline Complex<T> Complex<T>::operator + ( const Complex<T> & z ) const {
    return Complex<T>( r + z.r, i + z.i );
}

// (a + bi) + x -> (a + x) + bi

template <typename T> inline Complex<T> Complex<T>::operator + ( T x ) const {
    return Complex<T>( r + x, i );
}

// x + (a + bi) -> (a + x) + bi

template <typename T> Complex<T> operator + ( T x, const Complex<T> & z ) {
    return Complex<T>( x + z.real(), z.imag_i() );
}

template <typename T> Complex<T> operator + ( long x, const Complex<T> & z ) {
    return Complex<T>( static_cast<T>( x ) + z.real(), z.imag_i() );
}

// (a + bi) - (c + di) -> (a - c) + (b - d)i

template <typename T> inline Complex<T> Complex<T>::operator - ( const Complex<T> & z ) const {
    return Complex<T>( r - z.r, i - z.i );
}

// (a + bi) - x -> (a - x) + bi

template <typename T> inline Complex<T> Complex<T>::operator - ( T x ) const {
    return Complex<T>( r - x, i );
}

// x - (a + bi) -> (x + a) - bi

template <typename T> Complex<T> operator - ( T x, const Complex<T> & z ) {
    return Complex<T>( x - z.real(), -z.imag_i() );
}

template <typename T> Complex<T> operator - ( long x, const Complex<T> & z ) {
    return Complex<T>( static_cast<T>( x ) - z.real(), -z.imag_i() );
}

template <typename T> inline Complex<T> Complex<T>::operator * ( const Complex<T> & z ) const {
    return Complex<T>( r*z.r - i*z.i, r*z.i + i*z.r );
}

template <typename T> inline Complex<T> Complex<T>::operator * ( T x ) const {
    if ( Support<T>::is_inf( x ) && norm() > 0 ) {
        return Complex<T>(
            ( (  r == 0 )?Support<T>::sign(x):x )*r,
            ( (  i == 0 )?Support<T>::sign(x):x )*i
        );
    } else {
        return Complex<T>( x*r, x*i );
    }
}

template <typename T> Complex<T> operator * ( T x, const Complex<T> & z ) {
    return Complex<T>( x*z.real(), x*z.imag_i() );
}

template <typename T> Complex<T> operator * ( long x, const Complex<T> & z ) {
    return Complex<T>( static_cast<T>( x )*z.real(), static_cast<T>( x )*z.imag_i() );
}

template <typename T> inline Complex<T> Complex<T>::operator / ( const Complex<T> & z ) const {
    T denom = z.norm();
    return Complex<T>( (r*z.r + i*z.i)/denom, (-r*z.i + i*z.r)/denom );
}

template <typename T> inline Complex<T> Complex<T>::operator / ( T x ) const {
    if ( ( x == 0.0 ) && ( norm() > 0 ) ) {
        return Complex<T>(
            r / ( ( r == 0 )?Support<T>::sign( x ):x ),
            i / ( ( i == 0 )?Support<T>::sign( x ):x )
        );
    } else {
        return Complex<T>( r/x, i/x );
    }
}

template <typename T> Complex<T> operator / ( T x, const Complex<T> & z ) {
    T mltplr = x/z.norm();
    return Complex<T>( mltplr*z.real(), -mltplr*z.imag_i() );
}

template <typename T> Complex<T> operator / ( long x, const Complex<T> & z ) {
    T mltplr = static_cast<T>( x )/z.norm();
    return Complex<T>( mltplr*z.real(), -mltplr*z.imag_i() );
}

/******************************************
 * Unary Arithmetic Operators
 ******************************************/

// -(a + bi) -> -a - bi

template <typename T> inline Complex<T> Complex<T>::operator - () const {
    return Complex<T>( -r, -i );
}

/******************************************
 * Binary Boolean Operators
 ******************************************/

// (a + bi) == (c + di) -> (a == c) && (b == d)

template <typename T> inline bool Complex<T>::operator == ( const Complex<T> & z ) const {
    return ( r == z.r ) && ( i == z.i );
}

// (a + bi) == x -> (a == x) && (b == 0)

template <typename T> inline bool Complex<T>::operator == ( T x ) const {
    return ( r == x ) && ( i == 0.0 );
}

// (a + bi) != (c + di) -> (a != c) || (b != d)

template <typename T> inline bool Complex<T>::operator != ( const Complex<T> & z ) const {
    return ( r != z.r ) || ( i != z.i );
}

// (a + bi) != x -> (a != x) || (b != 0)

template <typename T> inline bool Complex<T>::operator != ( T x ) const {
    return ( r != x ) || ( i != 0.0 );
}

// x == (a + bi) -> (x == a) && (b == 0)

template <typename T> bool operator == ( T x, const Complex<T> & z ) {
    return ( z.real() == x ) && ( z.imag_i() == 0.0 );
}

template <typename T> bool operator == ( long x, const Complex<T> & z ) {
    return ( z.real() == static_cast<T>( x ) ) && ( z.imag_i() == 0.0 );
}

// x != (a + bi) -> (x != a) || (b != 0)

template <typename T> bool operator != ( T x, const Complex<T> & z ) {
    return ( z.real() != x ) || ( z.imag_i() != 0.0 );
}

template <typename T> bool operator != ( long x, const Complex<T> & z ) {
    return ( z.real() != static_cast<T>( x ) ) || ( z.imag_i() != 0.0 );
}

template<typename T> std::ostream & operator << ( std::ostream & out, const Complex<T> & z ) {
    Support<T>::print_real( z.real(), out );
    Support<T>::print_imaginary( z.imag_i(), Complex<T>::get_symbol(), out );

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

template <typename T> T real( const Complex<T> & z ) {
    return z.real();
}

template <typename T> T imag_i( const Complex<T> & z ) {
    return z.imag_i();
}

template <typename T> T csgn( const Complex<T> & z ) {
    return z.csgn();
}

template <typename T> T abs( const Complex<T> & z ) {
    return z.abs();
}

template <typename T> T norm( const Complex<T> & z ) {
    return z.norm();
}

template <typename T> T abs_imag( const Complex<T> & z ) {
    return z.abs_imag();
}

template <typename T> T norm_imag( const Complex<T> & z ) {
    return z.norm_imag();
}

/******************************************
 * Complex-valued Functions
 ******************************************/

template <typename T> Complex<T> imag( const Complex<T> & z ) {
    return z.imag();
}

template <typename T> Complex<T> conj( const Complex<T> & z ) {
    return z.conj();
}

template <typename T> Complex<T> signum( const Complex<T> & z ) {
    return z.signum();
}

template <typename T> Complex<T> sqr( const Complex<T> & z ) {
    return z.sqr();
}

template <typename T> Complex<T> sqrt( const Complex<T> & z ) {
    return z.sqrt();
}

/******************************************
 * Exponential and Logarithmic Functions
 ******************************************/

template <typename T> Complex<T> exp( const Complex<T> & z ) {
    return z.exp();
}

template <typename T> Complex<T> log( const Complex<T> & z ) {
    return z.log();
}

template <typename T> Complex<T> log10( const Complex<T> & z ) {
    return z.log10();
}

template <typename T> Complex<T> pow( const Complex<T> & z, const Complex<T> & w ) {
    return z.pow( w );
}

template <typename T> Complex<T> pow( const Complex<T> & z, T x ) {
    return z.pow( x );
}

template <typename T> Complex<T> inverse( const Complex<T> & z ) {
    return z.inverse();
}

/******************************************
 * Trigonometric and Hyperbolic Functions
 ******************************************/

template <typename T> Complex<T> sin( const Complex<T> & z ) {
    return z.sin();
}

template <typename T> Complex<T> cos( const Complex<T> & z ) {
    return z.cos();
}

template <typename T> Complex<T> tan( const Complex<T> & z ) {
    return z.tan();
}

template <typename T> Complex<T> sec( const Complex<T> & z ) {
    return z.sec();
}

template <typename T> Complex<T> csc( const Complex<T> & z ) {
    return z.csc();
}

template <typename T> Complex<T> cot( const Complex<T> & z ) {
    return z.cot();
}

template <typename T> Complex<T> sinh( const Complex<T> & z ) {
    return z.sinh();
}

template <typename T> Complex<T> cosh( const Complex<T> & z ) {
    return z.cosh();
}

template <typename T> Complex<T> tanh( const Complex<T> & z ) {
    return z.tanh();
}

template <typename T> Complex<T> sech( const Complex<T> & z ) {
    return z.sech();
}

template <typename T> Complex<T> csch( const Complex<T> & z ) {
    return z.csch();
}

template <typename T> Complex<T> coth( const Complex<T> & z ) {
    return z.coth();
}

template <typename T> Complex<T> asin( const Complex<T> & z ) {
    return z.asin();
}

template <typename T> Complex<T> acos( const Complex<T> & z ) {
    return z.acos();
}

template <typename T> Complex<T> atan( const Complex<T> & z ) {
    return z.atan();
}

template <typename T> Complex<T> asec( const Complex<T> & z ) {
    return z.asec();
}

template <typename T> Complex<T> acsc( const Complex<T> & z ) {
    return z.acsc();
}

template <typename T> Complex<T> acot( const Complex<T> & z ) {
    return z.acot();
}

template <typename T> Complex<T> asinh( const Complex<T> & z ) {
    return z.asinh();
}

template <typename T> Complex<T> acosh( const Complex<T> & z ) {
    return z.acosh();
}

template <typename T> Complex<T> atanh( const Complex<T> & z ) {
    return z.atanh();
}

template <typename T> Complex<T> asech( const Complex<T> & z ) {
    return z.asech();
}

template <typename T> Complex<T> acsch( const Complex<T> & z ) {
    return z.acsch();
}

template <typename T> Complex<T> acoth( const Complex<T> & z ) {
    return z.acoth();
}

template <typename T> Complex<T> bessel_J( int n, const Complex<T> & z ) {
    return z.bessel_J( n );
}

/******************************************
 * Integer Functions
 ******************************************/

template <typename T> Complex<T> floor( const Complex<T> & z ) {
    return z.floor();
}

template <typename T> Complex<T> ceil( const Complex<T> & z ) {
    return z.ceil();
}

/******************************************
 * Horner's Rule
 ******************************************/

template <typename T> Complex<T> horner( const Complex<T> & z, T * v, unsigned int n ) {
    return z.horner( v, n );
}

template <typename T> Complex<T> horner( const Complex<T> & z, T * v, T * c, unsigned int n ) {
    return z.horner( v, c, n );
}

template <typename T> Complex<T> horner( const Complex<T> & z, Complex<T> * v, unsigned int n ) {
    return z.horner( v, n );
}

template <typename T> Complex<T> horner( const Complex<T> & z, Complex<T> * v, T * c, unsigned int n ) {
    return z.horner( v, c, n );
}

template <typename T> Complex<T> horner( const Complex<T> & z, T * v, Complex<T> * c, unsigned int n ) {
    return z.horner( v, c, n );
}

template <typename T> Complex<T> horner( const Complex<T> & z, Complex<T> * v, Complex<T> * c, unsigned int n ) {
    return z.horner( v, c, n );
}

/**************************************************
 * ********************************************** *
 * *                                            * *
 * *    Double-precision Floating-point         * *
 * *    Instance of Template                    * *
 * *                                            * *
 * ********************************************** *
 **************************************************/

template class Complex<double>;

template <> const Complex<double> Complex<double>::ZERO = Complex<double>( 0, 0 );
template <> const Complex<double> Complex<double>::ONE  = Complex<double>( 1, 0 );
template <> const Complex<double> Complex<double>::I    = Complex<double>( 0, 1 );

template <> const Complex<double> Complex<double>::UNITS[2] = {
    Complex<double>::ONE,
    Complex<double>::I
};

template <> char Complex<double>::imaginary_symbol = 'i';

template std::ostream & operator << ( std::ostream & out, const Complex<double> & );

template Complex<double> operator + ( double, const Complex<double> & );
template Complex<double> operator + ( long, const Complex<double> & );
template Complex<double> operator - ( double, const Complex<double> & );
template Complex<double> operator - ( long, const Complex<double> & );
template Complex<double> operator * ( double, const Complex<double> & );
template Complex<double> operator * ( long, const Complex<double> & );
template Complex<double> operator / ( double, const Complex<double> & );
template Complex<double> operator / ( long, const Complex<double> & );

template bool operator == ( double, const Complex<double> & );
template bool operator == ( long, const Complex<double> & );
template bool operator != ( double, const Complex<double> & );
template bool operator != ( long, const Complex<double> & );

template double real( const Complex<double> & );
template double imag_i( const Complex<double> & );
template double csgn( const Complex<double> & );
template double abs( const Complex<double> & );
template double norm( const Complex<double> & );
template double abs_imag( const Complex<double> & );
template double norm_imag( const Complex<double> & );
template Complex<double> imag( const Complex<double> & );
template Complex<double> conj( const Complex<double> & );
template Complex<double> signum( const Complex<double> & );
template Complex<double> sqr( const Complex<double> & );
template Complex<double> sqrt( const Complex<double> & );
template Complex<double> exp( const Complex<double> & );
template Complex<double> log( const Complex<double> & );
template Complex<double> log10( const Complex<double> & );
template Complex<double> pow( const Complex<double> &, const Complex<double> & );
template Complex<double> pow( const Complex<double> &, double );
template Complex<double> inverse( const Complex<double> & );
template Complex<double> sin( const Complex<double> & );
template Complex<double> cos( const Complex<double> & );
template Complex<double> tan( const Complex<double> & );
template Complex<double> sec( const Complex<double> & );
template Complex<double> csc( const Complex<double> & );
template Complex<double> cot( const Complex<double> & );
template Complex<double> sinh( const Complex<double> & );
template Complex<double> cosh( const Complex<double> & );
template Complex<double> tanh( const Complex<double> & );
template Complex<double> sech( const Complex<double> & );
template Complex<double> csch( const Complex<double> & );
template Complex<double> coth( const Complex<double> & );
template Complex<double> asin( const Complex<double> & );
template Complex<double> acos( const Complex<double> & );
template Complex<double> atan( const Complex<double> & );
template Complex<double> asec( const Complex<double> & );
template Complex<double> acsc( const Complex<double> & );
template Complex<double> acot( const Complex<double> & );
template Complex<double> asinh( const Complex<double> & );
template Complex<double> acosh( const Complex<double> & );
template Complex<double> atanh( const Complex<double> & );
template Complex<double> asech( const Complex<double> & );
template Complex<double> acsch( const Complex<double> & );
template Complex<double> acoth( const Complex<double> & );
template Complex<double> bessel_J( int, const Complex<double> & );
template Complex<double> floor( const Complex<double> & );
template Complex<double> ceil( const Complex<double> & );
template Complex<double> horner( const Complex<double> &, double *, unsigned int );
template Complex<double> horner( const Complex<double> &, double *, double *, unsigned int );
template Complex<double> horner( const Complex<double> &, Complex<double> *, unsigned int );
template Complex<double> horner( const Complex<double> &, Complex<double> *, double *, unsigned int );
template Complex<double> horner( const Complex<double> &, double *, Complex<double> *, unsigned int );
template Complex<double> horner( const Complex<double> &, Complex<double> *, Complex<double> *, unsigned int );

/**************************************************
 * ********************************************** *
 * *                                            * *
 * *    Floating-point Instance of Template     * *
 * *                                            * *
 * ********************************************** *
 **************************************************/

template class Complex<float>;

template <> const Complex<float> Complex<float>::ZERO = Complex<float>( 0, 0 );
template <> const Complex<float> Complex<float>::ONE  = Complex<float>( 1, 0 );
template <> const Complex<float> Complex<float>::I    = Complex<float>( 0, 1 );

template <> const Complex<float> Complex<float>::UNITS[2] = {
    Complex<float>::ONE,
    Complex<float>::I
};

template <> char Complex<float>::imaginary_symbol = 'i';

template std::ostream & operator << ( std::ostream & out, const Complex<float> & );

template Complex<float> operator + ( float, const Complex<float> & );
template Complex<float> operator + ( long, const Complex<float> & );
template Complex<float> operator - ( float, const Complex<float> & );
template Complex<float> operator - ( long, const Complex<float> & );
template Complex<float> operator * ( float, const Complex<float> & );
template Complex<float> operator * ( long, const Complex<float> & );
template Complex<float> operator / ( float, const Complex<float> & );
template Complex<float> operator / ( long, const Complex<float> & );

template bool operator == ( float, const Complex<float> & );
template bool operator == ( long, const Complex<float> & );
template bool operator != ( float, const Complex<float> & );
template bool operator != ( long, const Complex<float> & );

template float real( const Complex<float> & );
template float imag_i( const Complex<float> & );
template float csgn( const Complex<float> & );
template float abs( const Complex<float> & );
template float norm( const Complex<float> & );
template float abs_imag( const Complex<float> & );
template float norm_imag( const Complex<float> & );
template Complex<float> imag( const Complex<float> & );
template Complex<float> conj( const Complex<float> & );
template Complex<float> signum( const Complex<float> & );
template Complex<float> sqr( const Complex<float> & );
template Complex<float> sqrt( const Complex<float> & );
template Complex<float> exp( const Complex<float> & );
template Complex<float> log( const Complex<float> & );
template Complex<float> log10( const Complex<float> & );
template Complex<float> pow( const Complex<float> &, const Complex<float> & );
template Complex<float> pow( const Complex<float> &, float );
template Complex<float> inverse( const Complex<float> & );
template Complex<float> sin( const Complex<float> & );
template Complex<float> cos( const Complex<float> & );
template Complex<float> tan( const Complex<float> & );
template Complex<float> sec( const Complex<float> & );
template Complex<float> csc( const Complex<float> & );
template Complex<float> cot( const Complex<float> & );
template Complex<float> sinh( const Complex<float> & );
template Complex<float> cosh( const Complex<float> & );
template Complex<float> tanh( const Complex<float> & );
template Complex<float> sech( const Complex<float> & );
template Complex<float> csch( const Complex<float> & );
template Complex<float> coth( const Complex<float> & );
template Complex<float> asin( const Complex<float> & );
template Complex<float> acos( const Complex<float> & );
template Complex<float> atan( const Complex<float> & );
template Complex<float> asec( const Complex<float> & );
template Complex<float> acsc( const Complex<float> & );
template Complex<float> acot( const Complex<float> & );
template Complex<float> asinh( const Complex<float> & );
template Complex<float> acosh( const Complex<float> & );
template Complex<float> atanh( const Complex<float> & );
template Complex<float> asech( const Complex<float> & );
template Complex<float> acsch( const Complex<float> & );
template Complex<float> acoth( const Complex<float> & );
template Complex<float> bessel_J( int, const Complex<float> & );
template Complex<float> floor( const Complex<float> & );
template Complex<float> ceil( const Complex<float> & );
template Complex<float> horner( const Complex<float> &, float *, unsigned int );
template Complex<float> horner( const Complex<float> &, float *, float *, unsigned int );
template Complex<float> horner( const Complex<float> &, Complex<float> *, unsigned int );
template Complex<float> horner( const Complex<float> &, Complex<float> *, float *, unsigned int );
template Complex<float> horner( const Complex<float> &, float *, Complex<float> *, unsigned int );
template Complex<float> horner( const Complex<float> &, Complex<float> *, Complex<float> *, unsigned int );

/************************************************************************
 * ******************************************************************** *
 * *                                                                  * *
 * *    Long Double-Precision Floating-point Instance of Template     * *
 * *                                                                  * *
 * ******************************************************************** *
 ************************************************************************/

template class Complex<long double>;

template <> const Complex<long double> Complex<long double>::ZERO = Complex<long double>( 0, 0 );
template <> const Complex<long double> Complex<long double>::ONE  = Complex<long double>( 1, 0 );
template <> const Complex<long double> Complex<long double>::I    = Complex<long double>( 0, 1 );

template <> const Complex<long double> Complex<long double>::UNITS[2] = {
    Complex<long double>::ONE,
    Complex<long double>::I
};

template <> char Complex<long double>::imaginary_symbol = 'i';

template std::ostream & operator << ( std::ostream & out, const Complex<long double> & );

template Complex<long double> operator + ( long double, const Complex<long double> & );
template Complex<long double> operator + ( long, const Complex<long double> & );
template Complex<long double> operator - ( long double, const Complex<long double> & );
template Complex<long double> operator - ( long, const Complex<long double> & );
template Complex<long double> operator * ( long double, const Complex<long double> & );
template Complex<long double> operator * ( long, const Complex<long double> & );
template Complex<long double> operator / ( long double, const Complex<long double> & );
template Complex<long double> operator / ( long, const Complex<long double> & );

template bool operator == ( long double, const Complex<long double> & );
template bool operator == ( long, const Complex<long double> & );
template bool operator != ( long double, const Complex<long double> & );
template bool operator != ( long, const Complex<long double> & );

template long double real( const Complex<long double> & );
template long double imag_i( const Complex<long double> & );
template long double csgn( const Complex<long double> & );
template long double abs( const Complex<long double> & );
template long double norm( const Complex<long double> & );
template long double abs_imag( const Complex<long double> & );
template long double norm_imag( const Complex<long double> & );
template Complex<long double> imag( const Complex<long double> & );
template Complex<long double> conj( const Complex<long double> & );
template Complex<long double> signum( const Complex<long double> & );
template Complex<long double> sqr( const Complex<long double> & );
template Complex<long double> sqrt( const Complex<long double> & );
template Complex<long double> exp( const Complex<long double> & );
template Complex<long double> log( const Complex<long double> & );
template Complex<long double> log10( const Complex<long double> & );
template Complex<long double> pow( const Complex<long double> &, const Complex<long double> & );
template Complex<long double> pow( const Complex<long double> &, long double );
template Complex<long double> inverse( const Complex<long double> & );
template Complex<long double> sin( const Complex<long double> & );
template Complex<long double> cos( const Complex<long double> & );
template Complex<long double> tan( const Complex<long double> & );
template Complex<long double> sec( const Complex<long double> & );
template Complex<long double> csc( const Complex<long double> & );
template Complex<long double> cot( const Complex<long double> & );
template Complex<long double> sinh( const Complex<long double> & );
template Complex<long double> cosh( const Complex<long double> & );
template Complex<long double> tanh( const Complex<long double> & );
template Complex<long double> sech( const Complex<long double> & );
template Complex<long double> csch( const Complex<long double> & );
template Complex<long double> coth( const Complex<long double> & );
template Complex<long double> asin( const Complex<long double> & );
template Complex<long double> acos( const Complex<long double> & );
template Complex<long double> atan( const Complex<long double> & );
template Complex<long double> asec( const Complex<long double> & );
template Complex<long double> acsc( const Complex<long double> & );
template Complex<long double> acot( const Complex<long double> & );
template Complex<long double> asinh( const Complex<long double> & );
template Complex<long double> acosh( const Complex<long double> & );
template Complex<long double> atanh( const Complex<long double> & );
template Complex<long double> asech( const Complex<long double> & );
template Complex<long double> acsch( const Complex<long double> & );
template Complex<long double> acoth( const Complex<long double> & );
template Complex<long double> bessel_J( int, const Complex<long double> & );
template Complex<long double> floor( const Complex<long double> & );
template Complex<long double> ceil( const Complex<long double> & );
template Complex<long double> horner( const Complex<long double> &, long double *, unsigned int );
template Complex<long double> horner( const Complex<long double> &, long double *, long double *, unsigned int );
template Complex<long double> horner( const Complex<long double> &, Complex<long double> *, unsigned int );
template Complex<long double> horner( const Complex<long double> &, Complex<long double> *, long double *, unsigned int );
template Complex<long double> horner( const Complex<long double> &, long double *, Complex<long double> *, unsigned int );
template Complex<long double> horner( const Complex<long double> &, Complex<long double> *, Complex<long double> *, unsigned int );
