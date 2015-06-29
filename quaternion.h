#ifndef QUATERNION_H
#define QUATERNION_H
/******************************************
 * C++ Quaternions
 * Version: 1.0.9
 * Author:  Douglas Wilhelm Harder
 * Date:    2008/03/03
 *
 * Copyright (c) 2006-2008 by Douglas Wilhelm Harder.
 * All rights reserved.
 ******************************************/

#include <iostream>
#include <cmath>
#include <string>

template <typename T> class Quaternion;

template <typename T> std::ostream & operator << ( std::ostream &, const Quaternion<T> & );

template <typename T> Quaternion<T> operator + ( T, const Quaternion<T> & );
template <typename T> Quaternion<T> operator + ( long, const Quaternion<T> & );
template <typename T> Quaternion<T> operator - ( T, const Quaternion<T> & );
template <typename T> Quaternion<T> operator - ( long, const Quaternion<T> & );
template <typename T> Quaternion<T> operator * ( T, const Quaternion<T> & );
template <typename T> Quaternion<T> operator * ( long, const Quaternion<T> & );
template <typename T> Quaternion<T> operator / ( T, const Quaternion<T> & );
template <typename T> Quaternion<T> operator / ( long, const Quaternion<T> & );

template <typename T> bool operator == ( T, const Quaternion<T> & );
template <typename T> bool operator == ( long, const Quaternion<T> & );
template <typename T> bool operator != ( T, const Quaternion<T> & );
template <typename T> bool operator != ( long, const Quaternion<T> & );

template <typename T> T real( const Quaternion<T> & );
template <typename T> T imag_i( const Quaternion<T> & );
template <typename T> T imag_j( const Quaternion<T> & );
template <typename T> T imag_k( const Quaternion<T> & );
template <typename T> T csgn( const Quaternion<T> & );
template <typename T> T abs( const Quaternion<T> & );
template <typename T> T norm( const Quaternion<T> & );
template <typename T> T abs_imag( const Quaternion<T> & );
template <typename T> T norm_imag( const Quaternion<T> & );
template <typename T> T arg( const Quaternion<T> & );
template <typename T> Quaternion<T> imag( const Quaternion<T> & );
template <typename T> Quaternion<T> conj( const Quaternion<T> & );
template <typename T> Quaternion<T> signum( const Quaternion<T> & );
template <typename T> Quaternion<T> sqr( const Quaternion<T> & );
template <typename T> Quaternion<T> sqrt( const Quaternion<T> & );
template <typename T> Quaternion<T> rotate( const Quaternion<T> &, const Quaternion<T> & );
template <typename T> Quaternion<T> exp( const Quaternion<T> & );
template <typename T> Quaternion<T> log( const Quaternion<T> & );
template <typename T> Quaternion<T> log10( const Quaternion<T> & );
template <typename T> Quaternion<T> pow( const Quaternion<T> &, const Quaternion<T> & );
template <typename T> Quaternion<T> pow( const Quaternion<T> &, T );
template <typename T> Quaternion<T> inverse( const Quaternion<T> & );
template <typename T> Quaternion<T> cross(const Quaternion<T> & );
template <typename T> Quaternion<T> sin( const Quaternion<T> & );
template <typename T> Quaternion<T> cos( const Quaternion<T> & );
template <typename T> Quaternion<T> tan( const Quaternion<T> & );
template <typename T> Quaternion<T> sec( const Quaternion<T> & );
template <typename T> Quaternion<T> csc( const Quaternion<T> & );
template <typename T> Quaternion<T> cot( const Quaternion<T> & );
template <typename T> Quaternion<T> sinh( const Quaternion<T> & );
template <typename T> Quaternion<T> cosh( const Quaternion<T> & );
template <typename T> Quaternion<T> tanh( const Quaternion<T> & );
template <typename T> Quaternion<T> sech( const Quaternion<T> & );
template <typename T> Quaternion<T> csch( const Quaternion<T> & );
template <typename T> Quaternion<T> coth( const Quaternion<T> & );
template <typename T> Quaternion<T> asin( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> acos( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> atan( const Quaternion<T> & );
template <typename T> Quaternion<T> asec( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> acsc( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> acot( const Quaternion<T> & );
template <typename T> Quaternion<T> asinh( const Quaternion<T> & );
template <typename T> Quaternion<T> acosh( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> atanh( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> asech( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> acsch( const Quaternion<T> & );
template <typename T> Quaternion<T> acoth( const Quaternion<T> &, const Quaternion<T> & = Quaternion<T>::I );
template <typename T> Quaternion<T> bessel_J( int, const Quaternion<T> & );
template <typename T> Quaternion<T> floor( const Quaternion<T> & );
template <typename T> Quaternion<T> ceil( const Quaternion<T> & );
template <typename T> Quaternion<T> horner( const Quaternion<T> &, T *, unsigned int );
template <typename T> Quaternion<T> horner( const Quaternion<T> &, T *, T *, unsigned int );

template <typename T = double> class Quaternion {
    private:
        T r, i, j, k;

        inline static Quaternion multiplier( T, T, const Quaternion & );
        inline static Quaternion make_inf( T, T );
        inline static Quaternion make_i( T, T );

    public:
        const static Quaternion ZERO;
        const static Quaternion ONE;
        const static Quaternion I;
        const static Quaternion J;
        const static Quaternion K;

        const static Quaternion UNITS[4];

        /******************************************
         * Constructor and Copy Constructor
         ******************************************/

        Quaternion( T, T, T, T );
        Quaternion( T = 0.0 );

        /******************************************
         * Assignment Operator
         ******************************************/

        const Quaternion & operator = ( const T & );

        /******************************************
         * Mutating Arithmetic Operators
         ******************************************/

        Quaternion & operator += ( const Quaternion & );
        Quaternion & operator -= ( const Quaternion & );
        Quaternion & operator *= ( const Quaternion & );
        Quaternion & operator /= ( const Quaternion & );

        Quaternion & operator += ( T );
        Quaternion & operator -= ( T );
        Quaternion & operator *= ( T );
        Quaternion & operator /= ( T );

        Quaternion & operator ++();
        Quaternion operator ++( int );
        Quaternion & operator --();
        Quaternion operator --( int );

        /******************************************
         * Real-valued Functions
         ******************************************/

        T real() const;
        T operator []( int ) const;
        T &operator []( int );
        T imag_i() const;
        T imag_j() const;
        T imag_k() const;
        T csgn() const;
        T abs() const;
        T norm() const;
        T abs_imag() const;
        T norm_imag() const;
        T arg() const;

        /******************************************
         * Quaternion-valued Functions
         ******************************************/

        Quaternion imag() const;
        Quaternion conj() const;
        Quaternion operator * () const;
        Quaternion signum() const;
        Quaternion sqr() const;
        Quaternion sqrt() const;
        Quaternion rotate( const Quaternion & ) const;

        /******************************************
         * Boolean-valued Functions
         ******************************************/

        bool is_imaginary() const;
        bool is_inf() const;
        bool is_nan() const;
        bool is_neg_inf() const;
        bool is_pos_inf() const;
        bool is_real() const;
        bool is_real_inf() const;
        bool is_zero() const;

        /******************************************
         * Exponential and Logarithmic Functions
         ******************************************/

        Quaternion exp() const;
        Quaternion log() const;
        Quaternion log10() const;
        Quaternion pow(const Quaternion & w) const;
        Quaternion pow(T x) const;
        Quaternion inverse() const;


        Quaternion cross(const Quaternion &) const;

        /******************************************
         * Trigonometric and Hyperbolic Functions
         ******************************************/

        Quaternion   sin() const;   Quaternion   cos() const;   Quaternion   tan() const;
        Quaternion   sec() const;   Quaternion   csc() const;   Quaternion   cot() const;
        Quaternion  sinh() const;   Quaternion  cosh() const;   Quaternion  tanh() const;
        Quaternion  sech() const;   Quaternion  csch() const;   Quaternion  coth() const;
        Quaternion  asin( const Quaternion & = I ) const;
        Quaternion  acos( const Quaternion & = I ) const;       Quaternion  atan() const;
        Quaternion  asec( const Quaternion & = I ) const;
        Quaternion  acsc( const Quaternion & = I ) const;       Quaternion  acot() const;

        Quaternion asinh() const;
        Quaternion acosh( const Quaternion & = I ) const;
        Quaternion atanh( const Quaternion & = I ) const;

        Quaternion asech( const Quaternion & = I ) const;
        Quaternion acsch() const;
        Quaternion acoth( const Quaternion & = I ) const;

        Quaternion bessel_J( int ) const;

        /******************************************
         * Integer Functions
         ******************************************/

        Quaternion ceil() const;
        Quaternion floor() const;

        /******************************************
         * Horner's Rule
         ******************************************/

        Quaternion horner( T *, unsigned int ) const;
        Quaternion horner( T *, T *, unsigned int ) const;

        /******************************************
         * Random Factories
         ******************************************/

        static Quaternion random();
        static Quaternion random_imag();
        static Quaternion random_real();

        /******************************************
         * Binary Arithmetic Operators
         ******************************************/

        Quaternion operator + ( const Quaternion & ) const;
        Quaternion operator + ( T ) const;

        Quaternion operator - ( const Quaternion & ) const;
        Quaternion operator - ( T ) const;

        Quaternion operator * ( const Quaternion & ) const;
        Quaternion operator * ( T ) const;

        Quaternion operator / ( const Quaternion & ) const;
        Quaternion operator / ( T ) const;

        /******************************************
         * Unary Arithmetic Operators
         ******************************************/

        Quaternion operator - () const;

        /******************************************
         * Binary Boolean Operators
         ******************************************/

        bool operator == ( const Quaternion & ) const;
        bool operator == ( T ) const;

        bool operator != ( const Quaternion & ) const;
        bool operator != ( T ) const;
};

#endif // QUATERNION_H

