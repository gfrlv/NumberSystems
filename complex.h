/******************************************
 * C++ Complex Numbers
 * Version: 1.0.9
 * Author:  Douglas Wilhelm Harder
 * Date:    2008/03/03
 *
 * Copyright (c) 2006-2008 by Douglas Wilhelm Harder.
 * All rights reserved.
 ******************************************/

#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include <cmath>
#include <string>

template <typename T> class Complex;

template <typename T> std::ostream & operator << ( std::ostream &, const Complex<T> & );

template <typename T> Complex<T> operator + ( T, const Complex<T> & );
template <typename T> Complex<T> operator + ( long, const Complex<T> & );
template <typename T> Complex<T> operator - ( T, const Complex<T> & );
template <typename T> Complex<T> operator - ( long, const Complex<T> & );
template <typename T> Complex<T> operator * ( T, const Complex<T> & );
template <typename T> Complex<T> operator * ( long, const Complex<T> & );
template <typename T> Complex<T> operator / ( T, const Complex<T> & );
template <typename T> Complex<T> operator / ( long, const Complex<T> & );

template <typename T> bool operator == ( T, const Complex<T> & );
template <typename T> bool operator == ( long, const Complex<T> & );
template <typename T> bool operator != ( T, const Complex<T> & );
template <typename T> bool operator != ( long, const Complex<T> & );

template <typename T> T real( const Complex<T> & );
template <typename T> T imag_i( const Complex<T> & );
template <typename T> T csgn( const Complex<T> & );
template <typename T> T abs( const Complex<T> & );
template <typename T> T norm( const Complex<T> & );
template <typename T> T abs_imag( const Complex<T> & );
template <typename T> T norm_imag( const Complex<T> & );
template <typename T> T arg( const Complex<T> & );
template <typename T> Complex<T> imag( const Complex<T> & );
template <typename T> Complex<T> conj( const Complex<T> & );
template <typename T> Complex<T> signum( const Complex<T> & );
template <typename T> Complex<T> sqr( const Complex<T> & );
template <typename T> Complex<T> sqrt( const Complex<T> & );
template <typename T> Complex<T> exp( const Complex<T> & );
template <typename T> Complex<T> log( const Complex<T> & );
template <typename T> Complex<T> log10( const Complex<T> & );
template <typename T> Complex<T> pow( const Complex<T> &, const Complex<T> & );
template <typename T> Complex<T> pow( const Complex<T> &, T );
template <typename T> Complex<T> inverse( const Complex<T> & );
template <typename T> Complex<T> sin( const Complex<T> & );
template <typename T> Complex<T> cos( const Complex<T> & );
template <typename T> Complex<T> tan( const Complex<T> & );
template <typename T> Complex<T> sec( const Complex<T> & );
template <typename T> Complex<T> csc( const Complex<T> & );
template <typename T> Complex<T> cot( const Complex<T> & );
template <typename T> Complex<T> sinh( const Complex<T> & );
template <typename T> Complex<T> cosh( const Complex<T> & );
template <typename T> Complex<T> tanh( const Complex<T> & );
template <typename T> Complex<T> sech( const Complex<T> & );
template <typename T> Complex<T> csch( const Complex<T> & );
template <typename T> Complex<T> coth( const Complex<T> & );
template <typename T> Complex<T> asin( const Complex<T> & );
template <typename T> Complex<T> acos( const Complex<T> & );
template <typename T> Complex<T> atan( const Complex<T> & );
template <typename T> Complex<T> asec( const Complex<T> & );
template <typename T> Complex<T> acsc( const Complex<T> & );
template <typename T> Complex<T> acot( const Complex<T> & );
template <typename T> Complex<T> asinh( const Complex<T> & );
template <typename T> Complex<T> acosh( const Complex<T> & );
template <typename T> Complex<T> atanh( const Complex<T> & );
template <typename T> Complex<T> asech( const Complex<T> & );
template <typename T> Complex<T> acsch( const Complex<T> & );
template <typename T> Complex<T> acoth( const Complex<T> & );
template <typename T> Complex<T> bessel_J( int, const Complex<T> & );
template <typename T> Complex<T> floor( const Complex<T> & );
template <typename T> Complex<T> ceil( const Complex<T> & );
template <typename T> Complex<T> horner( const Complex<T> &, T *, unsigned int );
template <typename T> Complex<T> horner( const Complex<T> &, T *, T *, unsigned int );
template <typename T> Complex<T> horner( const Complex<T> &, Complex<T> *, unsigned int );
template <typename T> Complex<T> horner( const Complex<T> &, Complex<T> *, T *, unsigned int );
template <typename T> Complex<T> horner( const Complex<T> &, T *, Complex<T> *, unsigned int );
template <typename T> Complex<T> horner( const Complex<T> &, Complex<T> *, Complex<T> *, unsigned int );

template <typename T = double> class Complex {
    private:
        T r, i;
        static char imaginary_symbol;

    public:
        static char use_symbol( char sym );
        static char get_symbol();

        const static Complex ZERO;
        const static Complex ONE;
        const static Complex I;

        const static Complex UNITS[2];

        /******************************************
         * Constructor and Copy Constructor
         ******************************************/

        Complex( T, T );
        Complex( T = 0.0 );

        /******************************************
         * Assignment Operator
         ******************************************/

        const Complex & operator = ( T );

        /******************************************
         * Mutating Arithmetic Operators
         ******************************************/

        Complex & operator += ( const Complex & );
        Complex & operator -= ( const Complex & );
        Complex & operator *= ( const Complex & );
        Complex & operator /= ( const Complex & );

        Complex & operator += ( T );
        Complex & operator -= ( T );
        Complex & operator *= ( T );
        Complex & operator /= ( T );

        Complex & operator ++();
        Complex operator ++( int );
        Complex & operator --();
        Complex operator --( int );

        /******************************************
         * Real-valued Functions
         ******************************************/

        T real() const;
        T operator []( int ) const;
        T& operator []( int );
        T imag_i() const;
        T csgn() const;
        T abs() const;
        T norm() const;
        T abs_imag() const;
        T norm_imag() const;
        T arg() const;

        /******************************************
         * Complex-valued Functions
         ******************************************/

        Complex imag() const;
        Complex conj() const;
        Complex operator * () const;
        Complex signum() const;
        Complex sqr() const;
        Complex sqrt() const;

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

        Complex exp() const;
        Complex log() const;
        Complex log10() const;
        Complex pow( const Complex & ) const;
        Complex pow( T ) const;
        Complex inverse() const;

        /******************************************
         * Trigonometric and Hyperbolic Functions
         ******************************************/

        Complex   sin() const;   Complex   cos() const;   Complex   tan() const;
        Complex   sec() const;   Complex   csc() const;   Complex   cot() const;
        Complex  sinh() const;   Complex  cosh() const;   Complex  tanh() const;
        Complex  sech() const;   Complex  csch() const;   Complex  coth() const;
        Complex  asin() const;   Complex  acos() const;   Complex  atan() const;
        Complex  asec() const;   Complex  acsc() const;   Complex  acot() const;
        Complex asinh() const;   Complex acosh() const;   Complex atanh() const;
        Complex asech() const;   Complex acsch() const;   Complex acoth() const;

        Complex bessel_J( int ) const;

        /******************************************
         * Integer Functions
         ******************************************/

        Complex ceil() const;
        Complex floor() const;

        /******************************************
         * Horner's Rule
         ******************************************/

        Complex horner( T *, unsigned int ) const;
        Complex horner( T *, T *, unsigned int ) const;
        Complex horner( Complex<T> *, unsigned int ) const;
        Complex horner( Complex<T> *, T *, unsigned int ) const;
        Complex horner( T *, Complex<T> *, unsigned int ) const;
        Complex horner( Complex<T> *, Complex<T> *, unsigned int ) const;

        /******************************************
         * Random Factories
         ******************************************/

        static Complex random();
        static Complex random_imag();
        static Complex random_real();

        /******************************************
         * Binary Arithmetic Operators
         ******************************************/

        Complex operator + ( const Complex & ) const;
        Complex operator + ( T ) const;

        Complex operator - ( const Complex & ) const;
        Complex operator - ( T ) const;

        Complex operator * ( const Complex & ) const;
        Complex operator * ( T ) const;

        Complex operator / ( const Complex & ) const;
        Complex operator / ( T ) const;

        /******************************************
         * Unary Arithmetic Operators
         ******************************************/

        Complex operator - () const;

        /******************************************
         * Binary Boolean Operators
         ******************************************/

        bool operator == ( const Complex & ) const;
        bool operator == ( T ) const;

        bool operator != ( const Complex & ) const;
        bool operator != ( T ) const;
};

#endif // COMPLEX_H
