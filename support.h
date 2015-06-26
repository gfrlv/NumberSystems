#ifndef SUPPORT_H
#define SUPPORT_H

/******************************************
 * C++ Support for the Number Systems Package
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
#include <cstdlib>

template <typename T>
class Support {
    public:
        static const T POS_INF;
        static const T NEG_INF;
        static const T NaN;
        static const T PI;
        static const T PI2;
        static const T GAMMA_TAB[171];

        static T Gamma( int );
        static T sec( T );
        static T csc( T );
        static T cot( T );
        static T sech( T );
        static T csch( T );
        static T coth( T );
        static T sign( T );
        static bool is_pos_inf( T );
        static bool is_neg_inf( T );
        static bool is_pos_zero( T );
        static bool is_neg_zero( T );
        static bool is_inf( T );
        static bool is_nan( T );

        static short bigendian;

        /******************************************
        * IO Stream Functions
        ******************************************/

        static void print_real( T, std::ostream & );
        static void print_imaginary( T, const std::string &, std::ostream & );
        static void print_imaginary( T, char, std::ostream & );
};


#endif // SUPPORT_H

