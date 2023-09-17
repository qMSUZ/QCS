/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2005 -- 2006 by Kamil Pawłowski                         *
 *                                                                         *
 *   Part of the Quantum Computing Simulator:                              *
 *   https://github.com/qMSUZ/QCS                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>

#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif

#ifdef PYTHON_SCRIPT
#include <Python.h>
#endif


#include "qcs.h"
#include "qcs_complex.h"
#include "qcs_misc.h"

const tf_qcs_complex z0 = {0, 0};
const tf_qcs_complex z1 = {1, 0};
const tf_qcs_complex zj = {0, 1};

DYNAMIC_LIB_DECORATION tf_qcs_complex *qcs_create_complex_float_arg(tf_qcs_real_number a, tf_qcs_real_number b)
{
    return qcs_create_complex( (tf_qcs_real_number)a, (tf_qcs_real_number)b );
}

DYNAMIC_LIB_DECORATION tf_qcs_complex *qcs_create_complex(tf_qcs_real_number a, tf_qcs_real_number b)
{
    tf_qcs_complex *n;

    n = (tf_qcs_complex*)malloc(sizeof(tf_qcs_complex));

    n->re=a;
    n->im=b;

    return n;
}

/* Complex addition */
DYNAMIC_LIB_DECORATION void qcs_complex_add(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    c_out->re = a_in->re + b_in->re;
    c_out->im = a_in->im + b_in->im;
}

/* Complex subtraction */
DYNAMIC_LIB_DECORATION void qcs_complex_sub(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    c_out->re = a_in->re - b_in->re;
    c_out->im = a_in->im - b_in->im;
}

/* Complex multiplication */
DYNAMIC_LIB_DECORATION void qcs_complex_mul(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    c_out->re = a_in->re * b_in->re - a_in->im * b_in->im;
    c_out->im = a_in->re * b_in->im + a_in->im * b_in->re;
}

/* Complex division */
DYNAMIC_LIB_DECORATION void qcs_complex_div(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    tf_qcs_real_number s = qcs_sqr(b_in->re) + qcs_sqr(b_in->im);
    tf_qcs_complex z;

    z.re = (a_in->re * b_in->re + a_in->im * b_in->im) / s;
    z.im = (a_in->im * b_in->re - a_in->re * b_in->im) / s;

    *c_out = z;
}

/*
  rounds real and imaginary part
  at specified power of 10
*/
DYNAMIC_LIB_DECORATION void qcs_complex_roundn(tf_qcs_complex *a_in, int n)
{
     tf_qcs_real_number factor;

     factor = powf(10, -n);
     a_in->re = roundf(a_in->re * factor) / factor;
     a_in->im = roundf(a_in->im * factor) / factor;
}

/* Complex & Real multiplication */
DYNAMIC_LIB_DECORATION void qcs_complex_mul_tf_qcs_real_number(tf_qcs_complex *a_in, tf_qcs_real_number b_in, tf_qcs_complex *c_out)
{
    c_out->re = a_in->re * b_in;
    c_out->im = a_in->im * b_in;
}

/* Real by Complex multiplication - the same as Complex by Real */
DYNAMIC_LIB_DECORATION void qcs_tf_qcs_real_number_mul_complex(tf_qcs_real_number a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    qcs_complex_mul_tf_qcs_real_number(b_in, a_in, c_out);
}

/* Complex by Real division */
DYNAMIC_LIB_DECORATION void qcs_complex_div_tf_qcs_real_number(tf_qcs_complex *a_in, tf_qcs_real_number b_in, tf_qcs_complex *c_out)
{
    c_out->re = a_in->re / b_in;
    c_out->im = a_in->im / b_in;
}

/* Real by Complex division */
DYNAMIC_LIB_DECORATION void qcs_tf_qcs_real_number_div_complex(tf_qcs_real_number a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out)
{
    tf_qcs_real_number s = qcs_sqr(b_in->re) + qcs_sqr(b_in->im);
    c_out->re = a_in * b_in->re / s;
    c_out->im = -a_in * b_in->im / s;
}

/* Inversion of Complex */
DYNAMIC_LIB_DECORATION void qcs_inv_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out)
{
    tf_qcs_real_number s = qcs_sqr(a_in->re) + qcs_sqr(a_in->im);
    b_out->re = a_in->re / s;
    b_out->im = -a_in->im / s;
}

/* Conjugate of Complex */
DYNAMIC_LIB_DECORATION void qcs_conj_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out)
{
    tf_qcs_complex tmp;

    tmp.re = a_in->re;
    tmp.im = -a_in->im;

    b_out->re = tmp.re;
    b_out->im = tmp.im;
}

/* For complex Z=X+i*Y, EXP(Z) = EXP(X)*(COS(Y)+i*SIN(Y)) */
DYNAMIC_LIB_DECORATION void qcs_exp_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out)
{

     b_out->re=exp(a_in->re)*cos(a_in->im);
     b_out->im=exp(a_in->re)*sin(a_in->im);

     if(fabs(b_out->re) < QCS_EPS) b_out->re=0;
     if(fabs(b_out->im) < QCS_EPS) b_out->im=0;
}

/* Modulo of Complex */
DYNAMIC_LIB_DECORATION void qcs_mod_complex(tf_qcs_complex *a_in, tf_qcs_real_number *b_out)
{
    *b_out = (tf_qcs_real_number)sqrt(qcs_sqr(a_in->re) + qcs_sqr(a_in->im));
}

/* square root of complex */
DYNAMIC_LIB_DECORATION void qcs_sqrt_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out)
{
     tf_qcs_complex c;
     tf_qcs_real_number x,y,w,r;

     if ( (a_in->re == 0.0) && (a_in->im == 0.0) )
     {
        b_out->re=0;
        b_out->im=0;
     }

     x=fabs(a_in->re);
     y=fabs(a_in->im);
     if (x >= y)
     {
        r=y/x;
        w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
     }
     else
     {
        r=x/y;
        w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
     }
     if (a_in->re >= 0.0)
     {
        c.re=w;
        c.im=a_in->im / (2.0*w);
     }
     else
     {
        c.im=(a_in->im >= 0) ? w : -w;
        c.re=a_in->im / (2.0*a_in->im);
     }
     b_out->re=c.re;
     b_out->im=c.im;
}

/* Print complex to file*/
DYNAMIC_LIB_DECORATION void qcs_print_complex_to_file(tf_qcs_complex *a_in, FILE *file_out)
{
    fprintf(file_out, "%2.9f + i%2.9f\t", a_in->re, a_in->im);
}

/* Print complex to stdout*/
DYNAMIC_LIB_DECORATION void qcs_print_complex(tf_qcs_complex *a_in)
{
#ifdef PYTHON_SCRIPT
     if(a_in->im>0)
         PySys_WriteStdout("%2.2f + %2.2fi ", a_in->re,  a_in->im);
     if(a_in->im<0)
         PySys_WriteStdout("%2.2f - %2.2fi ", a_in->re, -a_in->im);
     if(a_in->re!=0 && a_in->im == 0)
         PySys_WriteStdout("%2.2f ", a_in->re);
     if(a_in->re == 0 && a_in->im == 0)
         PySys_WriteStdout("0.00 ");
#else
     if(a_in->im>0)
         printf("%2.2f + %2.2fi ", a_in->re,  a_in->im);
     if(a_in->im<0)
         printf("%2.2f - %2.2fi ", a_in->re, -a_in->im);
     if(a_in->re!=0 && a_in->im == 0 )
         printf("%2.2f ", a_in->re);
     if(a_in->re == 0 && a_in->im == 0)
         printf("0.00 ");
#endif
}
