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

#ifndef __qcs_complex_h__
#define __qcs_complex_h__

#include <stdio.h>

#include "qcs_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    tf_qcs_real_number re;
    tf_qcs_real_number im;
} Complex;


typedef Complex tf_qcs_complex;
typedef Complex* pf_qcs_complex;

extern const tf_qcs_complex z0;
extern const tf_qcs_complex z1;
extern const tf_qcs_complex zj;

/* complex functions */

tf_qcs_complex *qcs_create_complex_float_arg(tf_qcs_real_number a, tf_qcs_real_number b);
tf_qcs_complex *qcs_create_complex(tf_qcs_real_number a, tf_qcs_real_number b);

void qcs_complex_add(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);
void qcs_complex_sub(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);
void qcs_complex_mul(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);
void qcs_complex_div(tf_qcs_complex *a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);

void qcs_complex_roundn(tf_qcs_complex *a_in, int n);

void qcs_complex_mul_tf_qcs_real_number(tf_qcs_complex *a_in, tf_qcs_real_number b_in, tf_qcs_complex *c_out);
void qcs_tf_qcs_real_number_mul_complex(tf_qcs_real_number a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);
void qcs_tf_qcs_real_number_div_complex(tf_qcs_real_number a_in, tf_qcs_complex *b_in, tf_qcs_complex *c_out);

void qcs_inv_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out);
void qcs_conj_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out);
void qcs_exp_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out);
void qcs_mod_complex(tf_qcs_complex *a_in, tf_qcs_real_number *b_out);
void qcs_sqrt_complex(tf_qcs_complex *a_in, tf_qcs_complex *b_out);

void qcs_print_complex_to_file(tf_qcs_complex *a_in, FILE *b_in);
void qcs_print_complex(tf_qcs_complex *a_in);

#ifdef __cplusplus
}
#endif

#endif /* __qcs_complex_h__ */
