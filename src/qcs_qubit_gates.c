/***************************************************************************
 *   Copyright (C) 2018, 2019 by Marek Sawerwain                            *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
 *   Part of the Quantum Computing Simulator:                              *
 *   http://code.google.com/p/qcs                                          *
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

#include <math.h>
#include <string.h>

#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif


#include "qcs.h"
#include "qcs_complex.h"
#include "qcs_matrix.h"
#include "qcs_gates.h"
#include "qcs_qubit.h"
#include "qcs_qubit_gates.h"

/**
   0 - zero matrix
   i - ident matriz
   n - not matrix

   h - hadamard matrix

   x - pauli x matrix
   y - pauli y matrix
   z - pauli z matrix

   S - square root matrix

   9x - x rot 90
   9y - y rot 90
   9z - z rot 90

   m9x - x rot 90
   m9y - y rot 90
   m9z - z rot 90

   c - cnot matrix
   T - Toffoli gate
   t - t gate
   s - shift gate
   r - rotate alpha gate
   g - rotate theta
   p - phase gate
   f - phase f gate

*/


static tf_qcs_matrix *zero_reset_gate=NULL;
static tf_qcs_matrix *id_gate=NULL;
static tf_qcs_matrix *not_gate=NULL;
static tf_qcs_matrix *pauli_x_gate=NULL;
static tf_qcs_matrix *pauli_y_gate=NULL;
static tf_qcs_matrix *pauli_z_gate=NULL;
static tf_qcs_matrix *hadamard_gate=NULL;
static tf_qcs_matrix *square_root_gate=NULL;
static tf_qcs_matrix *x_rot90_gate=NULL;
static tf_qcs_matrix *y_rot90_gate=NULL;
static tf_qcs_matrix *z_rot90_gate=NULL;
static tf_qcs_matrix *minus_x_rot90_gate=NULL;
static tf_qcs_matrix *minus_y_rot90_gate=NULL;
static tf_qcs_matrix *minus_z_rot90_gate=NULL;
static tf_qcs_matrix *rotate_alpha_gate=NULL;
static tf_qcs_matrix *rotate_theta_gate=NULL;
static tf_qcs_matrix *t_gate=NULL;
static tf_qcs_matrix *s_gate=NULL;
static tf_qcs_matrix *phase_gate=NULL;
static tf_qcs_matrix *phase_f_gate=NULL;
static tf_qcs_matrix *phase_m11_gate=NULL;
static tf_qcs_matrix *cnot_gate=NULL;
static tf_qcs_matrix *toffoli_gate=NULL;

DYNAMIC_LIB_DECORATION tf_qcs_matrix * give_qubit_matrix(char *arg)
{
//    printf("\ngive_matrix [%s]\n", arg);

    if (strcmp(arg,"0")==0) return get_zero_reset_gate();
    if (strcmp(arg,"i")==0) return get_id_gate();
    if (strcmp(arg,"n")==0) return get_not_gate();

    if (strcmp(arg,"h")==0) return get_hadamard_gate();

    if (strcmp(arg,"r")==0) return get_rotate_alpha_gate();

    if (strcmp(arg,"x")==0) return get_pauli_x_gate();
    if (strcmp(arg,"y")==0) return get_pauli_y_gate();
    if (strcmp(arg,"z")==0) return get_pauli_z_gate();

    if (strcmp(arg,"S")==0) return get_square_root_gate();

    if (strcmp(arg,"9x")==0) return get_x_rot90_gate();
    if (strcmp(arg,"9y")==0) return get_y_rot90_gate();
    if (strcmp(arg,"9z")==0) return get_z_rot90_gate();

    if (strcmp(arg,"m9x")==0) return get_minus_x_rot90_gate();
    if (strcmp(arg,"m9y")==0) return get_minus_y_rot90_gate();
    if (strcmp(arg,"m9z")==0) return get_minus_z_rot90_gate();

    if (strcmp(arg,"c")==0) return get_cnot_gate();

    if (strcmp(arg,"T")==0) return get_toffoli_gate();

    if (strcmp(arg,"t")==0) return get_t_gate();
    if (strcmp(arg,"s")==0) return get_s_gate();
    if (strcmp(arg,"p")==0) return get_phase_gate();
    if (strcmp(arg,"f")==0) return get_phase_f_gate();

    if (strcmp(arg,"g")==0) return get_rotate_theta_gate();

    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qubit_gate(enum _qcs_gate_type t)
{
    switch(t)
    {
        case QCS_ZERO_GATE:
            return zero_reset_gate;
        break;

        case QCS_ID_GATE:
            return id_gate;
        break;

        case QCS_NOT_GATE:
            return not_gate;
        break;

	    case QCS_PAULIX_GATE:
            return pauli_x_gate;
        break;

	    case QCS_PAULIY_GATE:
            return pauli_y_gate;
        break;

	    case QCS_PAULIZ_GATE:
            return pauli_z_gate;
        break;

	    case QCS_HADAMARD_GATE:
            return hadamard_gate;
        break;

	    case QCS_SQUARE_ROOT_GATE:
            return square_root_gate;
        break;

	    case QCS_X_ROT90_GATE:
            return x_rot90_gate;
        break;

	    case QCS_Y_ROT90_GATE:
            return y_rot90_gate;
        break;

	    case QCS_Z_ROT90_GATE:
            return z_rot90_gate;
        break;

	    case QCS_MX_ROT90_GATE:
            return minus_x_rot90_gate;
        break;

	    case QCS_MY_ROT90_GATE:
            return minus_y_rot90_gate;
        break;

	    case QCS_MZ_ROT90_GATE:
            return minus_z_rot90_gate;
        break;

	    case QCS_ROT_ALPHA_GATE:
            return rotate_alpha_gate;
        break;

	    case QCS_ROT_THETA_GATE:
            return rotate_theta_gate;
        break;

	    case QCS_T_GATE:
            return t_gate;
        break;

	    case QCS_S_GATE:
            return s_gate;
        break;

	    case QCS_PHASE_GATE:
            return phase_gate;
        break;

	    case QCS_PHASE_F_GATE:
            return phase_f_gate;
        break;

        default:
            return NULL;
    }
}

/* QUBIT GATES */

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_zero_reset_gate()
{
   return zero_reset_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_id_gate()
{
   return id_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_not_gate()
{
   return not_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_pauli_x_gate()
{
   return pauli_x_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_pauli_y_gate()
{
   return pauli_y_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_pauli_z_gate()
{
   return pauli_z_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_hadamard_gate()
{
   return hadamard_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_square_root_gate()
{
   return square_root_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_x_rot90_gate()
{
   return x_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_y_rot90_gate()
{
   return y_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_z_rot90_gate()
{
   return z_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_minus_x_rot90_gate()
{
   return minus_x_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_minus_y_rot90_gate()
{
   return minus_y_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_minus_z_rot90_gate()
{
   return minus_z_rot90_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_rotate_alpha_gate()
{
   return rotate_alpha_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_rotate_theta_gate()
{
   return rotate_theta_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_t_gate()
{
   return t_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_s_gate()
{
   return s_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_phase_gate()
{
   return phase_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_phase_f_gate()
{
   return phase_f_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_phase_m11_gate()
{
   return phase_m11_gate;
}


DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_cnot_gate()
{
   return cnot_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_toffoli_gate()
{
   return toffoli_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_swap_gate()
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_fredkin_gate()
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_rot_x_gate(tf_qcs_real_number theta)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_rot_y_gate(tf_qcs_real_number theta)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_rot_z_gate(tf_qcs_real_number theta)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_for_e_2q_gate_float_arg(float gamma)
{
    return qcs_create_matrix_for_e_2q_gate( gamma );
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_for_e_2q_gate(tf_qcs_real_number gamma)
{
        tf_qcs_matrix *m;

        tf_qcs_real_number gamma_div_2;

        gamma_div_2 = gamma / (tf_qcs_real_number)2.0;

        m=qcs_create_matrix(4,4);

        (m->m+0)->re=cos(gamma_div_2); (m->m+0)->im=0; (m->m+3)->re=0; (m->m+3)->im=sin(gamma_div_2);

        (m->m+5)->re=cos(gamma_div_2); (m->m+5)->im=0;      (m->m+6)->re=0; (m->m+6)->im=-sin(gamma_div_2);
        (m->m+9)->re=0.0; (m->m+9)->im=-sin(gamma_div_2);  (m->m+10)->re=cos(gamma_div_2); (m->m+10)->im=0;

        (m->m+12)->re=0.0; (m->m+12)->im=sin(gamma_div_2); (m->m+15)->re=cos(gamma_div_2); (m->m+15)->im=0;


        return m;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_for_ux_1q_gate_float_arg(float phi, float theta)
{
    return qcs_create_matrix_for_ux_1q_gate(phi, theta);
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_for_ux_1q_gate(tf_qcs_real_number phi, tf_qcs_real_number theta)
{
    tf_qcs_matrix *m;

    m=qcs_create_matrix(2, 2);

    tf_qcs_real_number vsin, vcos;
    tf_qcs_complex ein1, ein2, evalue1, evalue2, efin1, efin2;

    vsin = sin(theta / (tf_qcs_real_number)2.0f);
    vcos = cos(theta / (tf_qcs_real_number)2.0f);

    evalue1.re=0; evalue1.im=0;
    evalue2.re=0; evalue2.im=0;
    ein1.re=0; ein1.im =  phi;
    ein2.re=0; ein2.im = -phi;

    qcs_exp_complex( &ein1, &evalue1);
    qcs_exp_complex( &ein2, &evalue2);

    qcs_complex_mul_tf_qcs_real_number(&evalue1, vcos, &efin1);
    qcs_complex_mul_tf_qcs_real_number(&evalue2, vcos, &efin2);

    (m->m+0)->re =  efin1.re; (m->m+0)->im = efin1.im; (m->m+1)->re = vsin;     (m->m+1)->im = 0;
    (m->m+2)->re = -vsin;     (m->m+2)->im = 0;        (m->m+3)->re = efin2.re; (m->m+3)->im = efin2.im;

    return m;
}


/* build matrices */

DYNAMIC_LIB_DECORATION void qcs_build_qubit_gates_matrix()
{
        zero_reset_gate = qcs_create_matrix(2,2);
                id_gate = qcs_create_matrix(2,2);
               not_gate = qcs_create_matrix(2,2);
           pauli_x_gate = qcs_create_matrix(2,2);
           pauli_y_gate = qcs_create_matrix(2,2);
     pauli_z_gate=qcs_create_matrix(2,2);
     hadamard_gate=qcs_create_matrix(2,2);
     square_root_gate=qcs_create_matrix(2,2);
     x_rot90_gate=qcs_create_matrix(2,2);
     y_rot90_gate=qcs_create_matrix(2,2);
     z_rot90_gate=qcs_create_matrix(2,2);
     minus_x_rot90_gate=qcs_create_matrix(2,2);
     minus_y_rot90_gate=qcs_create_matrix(2,2);
     minus_z_rot90_gate=qcs_create_matrix(2,2);
     rotate_alpha_gate=qcs_create_matrix(2,2);
     rotate_theta_gate=qcs_create_matrix(2,2);
     t_gate=qcs_create_matrix(2,2);
     s_gate=qcs_create_matrix(2,2);
     phase_gate=qcs_create_matrix(2,2);
     phase_f_gate=qcs_create_matrix(2,2);
     phase_m11_gate=qcs_create_matrix(2,2);
     cnot_gate=qcs_create_matrix(4,4);
     toffoli_gate=qcs_create_matrix(8,8);

/* Zero matrix */
    (zero_reset_gate->m+0)->re=0; (zero_reset_gate->m+0)->im=0;     (zero_reset_gate->m+1)->re=0; (zero_reset_gate->m+1)->im=0;
    (zero_reset_gate->m+2)->re=0; (zero_reset_gate->m+2)->im=0;     (zero_reset_gate->m+3)->re=0; (zero_reset_gate->m+3)->im=0;

/* Ident matrix */
    (id_gate->m+0)->re=1; (id_gate->m+0)->im=0;     (id_gate->m+1)->re=0; (id_gate->m+1)->im=0;
    (id_gate->m+2)->re=0; (id_gate->m+2)->im=0;     (id_gate->m+3)->re=1; (id_gate->m+3)->im=0;

/* Not matrix */
    (not_gate->m+0)->re=0; (not_gate->m+0)->im=0;     (not_gate->m+1)->re=1; (not_gate->m+1)->im=0;
    (not_gate->m+2)->re=1; (not_gate->m+2)->im=0;     (not_gate->m+3)->re=0; (not_gate->m+3)->im=0;

/* Pauli X matrix */
    (pauli_x_gate->m+0)->re=0; (pauli_x_gate->m+0)->im=0;     (pauli_x_gate->m+1)->re=1; (pauli_x_gate->m+1)->im=0;
    (pauli_x_gate->m+2)->re=1; (pauli_x_gate->m+2)->im=0;     (pauli_x_gate->m+3)->re=0; (pauli_x_gate->m+3)->im=0;

/* Pauli Y matrix */
    (pauli_y_gate->m+0)->re=0; (pauli_y_gate->m+0)->im=0;     (pauli_y_gate->m+1)->re=0; (pauli_y_gate->m+1)->im=-1;
    (pauli_y_gate->m+2)->re=0; (pauli_y_gate->m+2)->im=1;     (pauli_y_gate->m+3)->re=0; (pauli_y_gate->m+3)->im=0;

/* Pauli Z matrix */
    (pauli_z_gate->m+0)->re=1; (pauli_z_gate->m+0)->im=0;     (pauli_z_gate->m+1)->re=0; (pauli_z_gate->m+1)->im=0;
    (pauli_z_gate->m+2)->re=0; (pauli_z_gate->m+2)->im=0;     (pauli_z_gate->m+3)->re=-1; (pauli_z_gate->m+3)->im=0;

/* Hadamard matrix */
    (hadamard_gate->m+0)->re=0.707107; (hadamard_gate->m+0)->im=0;     (hadamard_gate->m+1)->re=0.707107; (hadamard_gate->m+1)->im=0;
    (hadamard_gate->m+2)->re=0.707107; (hadamard_gate->m+2)->im=0;     (hadamard_gate->m+3)->re=-0.707107; (hadamard_gate->m+3)->im=0;

/* Square root matrix */
    (square_root_gate->m+0)->re=0.5; (square_root_gate->m+0)->im= 0.5;     (square_root_gate->m+1)->re=0.5; (square_root_gate->m+1)->im=-0.5;
    (square_root_gate->m+2)->re=0.5; (square_root_gate->m+2)->im=-0.5;     (square_root_gate->m+3)->re=0.5; (square_root_gate->m+3)->im= 0.5;

/*  X Rotate angle 90 gate matrix */
    (x_rot90_gate->m+0)->re=0.707107; (x_rot90_gate->m+0)->im=0;            (x_rot90_gate->m+1)->re=0;        (x_rot90_gate->m+1)->im=0.707107;
    (x_rot90_gate->m+2)->re=0;        (x_rot90_gate->m+2)->im=0.707107;     (x_rot90_gate->m+3)->re=0.707107; (x_rot90_gate->m+3)->im=0;

/*  Y Rotate angle 90 gate matrix */
    (y_rot90_gate->m+0)->re= 0.707107; (y_rot90_gate->m+0)->im=0;     (y_rot90_gate->m+1)->re=0.707107; (y_rot90_gate->m+1)->im=0;
    (y_rot90_gate->m+2)->re=-0.707107; (y_rot90_gate->m+2)->im=0;     (y_rot90_gate->m+3)->re=0.707107; (y_rot90_gate->m+3)->im=0;

/*  Z Rotate angle 90 gate matrix */
    (z_rot90_gate->m+0)->re=0.707107; (z_rot90_gate->m+0)->im=0.707107;    (z_rot90_gate->m+1)->re=0;        (z_rot90_gate->m+1)->im= 0;
    (z_rot90_gate->m+2)->re=0;        (z_rot90_gate->m+2)->im=0;           (z_rot90_gate->m+3)->re=0.707107; (z_rot90_gate->m+3)->im=-0.707107;

/*  Minus X Rotate angle 90 gate matrix */
    (minus_x_rot90_gate->m+0)->re=0.707107; (minus_x_rot90_gate->m+0)->im= 0;            (minus_x_rot90_gate->m+1)->re=0;        (minus_x_rot90_gate->m+1)->im=-0.707107;
    (minus_x_rot90_gate->m+2)->re=0;        (minus_x_rot90_gate->m+2)->im=-0.707107;     (minus_x_rot90_gate->m+3)->re=0.707107; (minus_x_rot90_gate->m+3)->im=0;

/*  Minus Y Rotate angle 90 gate matrix */
    (minus_y_rot90_gate->m+0)->re= 0.707107; (minus_y_rot90_gate->m+0)->im=0;     (minus_y_rot90_gate->m+1)->re=-0.707107; (minus_y_rot90_gate->m+1)->im=0;
    (minus_y_rot90_gate->m+2)->re= 0.707107; (minus_y_rot90_gate->m+2)->im=0;     (minus_y_rot90_gate->m+3)->re= 0.707107; (minus_y_rot90_gate->m+3)->im=0;

/*  Minus Z Rotate angle 90 gate matrix */
    (minus_z_rot90_gate->m+0)->re=0.707107; (minus_z_rot90_gate->m+0)->im=-0.707107;    (minus_z_rot90_gate->m+1)->re=0;        (minus_z_rot90_gate->m+1)->im=0;
    (minus_z_rot90_gate->m+2)->re=0;        (minus_z_rot90_gate->m+2)->im= 0;           (minus_z_rot90_gate->m+3)->re=0.707107; (minus_z_rot90_gate->m+3)->im=0.707107;

/* Rotate alpha gate matrix */
    (rotate_alpha_gate->m+0)->re=1; (rotate_alpha_gate->m+0)->im=0;     (rotate_alpha_gate->m+1)->re=0; (rotate_alpha_gate->m+1)->im=0;
    (rotate_alpha_gate->m+2)->re=0; (rotate_alpha_gate->m+2)->im=0;     (rotate_alpha_gate->m+3)->re=0; (rotate_alpha_gate->m+3)->im=0;

/* Rotate theta gate matrix */
    (rotate_alpha_gate->m+0)->re=1; (rotate_alpha_gate->m+0)->im=0;     (rotate_alpha_gate->m+1)->re=0; (rotate_alpha_gate->m+1)->im=0;
    (rotate_alpha_gate->m+2)->re=0; (rotate_alpha_gate->m+2)->im=0;     (rotate_alpha_gate->m+3)->re=1; (rotate_alpha_gate->m+3)->im=0;

/* T gate matrix */
    (t_gate->m+0)->re=1; (t_gate->m+0)->im=0;     (t_gate->m+1)->re=0;        (t_gate->m+1)->im=0;
    (t_gate->m+2)->re=0; (t_gate->m+2)->im=0;     (t_gate->m+3)->re=0.707107; (t_gate->m+3)->im=0.707107;

/* S gate matrix or V gate matrix*/
    (s_gate->m+0)->re=1; (s_gate->m+0)->im=0;     (s_gate->m+1)->re=0; (s_gate->m+1)->im=0;
    (s_gate->m+2)->re=0; (s_gate->m+2)->im=0;     (s_gate->m+3)->re=0; (s_gate->m+3)->im=1;

/* Phase gate matrix */
    (phase_gate->m+0)->re=1; (phase_gate->m+0)->im=0;     (phase_gate->m+1)->re=0; (phase_gate->m+1)->im=0;
    (phase_gate->m+2)->re=0; (phase_gate->m+2)->im=0;     (phase_gate->m+3)->re=0; (phase_gate->m+3)->im=1;

/* Phase F gate matrix */
    (phase_f_gate->m+0)->re=1; (phase_f_gate->m+0)->im=0;     (phase_f_gate->m+1)->re=0;  (phase_f_gate->m+1)->im=0;
    (phase_f_gate->m+2)->re=0; (phase_f_gate->m+2)->im=0;     (phase_f_gate->m+3)->re=-1; (phase_f_gate->m+3)->im=0;

/* Phase m11 gate matrix */
    (phase_m11_gate->m+0)->re=-1; (phase_m11_gate->m+0)->im=0;     (phase_m11_gate->m+1)->re= 0; (phase_m11_gate->m+1)->im=0;
    (phase_m11_gate->m+2)->re= 0; (phase_m11_gate->m+2)->im=0;     (phase_m11_gate->m+3)->re=-1; (phase_m11_gate->m+3)->im=0;


/* CNot gate matrix */
    (cnot_gate->m+ 0)->re=1; (cnot_gate->m+ 0)->im=0; (cnot_gate->m+ 1)->re=0; (cnot_gate->m+ 1)->im=0; (cnot_gate->m+ 2)->re=0; (cnot_gate->m+ 2)->im=0; (cnot_gate->m+ 3)->re=0; (cnot_gate->m+ 3)->im=0;
    (cnot_gate->m+ 4)->re=0; (cnot_gate->m+ 4)->im=0; (cnot_gate->m+ 5)->re=1; (cnot_gate->m+ 5)->im=0; (cnot_gate->m+ 6)->re=0; (cnot_gate->m+ 6)->im=0; (cnot_gate->m+ 7)->re=0; (cnot_gate->m+ 7)->im=0;
    (cnot_gate->m+ 8)->re=0; (cnot_gate->m+ 8)->im=0; (cnot_gate->m+ 9)->re=0; (cnot_gate->m+ 9)->im=0; (cnot_gate->m+10)->re=0; (cnot_gate->m+10)->im=0; (cnot_gate->m+11)->re=1; (cnot_gate->m+11)->im=0;
    (cnot_gate->m+12)->re=0; (cnot_gate->m+12)->im=0; (cnot_gate->m+13)->re=0; (cnot_gate->m+13)->im=0; (cnot_gate->m+14)->re=1; (cnot_gate->m+14)->im=0; (cnot_gate->m+15)->re=0; (cnot_gate->m+15)->im=0;

/* Toffoli gate matrix */
    (toffoli_gate->m+ 0)->re=1; (toffoli_gate->m+ 0)->im=0; (toffoli_gate->m+ 1)->re=0; (toffoli_gate->m+ 1)->im=0; (toffoli_gate->m+ 2)->re=0; (toffoli_gate->m+ 2)->im=0; (toffoli_gate->m+ 3)->re=0; (toffoli_gate->m+ 3)->im=0; (toffoli_gate->m+ 4)->re=0; (toffoli_gate->m+ 4)->im=0; (toffoli_gate->m+ 5)->re=0; (toffoli_gate->m+ 5)->im=0; (toffoli_gate->m+ 6)->re=0; (toffoli_gate->m+ 6)->im=0; (toffoli_gate->m+ 7)->re=0; (toffoli_gate->m+ 7)->im=0;
    (toffoli_gate->m+ 8)->re=0; (toffoli_gate->m+ 8)->im=0; (toffoli_gate->m+ 9)->re=1; (toffoli_gate->m+ 9)->im=0; (toffoli_gate->m+10)->re=0; (toffoli_gate->m+10)->im=0; (toffoli_gate->m+11)->re=0; (toffoli_gate->m+11)->im=0; (toffoli_gate->m+12)->re=0; (toffoli_gate->m+12)->im=0; (toffoli_gate->m+13)->re=0; (toffoli_gate->m+13)->im=0; (toffoli_gate->m+14)->re=0; (toffoli_gate->m+14)->im=0; (toffoli_gate->m+15)->re=0; (toffoli_gate->m+15)->im=0;
    (toffoli_gate->m+16)->re=0; (toffoli_gate->m+16)->im=0; (toffoli_gate->m+17)->re=0; (toffoli_gate->m+17)->im=0; (toffoli_gate->m+18)->re=1; (toffoli_gate->m+18)->im=0; (toffoli_gate->m+19)->re=0; (toffoli_gate->m+19)->im=0; (toffoli_gate->m+20)->re=0; (toffoli_gate->m+20)->im=0; (toffoli_gate->m+21)->re=0; (toffoli_gate->m+21)->im=0; (toffoli_gate->m+22)->re=0; (toffoli_gate->m+22)->im=0; (toffoli_gate->m+23)->re=0; (toffoli_gate->m+23)->im=0;
    (toffoli_gate->m+24)->re=0; (toffoli_gate->m+24)->im=0; (toffoli_gate->m+25)->re=0; (toffoli_gate->m+25)->im=0; (toffoli_gate->m+26)->re=0; (toffoli_gate->m+26)->im=0; (toffoli_gate->m+27)->re=1; (toffoli_gate->m+27)->im=0; (toffoli_gate->m+28)->re=0; (toffoli_gate->m+28)->im=0; (toffoli_gate->m+29)->re=0; (toffoli_gate->m+29)->im=0; (toffoli_gate->m+30)->re=0; (toffoli_gate->m+30)->im=0; (toffoli_gate->m+31)->re=0; (toffoli_gate->m+31)->im=0;
    (toffoli_gate->m+32)->re=0; (toffoli_gate->m+32)->im=0; (toffoli_gate->m+33)->re=0; (toffoli_gate->m+33)->im=0; (toffoli_gate->m+34)->re=0; (toffoli_gate->m+34)->im=0; (toffoli_gate->m+35)->re=0; (toffoli_gate->m+35)->im=0; (toffoli_gate->m+36)->re=1; (toffoli_gate->m+36)->im=0; (toffoli_gate->m+37)->re=0; (toffoli_gate->m+37)->im=0; (toffoli_gate->m+38)->re=0; (toffoli_gate->m+38)->im=0; (toffoli_gate->m+39)->re=0; (toffoli_gate->m+39)->im=0;
    (toffoli_gate->m+40)->re=0; (toffoli_gate->m+40)->im=0; (toffoli_gate->m+41)->re=0; (toffoli_gate->m+41)->im=0; (toffoli_gate->m+42)->re=0; (toffoli_gate->m+42)->im=0; (toffoli_gate->m+43)->re=0; (toffoli_gate->m+43)->im=0; (toffoli_gate->m+44)->re=0; (toffoli_gate->m+44)->im=0; (toffoli_gate->m+45)->re=1; (toffoli_gate->m+45)->im=0; (toffoli_gate->m+46)->re=0; (toffoli_gate->m+46)->im=0; (toffoli_gate->m+47)->re=0; (toffoli_gate->m+47)->im=0;
    (toffoli_gate->m+48)->re=0; (toffoli_gate->m+48)->im=0; (toffoli_gate->m+49)->re=0; (toffoli_gate->m+49)->im=0; (toffoli_gate->m+50)->re=0; (toffoli_gate->m+50)->im=0; (toffoli_gate->m+51)->re=0; (toffoli_gate->m+51)->im=0; (toffoli_gate->m+52)->re=0; (toffoli_gate->m+52)->im=0; (toffoli_gate->m+53)->re=0; (toffoli_gate->m+53)->im=0; (toffoli_gate->m+54)->re=0; (toffoli_gate->m+54)->im=0; (toffoli_gate->m+55)->re=1; (toffoli_gate->m+55)->im=0;
    (toffoli_gate->m+56)->re=0; (toffoli_gate->m+56)->im=0; (toffoli_gate->m+57)->re=0; (toffoli_gate->m+57)->im=0; (toffoli_gate->m+58)->re=0; (toffoli_gate->m+58)->im=0; (toffoli_gate->m+59)->re=0; (toffoli_gate->m+59)->im=0; (toffoli_gate->m+60)->re=0; (toffoli_gate->m+60)->im=0; (toffoli_gate->m+61)->re=0; (toffoli_gate->m+61)->im=0; (toffoli_gate->m+62)->re=1; (toffoli_gate->m+62)->im=0; (toffoli_gate->m+63)->re=0; (toffoli_gate->m+63)->im=0;

/* qcs_basic_sum_matrix(id_gate,not_gate, cnot_gate); */

/* cache construction */

   qudit_gates_cache = qcs_hash_table_create( 1024, NULL);
}


DYNAMIC_LIB_DECORATION void qcs_destroy_qubit_gates_matrix()
{

     qcs_hash_table_destroy( qudit_gates_cache );

     qcs_delete_matrix(zero_reset_gate);
     qcs_delete_matrix(id_gate);
     qcs_delete_matrix(not_gate);
     qcs_delete_matrix(pauli_x_gate);
     qcs_delete_matrix(pauli_y_gate);
     qcs_delete_matrix(pauli_z_gate);
     qcs_delete_matrix(hadamard_gate);
     qcs_delete_matrix(square_root_gate);
     qcs_delete_matrix(x_rot90_gate);
     qcs_delete_matrix(y_rot90_gate);
     qcs_delete_matrix(z_rot90_gate);
     qcs_delete_matrix(minus_x_rot90_gate);
     qcs_delete_matrix(minus_y_rot90_gate);
     qcs_delete_matrix(minus_z_rot90_gate);
     qcs_delete_matrix(rotate_alpha_gate);
     qcs_delete_matrix(rotate_theta_gate);
     qcs_delete_matrix(t_gate);
     qcs_delete_matrix(s_gate);
     qcs_delete_matrix(phase_gate);
     qcs_delete_matrix(phase_f_gate);
     qcs_delete_matrix(phase_m11_gate);
     qcs_delete_matrix(cnot_gate);
     qcs_delete_matrix(toffoli_gate);
}

/*****/

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_x()
{
    return qcs_clone_matrix( get_pauli_x_gate() );
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_y()
{
    return qcs_clone_matrix( get_pauli_y_gate() );
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_z()
{
    return qcs_clone_matrix( get_pauli_z_gate() );
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_hadamard()
{
    return qcs_clone_matrix( get_hadamard_gate() );
}

/*****/

DYNAMIC_LIB_DECORATION void qcs_1q_id_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{

//     out_qubit->a=in_qubit->a;
//     out_qubit->b=in_qubit->b;

     out_qubit->alpha.re=in_qubit->alpha.re;
     out_qubit->alpha.im=in_qubit->alpha.im;

     out_qubit->beta.re=in_qubit->beta.re;
     out_qubit->beta.im=in_qubit->beta.im;
}

DYNAMIC_LIB_DECORATION void qcs_1q_not_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
//     out_qubit->a=in_qubit->b;
//     out_qubit->b=in_qubit->a;

     out_qubit->beta.re=in_qubit->alpha.re;
     out_qubit->beta.im=in_qubit->alpha.im;

     out_qubit->alpha.re=in_qubit->beta.re;
     out_qubit->alpha.im=in_qubit->beta.im;

}

DYNAMIC_LIB_DECORATION void qcs_1q_pauli_x_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(pauli_x_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(pauli_x_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(pauli_x_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(pauli_x_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_pauli_y_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(pauli_y_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(pauli_y_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(pauli_y_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(pauli_y_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_pauli_z_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(pauli_z_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(pauli_z_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(pauli_z_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(pauli_z_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_hadamard_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(hadamard_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(hadamard_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(hadamard_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(hadamard_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_square_root_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(square_root_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(square_root_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(square_root_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(square_root_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}


DYNAMIC_LIB_DECORATION void qcs_1q_x_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(x_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(x_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(x_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(x_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_y_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(y_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(y_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(y_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(y_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_z_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(z_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(z_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(z_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(z_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}


DYNAMIC_LIB_DECORATION void qcs_1q_minus_x_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(minus_x_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(minus_x_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(minus_x_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(minus_x_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_minus_y_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(minus_y_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(minus_y_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(minus_y_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(minus_y_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_minus_z_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(minus_z_rot90_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(minus_z_rot90_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(minus_z_rot90_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(minus_z_rot90_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_rotate_alpha_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number k, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex a, t1, t2, t3, t4;

     a.re=0;
     a.im=(2.0*QCS_PI) / pow(2.0, k);

     qcs_exp_complex(&a, rotate_alpha_gate->m+3);

     qcs_complex_mul(rotate_alpha_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(rotate_alpha_gate->m+1, &in_qubit->beta,  &t2);

     qcs_complex_mul(rotate_alpha_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(rotate_alpha_gate->m+3, &in_qubit->beta,  &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_rotate_theta_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number theta, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex a, b, c, d, t1, t2, t3, t4;

     a.re= cos(theta); a.im=0; b.re=-sin(theta); b.im=0;
     c.re= sin(theta); c.im=0; d.re= cos(theta); d.im=0;

     (rotate_theta_gate->m+0)->re=a.re;
     (rotate_theta_gate->m+0)->im=a.im;

     (rotate_theta_gate->m+1)->re=b.re;
     (rotate_theta_gate->m+1)->im=b.im;

     (rotate_theta_gate->m+2)->re=c.re;
     (rotate_theta_gate->m+2)->im=c.im;

     (rotate_theta_gate->m+3)->re=d.re;
     (rotate_theta_gate->m+3)->im=d.im;

     qcs_complex_mul(rotate_theta_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(rotate_theta_gate->m+1, &in_qubit->beta,  &t2);

     qcs_complex_mul(rotate_theta_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(rotate_theta_gate->m+3, &in_qubit->beta,  &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}


DYNAMIC_LIB_DECORATION void qcs_1q_t_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(t_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(t_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(t_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(t_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_s_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(s_gate->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(s_gate->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(s_gate->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(s_gate->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_phase_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number k, tf_qcs_qubit *out_qubit)
{
    tf_qcs_complex a, t1,t2,t3,t4;

    a.re=0;
//    a.im=(2*MY_PI)/pow(2,k);
    a.im=k;
    qcs_exp_complex(&a, phase_gate->m+3);
    qcs_exp_complex(&a, phase_gate->m+0);

    qcs_complex_mul(phase_gate->m+0, &in_qubit->alpha, &t1);
    qcs_complex_mul(phase_gate->m+1, &in_qubit->beta, &t2);

    qcs_complex_mul(phase_gate->m+2, &in_qubit->alpha, &t3);
    qcs_complex_mul(phase_gate->m+3, &in_qubit->beta, &t4);

    qcs_complex_add(&t1, &t2, &out_qubit->alpha);
    qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_phase_f_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
    tf_qcs_complex a, t1,t2,t3,t4;

    qcs_complex_mul(phase_f_gate->m+0, &in_qubit->alpha, &t1);
    qcs_complex_mul(phase_f_gate->m+1, &in_qubit->beta, &t2);

    qcs_complex_mul(phase_f_gate->m+2, &in_qubit->alpha, &t3);
    qcs_complex_mul(phase_f_gate->m+3, &in_qubit->beta, &t4);

    qcs_complex_add(&t1, &t2, &out_qubit->alpha);
    qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION void qcs_1q_arbitrary_gate(tf_qcs_matrix *matrix, tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit)
{
     tf_qcs_complex t1,t2,t3,t4;

     qcs_complex_mul(matrix->m+0, &in_qubit->alpha, &t1);
     qcs_complex_mul(matrix->m+1, &in_qubit->beta, &t2);

     qcs_complex_mul(matrix->m+2, &in_qubit->alpha, &t3);
     qcs_complex_mul(matrix->m+3, &in_qubit->beta, &t4);

     qcs_complex_add(&t1, &t2, &out_qubit->alpha);
     qcs_complex_add(&t3, &t4, &out_qubit->beta);
}

DYNAMIC_LIB_DECORATION int qcs_2q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2)
{
/* qubit a znajduje siê w stanie |1> wiec zmienamy stan qubitu b */

     if(abs(in_qubit_1->beta.re)==1)
     {
        qcs_1q_not_gate(in_qubit_2, out_qubit_2);

//        printf("|1> beta re\n");

        return 0;

     }
     if(abs(in_qubit_1->beta.im)==1)
     {
        qcs_1q_not_gate(in_qubit_2, out_qubit_2);

//        printf("|1> beta im\n");

        return 0;
     }

/* qubit a znajduje siê w stanie |0> wiec nie zmienamy stanu qubitu b */

     if(in_qubit_1->alpha.re==1)
     {
        qcs_1q_id_gate(in_qubit_2, out_qubit_2);

//        printf("|0> alpha re\n");

        return 0;
     }
     if(in_qubit_1->alpha.im==1)
     {

//        printf("|0> alpha im\n");

        qcs_1q_id_gate(in_qubit_2, out_qubit_2);

        return 0;
     }

/* qubit a nie znajduje siê w stanie bazowym */
/* wprowadzenie splatania nie powinno wystepowac w trybie pqc  */

//     printf("|others>\n");
     qcs_1q_id_gate(in_qubit_2, out_qubit_2);

     return 1;
}

DYNAMIC_LIB_DECORATION int qcs_3q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3,
                     tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3)
{
     if(abs(in_qubit_1->beta.re)==1 && abs(in_qubit_2->beta.re)==1)
     {
        qcs_1q_not_gate(in_qubit_3, out_qubit_3);

//        printf("|1> beta re\n");

        return 0;

     }
     if(abs(in_qubit_1->beta.im)==1 && abs(in_qubit_2->beta.im)==1)
     {
        qcs_1q_not_gate(in_qubit_3, out_qubit_3);

//        printf("|1> beta im\n");

        return 0;
     }
     if(in_qubit_1->alpha.re==1 && in_qubit_2->alpha.re==1)
     {
        qcs_1q_id_gate(in_qubit_3, out_qubit_3);

//        printf("|0> alpha re\n");

        return 0;
     }
     if(in_qubit_1->alpha.im==1 && in_qubit_2->alpha.im==1)
     {

//        printf("|0> alpha im\n");

        qcs_1q_id_gate(in_qubit_3, out_qubit_3);

        return 0;
     }

/* qubit a nie znajduje siê w stanie bazowym */
/* wprowadzenie splatania nalezy dopracowac  */

//     printf("|others>\n");
     qcs_1q_id_gate(in_qubit_3, out_qubit_3);

     return 1;
}

DYNAMIC_LIB_DECORATION int qcs_4q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3, tf_qcs_qubit *in_qubit_4,
                     tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3, tf_qcs_qubit *out_qubit_4)
{
     if(abs(in_qubit_1->beta.re)==1 && abs(in_qubit_2->beta.re)==1 && abs(in_qubit_3->beta.re)==1)
     {
        qcs_1q_not_gate(in_qubit_4, out_qubit_4);

//        printf("|1> beta re\n");

        return 0;

     }
     if(abs(in_qubit_1->beta.im)==1 && abs(in_qubit_2->beta.im)==1 && abs(in_qubit_3->beta.im)==1)
     {
        qcs_1q_not_gate(in_qubit_4, out_qubit_4);

//        printf("|1> beta im\n");

        return 0;
     }
     if(abs(in_qubit_1->alpha.re)==1 && abs(in_qubit_2->alpha.re)==1 && abs(in_qubit_3->alpha.re)==1)
     {
        qcs_1q_id_gate(in_qubit_4, out_qubit_4);

//        printf("|0> alpha re\n");

        return 0;
     }
     if(abs(in_qubit_1->alpha.im)==1 && abs(in_qubit_2->alpha.im)==1 && abs(in_qubit_3->alpha.im)==1)
     {

//        printf("|0> alpha im\n");

        qcs_1q_id_gate(in_qubit_4, out_qubit_4);

        return 0;
     }

/* qubit a nie znajduje siê w stanie bazowym */
/* wprowadzenie splatania nalezy dopracowac  */

//     printf("|others>\n");
     qcs_1q_id_gate(in_qubit_4, out_qubit_4);

     return 1;
}

/*
void qcs_2q_swap_gate(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2)
{
}

void qcs_2q_cphase_gate(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2)
{
}

void qcs_2q_arbitrary_gate(tf_qcs_matrix4x4 *matrix, tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2)
{
}

void qcs_3q_toffoli_gate(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3)
{
}

void qcs_3q_fredkin_gate(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3)
{
}

void qcs_3q_arbitrary_gate(tf_qcs_matrix8x8 *matrix, tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3)
{
}
*/
/*

*/

DYNAMIC_LIB_DECORATION tf_qcs_matrix *make_matrix_for_one_qubit(char *gate_type, int n, int t)
{
    tf_qcs_matrix *u=NULL, *gate;

    int i, j, minimatrix=0, step=0;
    char bins[128];

    t++; // poprawka dotyczaca numeru qubitu

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    minimatrix=pow(2,n-1);
    step=pow(2,n-t)-1;
    gate=give_qubit_matrix(gate_type);

//     printf("smc=%d, s=%d\n", minimatrix, step);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;

        if (bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        (gate->m+0));
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        (gate->m+1));
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, (gate->m+2));
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, (gate->m+3));
            j++;
        }
    }

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* make_arbitrary_matrix_for_one_qubit(tf_qcs_matrix *gate, int n, int t)
{
    tf_qcs_matrix *u=NULL;

    int i, j, minimatrix=0, step=0;
    char bins[128];

    t++; // poprawka dotyczaca numeru qubitu

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    minimatrix=pow(2,n-1);
    step=pow(2,n-t)-1;

//     printf("smc=%d, s=%d\n", minimatrix, step);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;

        if (bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        (gate->m+0));
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        (gate->m+1));
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, (gate->m+2));
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, (gate->m+3));
            j++;
        }
    }

    return u;
}

DYNAMIC_LIB_DECORATION void add_one_qubit_gate_to_operation_matrix(tf_qcs_matrix *u, tf_qcs_matrix *gate, int n, int t)
{
    tf_qcs_complex t1, t2, t3, t4,
                   u1, u2, u3, u4,
                   o1, o2, o3, o4;

    int i, j, minimatrix=0, step=0;
    char bins[128];

    t++; // poprawka dotyczaca numeru qubitu

    minimatrix=pow(2,n-1);
    step=pow(2,n-t)-1;

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;

        if (bins[t-1]=='0' && j<=minimatrix)
        {

            t1 = *(gate->m+0);
            t2 = *(gate->m+1);
            t3 = *(gate->m+2);
            t4 = *(gate->m+3);

            u1 = *qcs_get_cell_at_matrix_complex(u, i,        i);
            u2 = *qcs_get_cell_at_matrix_complex(u, i,        i+1+step);
            u3 = *qcs_get_cell_at_matrix_complex(u, i+1+step, i);
            u4 = *qcs_get_cell_at_matrix_complex(u, i+1+step, i+1+step);

            o1.re=0; o1.im=0;
            o2.re=0; o2.im=0;
            o3.re=0; o3.im=0;
            o4.re=0; o4.im=0;

            qcs_complex_add(&u1, &t1, &o1);
            qcs_complex_add(&u2, &t2, &o2);
            qcs_complex_add(&u3, &t3, &o3);
            qcs_complex_add(&u4, &t4, &o4);

            qcs_set_cell_at_matrix_complex(u, i,        i,        &o1);
            qcs_set_cell_at_matrix_complex(u, i,        i+1+step, &o2);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i,        &o3);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &o4);

            j++;
        }
    }
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *crot45_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t)
{
    tf_qcs_complex tmp;
    tf_qcs_matrix *u=NULL;

    tmp.re=0;
    tmp.im=(2.0 * QCS_PI) / pow(2.0, 3);

    qcs_exp_complex(&tmp, get_rotate_alpha_gate()->m+3);

    u=two_qubit_syntesis_u_matrix_one_control_one_target(n, c1, t,
            *(get_rotate_alpha_gate()->m+0), *(get_rotate_alpha_gate()->m+1),
            *(get_rotate_alpha_gate()->m+2), *(get_rotate_alpha_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *crot90_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t)
{
    tf_qcs_complex tmp;
    tf_qcs_matrix *u=NULL;

    tmp.re=0;
    tmp.im=(2.0 * QCS_PI) / pow(2.0, 2);

    qcs_exp_complex(&tmp, get_rotate_alpha_gate()->m+3);

    u=two_qubit_syntesis_u_matrix_one_control_one_target(n, c1, t,
            *(get_rotate_alpha_gate()->m+0), *(get_rotate_alpha_gate()->m+1),
            *(get_rotate_alpha_gate()->m+2), *(get_rotate_alpha_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *crot_alpha_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t, tf_qcs_real_number k)
{
    tf_qcs_complex tmp;
    tf_qcs_matrix *u=NULL;

    tmp.re=0;
    tmp.im=(2.0 * QCS_PI) / pow(2.0, k);

    qcs_exp_complex(&tmp, get_rotate_alpha_gate()->m+3);

    u=two_qubit_syntesis_u_matrix_one_control_one_target(n, c1, t,
            *(get_rotate_alpha_gate()->m+0), *(get_rotate_alpha_gate()->m+1),
            *(get_rotate_alpha_gate()->m+2), *(get_rotate_alpha_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t)
{
    tf_qcs_matrix *u=NULL;

    u=two_qubit_syntesis_u_matrix_one_control_one_target(n, c1, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_two_control_one_target(int n, int c1, int c2, int t)
{
    tf_qcs_matrix *u=NULL;

    u=three_qubit_syntesis_u_matrix_two_control_one_target(n, c1, c2, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_three_control_one_target(int n, int c1, int c2, int c3, int t)
{
    tf_qcs_matrix *u=NULL;

    u=four_qubit_syntesis_u_matrix_three_control_one_target(n, c1, c2, c3, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_four_control_one_target(int n, int c1, int c2, int c3, int c4, int t)
{
    tf_qcs_matrix *u=NULL;

    u=five_qubit_syntesis_u_matrix_four_control_one_target(n, c1, c2, c3, c4, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_five_control_one_target(int n, int c1, int c2, int c3, int c4, int c5, int t)
{
    tf_qcs_matrix *u=NULL;

    u=six_qubit_syntesis_u_matrix_five_control_one_target(n, c1, c2, c3, c4, c5, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_one_control_one_target_zero_control(int n, int c1, int t)
{
    tf_qcs_matrix *u=NULL;

    u=two_qubit_syntesis_u_matrix_one_control_one_target_zero_control(n, c1, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_two_control_one_target_zero_control(int n, int c1, int c2, int t)
{
    tf_qcs_matrix *u=NULL;

    u=three_qubit_syntesis_u_matrix_two_control_one_target_zero_control(n, c1, c2, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_three_control_one_target_zero_control(int n, int c1, int c2, int c3, int t)
{
    tf_qcs_matrix *u=NULL;

    u=four_qubit_syntesis_u_matrix_three_control_one_target_zero_control(n, c1, c2, c3, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_four_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int t)
{
    tf_qcs_matrix *u=NULL;

    u=five_qubit_syntesis_u_matrix_four_control_one_target_zero_control(n, c1, c2, c3, c4, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_five_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int c5, int t)
{
    tf_qcs_matrix *u=NULL;

    u=six_qubit_syntesis_u_matrix_five_control_one_target_zero_control(n, c1, c2, c3, c4, c5, t,
            *(get_not_gate()->m+0), *(get_not_gate()->m+1),
            *(get_not_gate()->m+2), *(get_not_gate()->m+3));

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *make_arbitrary_matrix_for_two_qubit_gate(int n, int c1, int t, tf_qcs_matrix *u)
{
    tf_qcs_matrix *tmp;

    tmp = two_qubit_syntesis_u_matrix_one_control_one_target(n, c1, t, *(u->m+0), *(u->m+1), *(u->m+2), *(u->m+3));

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *make_arbitrary_matrix_for_three_qubit_gate(int n, int c1, int c2, int t, tf_qcs_matrix *u)
{
    tf_qcs_matrix *tmp;

    tmp = three_qubit_syntesis_u_matrix_two_control_one_target(n, c1, c2, t, *(u->m+0), *(u->m+1), *(u->m+2), *(u->m+3));

    return tmp;
}

/* gate with control state |1> */
DYNAMIC_LIB_DECORATION tf_qcs_matrix *two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
//         n=3,  // iloœæ linii w układzie
//         c=1,  // linia kontrolujaca
//         t=2,  // linia kontrolowana
    step, minimatrix;

    char bins[128];

    c1++; // dwie poprawki na numery qubitów
    t++;

//     printf("c1=%d t=%d\n", c1, t);

    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n)*0.25; // iloœæ ma³ych macierzy dla przypadku dwuqubitowej bramki cnot ()

//     printf("\nstep=%d\nminimatrix=%d\n", step,minimatrix);
//     printf("\ncontrol=%d\ntarget=%d\n", c, t);

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='1' && bins[t-1]=='0' && j<=minimatrix)
        {
            //printf("%*d %s*\n",n,i,&bins[0]);
            qcs_set_cell_at_matrix_complex(u, i,        i, &a);
            qcs_set_cell_at_matrix_complex(u, i,        i+1+step, &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
//        else
//            printf("%*d %s\n",n,i,&bins[0]);
    }
//     qcs_print_matrix_to_file_sqr(u, stdout);
//     qcs_delete_matrix(u);
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *three_qubit_syntesis_u_matrix_two_control_one_target(int n, int c1, int c2, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
//         n=4,  // iloœæ linii w uk³adzie
//        c1=1,  // pierwsza linia kontrolujaca
//        c2=2,  // druga linia kontrolujaca
//         t=3,  // linia kontrolowana
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    t++;


    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n)*0.125; // iloœæ ma³ych macierzy dla przypadku trzyqubitowej bramki cnot

//     printf("\nstep=%d\nminimatrix=%d\n", step,minimatrix);
//     printf("\ncontrol=%d\ntarget=%d\n", c, t);

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='1' && bins[c2-1]=='1' && bins[t-1]=='0' && j<=minimatrix)
        {
            //printf("%*d %s*\n",n,i,&bins[0]);
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
//        else
//            printf("%*d %s\n",n,i,&bins[0]);
    }
//     qcs_print_matrix_to_file_sqr(u, stdout);
//     qcs_delete_matrix(u);
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *four_qubit_syntesis_u_matrix_three_control_one_target(int n, int c1, int c2, int c3, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    t++;


    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-3); // iloœæ ma³ych macierzy dla przypadku czteroqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='1' && bins[c2-1]=='1' && bins[c3-1]=='1' &&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *five_qubit_syntesis_u_matrix_four_control_one_target(int n, int c1, int c2, int c3, int c4, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    c4++;
    t++;


    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-4); // iloœæ ma³ych macierzy dla przypadku piêcoqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='1' && bins[c2-1]=='1' && bins[c3-1]=='1' && bins[c4-1]=='1'&&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *six_qubit_syntesis_u_matrix_five_control_one_target(int n, int c1, int c2, int c3, int c4, int c5, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    c4++;
    c5++;
    t++;

    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-5); // iloœæ ma³ych macierzy dla przypadku szeœcioqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='1' && bins[c2-1]=='1' && bins[c3-1]=='1' && bins[c4-1]=='1' && bins[c5-1]=='1' &&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i,        i, &a);
            qcs_set_cell_at_matrix_complex(u, i,        i+1+step, &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

/* gate with control state |0> */

DYNAMIC_LIB_DECORATION tf_qcs_matrix *two_qubit_syntesis_u_matrix_one_control_one_target_zero_control(int n, int c1, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // dwie poprawki na numery qubitów
    t++;

    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n)*0.25; // iloœæ ma³ych macierzy dla przypadku dwuqubitowej bramki cnot ()

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='0' && bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *three_qubit_syntesis_u_matrix_two_control_one_target_zero_control(int n, int c1, int c2, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    t++;


    step=pow(2,n-t)-1; // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n)*0.125; // iloœæ ma³ych macierzy dla przypadku trzyqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='0' && bins[c2-1]=='0' && bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *four_qubit_syntesis_u_matrix_three_control_one_target_zero_control(int n, int c1, int c2, int c3, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    t++;


    step=pow(2,n-t)-1;     // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-3); // iloœæ ma³ych macierzy dla przypadku czteroqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='0' && bins[c2-1]=='0' && bins[c3-1]=='0' &&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *five_qubit_syntesis_u_matrix_four_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    c4++;
    t++;

    step=pow(2,n-t)-1;     // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-4); // iloœæ ma³ych macierzy dla przypadku piêcoqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='0' && bins[c2-1]=='0' && bins[c3-1]=='0' && bins[c4-1]=='0'&&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i, i,        &a);
            qcs_set_cell_at_matrix_complex(u, i, i+1+step,        &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *six_qubit_syntesis_u_matrix_five_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int c5, int t,
        tf_qcs_complex a, tf_qcs_complex b,
        tf_qcs_complex c, tf_qcs_complex d)
{
    tf_qcs_matrix *u=NULL;

    int i,j,
    step, minimatrix;

    char bins[128];

    c1++; // poprawki na numery qubitów
    c2++;
    c3++;
    c4++;
    c5++;
    t++;

    step=pow(2,n-t)-1;     // odstep pomiedzy elementami malych macierzy
    minimatrix=pow(2,n-5); // iloœæ ma³ych macierzy dla przypadku szeœcioqubitowej bramki cnot

    u=qcs_create_matrix(pow(2,n), pow(2,n));

    for (i=0;i<pow(2,n);i++) qcs_set_cell_at_matrix_direct(u,i,i,1,0);

    j=1;
    for (i=0;i<pow(2,n);i++)
    {
        qcs_dec2bin(i, n, &bins[0]);
        bins[n]=0;
        if (bins[c1-1]=='0' && bins[c2-1]=='0' && bins[c3-1]=='0' && bins[c4-1]=='0' && bins[c5-1]=='0' &&  bins[t-1]=='0' && j<=minimatrix)
        {
            qcs_set_cell_at_matrix_complex(u, i,        i, &a);
            qcs_set_cell_at_matrix_complex(u, i,        i+1+step, &b);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i, &c);
            qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, &d);
            j++;
        }
    }
    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *create_xy_spins_hamiltonian(int n)
{
    int i, sizemat;

    tf_qcs_complex j;
    pf_qcs_matrix h, t, t1, t2;

    sizemat = (int)pow(2, n) ;

    h = qcs_create_matrix( sizemat, sizemat );
    t = qcs_create_matrix( sizemat, sizemat );

    for ( i = 0; i < (n-1) ; i++)
    {
        // X_{i} X{i+1}

        t1 = make_arbitrary_matrix_for_one_qubit(get_pauli_x_gate(), n, i     );
        t2 = make_arbitrary_matrix_for_one_qubit(get_pauli_x_gate(), n, i + 1 );

        qcs_zero_matrix( t );
        qcs_mul_matrix( t1, t2, t);
        qcs_add_matrix(t, h, h);

        qcs_delete_matrix( t1 );
        qcs_delete_matrix( t2 );

        // Y_{i} Y{i+1}

        t1 = make_arbitrary_matrix_for_one_qubit(get_pauli_y_gate(), n, i     ) ;
        t2 = make_arbitrary_matrix_for_one_qubit(get_pauli_y_gate(), n, i + 1 );

        qcs_zero_matrix( t );
        qcs_mul_matrix( t1, t2, t);
        qcs_add_matrix(t, h, h);

        qcs_delete_matrix( t1 );
        qcs_delete_matrix( t2 );

    } // for ( i = 0; i < (n-1) ; i++)
    qcs_delete_matrix( t );

    j.re = 0.5;
    j.im = 0;
    qcs_mul_scalar_matrix( h, &j, h);

    return h;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *create_xy_spins_hamiltonian_with_jn(int n)
{
    int i, sizemat;

    tf_qcs_complex j;
    pf_qcs_matrix h, t, t1, t2;

    sizemat = (int)pow(2, n) ;

    h = qcs_create_matrix( sizemat, sizemat );
    t = qcs_create_matrix( sizemat, sizemat );

    for ( i = 0; i < (n-1) ; i++)
    {
        // X_{i} X{i+1}

        t1 = make_arbitrary_matrix_for_one_qubit(get_pauli_x_gate(), n, i     );
        t2 = make_arbitrary_matrix_for_one_qubit(get_pauli_x_gate(), n, i + 1 );

        qcs_zero_matrix( t );
        qcs_mul_matrix( t1, t2, t);
        j.re = sqrt((i+1)*(n-(i+1))) / 2.0f;
        j.im = 0;
        qcs_mul_scalar_matrix( t, &j, t);
        qcs_add_matrix(t, h, h);

        qcs_delete_matrix( t1 );
        qcs_delete_matrix( t2 );

        // Y_{i} Y{i+1}

        t1 = make_arbitrary_matrix_for_one_qubit(get_pauli_y_gate(), n, i     ) ;
        t2 = make_arbitrary_matrix_for_one_qubit(get_pauli_y_gate(), n, i + 1 );

        qcs_zero_matrix( t );
        qcs_mul_matrix( t1, t2, t);
        j.re = sqrt((i+1)*(n-(i+1))) / 2.0f;
        j.im = 0;
        qcs_mul_scalar_matrix( t, &j, t);
        qcs_add_matrix(t, h, h);

        qcs_delete_matrix( t1 );
        qcs_delete_matrix( t2 );

    } // for ( i = 0; i < (n-1) ; i++)

    qcs_delete_matrix( t );

    return h;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *create_xy_spins_hamiltonian_with_jn_for_qudit(int n,int freedom_level)
{
    int i, sizemat;

    tf_qcs_complex j;
    pf_qcs_matrix h, t, t1, t2;

    sizemat = (int)pow(freedom_level, n) ;

    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer(int n, tf_qcs_real_number t)
{
    tf_qcs_matrix *h;
    tf_qcs_complex tmp;

    h = create_xy_spins_hamiltonian_with_jn( n );

    tmp.re=0;
    tmp.im=-t ;

    qcs_scalar_mul_matrix(h, &tmp);

    qcs_exp_of_matrix( h );

    qcs_chop_matrix( h );

    return h;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer_float_arg(int n, float t)
{
        return qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer(n, t);
}
