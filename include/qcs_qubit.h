/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
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

#ifndef __qcs_qubit_h__
#define __qcs_qubit_h__

#include "qcs_matrix.h"
#include "qcs_complex.h"

typedef struct {
    Complex alpha, beta;
} Qubit;

typedef struct {
    int n;
    Qubit *arr;
} QubitArray;

#define STD_BASE 0
#define HADAMARD_BASE 1
#define M_BASE 2
#define PAULI_X 3
#define PAULI_Y 4
#define PAULI_Z 5

typedef struct {
    Complex v0_a, v0_b;
    Complex v1_a, v1_b;

    tf_qcs_matrix *gate_s0;
    tf_qcs_matrix *gate_s1;
} QubitBaseDesc;

typedef Qubit tf_qcs_qubit;
typedef Qubit* pf_qcs_qubit;

typedef QubitArray tf_qcs_qubit_array;
typedef QubitArray* pf_qcs_qubit_array;

typedef QubitBaseDesc tf_qcs_qubit_base_desc;
typedef QubitBaseDesc* pf_qcs_qubit_base_desc;

tf_qcs_qubit* qcs_new_qubit();
void qcs_delete_qubit(tf_qcs_qubit *q);

void qcs_set_base(tf_qcs_qubit_base_desc *base, tf_qcs_complex v0_a, tf_qcs_complex v0_b, tf_qcs_complex v1_a, tf_qcs_complex v1_b);

tf_qcs_qubit_base_desc *qcs_new_std_base();
tf_qcs_qubit_base_desc *qcs_new_hadamard_base();
tf_qcs_qubit_base_desc *qcs_new_pauli_x_base();
tf_qcs_qubit_base_desc *qcs_new_pauli_y_base();
tf_qcs_qubit_base_desc *qcs_new_pauli_z_base();
tf_qcs_qubit_base_desc *qcs_new_m_base(tf_qcs_real_number theta);
tf_qcs_qubit_base_desc *qcs_new_b_base(tf_qcs_real_number theta);

void qcs_make_std_base(tf_qcs_qubit_base_desc *base);
void qcs_make_hadamard_base(tf_qcs_qubit_base_desc *base);
void qcs_make_pauli_x_base(tf_qcs_qubit_base_desc *base);
void qcs_make_pauli_y_base(tf_qcs_qubit_base_desc *base);
void qcs_make_pauli_z_base(tf_qcs_qubit_base_desc *base);
void qcs_make_m_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta);
void qcs_make_b_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta);
//void qcs_update_m_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta); // to delete

void qcs_base_print(tf_qcs_qubit_base_desc *base);

void qcs_set_superposition_state(tf_qcs_qubit *q);
void qcs_set_random_real_state_qubit(tf_qcs_qubit *q);
void qcs_set_superposition_state_per(tf_qcs_qubit *q, tf_qcs_real_number per);
void qcs_set_superposition_state_theta_psi(tf_qcs_qubit *q, int _theta, int _psi);

void qcs_set_state0_qubit(tf_qcs_qubit *q);
void qcs_set_state1_qubit(tf_qcs_qubit *q);

int qcs_qubit_measure(tf_qcs_qubit *q_reg);

void qcs_printf_qubit(tf_qcs_qubit *q);
void qcs_printf_qubit_sqr(tf_qcs_qubit *q);

void qcs_printf_qubit_vec(tf_qcs_qubit *q);
void qcs_printf_qubit_vec_sqr(tf_qcs_qubit *q);

tf_qcs_qubit_array* qcs_new_qubit_array(int n);
void qcs_set_qubit_array_n(tf_qcs_qubit_array *qa, int i, tf_qcs_qubit *q);
tf_qcs_qubit *qcs_get_qubit_array_n(tf_qcs_qubit_array *qa, int n);
void qcs_delete_qubit_array(tf_qcs_qubit_array *qa);

#endif /* __qcs_qubit_h__ */
