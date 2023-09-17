/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2007 -- 2008 by Przemysław Ratajczak                    *
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

#ifndef __qcs_qudit_h__
#define __qcs_qudit_h__

#include "qcs_complex.h"
#include "qcs_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int freedom_level;
    Complex *alphas;
} Qudit;

typedef struct {
    int n;
    Qudit *arr;
} QuditArray;

typedef struct {
    int size;
    int freedom_level;
    Complex *vec_state;
} QCS_Qudit_State;

typedef QCS_Qudit_State ts_qcs_qudit_state;
typedef QCS_Qudit_State* ps_qcs_qudit_state;

typedef Qudit tf_qcs_qudit;
typedef Qudit* pf_qcs_qudit;

typedef QuditArray tf_qcs_qudit_array;
typedef QuditArray* pf_qcs_qudit_array;

QCS_Qudit_State *qcs_new_qudit_state(int size, int freedom_level);
void qcs_delete_qudit_state(ts_qcs_qudit_state *p);

void qcs_qudit_state_reset(ts_qcs_qudit_state *p);
void qcs_qudit_state_fill_zero(ts_qcs_qudit_state *p);
void qcs_qudit_set_state_dec(ts_qcs_qudit_state *p, int n);

/* optimised function to applied matrix to vector state -- qudit version */

void applied_1qudit_gate_to_quantum_reg(ts_qcs_qudit_state *q_state, int qn, tf_qcs_matrix *u);
void applied_1qudit_gate_to_quantum_reg_FASTER(ts_qcs_qudit_state *q_state, int t, tf_qcs_matrix *u);

/* gate with control state |a> for qudit d-level version where "a" fulfill assumption 0 <= a <= d - 1 */

/* gate with control state |1> -- qudit version */

void applied_2qudit_gate_to_quantum_reg(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u);

void applied_2qudit_gate_to_quantum_reg_d_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u);

void applied_2qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u);
void applied_3qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int t, tf_qcs_matrix *u);
void applied_4qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int t, tf_qcs_matrix *u);
void applied_5qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int t, tf_qcs_matrix *u);
void applied_6qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int c5, int t, tf_qcs_matrix *u);

/* gate with control state |0> -- qudit version */

void applied_2qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u);
void applied_3qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int t, tf_qcs_matrix *u);
void applied_4qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int t, tf_qcs_matrix *u);
void applied_5qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int t, tf_qcs_matrix *u);
void applied_6qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int c5, int t, tf_qcs_matrix *u);

void applied_cnot_gate_to_qudit_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t);
void applied_cnot_conj_gate_to_qudit_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t);

/* -------------------------------------------------- */

tf_qcs_qudit* qcs_new_qudit(int freedom_level);
void qcs_delete_qudit(tf_qcs_qudit *q);

void qcs_set_random_real_state_qudit(tf_qcs_qudit *q);
void qcs_set_state_n_qudit(tf_qcs_qudit *q, int s);

void qcs_printf_qudit(tf_qcs_qudit *q);
void qcs_printf_qudit_sqr(tf_qcs_qudit *q);


tf_qcs_qudit_array* qcs_new_qudit_array(int n, int freedom_level);
void qcs_set_qudit_array_n(tf_qcs_qudit_array *qa, int n, tf_qcs_qudit *q);
tf_qcs_qudit *qcs_get_qudit_array_n(tf_qcs_qudit_array *qa, int n);
void qcs_delete_qudit_array(tf_qcs_qudit_array *qa);

#ifdef __cplusplus
}
#endif

#endif /* __qcs_qudit_h__ */
