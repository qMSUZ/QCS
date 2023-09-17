/***************************************************************************
 *   Copyright (C) 2006 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                                                         *
 *   Part of the Quantum Computing Simulator:                              *
 *   http://code.google.com/p/qcs                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
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

#ifndef __qcs_symbolic_mode_h__
#define __qcs_symbolic_mode_h__

/* single qubit */

void* qcs_new_qubit_symbolic();
void qcs_delete_qubit_symbolic(void *q);

void qcs_set_symbolic_state(void *q, char *alpha_symbol_name, char *beta_symbol_name);
void qcs_set_superposition_state_symbolic(void *q);
void qcs_set_superposition_state_per_symbolic(void *q, float per);
void qcs_set_superposition_state_theta_psi_symbolic(void *q, int _theta, int _psi);

void qcs_set_state0_qubit_symbolic(void *q);
void qcs_set_state1_qubit_symbolic(void *q);

void qcs_printf_qubit_symbolic(void *q);
void qcs_printf_qubit_symbolic_sqr(void *q);

void qcs_printf_qubit_symbolic_vec(void *q);
void qcs_printf_qubit_symbolic_vec_sqr(void *q);

/* single qudit */

void* qcs_new_qudit_symbolic(int freedom_level);
void qcs_delete_qudit_symbolic(void *q);

void qcs_set_qudit_symbolic_state(void *q, char *alpha_symbol_names[]);
void qcs_set_superposition_state_qudit_symbolic(void *q);

void qcs_set_state_n_qudit_symbolic(void *q, int n);

void qcs_printf_qudit_symbolic(void *q);
void qcs_printf_qudit_symbolic_sqr(void *q);

void qcs_printf_qudit_symbolic_vec(void *q);
void qcs_printf_qudit_symbolic_vec_sqr(void *q);

/* qubit symbolic array */

void* qcs_new_qubit_symbolic_array(int n);
void qcs_set_qubit_symbolic_array_n(void *qa, int n, void *q);
void* qcs_get_qubit_symbolic_array_n(void *qa, int n);
void qcs_delete_qubit_symbolic_array(void *qa);

/* qudit symbolic array */

void* qcs_new_qudit_symbolic_array(int n);
void qcs_set_qudit_symbolic_array_n(void *qa, int n, void *q);
void* qcs_get_qudit_symbolic_array_n(void *qa, int n);
void qcs_delete_qudit_symbolic_array(void *qa);

/* quantum register */

void* qcs_new_qubit_symbolic_state(int size);
void qcs_delete_qubit_symbolic_state(void *p);

void qcs_qubit_symbolic_state_reset(void *p);
void qcs_qubit_symbolic_state_reset_to_full_symbolic_state(void *p);

void qcs_qubit_symbolic_set_state_dec(void *p, int n);
void qcs_qubit_symbolic_set_stat_from_qubit_symbolic_array(void *qr, void *qa);

void qcs_qubit_symbolic_state_substitution(void *q_state, char *symbol_name, float value);

void applied_1qubit_gate_to_symbolic_quantum_reg(void *q_state, int t, void *u);
void applied_2qudit_gate_to_symbolic_quantum_reg_one_control(void *q_state, int c1, int t, void *u);

void apply_matrix_to_symbolic_quantum_reg(void *q_state, void *u);

int qcs_qubit_symbolic_reg_measure(void *q_state, int k);
int qcs_qubit_symbolic_reg_measure_force(void *q_state, int k, int force_result);

char* qcs_qubit_symbolic_reg_get_n_amplitude_as_string(void *q, int n);

void qcs_quantum_symbolic_reg_print_dec_full( void *q_state);
void qcs_quantum_symbolic_reg_print_bin_full( void *q_state );

/* quantum qudit register */

void* qcs_new_qudit_symbolic_state(int size, int freedom_level);
void qcs_delete_qudit_symbolic_state(void *p);

void qcs_qudit_symbolic_state_reset(void *p);
void qcs_qudit_symbolic_state_reset_to_full_symbolic_state(void *p);

void qcs_qudit_symbolic_set_state_dec(void *p, int n);
void qcs_qudit_symbolic_set_stat_from_qudit_symbolic_array(void *qr, void *qa);

void qcs_qudit_symbolic_state_substitution(void *q_state, char *symbol_name, float value);

void applied_1qudit_gate_to_qudit_symbolic_quantum_reg(void *q_state, int t, void *u);
void applied_2qudit_gate_to_qudit_symbolic_quantum_reg_one_control(void *q_state, int c1, int t, void *u);

void apply_matrix_to_qudit_symbolic_quantum_reg(void *q_state, void *u);

int qcs_qudit_symbolic_reg_measure(void *q_state, int k);
int qcs_qudit_symbolic_reg_measure_force(void *q_state, int k, int force_result);

char* qcs_qudit_symbolic_reg_get_n_amplitude_as_string(void *q, int n);

void qcs_quantum_qudit_symbolic_reg_print_dec_full( void *q_state);
void qcs_quantum_qudit_symbolic_reg_print_bin_full( void *q_state );

/* single qubit */

/* tf_qcs_qubit_symbolic* qcs_new_qubit_symbolic();
void qcs_delete_qubit_symbolic(tf_qcs_qubit_symbolic *q);

void qcs_set_symbolic_state(tf_qcs_qubit_symbolic *q, char *alpha_symbol_name, char *beta_symbol_name);
void qcs_set_superposition_state_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_set_superposition_state_per_symbolic(tf_qcs_qubit_symbolic *q, float per);
void qcs_set_superposition_state_theta_psi_symbolic(tf_qcs_qubit_symbolic *q, int _theta, int _psi);

void qcs_set_state0_qubit_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_set_state1_qubit_symbolic(tf_qcs_qubit_symbolic *q);

void qcs_printf_qubit_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_printf_qubit_symbolic_sqr(tf_qcs_qubit_symbolic *q);

void qcs_printf_qubit_symbolic_vec(tf_qcs_qubit_symbolic *q);
void qcs_printf_qubit_symbolic_vec_sqr(tf_qcs_qubit_symbolic *q); */

/* qubit symbolic array */

/* tf_qcs_qubit_symbolic_array* qcs_new_qubit_symbolic_array(int n);
void qcs_set_qubit_symbolic_array_n(tf_qcs_qubit_symbolic_array *qa, int n, tf_qcs_qubit_symbolic *q);
tf_qcs_qubit_symbolic *qcs_get_qubit_symbolic_array_n(tf_qcs_qubit_symbolic_array *qa, int n);
void qcs_delete_qubit_symbolic_array(tf_qcs_qubit_symbolic_array *qa); */

/* quantum register */

/* ts_qcs_qubit_symbol_state *qcs_new_qubit_symbolic_state(int size);
void qcs_delete_qubit_symbolic_state(ts_qcs_qubit_symbol_state *p);

void qcs_qubit_symbolic_state_reset(ts_qcs_qubit_symbol_state *p);
void qcs_qubit_symbolic_state_reset_to_full_symbolic_state(ts_qcs_qubit_symbol_state *p);

void qcs_qubit_symbolic_set_state_dec(ts_qcs_qubit_symbol_state *p, int n);
void qcs_qubit_symbolic_set_stat_from_qubit_symbolic_array(ts_qcs_qubit_symbol_state *qr, tf_qcs_qubit_symbolic_array *qa);

void qcs_qubit_symbolic_state_substitution(ts_qcs_qubit_symbol_state *q_state, char *symbol_name, float value);

void applied_1qubit_gate_to_symbolic_quantum_reg(ts_qcs_qubit_symbol_state *q_state, int t, tf_qcs_matrix *u);
void applied_2qudit_gate_to_symbolic_quantum_reg_one_control(ts_qcs_qubit_symbol_state *q_state, int c1, int t, tf_qcs_matrix *u);

int qcs_qubit_symbolic_reg_measure(ts_qcs_qubit_symbol_state *q_state, int k);
int qcs_qubit_symbolic_reg_measure_force(ts_qcs_qubit_symbol_state *q_state, int k, int force_result);

void qcs_quantum_symbolic_reg_print_dec_full( ts_qcs_qubit_symbol_state *q_state);
void qcs_quantum_symbolic_reg_print_bin_full( ts_qcs_qubit_symbol_state *q_state ); */

#endif
