 /**************************************************************************
 *   Copyright (C) 2006 -- 2010 by Marek Sawerwain                         *
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

#ifndef __qcs_symbolic_mode_h__
#define __qcs_symbolic_mode_h__

#include <iostream>
#include <sstream>
#include <string>

#include "ginac/ginac.h"

using namespace std;
using namespace GiNaC;

typedef map<string, symbol> t_dictionary;

typedef struct {
	ex alpha;
	ex beta;
} QubitSymbolic;

typedef struct {
    int freedom_level;
	ex *alphas;
} QuditSymbolic;

typedef struct {
	int n;
	QubitSymbolic *arr;
} QubitSymbolicArray;

typedef struct {
	int n;
	QuditSymbolic *arr;
} QuditSymbolicArray;

typedef struct {
	int size;
	int error_level;
	ex *vec_state;
	t_dictionary *dictionary;
} QCS_QubitSymbol_State;

typedef struct {
	int size;
	int freedom_level;
	int error_level;
	ex *vec_state;
	t_dictionary *dictionary;
} QCS_QuditSymbol_State;

typedef QCS_QubitSymbol_State  ts_qcs_qubit_symbol_state;
typedef QCS_QubitSymbol_State* ps_qcs_qubit_symbol_state;

typedef QCS_QuditSymbol_State  ts_qcs_qudit_symbol_state;
typedef QCS_QuditSymbol_State* ps_qcs_qudit_symbol_state;

typedef QubitSymbolic tf_qcs_qubit_symbolic;
typedef QubitSymbolic* pf_qcs_qubit_symbolic;

typedef QuditSymbolic tf_qcs_qudit_symbolic;
typedef QuditSymbolic* pf_qcs_qudit_symbolic;

typedef QubitSymbolicArray tf_qcs_qubit_symbolic_array;
typedef QubitSymbolicArray* pf_qcs_qubit_symbolic_array;

typedef QuditSymbolicArray tf_qcs_qudit_symbolic_array;
typedef QuditSymbolicArray* pf_qcs_qudit_symbolic_array;


/* single qubit */

extern "C" {
const symbol& get_symbol(t_dictionary &dictionary, const string & s);

tf_qcs_qubit_symbolic* qcs_new_qubit_symbolic();
void qcs_delete_qubit_symbolic(tf_qcs_qubit_symbolic *q);

void qcs_set_symbolic_state(tf_qcs_qubit_symbolic *q, char *alpha_symbol_name, char *beta_symbol_name);
void qcs_set_superposition_state_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_set_superposition_state_per_symbolic(tf_qcs_qubit_symbolic *q, tf_qcs_real_number per);
void qcs_set_superposition_state_theta_psi_symbolic(tf_qcs_qubit_symbolic *q, int _theta, int _psi);

void qcs_set_state0_qubit_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_set_state1_qubit_symbolic(tf_qcs_qubit_symbolic *q);

void qcs_printf_qubit_symbolic(tf_qcs_qubit_symbolic *q);
void qcs_printf_qubit_symbolic_sqr(tf_qcs_qubit_symbolic *q);

void qcs_printf_qubit_symbolic_vec(tf_qcs_qubit_symbolic *q);
void qcs_printf_qubit_symbolic_vec_sqr(tf_qcs_qubit_symbolic *q);

/* single qudit */

tf_qcs_qudit_symbolic* qcs_new_qudit_symbolic(int freedom_level);
void qcs_delete_qudit_symbolic(tf_qcs_qudit_symbolic *q);

void qcs_set_qudit_symbolic_state(tf_qcs_qudit_symbolic *q, char *alpha_symbol_names[]);
void qcs_set_random_real_state_qudit_symbolic(tf_qcs_qudit_symbolic *q);
void qcs_set_superposition_state_qudit_symbolic(tf_qcs_qudit_symbolic *q);

void qcs_set_state_n_qudit_symbolic(tf_qcs_qudit_symbolic *q, int n);

void qcs_printf_qudit_symbolic(tf_qcs_qudit_symbolic *q);
void qcs_printf_qudit_symbolic_sqr(tf_qcs_qudit_symbolic *q);

void qcs_printf_qudit_symbolic_vec(tf_qcs_qudit_symbolic *q);
void qcs_printf_qudit_symbolic_vec_sqr(tf_qcs_qudit_symbolic *q);

/* qubit symbolic array */

tf_qcs_qubit_symbolic_array* qcs_new_qubit_symbolic_array(int n);
void qcs_set_qubit_symbolic_array_n(tf_qcs_qubit_symbolic_array *qa, int n, tf_qcs_qubit_symbolic *q);
tf_qcs_qubit_symbolic* qcs_get_qubit_symbolic_array_n(tf_qcs_qubit_symbolic_array *qa, int n);
void qcs_delete_qubit_symbolic_array(tf_qcs_qubit_symbolic_array *qa);

/* qudit symbolic array */

tf_qcs_qudit_symbolic_array* qcs_new_qudit_symbolic_array(int n);
void qcs_set_qudit_symbolic_array_n(tf_qcs_qudit_symbolic_array *qa, int n, tf_qcs_qudit_symbolic *q);
tf_qcs_qudit_symbolic* qcs_get_qudit_symbolic_array_n(tf_qcs_qudit_symbolic_array *qa, int n);
void qcs_delete_qudit_symbolic_array(tf_qcs_qudit_symbolic_array *qa);

/* quantum qubit register */

ts_qcs_qubit_symbol_state* qcs_new_qubit_symbolic_state(int size);
void qcs_delete_qubit_symbolic_state(ts_qcs_qubit_symbol_state *p);

int qcs_qubit_symbolic_state_get_size(ts_qcs_qubit_symbol_state *p);
int qcs_qubit_symbolic_state_get_freedom_level(ts_qcs_qubit_symbol_state *p);

void qcs_qubit_symbolic_state_reset(ts_qcs_qubit_symbol_state *p);
void qcs_qubit_symbolic_state_reset_to_full_symbolic_state(ts_qcs_qubit_symbol_state *p);

void qcs_qubit_symbolic_set_state_dec(ts_qcs_qubit_symbol_state *p, int n);
void qcs_qubit_symbolic_set_stat_from_qubit_symbolic_array(ts_qcs_qubit_symbol_state *qr, tf_qcs_qubit_symbolic_array *qa);

void qcs_qubit_symbolic_state_substitution(ts_qcs_qubit_symbol_state *q_state, char *symbol_name, tf_qcs_real_number value);
void qcs_qubit_symbolic_state_clear_amplitude(ts_qcs_qubit_symbol_state *q_state, int state);
void qcs_qubit_symbolic_state_clear_amplitude_for_given_qubit_and_state(ts_qcs_qubit_symbol_state *q_state, int n, int state);

void applied_1qubit_gate_to_symbolic_quantum_reg(ts_qcs_qubit_symbol_state *q_state, int t, tf_qcs_matrix *u);
void applied_2qubit_gate_to_symbolic_quantum_reg_one_control(ts_qcs_qubit_symbol_state *q_state, int c1, int t, tf_qcs_matrix *u);
void applied_2qubit_gate_to_symbolic_quantum_reg_zero_control(ts_qcs_qubit_symbol_state *q_state, int c1, int t, tf_qcs_matrix *u);

void apply_matrix_to_qubit_symbolic_quantum_reg(ts_qcs_qubit_symbol_state *q_state, tf_qcs_matrix *matrix);

int qcs_qubit_symbolic_reg_measure(ts_qcs_qubit_symbol_state *q_state, int k);
int qcs_qubit_symbolic_reg_measure_force(ts_qcs_qubit_symbol_state *q_state, int k, int force_result);

char* qcs_qubit_symbolic_reg_get_n_amplitude_as_string(ts_qcs_qubit_symbol_state *q, int n);

void qcs_quantum_symbolic_reg_print_dec_full( ts_qcs_qubit_symbol_state *q_state);
void qcs_quantum_symbolic_reg_print_bin_full( ts_qcs_qubit_symbol_state *q_state );

/* quantum qudit register */

ts_qcs_qudit_symbol_state* qcs_new_qudit_symbolic_state(int size, int freedom_level);
void qcs_delete_qudit_symbolic_state(ts_qcs_qudit_symbol_state *p);

int qcs_qudit_symbolic_state_get_size(ts_qcs_qudit_symbol_state *p);
int qcs_qudit_symbolic_state_get_freedom_level(ts_qcs_qudit_symbol_state *p);

void qcs_qudit_symbolic_state_reset(ts_qcs_qudit_symbol_state *p);
void qcs_qudit_symbolic_state_reset_to_full_symbolic_state(ts_qcs_qudit_symbol_state *p);

void qcs_qudit_symbolic_set_state_dec(ts_qcs_qudit_symbol_state *p, int n);
void qcs_qudit_symbolic_set_stat_from_qudit_symbolic_array(ts_qcs_qudit_symbol_state *qr, tf_qcs_qudit_symbolic_array *qa);

void qcs_qudit_symbolic_state_substitution(ts_qcs_qudit_symbol_state *q_state, char *symbol_name, tf_qcs_real_number value);
void qcs_qudit_symbolic_state_clear_amplitude(ts_qcs_qudit_symbol_state *q_state, int state);
void qcs_qudit_symbolic_state_clear_amplitude_for_given_qudit_and_state(ts_qcs_qudit_symbol_state *q_state, int n, int state);

void applied_1qudit_gate_to_qudit_symbolic_quantum_reg(ts_qcs_qudit_symbol_state *q_state, int t, tf_qcs_matrix *u);
void applied_2qudit_gate_to_qudit_symbolic_quantum_reg_one_control(ts_qcs_qudit_symbol_state *q_state, int c1, int t, tf_qcs_matrix *u);

void apply_matrix_to_qudit_symbolic_quantum_reg(ts_qcs_qudit_symbol_state *q_state, tf_qcs_matrix *matrix);

int qcs_qudit_symbolic_reg_measure(ts_qcs_qudit_symbol_state *q_state, int k);
int qcs_qudit_symbolic_reg_measure_force(ts_qcs_qudit_symbol_state *q_state, int k, int force_result);

char* qcs_qudit_symbolic_reg_get_n_amplitude_as_string(ts_qcs_qudit_symbol_state *q, int n);

void qcs_quantum_qudit_symbolic_reg_print_dec_full( ts_qcs_qudit_symbol_state *q_state );
void qcs_quantum_qudit_symbolic_reg_print_d_base_full( ts_qcs_qudit_symbol_state *q_state );
};

#endif
