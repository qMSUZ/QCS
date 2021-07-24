/***************************************************************************
 *   Copyright (C) 2018, 2019, 2021 by Marek Sawerwain                     *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
 *   Copyright (C) 2005 -- 2016 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
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

#ifndef  __qcs_quantum_register_h__
#define  __qcs_quantum_register_h__

#include "qcs_complex.h"
#include "qcs_matrix.h"

#define USE_NO_STATE_VECTOR    1001

#define USE_STATE_VECTOR_QUBIT 2001
#define USE_STATE_VECTOR_QUDIT 2002
#define USE_DENSITY_MATRIX     2003
#define USE_GRAPH_STATE_DESC   2004
#define USE_CHP_MODE		   2005	
#define USE_ONEWAY_MODEL	   2006	
#define USE_PQC_MODE	  	   2007	

#define USE_SYMBOLIC_STATE_VECTOR_QUBIT 3001
#define USE_SYMBOLIC_STATE_VECTOR_QUDIT 3002

typedef struct {

    int n,    // size of quantum register in qubits/qudits
	vec_state_size,
	fdl,  // freedomlevel, default value is 2 for qubits 
	mode, // mode of quantum register qubit,qudit
	el;   // error level

	Complex *vs; // vector state

} QuantumRegister;

typedef QuantumRegister  tf_qcs_quantum_register;
typedef QuantumRegister* pf_qcs_quantum_register;



tf_qcs_quantum_register* qcs_new_quantum_register(int size);
void qcs_delete_quantum_register(tf_qcs_quantum_register *q_reg);

void qcs_quantum_register_reset_error_level(tf_qcs_quantum_register *q_reg, int v);
void qcs_quantum_register_set_error_level(tf_qcs_quantum_register *q_reg, int v);
int qcs_quantum_register_get_error_level(tf_qcs_quantum_register *q_reg);

void qcs_quantum_register_reset(tf_qcs_quantum_register *q_reg);

void qcs_quantum_register_set_state_dec(tf_qcs_quantum_register *q_reg, int n);
void qcs_quantum_register_set_state_bin(tf_qcs_quantum_register *q_reg, char *state_desc);

void qcs_quantum_register_print_bin(tf_qcs_quantum_register *q_reg);
void qcs_quantum_register_print_bin_full(tf_qcs_quantum_register *q_reg);

void qcs_quantum_register_print_bin_with_prefix(tf_qcs_quantum_register *q_reg, char *prefix);

void qcs_quantum_register_print_dec(tf_qcs_quantum_register *q_reg);

void applied_1q_gate_to_quantum_register ( tf_qcs_quantum_register *q_reg, int t, tf_qcs_matrix *u);
void applied_2q_gate_to_quantum_register_one_control(tf_qcs_quantum_register *q_reg, int c1, int t, tf_qcs_matrix *u);

void qcs_quantum_register_pauli_x_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_pauli_y_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_pauli_z_gate(tf_qcs_quantum_register *q_reg, int i);

void qcs_quantum_register_had_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_had_n_conj_gate(tf_qcs_quantum_register *q_reg, int i);

void qcs_quantum_register_square_root_not_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_x_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_y_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_z_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_mx_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_my_rot90_n(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_mz_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_rotate_alpha_n_gate(tf_qcs_quantum_register *q_reg, tf_qcs_real_number k, int i);
void qcs_quantum_register_rotate_theta_n_gate(tf_qcs_quantum_register *q_reg, tf_qcs_real_number theta, int i);
void qcs_quantum_register_t_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_v_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_s_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_s_adj_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_phase_n_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_phase_f_n_gate(tf_qcs_quantum_register *q_reg, int i);

void qcs_qubit_arbitrary_one_qubit_gate(tf_qcs_quantum_register *q_reg, tf_qcs_matrix *gate, int i);
void qcs_arbitrary_single_gate(tf_qcs_quantum_register *q_reg, tf_qcs_matrix *m, int i);


void qcs_quantum_register_cnot(tf_qcs_quantum_register *q_reg, int c1, int t);
void qcs_quantum_register_cnot_conj(tf_qcs_quantum_register *q_reg, int c1, int t);


int qcs_quantum_register_measure_one_qubit(tf_qcs_quantum_register *q_reg, int k);
int qcs_quantum_register_measure_one_qubit_in_std_base(tf_qcs_quantum_register *q_reg, int k);
int qcs_quantum_register_measure_one_qubit_in_std_base_force(tf_qcs_quantum_register *q_reg, int k, int force_result);

void qcs_quantum_register_probe_one_qubit_in_std_base(tf_qcs_quantum_register *q_reg, int k, tf_qcs_real_number *out_value_0, tf_qcs_real_number *out_value_1);


#endif /* __qcs_quantum_register_h__ */
