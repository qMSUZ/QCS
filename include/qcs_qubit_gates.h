/***************************************************************************
 *   Copyright (C) 2005 -- 2013 by Marek Sawerwain                         *
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

#ifndef __qcs_qubit_gates_h__
#define __qcs_qubit_gates_h__

#include <math.h>

#include "qcs_matrix.h"
#include "qcs_qubit.h"

#ifdef __cplusplus
extern "C" {
#endif

/* qubit gate generators */

tf_qcs_matrix *get_zero_reset_gate();
tf_qcs_matrix *get_id_gate();
tf_qcs_matrix *get_not_gate();
tf_qcs_matrix *get_pauli_x_gate();
tf_qcs_matrix *get_pauli_y_gate();
tf_qcs_matrix *get_pauli_z_gate();
tf_qcs_matrix *get_hadamard_gate();
tf_qcs_matrix *get_square_root_gate();
tf_qcs_matrix *get_x_rot90_gate();
tf_qcs_matrix *get_y_rot90_gate();
tf_qcs_matrix *get_z_rot90_gate();
tf_qcs_matrix *get_minus_x_rot90_gate();
tf_qcs_matrix *get_minus_y_rot90_gate();
tf_qcs_matrix *get_minus_z_rot90_gate();
tf_qcs_matrix *get_rotate_alpha_gate();
tf_qcs_matrix *get_rotate_theta_gate();
tf_qcs_matrix *get_t_gate();
tf_qcs_matrix *get_s_gate();
tf_qcs_matrix *get_phase_gate();
tf_qcs_matrix *get_phase_f_gate();
tf_qcs_matrix *get_phase_m11_gate();
tf_qcs_matrix *get_cnot_gate();
tf_qcs_matrix *get_toffoli_gate();

tf_qcs_matrix *get_swap_gate();
tf_qcs_matrix *get_fredkin_gate();

tf_qcs_matrix *qcs_rot_x_gate(tf_qcs_real_number theta);
tf_qcs_matrix *qcs_rot_y_gate(tf_qcs_real_number theta);
tf_qcs_matrix *qcs_rot_z_gate(tf_qcs_real_number theta);

tf_qcs_matrix *qcs_create_matrix_for_e_2q_gate_float_arg(float gamma);
tf_qcs_matrix *qcs_create_matrix_for_e_2q_gate(tf_qcs_real_number gamma);

tf_qcs_matrix *qcs_create_matrix_for_ux_1q_gate_float_arg(float phi, float theta);
tf_qcs_matrix *qcs_create_matrix_for_ux_1q_gate(tf_qcs_real_number phi, tf_qcs_real_number theta);

tf_qcs_matrix* give_qubit_matrix(char *arg);

tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_x();
tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_y();
tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_z();
tf_qcs_matrix* qcs_create_matrix_for_qubit_gate_hadamard();


void qcs_1q_id_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_not_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_pauli_x_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_pauli_y_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_pauli_z_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_hadamard_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_square_root_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_x_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_y_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_z_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_minus_x_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_minus_y_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_minus_z_rotate90_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_rotate_alpha_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number k, tf_qcs_qubit *out_qubit);
void qcs_1q_rotate_theta_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number theta, tf_qcs_qubit *out_qubit);

void qcs_1q_t_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);
void qcs_1q_s_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_phase_gate(tf_qcs_qubit *in_qubit, tf_qcs_real_number k, tf_qcs_qubit *out_qubit);
void qcs_1q_phase_f_gate(tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

void qcs_1q_arbitrary_gate(tf_qcs_matrix *matrix, tf_qcs_qubit *in_qubit, tf_qcs_qubit *out_qubit);

/*
	two and more qubit gates for pqc mode
*/

int qcs_2q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2);
int qcs_3q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3,
                        tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3);
int qcs_4q_cnot_gate_pqc(tf_qcs_qubit *in_qubit_1, tf_qcs_qubit *in_qubit_2, tf_qcs_qubit *in_qubit_3, tf_qcs_qubit *in_qubit_4,
                    tf_qcs_qubit *out_qubit_1, tf_qcs_qubit *out_qubit_2, tf_qcs_qubit *out_qubit_3, tf_qcs_qubit *out_qubit_4);


tf_qcs_matrix *make_matrix_for_one_qubit(char *gate_type, int n, int t);
tf_qcs_matrix *make_arbitrary_matrix_for_one_qubit(tf_qcs_matrix *gate, int n, int t);

void add_one_qubit_gate_to_operation_matrix(tf_qcs_matrix *u, tf_qcs_matrix *gate, int n, int t);

tf_qcs_matrix *crot45_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t);
tf_qcs_matrix *crot90_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t);
tf_qcs_matrix *crot_alpha_two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t, tf_qcs_real_number k);

tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_two_control_one_target(int n, int c1, int c2, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_three_control_one_target(int n, int c1, int c2, int c3, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_four_control_one_target(int n, int c1, int c2, int c3, int c4, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_five_control_one_target(int n, int c1, int c2, int c3, int c4, int c5, int t);

tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_one_control_one_target_zero_control(int n, int c1, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_two_control_one_target_zero_control(int n, int c1, int c2, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_three_control_one_target_zero_control(int n, int c1, int c2, int c3, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_four_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int t);
tf_qcs_matrix *cnot_qubit_syntesis_u_matrix_five_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int c5, int t);

/* gate with control state |1> */

tf_qcs_matrix *make_arbitrary_matrix_for_two_qubit_gate(int n, int c1, int t, tf_qcs_matrix *u);
tf_qcs_matrix *make_arbitrary_matrix_for_three_qubit_gate(int n, int c1, int c2, int t, tf_qcs_matrix *u);

tf_qcs_matrix *two_qubit_syntesis_u_matrix_one_control_one_target(int n, int c1, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *three_qubit_syntesis_u_matrix_two_control_one_target(int n, int c1, int c2, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *four_qubit_syntesis_u_matrix_three_control_one_target(int n, int c1, int c2, int c3, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *five_qubit_syntesis_u_matrix_four_control_one_target(int n, int c1, int c2, int c3, int c4, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *six_qubit_syntesis_u_matrix_five_control_one_target(int n, int c1, int c2, int c3, int c4, int c5, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

/* gate with control state |0> */

tf_qcs_matrix *two_qubit_syntesis_u_matrix_one_control_one_target_zero_control(int n, int c1, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *three_qubit_syntesis_u_matrix_two_control_one_target_zero_control(int n, int c1, int c2, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *four_qubit_syntesis_u_matrix_three_control_one_target_zero_control(int n, int c1, int c2, int c3, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *five_qubit_syntesis_u_matrix_four_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

tf_qcs_matrix *six_qubit_syntesis_u_matrix_five_control_one_target_zero_control(int n, int c1, int c2, int c3, int c4, int c5, int t,
              tf_qcs_complex a, tf_qcs_complex b,
              tf_qcs_complex c, tf_qcs_complex d);

/* spins hamiltonian */

tf_qcs_matrix *create_xy_spins_hamiltonian(int n);
tf_qcs_matrix *create_xy_spins_hamiltonian_with_jn(int n);
tf_qcs_matrix *create_xy_spins_hamiltonian_with_jn_for_qudit(int n,int freedom_level);

tf_qcs_matrix *qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer(int n, tf_qcs_real_number t);
tf_qcs_matrix *qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer_float_arg(int n, float t);

void qcs_build_qubit_gates_matrix();
void qcs_destroy_qubit_gates_matrix();

#ifdef __cplusplus
}
#endif

#endif /* __qcs_qubit_gates_h__ */

