/***************************************************************************
 *   Copyright (C) 2005 -- 2013 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2007 -- 2008 by Przemyslaw Ratajczak                    *
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

#ifndef __qcs_qudit_gates_h__
#define __qcs_qudit_gates_h__

#ifdef __cplusplus
extern "C" {
#endif

/* qudit gate generators */

tf_qcs_matrix *get_qudit_zero_reset_gate(int freedom_level);
tf_qcs_matrix *get_qudit_id_gate(int freedom_level);
tf_qcs_matrix *get_qudit_not_gate(int freedom_level);
tf_qcs_matrix *get_qudit_pauli_x_gate(int freedom_level);
tf_qcs_matrix *get_qudit_pauli_y_gate(int freedom_level);
tf_qcs_matrix *get_qudit_pauli_z_gate(int freedom_level);
tf_qcs_matrix *get_qudit_hadamard_gate(int freedom_level);
tf_qcs_matrix *get_qudit_square_root_gate(int freedom_level);
tf_qcs_matrix *get_qudit_x_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_y_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_z_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_minus_x_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_minus_y_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_minus_z_rot90_gate(int freedom_level);
tf_qcs_matrix *get_qudit_rotate_alpha_gate(int freedom_level);
tf_qcs_matrix *get_qudit_rotate_theta_gate(int freedom_level);
tf_qcs_matrix *get_qudit_t_gate(int freedom_level);
tf_qcs_matrix *get_qudit_s_gate(int freedom_level);
tf_qcs_matrix *get_qudit_phase_gate(int freedom_level);
tf_qcs_matrix *get_qudit_phase_f_gate(int freedom_level);
tf_qcs_matrix *get_qudit_cnot_gate(int freedom_level);
tf_qcs_matrix *get_qudit_cnot_swaped_gate(int freedom_level); //control & target qubits swaped
tf_qcs_matrix *get_qudit_toffoli_gate(int freedom_level);
tf_qcs_matrix *get_qudit_swap_gate(int freedom_level);

/* different implementation of qudit gates */

tf_qcs_matrix *get_qudit_hadamard_d_gate(int freedom_level);
tf_qcs_matrix *get_qudit_pauli_x_d_gate(int freedom_level);
tf_qcs_matrix *get_qudit_pauli_z_d_gate(int freedom_level);

tf_qcs_matrix *get_qudit_power_of_x_gate(int freedom_level, int power);
tf_qcs_matrix *get_qudit_power_of_z_gate(int freedom_level, int power);
tf_qcs_matrix *get_qudit_power_of_u_gate(tf_qcs_matrix *u, int power);

tf_qcs_matrix *get_qudit_pauli_operator_gate(int freedom_level, int r, int s);


tf_qcs_matrix *make_arbitrary_matrix_for_one_qudit(tf_qcs_matrix *gate, int n, int fd, int t);


/* SU(d) generators */
tf_qcs_matrix *get_general_e_matrix(int fd, int k, int j); /*1 <= k,j <= d*/

tf_qcs_matrix *get_general_theta_matrix(int fd, int k, int j); /* 1 <= k < j <= d */
tf_qcs_matrix *get_general_beta_matrix(int fd, int k, int j); /* 1 <= k < j <= d */
tf_qcs_matrix *get_general_eta_matrix(int fd, int r); /* 1 <= r < d*/


tf_qcs_matrix *qudit_d_level_teleportation_correction_matrix(int n, int m, int d);

tf_qcs_matrix *cnot_qudit_syntesis_one_control_one_target(int n, int d, int c1, int t);
tf_qcs_matrix *cnot_qudit_syntesis_two_control_one_target(int n, int d, int c1, int c2, int t);
tf_qcs_matrix *cnot_qudit_syntesis_three_control_one_target(int n, int d, int c1, int c2, int c3, int t);

/*
	one qubit and qudit gates gates
*/

void qcs_build_qubit_gates_matrix();
void qcs_destroy_qubit_gates_matrix();

#ifdef __cplusplus
}
#endif

#endif /* __qcs_qudit_gates_h__ */
