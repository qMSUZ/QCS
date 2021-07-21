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

#define USE_STATE_VECTOR_QUBIT 2001
#define USE_STATE_VECTOR_QUDIT 2002
#define USE_DENSITY_MATRIX     2003
#define USE_GRAPH_STATE_DESC   2004
#define USE_CHP_MODE		   2005	
#define USE_ONEWAY_MODEL	   2006	

#define USE_SYMBOLIC_STATE_VECTOR_QUBIT 3000
#define USE_SYMBOLIC_STATE_VECTOR_QUDIT 3001

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

void qcs_quantum_register_print_dec(tf_qcs_quantum_register *q_reg);

void qcs_quantum_register_pauli_x_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_pauli_y_gate(tf_qcs_quantum_register *q_reg, int i);
void qcs_quantum_register_pauli_z_gate(tf_qcs_quantum_register *q_reg, int i);


#endif /* __qcs_quantum_register_h__ */
