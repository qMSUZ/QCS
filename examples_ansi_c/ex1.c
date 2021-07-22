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


#include <stdio.h>

#include "qcs.h"


void qcs_qubit_test1()
{
	
 pf_qcs_qubit q1, q2;


 q1=qcs_new_qubit();
 q2=qcs_new_qubit();

 qcs_set_state1_qubit(q1);
 qcs_set_state1_qubit(q2);

 printf("---\n");

 qcs_printf_qubit(q1);
 qcs_printf_qubit(q2);

 printf("---\n");

 qcs_1q_rotate_alpha_gate(q1, 2, q2);

 qcs_print_matrix_to_file(get_rotate_alpha_gate(), stdout);

 printf("---\n");

 qcs_printf_qubit(q1);
 qcs_printf_qubit(q2);
}

void qcs_quantum_reg_test1()
{
	int i;
	pf_qcs_quantum_register r;

	r=qcs_new_quantum_register(2);

	qcs_quantum_register_reset(r);

	for(i=0;i<4;i++)
	{
		qcs_quantum_register_set_state_dec(r, i);

		qcs_quantum_register_print_dec(r);
		qcs_quantum_register_print_bin(r);
	}

	qcs_delete_quantum_register(r);
}

void qcs_quantum_reg_test2()
{
 	int i;
 	pf_qcs_quantum_register r;

 	r=qcs_new_quantum_register(3);

	for(i=0;i<3;i++)
	{
		qcs_quantum_register_reset(r);
		printf("Initial state after reset:\n");
		qcs_quantum_register_print_bin(r);

		printf("After hadamard gate (qubit %d)):\n", i);
		qcs_quantum_register_had_n_gate(r, i);
		qcs_quantum_register_print_bin(r);
	}

 	qcs_delete_quantum_register(r);
}

void qcs_quantum_reg_test3()
{
 	int i;
 	pf_qcs_quantum_register r;

 	r=qcs_new_quantum_register(2);

	qcs_quantum_register_reset(r);
	qcs_quantum_register_had_n_gate(r, 0);
	qcs_quantum_register_cnot(r, 0, 1);

	printf("Bell state: 1/(√2) |00> + |11>:\n");
	qcs_quantum_register_print_bin_with_prefix(r,"    "); // four additionall spaces

 	qcs_delete_quantum_register(r);
}


int main( int argc, char *argv[] )
{
	qcs_core_library_initialization ( );

#if 0
	printf(": qcs_qubit_test1\n");
	qcs_qubit_test1();

	printf(": qcs_quantum_reg_test1\n");
	qcs_quantum_reg_test1();

	printf(": qcs_quantum_reg_test2\n");
	qcs_quantum_reg_test2();
#endif

	printf(": qcs_quantum_reg_test3\n");
	qcs_quantum_reg_test3();

	qcs_core_library_deinitialization();

	return 0;
}