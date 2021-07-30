/***************************************************************************
 *   Copyright (C) 2021 by Marek Sawerwain                                 *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
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

void spectral_decomposition_test()
{
    tf_qcs_quantum_register *r = NULL;

    tf_qcs_matrix *denmat = NULL,
				  *eigenvectors = NULL,
	 			  *eigenvalues = NULL;

    r = qcs_new_quantum_register ( 2 );
    qcs_quantum_register_reset ( r );

    qcs_quantum_register_had_n_gate(r, 0);
    qcs_quantum_register_cnot(r, 0, 1);
	printf("Bell state: 1/(√2) |00> + |11>:\n");
	qcs_quantum_register_print_bin_with_prefix(r,"    ");

    denmat = qcs_quantum_register_generate_density_matrix(r);
    eigenvectors = qcs_create_matrix(denmat->rows, denmat->cols);
    eigenvalues = qcs_create_matrix(denmat->rows, 1);

	printf("Density matrix\n");
	qcs_print_matrix( denmat );
	
	qcs_spectral_decompose_of_matrix(denmat, eigenvalues, eigenvectors);

	printf("Eigenvectors\n");
	qcs_print_matrix( eigenvectors );
	
	printf("Eigenvalues\n");
	qcs_print_matrix( eigenvalues );

	qcs_delete_matrix( eigenvalues );
	qcs_delete_matrix( eigenvectors );
	qcs_delete_matrix( denmat );

    qcs_delete_quantum_register ( r ) ;
}

int main( int argc, char *argv[] )
{
	qcs_core_library_initialization ( );


	printf(": qcs spectral decomposition test\n");
    spectral_decomposition_test();


	qcs_core_library_deinitialization();

	return 0;
}