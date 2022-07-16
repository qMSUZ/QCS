/***************************************************************************
 *   Copyright (C) 2022 by Marek Sawerwain                                 *
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

void matrix_test_no_1()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_complex num;

    num.re=2;
    num.im=4;

    printf(".: matrix_test_no_1 :.\n\n");

    a=qcs_create_matrix(3, 3);
    b=qcs_create_matrix(3, 3);
    c=qcs_create_matrix(3, 3);

    qcs_set_cell_at_matrix_complex(a, 0, 0, &num);
    qcs_set_cell_at_matrix_direct(a, 0, 1, 2, 0);
    qcs_set_cell_at_matrix_direct(a, 0, 2, 3, 0);

    qcs_set_cell_at_matrix_complex(b, 0, 1, &num);
    qcs_set_cell_at_matrix_direct(b, 0, 0, 3, 0);
    qcs_set_cell_at_matrix_direct(b, 0, 2, 3, 0);

    qcs_mul_matrix(a, b, c);

    qcs_print_matrix_to_file(a, stdout);
    printf("\n");
    qcs_print_matrix_to_file(b, stdout);
    printf("\n");
    qcs_print_matrix_to_file(c, stdout);

    printf("in matlab/octave format\n");
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");
    qcs_print_matrix_in_matlab_format_with_prefix(b, "b =");
    printf("\n");
    qcs_print_matrix_in_matlab_format_with_prefix(c, "c =");
    printf("\n");

    qcs_delete_matrix(a);
    qcs_delete_matrix(b);
    qcs_delete_matrix(c);
}

void vector_test_no_1()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_complex num;
    tf_qcs_real_number norm_val=(tf_qcs_real_number)0.0;

    num.re=2;
    num.im=4;

    printf(".: vector_test_no_1 :.\n\n");

    //a=qcs_create_matrix(4, 1);
    a=qcs_create_vector_column(4);
    qcs_set_cell_at_matrix_complex(a, 0, 0, &num);
    qcs_set_cell_at_matrix_direct(a, 1, 0, 2, 0);
    qcs_set_cell_at_matrix_direct(a, 2, 0, 3, 0);
    qcs_set_cell_at_matrix_direct(a, 3, 0, 4, 0);

    
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");

    qcs_calculate_norm1_of_column_vector(a, &norm_val);
    printf("norm 1 of a = %lf\n", (double)norm_val);
    
    qcs_delete_matrix(a);
}

void vector_test_no_2()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_complex num;
    tf_qcs_real_number norm_val=(tf_qcs_real_number)0.0;

    num.re=2;
    num.im=4;

    printf(".: vector_test_no_2 :.\n\n");

    a=qcs_create_random_real_vector_column(4);
    
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");

    qcs_calculate_norm1_of_column_vector(a, &norm_val);
    printf("norm 1 of a = %lf\n", (double)norm_val);

    qcs_calculate_norm2_of_column_vector(a, &norm_val);
    printf("norm 2 of a = %lf\n", (double)norm_val);


    qcs_delete_matrix(a);
}

void vector_test_no_2b()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_complex num;
    tf_qcs_real_number norm_val=(tf_qcs_real_number)0.0;

    num.re=2;
    num.im=4;

    printf(".: vector_test_no_2b :.\n\n");

    a=qcs_create_random_vector_column(4);
    
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");

    qcs_calculate_norm1_of_column_vector(a, &norm_val);
    printf("norm 1 of a = %lf\n", (double)norm_val);
    
    qcs_delete_matrix(a);
}

void vector_test_no_2c()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_real_number norm_val=(tf_qcs_real_number)0.0;

    printf(".: vector_test_no_2c :.\n\n");

    a=qcs_create_random_vector_column_normalised(4);
    
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");

    qcs_calculate_norm1_of_column_vector(a, &norm_val);
    printf("norm 1 of a = %lf\n", (double)norm_val);
    
    qcs_delete_matrix(a);
}


void vector_test_no_3()
{
    tf_qcs_matrix *a=NULL, *b=NULL, *c=NULL;
    tf_qcs_real_number norm_val=(tf_qcs_real_number)0.0;

    printf(".: vector_test_no_3 :.\n\n");

    a=qcs_create_random_real_vector_column_normalised(4);
    
    qcs_print_matrix_in_matlab_format_with_prefix(a, "a =");
    printf("\n");

    qcs_calculate_norm1_of_column_vector(a, &norm_val);
    printf("norm 1 of a = %lf\n", (double)norm_val);
    
    qcs_delete_matrix(a);
}

int main( int argc, char *argv[] )
{
	qcs_core_library_initialization();

    //matrix_test_no_1();

    //vector_test_no_1();
    vector_test_no_2();
    //vector_test_no_2b();
    //vector_test_no_2c();
    //vector_test_no_3();

	qcs_core_library_deinitialization();

	return 0;
}