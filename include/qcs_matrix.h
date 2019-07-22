/***************************************************************************
 *   Copyright (C) 2019 by Marek Sawerwain                                 *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
 *   Copyright (C) 2005 -- 2012 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2007 -- 2008 by Przemysław Ratajczak                    *
 *   Copyright (C) 2005 -- 2006 by Kamil Pawłowski                         *
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

#ifndef __qcs_matrix_h__
#define __qcs_matrix_h__

#include "qcs_misc.h"
#include "qcs_complex.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    tf_qcs_complex *m;
    unsigned int rows, cols;
    unsigned int q, freedom_level;
} Matrix;

typedef struct {
    tf_qcs_complex *m;
    unsigned int w1, w2, w3, w4;
    unsigned int x1, x2, x3, x4;
    unsigned int q, freedom_level;
} Matrix4x4;

typedef struct {
    int *m;
    unsigned int rows, cols;
} IMatrix;

typedef struct {
    int *m;
    unsigned int w1, w2, w3, w4;
    unsigned int x1, x2, x3, x4;
} IMatrix4x4;

typedef Matrix tf_qcs_matrix;
typedef Matrix* pf_qcs_matrix;

typedef Matrix4x4 tf_qcs_matrix4x4;
typedef Matrix4x4* pf_qcs_matrix4x4;

typedef IMatrix ti_qcs_matrix;
typedef IMatrix* pi_qcs_matrix;

typedef IMatrix4x4 ti_qcs_matrix4x4;
typedef IMatrix4x4* pi_qcs_matrix4x4;

typedef struct schmidt_decomposition {
    int is_entangled, schmidt_coeff_number;
    tf_qcs_matrix *base1;
    tf_qcs_matrix *base2;
    tf_qcs_matrix *schmidt_coeff;
} SchmidtDecomposition;

typedef SchmidtDecomposition tf_qcs_schmidt_decomposition;
typedef SchmidtDecomposition* pf_qcs_schmidt_decomposition;

typedef struct spectral_decomposition {
    tf_qcs_matrix *eigenvectors;
    tf_qcs_matrix *eigenvalues;
} SpectralDecomposition;

typedef SpectralDecomposition tf_qcs_spectral_decomposition;
typedef SpectralDecomposition* pf_qcs_spectral_decomposition;

typedef struct svd_decomposition {
    tf_qcs_matrix *U;
    tf_qcs_matrix *S;
    tf_qcs_matrix *V;
} SVDDecomposition;

typedef SVDDecomposition tf_qcs_svd_decomposition;
typedef SVDDecomposition* pf_qcs_svd_decomposition;

/* matrix of any size which contains complex numbers */

tf_qcs_matrix *qcs_create_matrix(int rows, int cols);
tf_qcs_matrix *qcs_create_eye_matrix(int size);
tf_qcs_matrix *qcs_create_matrix_for_quantum_state(int q, int freedom_level);
tf_qcs_matrix *qcs_create_empty_matrix();

tf_qcs_matrix *qcs_create_unitary_random_real_matrix(int size);
tf_qcs_matrix *qcs_create_unitary_random_matrix(int size);

tf_qcs_matrix *qcs_create_unitary_random_real_matrix_by_qr(int size);
tf_qcs_matrix *qcs_create_unitary_random_matrix_by_qr(int size);

tf_qcs_matrix *qcs_create_density_random_real_matrix(int size);
tf_qcs_matrix *qcs_create_density_random_matrix(int size);

tf_qcs_matrix *qcs_create_hermitian_random_real_matrix(int size);
tf_qcs_matrix *qcs_create_hermitian_random_matrix(int size);

tf_qcs_matrix *qcs_create_matrix_arange_operation(tf_qcs_real_number _start, tf_qcs_real_number _end, tf_qcs_real_number _step);
tf_qcs_matrix *qcs_create_matrix_linspace_operation_with_endpoint(tf_qcs_real_number _start, tf_qcs_real_number _end, int _n_count);
tf_qcs_matrix *qcs_create_matrix_linspace_operation_without_endpoint(tf_qcs_real_number _start, tf_qcs_real_number _end, int _n_count);

tf_qcs_matrix *qcs_create_matrix_arange_operation_with_float_args(float _start, float _end, float _step);
tf_qcs_matrix *qcs_create_matrix_linspace_operation_with_endpoint_with_float_args(float _start, float _end, int _n_count);
tf_qcs_matrix *qcs_create_matrix_linspace_operation_without_endpoint_with_float_args(float _start, float _end, int _n_count);

tf_qcs_matrix *qcs_create_create_operator(int N);
tf_qcs_matrix *qcs_create_destroy_operator(int N);

void qcs_delete_matrix(tf_qcs_matrix *a_out);
void qcs_copy_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_out);
tf_qcs_matrix *qcs_clone_matrix(tf_qcs_matrix *a_in);
void qcs_eye_matrix(tf_qcs_matrix *a_in);
void qcs_eye_matrix_with_param(tf_qcs_matrix *a_in, tf_qcs_real_number p);
void qcs_zero_matrix(tf_qcs_matrix *a_in);
void qcs_roundn_matrix(tf_qcs_matrix *a_in, int n);
void qcs_chop_matrix(tf_qcs_matrix *a_in);
void qcs_fill_matrix(tf_qcs_matrix *a_in);
void qcs_fill_random_real_values(tf_qcs_matrix *a_in);
void qcs_fill_random_complex_values(tf_qcs_matrix *a_in);

void qcs_make_projector_from_vector(tf_qcs_matrix *vec_in, tf_qcs_matrix *prj_out);

void qcs_normalise_real_matrix(tf_qcs_matrix *a_in);

tf_qcs_matrix *qcs_get_column_from_matrix(tf_qcs_matrix *a_in, int c);
tf_qcs_matrix *qcs_get_row_from_matrix(tf_qcs_matrix *a_in, int r);

void qcs_get_column_from_matrix_into_vector(tf_qcs_matrix *a_in, int c, tf_qcs_matrix *vec_out);

void qcs_set_cell_at_matrix_complex(tf_qcs_matrix *a_in, int r, int c, tf_qcs_complex *b_in);
void qcs_set_cell_at_matrix_direct(tf_qcs_matrix *a_in, int r, int c, tf_qcs_real_number re, tf_qcs_real_number im);

tf_qcs_complex* qcs_get_cell_at_matrix_complex(tf_qcs_matrix *a_in, int r, int c);

tf_qcs_complex qcs_calc_dot_of_vector(tf_qcs_matrix *vec_in);
tf_qcs_complex qcs_calc_dot_of_two_vector(tf_qcs_matrix *vec_in_1, tf_qcs_matrix *vec_in_2);

int qcs_non_zero_elements_of_matrix(tf_qcs_matrix *a_in);

void qcs_det_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_out);
void qcs_inv_matrix(tf_qcs_matrix *a_in);

void qcs_inv_2x2_matrix(tf_qcs_matrix *a_in);

void qcs_transpose_matrix(tf_qcs_matrix *a_in);
void qcs_partial_transpose_matrix_qudit(tf_qcs_matrix *a_in, int i);

// realize D * D' operations
tf_qcs_matrix *qcs_calculate_d_dot_dt_matrix(tf_qcs_matrix *a_in);

void qcs_matrix_realignment(tf_qcs_matrix *a_in, tf_qcs_matrix *out_in);
void qcs_matrix_realignment2(tf_qcs_matrix *a_in, tf_qcs_matrix *out_in);

void qcs_swap_block_matrix(tf_qcs_matrix *a_in, int r1, int c1, int r2, int c2, int size_r, int size_c);

void qcs_add_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out);
void qcs_sub_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out);
void qcs_mul_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out);

// A = A + cB
void qcs_add_isM_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_in);
// A = A - cB
void qcs_sub_isM_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_in);

void qcs_add_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out);
void qcs_sub_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out);
void qcs_mul_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out);
void qcs_div_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out);

void qcs_scalar_mul_matrix(tf_qcs_matrix *a_mat, tf_qcs_complex *b_in);

void qcs_basic_sum_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out);
void qcs_tensor_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out);

tf_qcs_matrix *qcs_basic_sum(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in);
tf_qcs_matrix *qcs_tensor_product(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in);

void qcs_partial_transpose_matrix(tf_qcs_matrix *a_in, int n);
void qcs_partial_transpose_matrix_qudit(tf_qcs_matrix *a_in, int i);

tf_qcs_matrix* qcs_create_partial_trace_matrix_1_qubit(tf_qcs_matrix *a_in, int n);
tf_qcs_matrix* qcs_create_partial_trace_matrix_n_qubit(tf_qcs_matrix *a_in, int _from, int _to);
tf_qcs_matrix* qcs_create_partial_trace_matrix_1_qudit(tf_qcs_matrix *a_in, int i);

void qcs_spectral_decompose_of_matrix(tf_qcs_matrix *a_mat, tf_qcs_matrix *eigenvalues, tf_qcs_matrix *eigenvectors);
void qcs_svd_decompose_of_matrix(tf_qcs_matrix *state, tf_qcs_matrix *out_coeff, tf_qcs_matrix *out_base1, tf_qcs_matrix *out_base2);

tf_qcs_matrix* qcs_square_root_of_operator_matrix(tf_qcs_matrix *a_mat);
void qcs_square_root_of_operator_matrix_self(tf_qcs_matrix *a_mat);

tf_qcs_matrix* qcs_norm_of_matrix(tf_qcs_matrix *a_mat);
void qcs_norm_of_matrix_self(tf_qcs_matrix *a_mat);

void qcs_exp_of_matrix(tf_qcs_matrix *a_in);

tf_qcs_complex qcs_infinity_norm_of_matrix(tf_qcs_matrix *a);

tf_qcs_real_number qcs_entropy_of_matrix(tf_qcs_matrix *a);
tf_qcs_real_number qcs_linear_entropy_of_matrix(tf_qcs_matrix *a);
tf_qcs_real_number qcs_negativity_of_matrix(tf_qcs_matrix *a);

tf_qcs_complex qcs_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);
tf_qcs_complex qcs_square_of_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);
tf_qcs_complex qcs_super_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);
tf_qcs_complex qcs_trace_distance(tf_qcs_matrix *rho, tf_qcs_matrix *sigma);
tf_qcs_complex qcs_hilbert_schmidt_distance(tf_qcs_matrix *rho, tf_qcs_matrix *sigma);
tf_qcs_complex qcs_bures_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);
tf_qcs_complex qcs_angle_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);
tf_qcs_complex qcs_sine_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho);

tf_qcs_complex* qcs_eigenvalues_of_matrix(tf_qcs_matrix *a_in);

void qcs_add_noise_to_matrix(tf_qcs_matrix *a, tf_qcs_real_number n);

tf_qcs_complex qcs_trace_matrix(tf_qcs_matrix *a_in);
tf_qcs_complex qcs_trace_square_matrix(tf_qcs_matrix *a_in);
tf_qcs_complex qcs_trace_square_matrix_fast(tf_qcs_matrix *a_in);
tf_qcs_complex qcs_calculate_sum_of_matrix(tf_qcs_matrix *a_in);

void qcs_print_eigenvalue_of_matrix(tf_qcs_matrix *a);
void qcs_make_eigenvectors_of_matrix(tf_qcs_matrix *a);

void qcs_print_matrix(tf_qcs_matrix *a_in);
void qcs_print_matrix_only_real_part(tf_qcs_matrix *a_in);
void qcs_print_matrix_dot(tf_qcs_matrix *a_in, int bits, int base);

void qcs_print_matrix_in_matlab_format(tf_qcs_matrix *a_in);
void qcs_print_matrix_in_mathematica_format(tf_qcs_matrix *a_in);

void qcs_print_matrix_to_file(tf_qcs_matrix *a_in, FILE *b_out);
void qcs_print_matrix_to_file_sqr(tf_qcs_matrix *a_in, FILE *b_out);
void qcs_print_matrix_to_file_sqr_01(tf_qcs_matrix *a_in, FILE *b_out);
void qcs_print_matrix_to_file_sqr_mathematica(tf_qcs_matrix *a_in, FILE *b_out);
void qcs_print_matrix_to_file_sqr_matlab(tf_qcs_matrix *a_in, FILE *b_out);

void qcs_save_matrix_as_flat_text_file_for_gnuplot(tf_qcs_matrix *a_in, char *fname_out);

tf_qcs_schmidt_decomposition* qcs_create_schmidt_decomposition_empty();
tf_qcs_schmidt_decomposition* qcs_create_schmidt_decomposition(tf_qcs_matrix *b1, tf_qcs_matrix *b2, tf_qcs_matrix *sc);
void qcs_delete_schmidt_decomposition( tf_qcs_schmidt_decomposition *sd);

int qcs_compare_schmidt_decompositions(tf_qcs_schmidt_decomposition *a, tf_qcs_schmidt_decomposition *b, int test_type);

tf_qcs_spectral_decomposition* qcs_create_spectral_decomposition();
void qcs_delete_spectral_decomposition( tf_qcs_spectral_decomposition *sd);

tf_qcs_svd_decomposition* qcs_create_svd_decomposition();
void qcs_delete_svd_decomposition( tf_qcs_svd_decomposition *sd);

int qcs_matrix_ispure(tf_qcs_matrix *a);

tf_qcs_matrix* qcs_mixing_two_density_matrix(tf_qcs_real_number p1, tf_qcs_matrix *m1, tf_qcs_real_number p2, tf_qcs_matrix* m2);

tf_qcs_matrix* qcs_create_horodecky_9x9_state(tf_qcs_real_number a);
tf_qcs_matrix* qcs_create_horodecky_9x9_state_with_param(tf_qcs_real_number a, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_horodecki_4x4_state(tf_qcs_real_number p, tf_qcs_real_number a, tf_qcs_real_number b);
tf_qcs_matrix* qcs_create_ha_9x9_state(tf_qcs_real_number g);
tf_qcs_matrix* qcs_create_ha_9x9_state_with_param(tf_qcs_real_number g, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_w0_9x9_matrix();

tf_qcs_matrix* qcs_create_maximally_mixed_state(int x);
tf_qcs_matrix* qcs_create_maximally_mixed_state_with_param(int x, tf_qcs_real_number p);

void qcs_update_horodecky_9x9_state(tf_qcs_matrix *tmp, tf_qcs_real_number a);
void qcs_update_horodecki_4x4_state(tf_qcs_matrix *m, tf_qcs_real_number p, tf_qcs_real_number a, tf_qcs_real_number b);
void qcs_update_ha_9x9_state(tf_qcs_matrix *m, tf_qcs_real_number g);

/* entangled detection */

tf_qcs_real_number qcs_witnesses_application_to_matrix(tf_qcs_matrix *w, tf_qcs_matrix *m);
tf_qcs_real_number qcs_ppt_criterion(tf_qcs_matrix *m);
tf_qcs_real_number qcs_ccnr_criterion(tf_qcs_matrix *m);

/* depolarizing channel Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_depolarizing_channel_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_depolarizing_channel_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* amplitude damping Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_amplitude_damping_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_amplitude_damping_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* phase damping Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_phase_damping_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_phase_damping_v2_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_phase_damping_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* phase flip Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_phase_flip_operator_for_qubit(int m, tf_qcs_real_number p);
//tf_qcs_matrix* qcs_create_phase_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* phase flip Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_bit_flip_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_bit_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* bit-phase flip Kraus operators for qubit and qudit */

tf_qcs_matrix* qcs_create_bit_phase_flip_operator_for_qubit(int m, tf_qcs_real_number p);
tf_qcs_matrix* qcs_create_bit_phase_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p);

/* fully correlated phase flip for qubit */

//tf_qcs_matrix* qcs_create_fully_correlated_phase_flip_operator_for_qubit_register(int m, int q, tf_qcs_real_number p);

/* four dimensional matrix of any size which contains complex numbers */

tf_qcs_matrix4x4 *qcs_create_matrix4x4(int x1, int x2, int x3, int x4);
void qcs_delete_matrix4x4(tf_qcs_matrix4x4 *m);

void qcs_set_cell_at_matrix4x4(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, tf_qcs_complex *m);
void qcs_set_cell_at_matrix4x4_direct(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, tf_qcs_real_number re, tf_qcs_real_number im);
tf_qcs_complex* qcs_get_cell_at_matrix4x4(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4);

void qcs_zero_matrix4x4(tf_qcs_matrix4x4 *a_in);

void qcs_print_matrix4x4(tf_qcs_matrix4x4 *a_in);

/* two dimensional matrix of any size which contains integer numbers */

ti_qcs_matrix *qcs_create_int_matrix(int rows, int cols);
void qcs_delete_int_matrix(ti_qcs_matrix *m);

void qcs_set_cell_at_int_matrix(ti_qcs_matrix *a_in, int r, int c, int v);
int qcs_get_cell_at_int_matrix(ti_qcs_matrix *a_in, int r, int c);

void qcs_print_int_matrix(ti_qcs_matrix *a_in);

/* four dimensional matrix of any size which contains integer numbers */

ti_qcs_matrix4x4 *qcs_create_int_matrix4x4(int x1, int x2, int x3, int x4);
void qcs_delete_int_matrix4x4(ti_qcs_matrix4x4 *a_in);

void qcs_set_cell_at_int_matrix4x4(ti_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, int v);
int qcs_get_cell_at_int_matrix4x4(ti_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4);

void qcs_zero_int_matrix4x4(tf_qcs_matrix4x4 *a_in);

void qcs_print_int_matrix4x4(ti_qcs_matrix4x4 *a_in);


#ifdef __cplusplus
};
#endif


#endif /* __qcs_matrix_h__ */
