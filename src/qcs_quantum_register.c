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

#include <string.h>
#include <math.h>

#ifdef PYTHON_SCRIPT
#include <Python.h>
#endif

#include "qcs_quantum_register.h"
#include "qcs_qubit_gates.h"

#ifdef PYTHON_SCRIPT
#define _PRINT  PySys_WriteStdout
#else
#define _PRINT  printf
#endif


DYNAMIC_LIB_DECORATION tf_qcs_quantum_register* qcs_new_quantum_register(int size)
{
    int i;
    tf_qcs_quantum_register* tmp;

    tmp = (tf_qcs_quantum_register*)malloc( sizeof(tf_qcs_quantum_register) );

    tmp->n = size;
    tmp->vec_state_size = 1 << tmp->n;

    tmp->vs = (tf_qcs_complex*)malloc( sizeof(tf_qcs_complex) * (tmp->vec_state_size) );

    for ( i = 0 ; i < (  tmp->vec_state_size ) ; i++)
    {
        tmp->vs[i].re = (tf_qcs_real_number)0.0;
        tmp->vs[i].im = (tf_qcs_real_number)0.0;
    }

    tmp->mode = USE_STATE_VECTOR_QUBIT;
    tmp->fdl = 2;
    tmp->el = 0;

    tmp->dim_of_qs = NULL;

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_quantum_register(tf_qcs_quantum_register *q_reg)
{
    int i;

    if (q_reg->mode == USE_STATE_VECTOR_QUBIT)
    {

        if (q_reg->vs != NULL)
        {
            free((void*)(q_reg->vs));
        }
        q_reg->vs=NULL;
    }

    if (q_reg->mode == USE_CHP_MODE)
    {
    }

    if (q_reg->mode == USE_ONEWAY_MODEL)
    {
    }

    if(q_reg->mode == USE_DENSITY_MATRIX)
    {
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
    }


    if ( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
    }

    if ( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT )
    {
    }

    free((void*)q_reg);
    q_reg=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_reset(tf_qcs_quantum_register *q_reg)
{
    int i;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        for ( i = 1 ; i < (  q_reg->vec_state_size ) ; i++)
        {
            q_reg->vs[i].re = (tf_qcs_real_number)0.0;
            q_reg->vs[i].im = (tf_qcs_real_number)0.0;
        }
        q_reg->vs[0].re = (tf_qcs_real_number)1.0;
        q_reg->vs[0].im = (tf_qcs_real_number)0.0;

        q_reg->el = 0;
    }

    
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
	    // empty code
    }

    if (q_reg->mode == USE_DENSITY_MATRIX)
    {
	    // empty code
    }

    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
    {
	    // empty code
    }

    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
	    // empty code
    }
    

    q_reg->el = 0;
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_reset_error_level(tf_qcs_quantum_register *q_reg, int v)
{
    q_reg->el=0;
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_set_error_level(tf_qcs_quantum_register *q_reg, int v)
{
    q_reg->el=v;
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_get_error_level(tf_qcs_quantum_register *q_reg)
{
    return q_reg->el;
}


DYNAMIC_LIB_DECORATION void qcs_quantum_register_set_state_dec(tf_qcs_quantum_register *q_reg, int n)
{
    int i;
    char t[128];

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        t[0]=0;
        qcs_dec2bin(n, q_reg->n, &t[0]);

        if (q_reg->vs!=NULL)
        {
            for (i=0;i<(q_reg->vec_state_size);i++)
            {
                q_reg->vs[i].re=0;
                q_reg->vs[i].im=0;
            }
            q_reg->vs[n].re=1;
            q_reg->vs[n].im=0;
        }
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
	    // empty code
    }

    if(q_reg->mode == USE_DENSITY_MATRIX)
    {
	    // empty code
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_set_state_bin(tf_qcs_quantum_register *q_reg, char *state_desc)
{
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        int v;

        v=qcs_bin2dec(state_desc);
        qcs_quantum_register_set_state_dec(q_reg, v);
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        //int v;

        //v = qcs_base_d2dec(state_desc, q_reg->qudit_state->freedom_level);
        //qcs_quantum_register_set_state_dec(q_reg, v);
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_fill_zero(tf_qcs_quantum_register *q_reg)
{
    int i;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        for (i=0;i<(q_reg->vec_state_size);i++)
        {
            q_reg->vs[i].re=0.0;
            q_reg->vs[i].im=0.0;
        }
    }

/*
    if (q_reg->mode==USE_STATE_VECTOR_QUDIT)
    {
        qcs_qudit_state_fill_zero( q_reg->qudit_state );
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_set_ghz_state(tf_qcs_quantum_register *q_reg)
{
    qcs_quantum_register_fill_zero( q_reg );

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        q_reg->vs[0].re=(tf_qcs_real_number)(1.0 / sqrt(2.0));
        q_reg->vs[0].im=(tf_qcs_real_number)0.0;
        q_reg->vs[(q_reg->vec_state_size)-1].re=(tf_qcs_real_number)(1.0 / sqrt(2.0));
        q_reg->vs[(q_reg->vec_state_size)-1].im=(tf_qcs_real_number)0.0;
    }

}



/**
 * 
 * Unitary gates functions
 * for single qubit/qudits 
 * 
**/


static inline void oper_on_rows(tf_qcs_quantum_register *q_reg, int r1, int r2, tf_qcs_matrix *u)
{
    tf_qcs_complex tmp1, tmp2, o_mul_tmp1, o_mul_tmp2, o_sum_tmp1, o_sum_tmp2 ;

    tmp1.re=0;
    tmp1.im=0;
    tmp2.re=0;
    tmp2.im=0;

    o_sum_tmp1.re=0;
    o_sum_tmp1.im=0;
    o_sum_tmp2.re=0;
    o_sum_tmp2.im=0;

    o_mul_tmp1.re=0;
    o_mul_tmp1.im=0;
    o_mul_tmp2.re=0;
    o_mul_tmp2.im=0;

    tmp1 = *( q_reg->vs + r1 );
    tmp2 = *( q_reg->vs + r2 );

    qcs_complex_mul((u->m + 0), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((u->m + 1), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp1);

    qcs_complex_mul((u->m + 2), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((u->m + 3), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp2);

    (*( q_reg->vs + r1 )).re = o_sum_tmp1.re;
    (*( q_reg->vs + r1 )).im = o_sum_tmp1.im;
    (*( q_reg->vs + r2 )).re = o_sum_tmp2.re;
    (*( q_reg->vs + r2 )).im = o_sum_tmp2.im;
}

static int compare_int ( const void *arg1, const void *arg2 )
{
    if ( *((int*)arg1) <  *((int*)arg2) ) return -1;
    if ( *((int*)arg1) == *((int*)arg2) ) return  0;
    if ( *((int*)arg1) >  *((int*)arg2) ) return  1;
}

DYNAMIC_LIB_DECORATION void applied_1q_gate_to_quantum_register ( tf_qcs_quantum_register *q_reg, int t, tf_qcs_matrix *u)
{
    int irow=0, ip=0, i=0, n=0, m=0, step=0, vstep=0, p=0, r1=0, r2=0 ;

    char r1_bin[127], r2_bin[127];

    if ( u == NULL )
    {
        q_reg->el = -1;
        return ;
    }

    if ( u->rows != 2 && u->cols != 2 )
    {
        q_reg->el = -1;
        return ;
    }

    n=q_reg->n;
    t++;

    m = pow( 2, n - 1 );      // mini matrices
    step = pow( 2, n - t ) - 1;  // step between matrices
    p = pow( 2, t - 1 );      // block count
    vstep = pow( 2, n ) / p ;     // step in vector state

    //printf("m=%d step=%d p=%d vstep=%d\n", m, step, p, vstep);

    irow = 0;
    for ( ip = 0 ; ip < p ; ip++ )
    {
        for ( i = 0 ; i < step + 1 ; i++ )
        {

            r1=irow + i;
            r2=irow + i + step + 1;
/*
            memset(&r1_bin[0], 0, 127);
            memset(&r2_bin[0], 0, 127);

            qcs_dec2bin( r1, q_reg->size, &r1_bin[0]);
            qcs_dec2bin( r2, q_reg->size, &r2_bin[0]);

            printf("%-> %.4d %.4d [%.4s, %.4s]\n", r1, r2, &r1_bin[0], &r2_bin[0]);
*/
            oper_on_rows(q_reg, r1, r2, u);

        }
        irow = irow + vstep;
    }
}

DYNAMIC_LIB_DECORATION void applied_2q_gate_to_quantum_register_one_control(tf_qcs_quantum_register *q_reg, int c1, int t, tf_qcs_matrix *u)
{
    int irow=0, i=0, n=0, m=0, step=0, r1=0, r2=0;
    char irow_addr[128];

    if ( u == NULL )
    {
        q_reg->el = -1;
        return ;
    }

    if ( u->rows != 2 && u->cols != 2 )
    {
        q_reg->el = -1;
        return ;
    }

    n=q_reg->n;
    c1++;
    t++;

    step=pow(2, n - t) - 1; // odstep pomiedzy elementami malych macierzy
    m=pow(2, n - 2);        // ilość małych macierzy dla przypadku dwuqubitowej bramki

    for ( i = 0 ; i < m ; i++ )
    {
        memset(&irow_addr[0], 0, 128);
        qcs_dec2bin(i, n - 2, &irow_addr[0]);

        // odpowiednia kolejność dodawania bitów
        if ( c1 < t )
        {
            qcs_str_insert_at(c1 - 1, '1', &irow_addr[0]);
            qcs_str_insert_at(t  - 1, '0', &irow_addr[0]);
        }
        if ( c1 > t )
        {
            qcs_str_insert_at(t  - 1, '0', &irow_addr[0]);
            qcs_str_insert_at(c1 - 1, '1', &irow_addr[0]);
        }

        irow = qcs_bin2dec(irow_addr);

        r1 = irow ;
        r2 = irow + step + 1;

        oper_on_rows(q_reg, r1, r2, u);
        // printf("[%d] [%s] irow = %d %d [c1=%d,t=%d]\n", i, &irow_addr[0], irow, irow+step+1, c1-1, t-1);
    }

    // printf("m=%d step=%d \n", m, step);
}


DYNAMIC_LIB_DECORATION void qcs_quantum_register_pauli_x_gate(tf_qcs_quantum_register *q_reg, int i)
{
    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        applied_1q_gate_to_quantum_register( q_reg, i, get_pauli_x_gate() );
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_pauli_y_gate(tf_qcs_quantum_register *q_reg, int i)
{
    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        applied_1q_gate_to_quantum_register( q_reg, i, get_pauli_y_gate() );
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_pauli_z_gate(tf_qcs_quantum_register *q_reg, int i)
{
    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        applied_1q_gate_to_quantum_register( q_reg, i, get_pauli_z_gate() );
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_had_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    //char t[128];
    //int j;
    //tf_qcs_qubit a;
    //tf_qcs_matrix *m;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_1q_gate_to_quantum_register(q_reg, i, get_hadamard_gate());
    }

    if ( q_reg->mode == USE_DENSITY_MATRIX)
    {
        //qcs_density_matrix_apply_arbitrary_unitary_operator_qubit_fast(q_reg->dms, i, get_hadamard_gate());
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        //applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_hadamard_gate( q_reg->freedom_level) );
        //applied_1qudit_gate_to_quantum_reg_FASTER( q_reg->qudit_state, i, get_qudit_hadamard_gate( q_reg->freedom_level) );
    }

    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
    {
        //applied_1qubit_gate_to_symbolic_quantum_reg( q_reg->qubit_symbolic_state, i, get_hadamard_gate());
    }
    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
        //applied_1qudit_gate_to_qudit_symbolic_quantum_reg( q_reg->qudit_symbolic_state, i, get_qudit_hadamard_gate( q_reg->freedom_level) );
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_had_n_conj_gate(tf_qcs_quantum_register *q_reg, int i)
{

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_1q_gate_to_quantum_register(q_reg, i, get_hadamard_gate());
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        //tf_qcs_matrix *m=qcs_clone_matrix(get_qudit_hadamard_gate( q_reg->freedom_level));
        //qcs_transpose_matrix( m );
        //applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i,  m);
        //qcs_delete_matrix( m );
    }

    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
    {
    }
    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
    }

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_had_gate_for_whole_register(tf_qcs_quantum_register *q_reg)
{
    int i;

    for ( i=0 ; i<q_reg->n ; i++ )
    {
        qcs_quantum_register_had_n_gate(q_reg, i);
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_square_root_not_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_square_root_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }

    if (q_reg->mode==USE_CHP_MODE)
    {
    }
*/
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_1q_gate_to_quantum_register(q_reg, i, get_square_root_gate());
    }
/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_square_root_gate(q_reg->freedom_level) );
    }

    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
    {
        //applied_1qubit_gate_to_symbolic_quantum_reg( q_reg->qubit_symbolic_state, i, get_square_root_gate());
    }
    if (q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
        //applied_1qudit_gate_to_qudit_symbolic_quantum_reg( q_reg->qudit_symbolic_state, i, get_qudit_square_root_gate(q_reg->freedom_level) );
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_x_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_x_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<=q_reg->size;j++) t[j]='i';
//           t[i]='9';
//           t[i+1]='x';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='x';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);
//           m=make_matrix_for_one_qubit("9x", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register(q_reg, i, get_x_rot90_gate());
    }
/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_x_rot90_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_y_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;
/*
    if (q_reg->mode==NO_STATE_VECTOR)
    {
        qcs_1q_y_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<=q_reg->size;j++) t[j]='i';
//           t[i]='9';
//           t[i+1]='y';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='y';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("9y", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i,  get_y_rot90_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_y_rot90_gate(q_reg->freedom_level) );
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_z_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_z_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/    
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<=q_reg->size;j++) t[j]='i';
//           t[i]='9';
//           t[i+1]='z';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='z';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("9z", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i,  get_z_rot90_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_z_rot90_gate(q_reg->freedom_level) );
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_mx_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_minus_x_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/    
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<=q_reg->size+1;j++) t[j]='i';
//           t[i]='m';
//           t[i+1]='9';
//           t[i+2]='x';
//           t[q_reg->size - i - 2]='m';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='x';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);
//           m=make_matrix_for_one_qubit("m9x", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

//           qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_minus_x_rot90_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_minus_x_rot90_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_my_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_minus_y_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/
   
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<=q_reg->size+1;j++) t[j]='i';
//           t[i]='m';
//           t[i+1]='9';
//           t[i+2]='y';
//           t[q_reg->size - i - 2]='m';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='y';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("m9y", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

//           qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_minus_y_rot90_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_minus_y_rot90_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_mz_rot90_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_minus_z_rotate90_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        memset(&t[0],0,sizeof(t));

//           for(j=0;j<=q_reg->size+1;j++) t[j]='i';
//           t[i]='m';
//           t[i+1]='9';
//           t[i+2]='z';
//           t[q_reg->size - i - 2]='m';
//           t[q_reg->size - i - 1]='9';
//           t[q_reg->size - i ]='z';
//           printf("--[%s]--\n", &t[0]);
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("m9z", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_minus_z_rot90_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_minus_z_rot90_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_x_rot_n_gate(tf_qcs_quantum_register *q_reg, int i, tf_qcs_real_number theta)
{
    tf_qcs_matrix *m;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        m = qcs_get_rot_x_gate( theta );
        applied_1q_gate_to_quantum_register( q_reg, i, m);
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_y_rot_n_gate(tf_qcs_quantum_register *q_reg, int i, tf_qcs_real_number theta)
{
    tf_qcs_matrix *m;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        m = qcs_get_rot_y_gate( theta );
        applied_1q_gate_to_quantum_register( q_reg, i, m);
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_z_rot_n_gate(tf_qcs_quantum_register *q_reg, int i, tf_qcs_real_number theta)
{
    tf_qcs_matrix *m;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        m = qcs_get_rot_z_gate( theta );
        applied_1q_gate_to_quantum_register( q_reg, i, m);
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_rotate_alpha_n_gate(tf_qcs_quantum_register *q_reg, int i, tf_qcs_real_number alpha)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_complex c;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_rotate_alpha_gate(q_reg->qubits[i], k, &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {

        c.re=0;
        c.im=(2.0*QCS_PI) / pow(2.0, alpha);
        qcs_exp_complex(&c, get_rotate_alpha_gate()->m+3);

//           m=qcs_build_matrix_from_tensor(t);
//           m=make_matrix_for_one_qubit("r", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_rotate_alpha_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_rotate_alpha_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_rotate_theta_n_gate(tf_qcs_quantum_register *q_reg, int i, tf_qcs_real_number theta)
{
    char t[128];
    int j;

    tf_qcs_qubit aq;
    tf_qcs_complex a,b,c,d;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_rotate_theta_gate(q_reg->qubits[i], theta, &aq);

        q_reg->qubits[i]->alpha=aq.alpha;
        q_reg->qubits[i]->beta=aq.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        a.re= cos(theta);
        a.im=0;
        b.re=-sin(theta);
        b.im=0;
        c.re= sin(theta);
        c.im=0;
        d.re=cos(theta);
        d.im=0;

        (get_rotate_theta_gate()->m+0)->re=a.re;
        (get_rotate_theta_gate()->m+0)->im=a.im;

        (get_rotate_theta_gate()->m+1)->re=b.re;
        (get_rotate_theta_gate()->m+1)->im=b.im;

        (get_rotate_theta_gate()->m+2)->re=c.re;
        (get_rotate_theta_gate()->m+2)->im=c.im;

        (get_rotate_theta_gate()->m+3)->re=d.re;
        (get_rotate_theta_gate()->m+3)->im=d.im;

//           m=make_matrix_for_one_qubit("g", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        applied_1q_gate_to_quantum_register( q_reg, i, get_rotate_theta_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_rotate_theta_gate(q_reg->freedom_level) );
    }

    if( q_reg->mode == USE_DENSITY_MATRIX )
    {
        a.re= cos(theta);
        a.im=0;
        b.re=-sin(theta);
        b.im=0;
        c.re= sin(theta);
        c.im=0;
        d.re=cos(theta);
        d.im=0;

        (get_rotate_theta_gate()->m+0)->re=a.re;
        (get_rotate_theta_gate()->m+0)->im=a.im;

        (get_rotate_theta_gate()->m+1)->re=b.re;
        (get_rotate_theta_gate()->m+1)->im=b.im;

        (get_rotate_theta_gate()->m+2)->re=c.re;
        (get_rotate_theta_gate()->m+2)->im=c.im;

        (get_rotate_theta_gate()->m+3)->re=d.re;
        (get_rotate_theta_gate()->m+3)->im=d.im;

        qcs_density_matrix_apply_arbitrary_unitary_operator_qubit_fast(q_reg->dms, i, get_rotate_theta_gate() );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_t_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if ( q_reg->mode == USE_NO_STATE_VECTOR )
    {
        qcs_1q_t_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<q_reg->size;j++) t[j]='i';
//           t[i]='t';
//           t[q_reg->size - i - 1]='h';
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("t", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_t_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_t_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_v_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    qcs_quantum_register_s_n_gate(q_reg, i);
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_s_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_s_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<q_reg->size;j++) t[j]='i';
//           t[i]='s';
//           t[q_reg->size - i - 1]='h';
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("s", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        applied_1q_gate_to_quantum_register( q_reg, i, get_s_gate() );
    }
/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_s_gate(q_reg->freedom_level) );
    }
*/    
}

void qcs_quantum_register_s_adj_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;

/*
    if (q_reg->mode==NO_STATE_VECTOR)
    {
        qcs_1q_s_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha=a.alpha;
        q_reg->qubits[i]->beta=a.beta;
    }
*/
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
//           memset(&t[0],0,sizeof(t));
//           for(j=0;j<q_reg->size;j++) t[j]='i';
//           t[i]='s';
//           t[q_reg->size - i - 1]='h';
//           m=qcs_build_matrix_from_tensor(t);

//           m=make_matrix_for_one_qubit("s", q_reg->size, i);
//           qcs_quantum_reg_mult_by_matrix(q_reg, m);
//           qcs_delete_matrix(m);

        //qcs_quantum_reg_rebuild_state_vector(q_reg);
        m = qcs_create_matrix(2, 2);
        qcs_copy_matrix(get_s_gate(), m);
        qcs_transpose_matrix(m);
        applied_1q_gate_to_quantum_register( q_reg, i, m );
        qcs_delete_matrix( m );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        m = get_qudit_s_gate(q_reg->freedom_level);
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, m );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_phase_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];

    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;
    tf_qcs_complex c;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        //qcs_1q_phase_gate(q_reg->qubits[i], &a);

        //q_reg->qubits[i]->alpha = a.alpha;
        //q_reg->qubits[i]->beta  = a.beta;
    }


    if (q_reg->mode==USE_CHP_MODE)
    {
        qcs_chp_phase(q_reg->chp_state, i);
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_1q_gate_to_quantum_register( q_reg, i, get_phase_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_phase_gate(q_reg->freedom_level) );
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_phase_f_n_gate(tf_qcs_quantum_register *q_reg, int i)
{
    char t[128];
    int j;
    tf_qcs_qubit a;
    tf_qcs_matrix *m;
    tf_qcs_complex c;

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
        qcs_1q_phase_f_gate(q_reg->qubits[i], &a);

        q_reg->qubits[i]->alpha = a.alpha;
        q_reg->qubits[i]->beta  = a.beta;
    }

    if (q_reg->mode==USE_CHP_MODE)
    {
        qcs_chp_phase(q_reg->chp_state, i);
    }
*/
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_1q_gate_to_quantum_register( q_reg, i, get_phase_f_gate() );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, get_qudit_phase_f_gate(q_reg->freedom_level) );
    }
*/

}

DYNAMIC_LIB_DECORATION void qcs_qubit_arbitrary_one_qubit_gate(tf_qcs_quantum_register *q_reg, tf_qcs_matrix *gate, int i)
{
    tf_qcs_matrix *m;

    if (gate == NULL)
    {
        q_reg->el = -1;
        return ;
    }

/*
    if (q_reg->mode==USE_NO_STATE_VECTOR)
    {
    }

    if (q_reg->mode==CHP_MODE)
    {
    }
*/

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        if (gate->rows != 2 && gate->cols != 2)
        {
            q_reg->el = -1;
            return ;
        }

        applied_1q_gate_to_quantum_register( q_reg, i, gate );
    }

/*
    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        applied_1qudit_gate_to_quantum_reg( q_reg->qudit_state, i, gate);
    }
*/
}

DYNAMIC_LIB_DECORATION void qcs_arbitrary_single_gate(tf_qcs_quantum_register *q_reg, tf_qcs_matrix *m, int i)
{
    qcs_qubit_arbitrary_one_qubit_gate(q_reg, m, i);
}



/**
 * 
 * two qubit gates
 * 
 *
**/

DYNAMIC_LIB_DECORATION void qcs_quantum_register_cnot(tf_qcs_quantum_register *q_reg, int c1, int t)
{
    tf_qcs_qubit qa, qb;
    tf_qcs_matrix *m;

    if (q_reg->mode==USE_CHP_MODE)
    {
    }

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_2q_gate_to_quantum_register_one_control(q_reg, c1, t, get_not_gate());
    }

    if (q_reg->mode==USE_PQC_MODE)
    {
    }

    if ( q_reg->mode == USE_DENSITY_MATRIX )
    {
        //qcs_density_matrix_apply_arbitrary_two_qubit_operator_one_control(q_reg->dms, c1, t, get_not_gate());
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        //applied_2qudit_gate_to_quantum_reg(q_reg->qudit_state, c1, t, get_qudit_not_gate(q_reg->freedom_level) );
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
    {
        //applied_2qubit_gate_to_symbolic_quantum_reg_one_control(q_reg->qubit_symbolic_state, c1, t, get_not_gate());
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
        //applied_2qudit_gate_to_qudit_symbolic_quantum_reg_one_control(q_reg->qudit_symbolic_state, c1, t, get_qudit_not_gate(q_reg->freedom_level) );

        //tf_qcs_matrix *m;
        //m = cnot_qudit_syntesis_one_control_one_target( qcs_qudit_symbolic_state_get_size(q_reg->qudit_symbolic_state) , qcs_qudit_symbolic_state_get_freedom_level(q_reg->qudit_symbolic_state), c1, t);
        //apply_matrix_to_qudit_symbolic_quantum_reg(q_reg->qudit_symbolic_state, m);
        //qcs_delete_matrix( m );
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_cnot_conj(tf_qcs_quantum_register *q_reg, int c1, int t)
{
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        applied_2q_gate_to_quantum_register_one_control(q_reg, c1, t, get_not_gate());
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        //applied_cnot_conj_gate_to_qudit_quantum_reg_one_control(q_reg->qudit_state, c1, t);
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
    {
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_cnot_by_elem_permutation(tf_qcs_quantum_register *q_reg, int c1, int t)
{
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
    }
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_toffoli_by_elem_permutation(tf_qcs_quantum_register *q_reg, int c1, int t)
{
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
    }
}

/**
 * 
 * multiqubit gates
 * 
 * 
**/

// TODO (marek#1#): finalise implementation
DYNAMIC_LIB_DECORATION void qcs_quantum_register_swap_gate(tf_qcs_quantum_register *q_reg, int a, int b)
{
}

// TODO (marek#1#): finalise implementation
DYNAMIC_LIB_DECORATION void qcs_quantum_reg_fredkin_gate(tf_qcs_quantum_register *q_reg, int a, int b, int c)
{
}

/**
 * 
 * Measurement functions
 * for qubit/qudits and quantum registers
 * 
**/


DYNAMIC_LIB_DECORATION void qcs_quantum_register_probe_one_qubit(tf_qcs_quantum_register *q_reg, int k, tf_qcs_real_number *out_value_0, tf_qcs_real_number *out_value_1)
{
    qcs_quantum_register_probe_one_qubit_in_std_base(q_reg, k, out_value_0, out_value_1);
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_probe_one_qubit_in_std_base(tf_qcs_quantum_register *q_reg, int t, tf_qcs_real_number *out_value_0, tf_qcs_real_number *out_value_1)
{
    char bin[127];
    int m_result = -1;
    int i, j, q_mask, n, max_i;

    tf_qcs_real_number fmod_a, fmod_b;
    tf_qcs_complex tmp, d_tmp;

    tf_qcs_real_number  fmod = 0, probe_1 = 0, probe_0 = 0, probe ;

    n=q_reg->n;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        // check the number of qubit
        if ( t < 0 || t >= n )
        {
            q_reg->el = ERROR_BAD_QUBIT_NUMBER;

            return;
        }

        max_i = 1 << q_reg->n ; // max_i=(int)pow(2, n)-1;

        memset(&bin[0], 0, 127);

        for (j=0;j<n;j++) bin[j]='0';
        bin[t]='1';

        q_mask=qcs_bin2dec(&bin[0]);

        // collection of probability value

        for (i = 0 ; i < max_i ; i++)
        {
            if ((i & q_mask) == q_mask) // for qubit in state 1
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs+i, &fmod);
                probe_1 = probe_1 + ( fmod * fmod );
            }

            if (((~i) & q_mask) == q_mask) // for qubit in state 0
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs+i, &fmod);
                probe_0 = probe_0 + ( fmod * fmod );
            }
        } // for (i = 0 ; i < max_i ; i++)
    }

    *out_value_0 = probe_0;
    *out_value_1 = probe_1;
}

DYNAMIC_LIB_DECORATION void qcs_quantum_reg_probe_one_qubit_in_base(tf_qcs_quantum_register *q_reg, int k, tf_qcs_qubit_base_desc *q_base, tf_qcs_real_number *out_value_0, tf_qcs_real_number *out_value_1)
{
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit(tf_qcs_quantum_register *q_reg, int k)
{
    if( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        return qcs_quantum_register_measure_one_qubit_in_std_base(q_reg, k);
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
        //return qcs_qubit_symbolic_reg_measure(q_reg->qubit_symbolic_state, k);
        return -1;
    }
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit_force(tf_qcs_quantum_register *q_reg, int k, int force_result)
{
    if( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        return qcs_quantum_register_measure_one_qubit_in_std_base_force(q_reg, k, force_result);
    }
    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
        //return qcs_qubit_symbolic_reg_measure_force(q_reg->qubit_symbolic_state, k, force_result);
        return -1;
    }
}

static void remove_states_1(tf_qcs_quantum_register *q_reg, int max_i, int q_mask)
{
    int i;

    for (i = 0 ; i < max_i ; i++)
    {
        if ((i & q_mask) == q_mask) // for qubit in state 1
        {
            q_reg->vs[i].re=(tf_qcs_real_number)0.0;
            q_reg->vs[i].im=(tf_qcs_real_number)0.0;
        }
    }
}

static void remove_states_0(tf_qcs_quantum_register *q_reg, int max_i, int q_mask)
{
    int i;

    for (i = 0 ; i < max_i ; i++)
    {
        if (((~i) & q_mask) == q_mask) // for qubit in state 0
        {
            q_reg->vs[i].re=0;
            q_reg->vs[i].im=0;
        }
    }
}

static void oper_on_rows_collect_prob(tf_qcs_quantum_register *q_reg, int r1, int r2, tf_qcs_matrix *s0, tf_qcs_matrix *s1, tf_qcs_real_number *pr_0, tf_qcs_real_number *pr_1)
{
    tf_qcs_real_number fmod;

    tf_qcs_complex
    tmp1, tmp2,
    o_mul_tmp1, o_mul_tmp2,
    o_sum_tmp1, o_sum_tmp2,
    o_sum_final ;

    /* for matrix s0 */

    tmp1.re=0;
    tmp1.im=0;
    tmp2.re=0;
    tmp2.im=0;

    o_sum_tmp1.re=0;
    o_sum_tmp1.im=0;
    o_sum_tmp2.re=0;
    o_sum_tmp2.im=0;

    o_mul_tmp1.re=0;
    o_mul_tmp1.im=0;
    o_mul_tmp2.re=0;
    o_mul_tmp2.im=0;

    tmp1 = *( q_reg->vs + r1 );
    tmp2 = *( q_reg->vs + r2 );

    qcs_complex_mul((s0->m + 0), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((s0->m + 1), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp1);

    qcs_complex_mul((s0->m + 2), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((s0->m + 3), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp2);

    qcs_complex_mul( ( q_reg->vs + r1 ), &o_sum_tmp1, &tmp1 );
    qcs_complex_mul( ( q_reg->vs + r2 ), &o_sum_tmp2, &tmp2 );
    qcs_complex_add( &tmp1, &tmp2, &o_sum_final);

    qcs_mod_complex(&o_sum_final, &fmod);
    *pr_0 = *pr_0 + ( fmod );

    /* for matrix s1 */

    tmp1.re=0;
    tmp1.im=0;
    tmp2.re=0;
    tmp2.im=0;

    o_sum_tmp1.re=0;
    o_sum_tmp1.im=0;
    o_sum_tmp2.re=0;
    o_sum_tmp2.im=0;

    o_mul_tmp1.re=0;
    o_mul_tmp1.im=0;
    o_mul_tmp2.re=0;
    o_mul_tmp2.im=0;

    tmp1 = *( q_reg->vs + r1 );
    tmp2 = *( q_reg->vs + r2 );

    qcs_complex_mul((s1->m + 0), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((s1->m + 1), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp1);

    qcs_complex_mul((s1->m + 2), &tmp1, &o_mul_tmp1);
    qcs_complex_mul((s1->m + 3), &tmp2, &o_mul_tmp2);
    qcs_complex_add(&o_mul_tmp1, &o_mul_tmp2, &o_sum_tmp2);


    qcs_complex_mul( ( q_reg->vs + r1 ), &o_sum_tmp1, &tmp1 );
    qcs_complex_mul( ( q_reg->vs + r2 ), &o_sum_tmp2, &tmp2 );
    qcs_complex_add( &tmp1, &tmp2, &o_sum_final);

    qcs_mod_complex(&o_sum_final, &fmod);
    *pr_1 = *pr_1 + ( fmod );
}


static void collect_prob_in_base(tf_qcs_quantum_register *q_reg, int t, tf_qcs_matrix *s0, tf_qcs_matrix *s1, tf_qcs_real_number *pr_0, tf_qcs_real_number *pr_1)
{
    int irow=0, ip=0, i=0, n=0, m=0, step=0, vstep=0, p=0, r1=0, r2=0 ;

    n=q_reg->n;
    t++;

    m = pow( 2, n - 1 );         // mini matrices
    step = pow( 2, n - t ) - 1;  // step between matrices
    p = pow( 2, t - 1 );         // block count
    vstep = pow( 2, n ) / p ;    // step in vector state

    // printf("m=%d step=%d p=%d vstep=%d\n", m, step, p, vstep);

    irow = 0;
    for ( ip = 0 ; ip < p ; ip++ )
    {
        for ( i = 0 ; i < step + 1 ; i++ )
        {
            // printf("%-> %d %d\n", irow + i, irow + i + step + 1);
            r1=irow + i;
            r2=irow + i + step + 1;

            oper_on_rows_collect_prob(q_reg, r1, r2, s0, s1, pr_0, pr_1);

        }
        irow = irow + vstep;
    }
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit_in_std_base(tf_qcs_quantum_register *q_reg, int k)
{
    int m_result = -1;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        char bin[127];
        int i, j, q_mask, n, max_i;

        tf_qcs_real_number fmod_a, fmod_b;
        tf_qcs_complex tmp, d_tmp;

        tf_qcs_real_number fmod = 0, probe_1 = 0, probe_0 = 0, probe ;

        n=q_reg->n;

        // check the number of qubit
        if ( k < 0 || k >= n )
        {
            return -1;
        }

        max_i = 1 << q_reg->n ; // max_i=(int)pow(2, n)-1;

#ifdef __qcs_core_library_debug_mode__
        printf("i_max=%d\n", max_i);
#endif

        memset(&bin[0], 0, 127);

        for (j=0;j<n;j++) bin[j]='0';
        bin[k]='1';

        q_mask=qcs_bin2dec(&bin[0]);

#ifdef __qcs_core_library_debug_mode__
        printf("q_mask=%d [mask %sb]\n", q_mask, &bin[0]);
#endif

        // collection of probability value

        for (i = 0 ; i < max_i ; i++)
        {
            if ((i & q_mask) == q_mask) // for qubit in state 1
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs + i, &fmod);
                probe_1 = probe_1 + ( fmod * fmod );
            }
            if (((~i) & q_mask) == q_mask) // for qubit in state 0
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs + i, &fmod);
                probe_0 = probe_0 + ( fmod * fmod );
            }
        }


        // removing amplitudes
#ifdef __qcs_core_library_debug_mode__
        printf("p1=%f p0=%f sum ( %f )\n", probe_1, probe_0, probe_1 + probe_0);
#endif
        if (probe_1 == probe_0)
        {
            if (rand() % 2 == 0) // measuring qubit in state 1
            {
#ifdef __qcs_core_library_debug_mode__
                printf("measuring qubit in state 1\n");
#endif
                m_result = 1;
                remove_states_0(q_reg, max_i, q_mask);
                probe=probe_1;
            }
            else // measuring qubit in state 0
            {
#ifdef __qcs_core_library_debug_mode__
                printf("measuring qubit in state 0\n");
#endif
                m_result = 0;
                remove_states_1(q_reg, max_i, q_mask);
                probe=probe_0;
            }
        }

        if ( probe_1 > probe_0) // measuring qubit in state 1
        {
#ifdef __qcs_core_library_debug_mode__
            printf("measuring qubit in state 1\n");
#endif
            m_result = 1;
            remove_states_0(q_reg, max_i, q_mask);
            probe=probe_1;
        }
        if ( probe_1 < probe_0) // measuring qubit in state 0
        {
#ifdef __qcs_core_library_debug_mode__
            printf("measuring qubit in state 0\n");
#endif
            m_result = 0;
            remove_states_1(q_reg, max_i, q_mask);
            probe=probe_0;
        }

        // state normalization
        for ( i=0 ; i<max_i ; i++ )
        {
            qcs_mod_complex(q_reg->vs + i, &fmod);

            if (fmod!=0)
            {
                tmp.re=1/sqrt(probe);
                tmp.im=0;
                qcs_complex_mul(&tmp, q_reg->vs + i, &d_tmp);
                q_reg->vs[i].re=d_tmp.re;
                q_reg->vs[i].im=d_tmp.im;

            }
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT )

    return m_result;
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit_in_std_base_force(tf_qcs_quantum_register *q_reg, int k, int force_result)
{
    char bin[127];
    int m_result = -1;
    int i, j, q_mask, n, max_i;

    tf_qcs_real_number fmod_a, fmod_b;
    tf_qcs_complex tmp, d_tmp;

    tf_qcs_real_number  fmod = 0, probe_1 = 0, probe_0 = 0, probe ;

    n=q_reg->n;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        // check the number of qubit
        if ( k < 0 || k >= n )
        {
            return -1;
        }

        if (force_result != 0 && force_result !=1)
        {
            return -2;
        }

        max_i = 1 << q_reg->n ; // max_i=(int)pow(2, n)-1;

#ifdef __qcs_core_library_debug_mode__
        printf("i_max=%d\n", max_i);
#endif

        memset(&bin[0], 0, 127);

        for (j=0;j<n;j++) bin[j]='0';
        bin[k]='1';

        q_mask=qcs_bin2dec(&bin[0]);

#ifdef __qcs_core_library_debug_mode__
        printf("q_mask=%d [mask %sb]\n", q_mask, &bin[0]);
#endif

        // collection of probability value

        for (i = 0 ; i < max_i ; i++)
        {
            if ((i & q_mask) == q_mask) // for qubit in state 1
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs + i, &fmod);
                probe_1 = probe_1 + ( fmod * fmod );
            }
            if (((~i) & q_mask) == q_mask) // for qubit in state 0
            {
                fmod=0;
                qcs_mod_complex(q_reg->vs + i, &fmod);
                probe_0 = probe_0 + ( fmod * fmod );
            }
        }


        if (force_result == 0)
        {
            m_result = 0;
            remove_states_1(q_reg, max_i, q_mask);
            probe=probe_0;
        }

        if (force_result == 1)
        {
            m_result = 1;
            remove_states_0(q_reg, max_i, q_mask);
            probe=probe_1;
        }

        // state normalization
        for (i=0;i<max_i;i++)
        {
            qcs_mod_complex(q_reg->vs + i, &fmod);

            if (fmod!=0)
            {
                tmp.re=1/sqrt(probe);
                tmp.im=0;
                qcs_complex_mul(&tmp, q_reg->vs + i, &d_tmp);
                q_reg->vs[i].re=d_tmp.re;
                q_reg->vs[i].im=d_tmp.im;

            }
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT )

    return m_result;
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit_in_base(tf_qcs_quantum_register *q_reg, int k, tf_qcs_qubit_base_desc *q_base)
{
#ifdef __qcs_core_library_debug_mode__
    char bin[127];
    int q_mask;
#endif
    int m_result = -1;
    int i, j, n, max_i;

    tf_qcs_real_number fmod_a, fmod_b;
    tf_qcs_complex tmp, d_tmp;

    tf_qcs_real_number  fmod = 0, probe_1 = 0, probe_0 = 0, probe ;

    n=q_reg->n;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        // check the number of qubit
        if ( k < 0 || k >= n )
        {
            return -1;
        }

        max_i = 1 << q_reg->n ; // max_i=(int)pow(2, n)-1;

#ifdef __qcs_core_library_debug_mode__
        printf("i_max=%d\n", max_i);
#endif


#ifdef __qcs_core_library_debug_mode__
        memset(&bin[0], 0, 127);

        for (j=0;j<n;j++) bin[j]='0';
        bin[k]='1';

        q_mask=qcs_bin2dec(&bin[0]);
        printf("q_mask=%d [mask %sb]\n", q_mask, &bin[0]);
#endif

        // collection of probability value
        collect_prob_in_base( q_reg, k, q_base->gate_s0, q_base->gate_s1, &probe_0, &probe_1);

        // removing amplitudes
#ifdef __qcs_core_library_debug_mode__
        printf("p1=%f p0=%f sum ( %f )\n", probe_1, probe_0, probe_1 + probe_0);
#endif
        if ( probe_1 == probe_0 )
        {
            if ( rand() % 2 == 0 ) // measuring qubit in state 1
            {
#ifdef __qcs_core_library_debug_mode__
                printf("measuring qubit in base state 1\n");
#endif
                //remove_states_0(q_reg, max_i, q_mask);
                applied_1q_gate_to_quantum_register( q_reg, k, q_base->gate_s1 );
                
                probe=probe_1;
                m_result = 1;
            }
            else // measuring qubit in state 0
            {
#ifdef __qcs_core_library_debug_mode__
                printf("measuring qubit in base state 0\n");
#endif
                //remove_states_1(q_reg, max_i, q_mask);
                applied_1q_gate_to_quantum_register( q_reg, k, q_base->gate_s0 );
                probe=probe_0;
                m_result = 0;
            }
        }

        if ( probe_1 > probe_0) // measuring qubit in state 1
        {
#ifdef __qcs_core_library_debug_mode__
            printf("measuring qubit in base state 1\n");
#endif
            //remove_states_0(q_reg, max_i, q_mask);
            applied_1q_gate_to_quantum_register( q_reg, k, q_base->gate_s1 );
            probe=probe_1;
            m_result = 1;
        }
        if ( probe_1 < probe_0) // measuring qubit in state 0
        {
#ifdef __qcs_core_library_debug_mode__
            printf("measuring qubit in base state 0\n");
#endif
            //remove_states_1(q_reg, max_i, q_mask);
            applied_1q_gate_to_quantum_register( q_reg, k, q_base->gate_s0 );
            probe=probe_0;
            m_result = 0;
        }
        // state normalization
        for (i=0;i<max_i;i++)
        {
            qcs_mod_complex(q_reg->vs + i, &fmod);

            if (fmod!=0)
            {
                tmp.re=1/sqrt(probe);
                tmp.im=0;
                qcs_complex_mul(&tmp, q_reg->vs + i, &d_tmp);
                q_reg->vs[i].re=d_tmp.re;
                q_reg->vs[i].im=d_tmp.im;

            }
        }

    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT )


    return m_result;
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_one_qubit_with_observable(tf_qcs_quantum_register *q_reg, int k, tf_qcs_matrix *m)
{

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
    }

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
    }


    return -1; // unknown simulation mode
}

DYNAMIC_LIB_DECORATION int qcs_quantum_register_measure_from_to(tf_qcs_quantum_register *q_reg, int q_from, int q_to)
{
    int idx, i, q, qubits, qubit_value;
    //char t1[128]; 
    char t2[128];

    qubits=q_to - q_from + 1;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {

        // we check whether qubits numbers are correct
        if ( q_from > q_to )
            return -1;
        if ( q_from < 0 || q_from >= q_reg->n )
            return -1;
        if ( q_to < 0 || q_to >= q_reg->n )
            return -1;

        idx=0;
        memset(t2, 0, sizeof(t2));
        for (q=q_from;q<=q_to;q++)
        {
                qubit_value = qcs_quantum_register_measure_one_qubit_in_std_base(q_reg, q);
                if (qubit_value==0) t2[idx]='0';
                if (qubit_value==1) t2[idx]='1';
                idx++;
        }
        t2[idx]=0;
        return qcs_bin2dec(t2);
    }

/*
    if( q_reg->mode == USE_STATE_VECTOR_QUDIT )
    {
        //  we check whether qudits numbers are correct
        if ( q_from > q_to )
            return -1;
        if ( q_from < 0 || q_from >= q_reg->n )
            return -1;
        if ( q_to < 0 || q_to >= q_reg->n )
            return -1;
    }
*/

    return -1;

}


DYNAMIC_LIB_DECORATION void qcs_quantum_register_print_bin(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;
    char msg[512];

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT)
    {
        memset( msg, 0, sizeof(msg) );
        max_value = 1 << q_reg->n;

        for ( i=0; i<max_value; i++ )
        {
            qcs_dec2bin( i, q_reg->n, msg );
            if (q_reg->vs[i].re != 0 || q_reg->vs[i].im != 0)
            {

//#ifdef PYTHON_SCRIPT
//                PySys_WriteStdout( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
//#else
//                printf( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
//#endif
                _PRINT( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
            }
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT)


    if ( q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        // qcs_quantum_reg_print_d_base( q_reg );
    }


/*
    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
        //qcs_quantum_symbolic_reg_print_bin_full( q_reg->qubit_symbolic_state );
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT )
    {
        //qcs_quantum_qudit_symbolic_reg_print_d_base_full( q_reg->qudit_symbolic_state );
    }
*/
}

void qcs_quantum_register_print_bin_in_matlab_format(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;
    char msg[512];

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT)
    {
        memset( msg, 0, sizeof(msg) );
        max_value = 1 << q_reg->n;

        _PRINT( "[ ... \n" );
        for ( i=0; i<max_value; i++ )
        {
            qcs_dec2bin( i, q_reg->n, msg );
            if (q_reg->vs[i].re != 0 || q_reg->vs[i].im != 0)
            {
                _PRINT( "(%2.6f + %2.6fi) ; ... %% |%s> \n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
            }
        }
        _PRINT( "]\n" );
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT)


/*
    if ( q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
    }
*/

/*
    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
    }

    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT )
    {
    }
*/   
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_print_bin_sqr(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;
    char msg[512];
    tf_qcs_real_number ma, sum = 0;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        memset(msg, 0, sizeof(msg));
        max_value=1 << q_reg->n;

        for (i=0;i<max_value;i++)
        {
            qcs_mod_complex(&q_reg->vs[i], &ma);
            qcs_dec2bin(i, q_reg->n, msg);
            if ( q_reg->vs[i].re!=0 || q_reg->vs[i].im!=0 )
            {
//#ifdef PYTHON_SCRIPT
//                PySys_WriteStdout("%2.6f |%s>\n", ma * ma, msg);
//#else
//                printf("%2.6f |%s>\n", ma * ma, msg);
//#endif
                _PRINT("%2.6f |%s>\n", ma * ma, msg);
                sum+=ma*ma;
            }
        }
//#ifdef PYTHON_SCRIPT
//                PySys_WriteStdout("sum=%2.6f\n", sum);
//#else
//                printf("amplitude sum=%2.6f\n", sum);
//#endif
        _PRINT("amplitude sum=%2.6f\n", sum);

    } // if (q_reg->mode==USE_STATE_VECTOR_QUBIT)

/*
    if(q_reg->mode==USE_STATE_VECTOR_QUDIT)
    {

        memset ( msg, 0, sizeof(msg) );
        max_value = (int)pow(q_reg->freedom_level, q_reg->size);
        for ( i = 0 ; i < max_value ; i++ )
        {
            qcs_dec2base_d(i, q_reg->size, q_reg->freedom_level, msg);
            if (q_reg->qudit_state->vec_state[i].re!=0 || q_reg->qudit_state->vec_state[i].im!=0)
            {
                qcs_mod_complex(&q_reg->qudit_state->vec_state[i], &ma);
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout("%2.6f |%s>\n", ma * ma, msg);
#else
                printf("%2.6f  |%s>\n", ma * ma, msg);
#endif
                sum+=ma*ma;
            }
        }

#ifdef PYTHON_SCRIPT
                PySys_WriteStdout("sum=%2.6f\n", sum);
#else
                printf("sum=%2.6f\n", sum);
#endif
    } // if (q_reg->mode==USE_STATE_VECTOR_QUDIT)
*/
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_print_bin_full(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;
    char msg[512];

/*
    if ( q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        qcs_quantum_reg_print_d_base_full( q_reg );
    }
*/

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        memset( msg, 0, sizeof( msg ) );
        max_value = 1 << q_reg->n;

        for ( i=0; i<max_value; i++ )
        {
            qcs_dec2bin( i, q_reg->n, msg );
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
#else
            printf( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
#endif
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT)

/*
    if( q_reg->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    {
        //qcs_quantum_symbolic_reg_print_bin_full( q_reg->qubit_symbolic_state );
    }
*/
}

void qcs_quantum_register_print_bin_full_sqr(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;
    char msg[512];
    tf_qcs_real_number ma, sum = 0;

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT )
    {
        memset( msg, 0, sizeof( msg ) );
        max_value = 1 << q_reg->n;

        for ( i=0; i<max_value; i++ )
        {
            qcs_mod_complex(&q_reg->vs[i], &ma);
            qcs_dec2bin( i, q_reg->n, msg );
            
            _PRINT("%2.6f |%s>\n", ma * ma, msg);
            sum+=ma*ma;
        }

        _PRINT("amplitude sum=%2.6f\n", sum);
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT)

}

void qcs_quantum_register_print_bin_with_prefix(tf_qcs_quantum_register *q_reg, char *prefix)
{
    int i, max_value;
    char msg[512];

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT)
    {
        memset( msg, 0, sizeof(msg) );
        max_value = 1 << q_reg->n;

        for ( i=0; i<max_value; i++ )
        {
            qcs_dec2bin( i, q_reg->n, msg );
            if (q_reg->vs[i].re != 0 || q_reg->vs[i].im != 0)
            {
//#ifdef PYTHON_SCRIPT
//                PySys_WriteStdout( "%s%2.6f + %2.6fi |%s>\n", prefix, q_reg->vs[i].re, q_reg->vs[i].im, msg );
//#else
//                printf( "%s%2.6f + %2.6fi |%s>\n", prefix, q_reg->vs[i].re, q_reg->vs[i].im, msg );
//#endif
                _PRINT( "%s%2.6f + %2.6fi |%s>\n", prefix, q_reg->vs[i].re, q_reg->vs[i].im, msg );
            }
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT)    
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_print_dec(tf_qcs_quantum_register *q_reg)
{
    int i, max_value;

    if ( q_reg->vs == NULL )
    {
#ifdef PYTHON_SCRIPT
        PySys_WriteStdout("no state vec\n");
#else
        printf("no state vec\n");
#endif
    }

    if ( q_reg->mode == USE_STATE_VECTOR_QUBIT)
    {
        max_value = 1 << q_reg->n;
        for ( i = 0 ; i < max_value ; i++ )
        {
            if (q_reg->vs[i].re!=0 || q_reg->vs[i].im!=0)
            {
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout("%2.6f + %2.6fi |%d>\n", q_reg->vs[i].re, q_reg->vs[i].im, i);
#else
                printf("%2.6f + %2.6fi |%d>\n", q_reg->vs[i].re, q_reg->vs[i].im, i);
#endif
            }
        }
    } // if( q_reg->mode == USE_STATE_VECTOR_QUBIT_QUBIT)
}

DYNAMIC_LIB_DECORATION void qcs_quantum_register_get_probability_amplitude(tf_qcs_quantum_register *q_reg, int idx, tf_qcs_real_number *out_value_re, tf_qcs_real_number *out_value_im)
{
    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {

        if (idx <= 1 << q_reg->n )
        {
            *out_value_re =  (tf_qcs_real_number) q_reg->vs[idx].re;
            *out_value_im =  (tf_qcs_real_number) q_reg->vs[idx].im;
        }

    } // if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_quantum_register_generate_density_matrix(tf_qcs_quantum_register *q_reg)
{
    int x, y, max_value;
    tf_qcs_complex a_tmp, b_tmp, out_tmp;
    tf_qcs_matrix *d=NULL;


    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        max_value=1 << q_reg->n;
        d = qcs_create_matrix(max_value, max_value);
        d->q = q_reg->n;
        d->freedom_level = 2;

        for (y=0;y<max_value;y++)
        {
            a_tmp.re = q_reg->vs[y].re;
            a_tmp.im = q_reg->vs[y].im;
            for (x=0;x<max_value;x++)
            {
                qcs_conj_complex(&q_reg->vs[x], &b_tmp);
                qcs_complex_mul(&a_tmp, &b_tmp, &out_tmp);

                qcs_set_cell_at_matrix_complex(d, x, y, &out_tmp);
            }
        }
    } // if (q_reg->mode==USE_STATE_VECTOR_QUBIT)

    if (q_reg->mode == USE_STATE_VECTOR_QUDIT)
    {
        /*
        max_value=(int)powf(q_reg->qudit_state->freedom_level, q_reg->qudit_state->size);

        d = qcs_create_matrix(max_value,max_value);
        d->q = q_reg->qudit_state->size;
        d->freedom_level = q_reg->qudit_state->freedom_level;

        for (y=0;y<max_value;y++)
        {
            a_tmp.re = q_reg->qudit_state->vec_state[y].re;
            a_tmp.im = q_reg->qudit_state->vec_state[y].im;
            for (x=0;x<max_value;x++)
            {
                qcs_conj_complex(&q_reg->qudit_state->vec_state[x], &b_tmp);
                qcs_complex_mul(&a_tmp, &b_tmp, &out_tmp);

                qcs_set_cell_at_matrix_complex(d, x, y, &out_tmp);
            }
        }
        */
    } // if (q_reg->mode == USE_STATE_VECTOR_QUDIT)

    if ( q_reg->mode == USE_DENSITY_MATRIX )
    {
        /*
        d = qcs_clone_matrix(q_reg->dms->state);
        */
    }

    return d;
}