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

DYNAMIC_LIB_DECORATION tf_qcs_quantum_register* qcs_new_quantum_register(int size)
{
    int i, vec_state_size;
    tf_qcs_quantum_register* tmp;

    tmp = (tf_qcs_quantum_register*)malloc( sizeof(tf_qcs_quantum_register) );

    tmp->n = size;
    tmp->vec_state_size = 1 << tmp->n;

    tmp->vs = (Complex*)malloc( sizeof(Complex) * (vec_state_size) );

    for ( i = 0 ; i < (  tmp->vec_state_size ) ; i++)
    {
        tmp->vs[i].re = 0;
        tmp->vs[i].im = 0;
    }

    tmp->mode = USE_STATE_VECTOR_QUBIT;
    tmp->fdl = 2;
    tmp->el = 0;

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
    int i, vec_state_size;

    q_reg->vec_state_size = 1 << q_reg->n;

    if (q_reg->mode==USE_STATE_VECTOR_QUBIT)
    {
        for (i=1;i<q_reg->vec_state_size;i++)
        {
            q_reg->vs[i].re=0;
            q_reg->vs[i].im=0;
        }
        q_reg->vs[0].re=1;
        q_reg->vs[0].im=0;
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
            for (i=0;i<(1 << q_reg->vec_state_size);i++)
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
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
#else
                printf( "%2.6f + %2.6fi |%s>\n", q_reg->vs[i].re, q_reg->vs[i].im, msg );
#endif
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
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout( "%s%2.6f + %2.6fi |%s>\n", prefix, q_reg->vs[i].re, q_reg->vs[i].im, msg );
#else
                printf( "%s%2.6f + %2.6fi |%s>\n", prefix, q_reg->vs[i].re, q_reg->vs[i].im, msg );
#endif
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
