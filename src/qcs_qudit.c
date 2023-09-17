/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2007 -- 2008 by Przemysław Ratajczak                    *
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "qcs.h"
#include "qcs_complex.h"
#include "qcs_matrix.h"
#include "qcs_qudit.h"

DYNAMIC_LIB_DECORATION QCS_Qudit_State *qcs_new_qudit_state(int size, int freedom_level)
{
    int i = 0;
    QCS_Qudit_State* tmp;

    tmp = (QCS_Qudit_State*)malloc( sizeof(QCS_Qudit_State) );

    tmp->size = size;
    tmp->freedom_level = freedom_level;

    tmp->vec_state = (Complex*)malloc( sizeof(Complex) * (pow(freedom_level, size)));

    for( i = 0 ; i < ( pow(freedom_level, size) ) ; i++)
    {
        tmp->vec_state[i].re = 0;
        tmp->vec_state[i].im = 0;
    }
    int a=0;//pow(freedom_level, size)-1;stan |n>
    tmp->vec_state[a].re = 1; //stan |0>
    tmp->vec_state[a].im = 0;

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_qudit_state( ts_qcs_qudit_state *p )
{
       if ((p)->vec_state != NULL)
       {
             free((void*)(p)->vec_state);
             (p)->vec_state = NULL;
       }
       free((void*)p);
}

DYNAMIC_LIB_DECORATION void qcs_qudit_state_reset(ts_qcs_qudit_state *p)
{
    int i = 0;
    if ((p)->vec_state != NULL)
    {
        for( i = 0 ; i < ( pow(p->freedom_level, p->size) ) ; i++)
        {
            p->vec_state[i].re = 0;
            p->vec_state[i].im = 0;
        }
    }
    p->vec_state[0].re = 1;
    p->vec_state[0].im = 0;
}

DYNAMIC_LIB_DECORATION void qcs_qudit_state_fill_zero(ts_qcs_qudit_state *p)
{
    int i = 0;

    if ((p)->vec_state != NULL)
    {
        for( i = 0 ; i < ( pow(p->freedom_level, p->size) ) ; i++)
        {
            p->vec_state[i].re = 0;
            p->vec_state[i].im = 0;
        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_qudit_set_state_dec(ts_qcs_qudit_state *p, int n)
{
    int i = 0;
    if ((p)->vec_state != NULL)
    {
        for( i = 0 ; i < ( pow(p->freedom_level, p->size) ) ; i++)
        {
            p->vec_state[i].re = 0;
            p->vec_state[i].im = 0;
        }
    }
    p->vec_state[n].re = 1;
    p->vec_state[n].im = 0;
}

DYNAMIC_LIB_DECORATION void applied_1qudit_gate_to_quantum_reg(ts_qcs_qudit_state *q_state, int qn, tf_qcs_matrix *u)
{
    qn = qn + 1; /// IMPORTANT!!! operates on numbers of target qudits in [1..n], not in [0..d-1] !!! ///

    // 1 step: creating matrix (prev) by multiplication (NO OPTIMALIZATION DONE)
    int i = 1, j = 1;
    int freedom_lev = q_state->freedom_level;
    int prev_size = freedom_lev;  //number of cols&rows in prev matrix
    int output_size = freedom_lev;  //number of cols&rows in output matrix
    tf_qcs_matrix *prev = NULL, *output=NULL;
    tf_qcs_matrix *id_gate = get_qudit_id_gate( freedom_lev );

    prev = qcs_create_matrix(prev_size, prev_size);
    output = qcs_create_matrix(output_size, output_size);

    if (qn == 1) qcs_copy_matrix(u, prev);
        else qcs_copy_matrix(id_gate, prev);
    for (i=1; i<qn-1; i++)
    {
        qcs_delete_matrix(output);
        output_size = output_size * freedom_lev;
        output = qcs_create_matrix(output_size, output_size);
        qcs_tensor_matrix(prev, id_gate, output);

        qcs_delete_matrix(prev);
        prev_size = output_size;
        prev = qcs_create_matrix(prev_size, prev_size);
        qcs_copy_matrix(output,prev);
    }
    if (qn != 1)
    {
        qcs_delete_matrix(output);
        output_size = output_size * freedom_lev;
        output = qcs_create_matrix(output_size, output_size);
        qcs_tensor_matrix(prev, u, output);

        qcs_delete_matrix(prev);
        prev_size = output_size;
        prev = qcs_create_matrix(prev_size, prev_size);
        qcs_copy_matrix(output, prev);
    }
    for (i=qn; i<q_state->size; i++)
    {
        qcs_delete_matrix(output);
        output_size = output_size * freedom_lev;
        output = qcs_create_matrix(output_size, output_size);
        qcs_tensor_matrix(prev, id_gate, output);

       // if (i != q_state->size - 1) //skip last copying
       // {
            qcs_delete_matrix(prev);
            prev_size = output_size;
            prev = qcs_create_matrix(prev_size, prev_size);
            qcs_copy_matrix(output,prev);
       // }
    }

    // 2 step: multiplcation vector state by previously created matrix(prev)
//FILE *f;
// f=fopen("c:\\had0d3s3.txt","w");
// qcs_print_matrix_to_file_sqr_matlab(prev, f);

     tf_qcs_complex sum, result, tmp_result;

     Complex *tmp_vec_state=(Complex*)malloc(sizeof(Complex)*(pow(freedom_lev, q_state->size)));

     for(i = 0 ; i < (pow(freedom_lev, q_state->size)) ; i++)
     {
       tmp_vec_state[i].re = 0;
       tmp_vec_state[i].im = 0;
     }

    for(i = 0; i<prev->rows; i++)
     {
        sum.re = 0; sum.im = 0;
        for(j = 0; j<prev->cols; j++)
        {
           qcs_complex_mul(prev->m + i * prev->cols + j, q_state->vec_state + j, &result);
           qcs_complex_add(&sum, &result, &tmp_result);
           sum.re=tmp_result.re;
           sum.im=tmp_result.im;
        }
        tmp_vec_state[i].re = sum.re;
        tmp_vec_state[i].im = sum.im;
     }
     for(i = 0 ; i < (pow(freedom_lev, q_state->size)) ; i++)
     {
       q_state->vec_state[i].re = tmp_vec_state[i].re;
       q_state->vec_state[i].im = tmp_vec_state[i].im;
     }
}

inline void oper_on_rows_qudit(ts_qcs_qudit_state *q_state, long* r_tab, tf_qcs_matrix *u)
{
    int d = q_state->freedom_level;
    long ind = 0;
    long ind2 = 0;
    tf_qcs_complex temp[d];
    tf_qcs_complex o_sum;
    tf_qcs_complex result;

    result.re=0; result.im=0;

    for (ind = 0; ind < d; ind++)
        temp[ind] = *(q_state->vec_state + r_tab[ind]);

    for (ind = 0; ind < d; ind++)
    {
        o_sum.re=0; o_sum.im=0;
        for (ind2 = 0; ind2 < d; ind2++)
        {

            qcs_complex_mul((u->m + ind*d + ind2), &temp[ind2], &result);
            qcs_complex_add(&o_sum, &result, &o_sum); //     sum = sum+res;
            result.re=0; result.im=0;
        }
        (*(q_state->vec_state + r_tab[ind])).re = o_sum.re;
        (*(q_state->vec_state + r_tab[ind])).im = o_sum.im;
    }
}


DYNAMIC_LIB_DECORATION void applied_1qudit_gate_to_quantum_reg_FASTER(ts_qcs_qudit_state *q_state, int t, tf_qcs_matrix *u)
{

    long irow=0, ip=0, i=0, n=0, m=0, step=0, vstep=0, p=0, r1=0, r2=0 ;
    int d = q_state->freedom_level;
    long r_tab[d];
    long ind = 0;
    int tt=t;

    if( u == NULL )
    {
       // q_reg->error_level = -1;
        return ;
    }

    if( u->rows != d && u->cols != d )
    {
       // q_reg->error_level = -1;
        return ;
    }

    n = q_state->size;
    tt++;    /// !!! operates on target no in [1..d], not in [0..d-1]

    //    m = pow( d, n - 1 );      //ok mini matrices
        step = pow( d, n - tt ) - 1;  //ok step between matrices// nbr of zeros betwen nonzero elements
        p = pow( d, tt - 1 );      //ok block count
        vstep = pow( d, n ) / p ;     // step in vector state

#ifdef __qcs_core_library_debug_mode__
    printf("m=%d step=%d p=%d vstep=%d\n", m, step, p, vstep);
#endif

    irow = 0; //presuniecie na kolejne wiersze
    for( ip = 0 ; ip < p ; ip++ )
    {
        for( i = 0 ; i < step + 1 ; i++ )
        {
            for (ind = 0; ind < d; ind++)
            {
                r_tab[ind] = irow + ind*(step + 1) + i;
            }
            oper_on_rows_qudit(q_state, &r_tab[0], u);
        }
        irow = irow + vstep;
    }
}

DYNAMIC_LIB_DECORATION void applied_2qudit_gate_to_quantum_reg(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u)
{
    long irow=0, i=0, n=0, m=0, step=0, r1=0, r2=0;
    int d,fd = q_state->freedom_level;
    long r_tab[fd];
    char irow_addr[128];
    long ind = 0;
    tf_qcs_matrix *a;

    if ( u == NULL )
    {
        //q_reg->error_level = -1;
        printf("u is null\n");
        return ;
    }

    if ( u->rows != fd && u->cols != fd )
    {
        //q_reg->error_level = -1;
        printf("bad dimensions of u\n");
        return ;
    }

    n=q_state->size;
    c1++;
    t++;

    step=pow(fd, n - t) - 1; // odstep pomiedzy elementami malych macierzy
    m=pow(fd, n - 2);        // ilo�� ma�ych macierzy dla przypadku dwuquditowej bramki

    for( d=0; d < fd-1;d++)
    {
        for ( i = 0 ; i < m ; i++ )
        {
            memset(&irow_addr[0], 0, 128);
            //qcs_dec2bin(i, n - 2, &irow_addr[0]);
            qcs_dec2base_d(i, n-2, fd, &irow_addr[0]);
            // odpowiednia kolejno�� dodawania dit�w
            if ( c1 < t )
            {
                qcs_str_insert_at(c1 - 1, '1'+d, &irow_addr[0]);
                qcs_str_insert_at(t  - 1, '0', &irow_addr[0]);
            }
            if ( c1 > t )
            {
                qcs_str_insert_at(t  - 1, '0', &irow_addr[0]);
                qcs_str_insert_at(c1 - 1, '1'+d, &irow_addr[0]);
            }

            //irow = qcs_bin2dec(irow_addr);
            irow = qcs_base_d2dec(irow_addr, fd);

            //r1 = irow ;
            //r2 = irow + step + 1;

            for (ind = 0; ind < fd; ind++)
            {
                r_tab[ind] = irow + ind*(step + 1);
            }
            a = get_qudit_power_of_u_gate(u, d+1);
            //qcs_print_matrix(a);
            oper_on_rows_qudit(q_state, &r_tab[0], a);
            /*printf("u^%d ", d+1);
            printf("dec[%d] bin[%s] irow = (", i, &irow_addr[0] );
            for (ind = 0; ind < fd; ind++)
            {
                printf("%d ", r_tab[ind]);
            }
            printf(") [c1=%d,t=%d]\n",  c1-1, t-1);*/
        } // for ( i = 0 ; i < m ; i++ )
    } // for(d=0;d<fd;d++)

    //printf("m=%d step=%d \n", m, step);
}

DYNAMIC_LIB_DECORATION void applied_2qudit_gate_to_quantum_reg_d_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u)
{
}

/* gate with control state |1> -- qudit version */

DYNAMIC_LIB_DECORATION void applied_2qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u)
{
    long irow=0, i=0, n=0, m=0, step=0, r1=0, r2=0;
    int d = q_state->freedom_level;
    long r_tab[d];
    char irow_addr[128];
    long ind = 0;

    if ( u == NULL )
    {
        //q_reg->error_level = -1;
        printf("u is null\n");
        return ;
    }

    if ( u->rows != d && u->cols != d )
    {
        //q_reg->error_level = -1;
        printf("bad dimensions of u\n");
        return ;
    }

    n=q_state->size;
    c1++;
    t++;

    step=pow(d, n - t) - 1; // odstep pomiedzy elementami malych macierzy
    m=pow(d, n - 2);        // ilo�� ma�ych macierzy dla przypadku dwuquditowej bramki

    for ( i = 0 ; i < m ; i++ )
    {
        memset(&irow_addr[0], 0, 128);
        //qcs_dec2bin(i, n - 2, &irow_addr[0]);
        qcs_dec2base_d(i, n-2, d, &irow_addr[0]);
        // odpowiednia kolejno�� dodawania bit�w
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

        //irow = qcs_bin2dec(irow_addr);
        irow = qcs_base_d2dec(irow_addr, d);

        //r1 = irow ;
        //r2 = irow + step + 1;

        for (ind = 0; ind < d; ind++)
        {
            r_tab[ind] = irow + ind*(step + 1) + i;
        }
        oper_on_rows_qudit(q_state, &r_tab[0], u);
        //printf("dec[%d] bin[%s] irow = %d %d [c1=%d,t=%d]\n", i, &irow_addr[0], irow, irow+step+1, c1-1, t-1);
    }

//    printf("m=%d step=%d \n", m, step);

}

// TODO (marek#1#): to finish applied_3qudit_gate_to_quantum_reg_one_control
DYNAMIC_LIB_DECORATION void applied_3qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int t, tf_qcs_matrix *u)
{
}

// TODO (marek#1#): to finish applied_4qudit_gate_to_quantum_reg_one_control
DYNAMIC_LIB_DECORATION void applied_4qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int t, tf_qcs_matrix *u)
{
}

// TODO (marek#1#): to finish applied_5qudit_gate_to_quantum_reg_one_control
DYNAMIC_LIB_DECORATION void applied_5qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int t, tf_qcs_matrix *u)
{
}

// TODO (marek#1#): to finish applied_6qudit_gate_to_quantum_reg_one_control
DYNAMIC_LIB_DECORATION void applied_6qudit_gate_to_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int c5, int t, tf_qcs_matrix *u)
{
}


/* gate with control state |0> -- qudit version */

DYNAMIC_LIB_DECORATION void applied_2qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int t, tf_qcs_matrix *u)
{
}

DYNAMIC_LIB_DECORATION void applied_3qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int t, tf_qcs_matrix *u)
{
}

DYNAMIC_LIB_DECORATION void applied_4qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int t, tf_qcs_matrix *u)
{
}

DYNAMIC_LIB_DECORATION void applied_5qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int t, tf_qcs_matrix *u)
{
}

DYNAMIC_LIB_DECORATION void applied_6qudit_gate_to_quantum_reg_zero_control(ts_qcs_qudit_state *q_state, int c1, int c2, int c3, int c4, int c5, int t, tf_qcs_matrix *u)
{
}

DYNAMIC_LIB_DECORATION void applied_cnot_gate_to_qudit_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t)
{
    if (c1==t)
    {
        printf("\n CNOT GENERATION : Control = target !! \n");
        return;
    }

    tf_qcs_matrix *cnot = NULL;
    cnot = cnot_qudit_syntesis_one_control_one_target(q_state->size, q_state->freedom_level, c1, t);  //!! TODO : cash !!//
//qcs_print_matrix(cnot);

    // mult by register
     tf_qcs_complex sum, result, tmp_result;

     Complex *tmp_vec_state=(Complex*)malloc(sizeof(Complex)*(pow(q_state->freedom_level, q_state->size)));
     int i = 0, j = 0;
     for(i = 0 ; i < (pow(q_state->freedom_level, q_state->size)) ; i++)
     {
       tmp_vec_state[i].re = 0;
       tmp_vec_state[i].im = 0;
     }

    for(i = 0; i<cnot->rows; i++)
     {
        sum.re = 0; sum.im = 0;
        for(j = 0; j<cnot->cols; j++)
        {
           qcs_complex_mul(cnot->m + i * cnot->cols + j, q_state->vec_state + j, &result);
           qcs_complex_add(&sum, &result, &tmp_result);
           sum.re=tmp_result.re;
           sum.im=tmp_result.im;
        }
        tmp_vec_state[i].re = sum.re;
        tmp_vec_state[i].im = sum.im;
     }
     for(i = 0 ; i < (pow(q_state->freedom_level, q_state->size)) ; i++)
     {
       q_state->vec_state[i].re = tmp_vec_state[i].re;
       q_state->vec_state[i].im = tmp_vec_state[i].im;
     }
}


DYNAMIC_LIB_DECORATION void applied_cnot_conj_gate_to_qudit_quantum_reg_one_control(ts_qcs_qudit_state *q_state, int c1, int t)
{
    if (c1==t)
    {
        printf("\n CNOT GENERATION : Control = target !! \n");
        return;
    }

    tf_qcs_matrix *cnot = NULL;
    cnot = cnot_qudit_syntesis_one_control_one_target(q_state->size, q_state->freedom_level, c1, t);  //!! TODO : cash !!//
    qcs_transpose_matrix(cnot);
//qcs_print_matrix(cnot);

    // mult by register
     tf_qcs_complex sum, result, tmp_result;

     Complex
    *tmp_vec_state=(Complex*)malloc(sizeof(Complex)*(pow(q_state->freedom_level, q_state->size)));
     int i = 0, j = 0;
     for(i = 0 ; i < (pow(q_state->freedom_level, q_state->size)) ; i++)
     {
       tmp_vec_state[i].re = 0;
       tmp_vec_state[i].im = 0;
     }

    for(i = 0; i<cnot->rows; i++)
     {
        sum.re = 0; sum.im = 0;
        for(j = 0; j<cnot->cols; j++)
        {
           qcs_complex_mul(cnot->m + i * cnot->cols + j, q_state->vec_state + j, &result);
           qcs_complex_add(&sum, &result, &tmp_result);
           sum.re=tmp_result.re;
           sum.im=tmp_result.im;
        }
        tmp_vec_state[i].re = sum.re;
        tmp_vec_state[i].im = sum.im;
     }
     for(i = 0 ; i < (pow(q_state->freedom_level, q_state->size)) ; i++)
     {
       q_state->vec_state[i].re = tmp_vec_state[i].re;
       q_state->vec_state[i].im = tmp_vec_state[i].im;
     }
}


DYNAMIC_LIB_DECORATION tf_qcs_qudit* qcs_new_qudit(int freedom_level)
{
        tf_qcs_qudit *q;

        q=(tf_qcs_qudit *)malloc(sizeof(tf_qcs_qudit));

        q->freedom_level = freedom_level;
        q->alphas = (Complex*) malloc(sizeof(Complex) * freedom_level);

        memset(q->alphas, 0, sizeof(Complex) * freedom_level);

        return q;
}

DYNAMIC_LIB_DECORATION void qcs_delete_qudit(tf_qcs_qudit *q)
{
    if(q==NULL) return ;

    if(q->alphas!=NULL)
    {
        free(q->alphas);
        q->alphas = NULL;
    }

    free(q);
    q=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_set_random_real_state_qudit(tf_qcs_qudit *q)
{
    int i;
    double sum = 0, tmp_r;

    for(i=0;i<q->freedom_level;i++)
    {
        tmp_r = qcs_randu();
        sum = sum + ( tmp_r * tmp_r );
        //printf("%lf\n", tmp_r);

        q->alphas[i].re=(tf_qcs_real_number)tmp_r;
        q->alphas[i].im=0;
    }
    for(i=0;i<q->freedom_level;i++)
    {
        q->alphas[i].re= q->alphas[i].re / sqrtf(sum);
    }
}

DYNAMIC_LIB_DECORATION void qcs_set_state_n_qudit(tf_qcs_qudit *q, int s)
{
    if(q!=NULL)
    {
        memset(q->alphas, 0, sizeof(Complex) * q->freedom_level);
        q->alphas[s].re=1;
    }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qudit(tf_qcs_qudit *q)
{
    int i;

    if(q!=NULL)
    {
        for(i=0;i<q->freedom_level;i++)
        {
#ifdef PYTHON_SCRIPT
     		PySys_WriteStdout("(%2.6f + %2.6fi) |{%d}_dec> \n",   q->alphas[i].re, q->alphas[i].im, i);
#else
     		printf("(%2.6f + %2.6fi) |{%d}_dec> \n", q->alphas[i].re, q->alphas[i].im, i);
#endif
        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qudit_sqr(tf_qcs_qudit *q)
{
    int i;
    tf_qcs_real_number me, sum = 0;
    if(q!=NULL)
    {
        for(i=0;i<q->freedom_level;i++)
        {
            qcs_mod_complex(&q->alphas[i], &me);
            sum+=(me*me);
#ifdef PYTHON_SCRIPT
     		PySys_WriteStdout("(%-f |{%d}_dec> \n",   me * me, i);
#else
     		printf("(%-f |{%d}_dec> \n",   me * me, i);
#endif
        }
        printf("sum=%f\n", sum);
    }
}

DYNAMIC_LIB_DECORATION tf_qcs_qudit_array* qcs_new_qudit_array(int n, int freedom_level)
{
    tf_qcs_qudit_array* t;
    int i, j;

    t = (tf_qcs_qudit_array*)malloc(sizeof(tf_qcs_qudit_array));
    t->arr = (tf_qcs_qudit*)malloc(n*sizeof(tf_qcs_qudit));
    t->n = n;

    for( i=0 ; i < n ; i++ )
    {
        (t->arr+i)->alphas = (Complex*) malloc(sizeof(Complex) * freedom_level);
        (t->arr+i)->freedom_level = freedom_level;
        for(j=0;j<freedom_level;j++)
        {
            (t->arr+i)->alphas[j].re=0;
            (t->arr+i)->alphas[j].im=0;
        }
    }

    return t;
}

DYNAMIC_LIB_DECORATION void qcs_set_qudit_array_n(tf_qcs_qudit_array *qa, int n, tf_qcs_qudit *q)
{
    int k;

    for(k=0;k<q->freedom_level;k++)
    {
        (qa->arr + n)->alphas[k].re=q->alphas[k].re;
        (qa->arr + n)->alphas[k].im=q->alphas[k].im;
    }
}

DYNAMIC_LIB_DECORATION tf_qcs_qudit *qcs_get_qudit_array_n(tf_qcs_qudit_array *qa, int n)
{
    return (qa->arr + n);
}

DYNAMIC_LIB_DECORATION void qcs_delete_qudit_array(tf_qcs_qudit_array *qa)
{
     free(qa->arr);
     qa->n=0;
     free(qa);
     qa=NULL;
}
