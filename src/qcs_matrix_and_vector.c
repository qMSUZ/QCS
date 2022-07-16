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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <complex.h>


#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif

#ifdef PYTHON_SCRIPT
#include <Python.h>
#endif

#ifdef PYTHON_SCRIPT
#define _PRINT  PySys_WriteStdout
#else
#define _PRINT  printf
#endif

#include "qcs.h"
#include "qcs_misc.h"
#include "qcs_complex.h"
#include "qcs_matrix_and_vector.h"

/* matrix any size */
int compare_tf_qcs_complex(const void *_a, const void *_b)
{
    const tf_qcs_complex *a = (tf_qcs_complex *)_a;
    const tf_qcs_complex *b = (tf_qcs_complex *)_b;

    if(a->re > b->re) return 1;
    if(a->re < b->re) return -1;

    return 0;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix(int rows, int cols)
{
    tf_qcs_matrix *t = NULL;

    t = (tf_qcs_matrix*)malloc(sizeof(tf_qcs_matrix));

    t->rows = rows;
    t->cols = cols;

    t->m = (tf_qcs_complex*)malloc(sizeof(tf_qcs_complex)*rows*cols);

    t->q = -1;
    t->freedom_level = -1;

    memset(t->m, 0, sizeof(tf_qcs_complex)*rows*cols);

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_eye_matrix(int size)
{
    tf_qcs_matrix *t;

    t = qcs_create_matrix(size, size);

    qcs_eye_matrix( t );

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_for_quantum_state(int q, int freedom_level)
{
    tf_qcs_matrix *t;

    t = qcs_create_matrix(powf(freedom_level, q), powf(freedom_level, q));

    t->q = q;
    t->freedom_level = freedom_level;

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_empty_matrix()
{
    tf_qcs_matrix *t;

    t = (tf_qcs_matrix*)malloc(sizeof(tf_qcs_matrix));

    t->rows = 0 ;
    t->cols = 0 ;

    t->m = NULL;

    return t;
}

tf_qcs_matrix *qcs_create_vector_column(int rows)
{
    tf_qcs_matrix *t;

    t = qcs_create_matrix( rows, 1 );

    return t;
}

tf_qcs_matrix *qcs_create_random_vector_column(int rows)
{
    tf_qcs_matrix *t;
    int i;

    t = qcs_create_matrix( rows, 1 );

    for(i=0;i<rows;i++) 
    {
        qcs_set_cell_at_matrix_direct(t, i, 0, qcs_randn(), qcs_randn());
    }

    return t;
}

tf_qcs_matrix *qcs_create_random_vector_column_normalised(int rows)
{
    tf_qcs_matrix *t;
    tf_qcs_real_number total_val;
    tf_qcs_complex v;
    int i;

    t = qcs_create_matrix( rows, 1 );

    for(i=0;i<rows;i++) 
    {
        qcs_set_cell_at_matrix_direct(t, i, 0, qcs_randn(), qcs_randn());
    }

    qcs_calculate_norm1_of_column_vector(t, &total_val);

    for(i=0;i<rows;i++)
    {
        v = *qcs_get_cell_at_matrix_complex(t, i, 0);
        
        v.re = v.re / total_val;
        v.im = v.im / total_val;

        qcs_set_cell_at_matrix_complex(t, i, 0, &v);
    }

    return t;

}

tf_qcs_matrix *qcs_create_random_real_vector_column(int rows)
{
    tf_qcs_matrix *t;
    int i;

    t = qcs_create_matrix( rows, 1 );

    for(i=0;i<rows;i++) 
    {
        qcs_set_cell_at_matrix_direct(t, i, 0, qcs_randn(), 0);
    }

    return t;
}

tf_qcs_matrix *qcs_create_random_real_vector_column_normalised(int rows)
{
    tf_qcs_matrix *t;
    tf_qcs_real_number total_val=(tf_qcs_real_number)0.0, v;
    int i;

    t = qcs_create_matrix( rows, 1 );

    for(i=0;i<rows;i++) 
    {
        qcs_set_cell_at_matrix_direct(t, i, 0, qcs_randn(), 0);
    }

    qcs_calculate_norm1_of_column_vector(t, &total_val);

    for(i=0;i<rows;i++)
    {
        v=qcs_get_cell_at_matrix_complex(t, i, 0)->re;
        v=v / total_val;
        qcs_set_cell_at_matrix_direct(t, i, 0, v, 0);
    }

    return t;    
}

/*
U=zeros(d^N,d^N);
% Create a random d^N x d^N unitary
% from d^N orthogonal vectors
for k=1:d^N
    %vv=randn(d^N,1)+i*randn(d^N,1);
	vv=randn(d^N,1);
    for m=1:k-1
        vv=vv-U(:,m)*(U(:,m)'*vv);
    end %for
    U(:,k)=vv/sqrt(vv'*vv);
end %for
*/

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_unitary_random_real_matrix(int size)
{
    tf_qcs_matrix *t;
    tf_qcs_real_number *vv, r;
    int k,i,m;

    vv = malloc(sizeof(tf_qcs_real_number)*size);
    t = qcs_create_matrix(size, size);


    for(k=0;k<size;k++)
    {
        for(i=0;i<size;i++) { vv[i]=(qcs_randn()); }
        //for(i=0;i<size;i++) { printf("vv[%d]=%f ", i,vv[i]); } printf("\n");
        //if(k==0) {vv[0]=0.42408; vv[1]=0.98382;}
        //if(k==1) {vv[0]=0.19568; vv[1]=0.32126;}

        for(m=0;m<=(k-1);m++)
        {
            r=0; for(i=0;i<size;i++) { r+=qcs_get_cell_at_matrix_complex(t,i,m)->re * vv[i]; }
            //printf("m: r=%f\n", r);

            for(i=0;i<size;i++) { vv[i]=vv[i] - (qcs_get_cell_at_matrix_complex(t,i,m)->re * r); };
        }
        r=0; for(i=0;i<size;i++) { r=r+(vv[i] * vv[i]); } r=sqrt(r); //printf("k[%d]: r=%f\n", k, r);

        for(i=0;i<size;i++) qcs_set_cell_at_matrix_direct(t, i, k, vv[i] / r, 0);
    }

    free( vv );

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_unitary_random_matrix(int size)
{
    tf_qcs_matrix *u, *vv, *vv_out, *um;
    tf_qcs_complex dv, sqrt_dv, out;
    int i, k, m;

    u = qcs_create_matrix(size, size);

    vv = qcs_create_matrix(size, 1);
    um = qcs_create_matrix(size, 1);
    vv_out = qcs_create_matrix(size, 1);

    for(k=1; k<=size; k++)
    {
        qcs_fill_random_complex_values(vv);

        for(m=1; m<=k-1; m++)
        {
            qcs_get_column_from_matrix_into_vector(u, m-1, um);

            dv = qcs_calc_dot_of_two_vector( um, vv );

            for(i=0;i<size;i++)
            {
                qcs_complex_mul(um->m+i, &dv, vv_out->m+i);
                qcs_complex_sub(vv->m+i, vv_out->m+i, vv->m+i);
            }
        }
        dv = qcs_calc_dot_of_vector( vv );
        qcs_sqrt_complex( &dv, &sqrt_dv );

        for( i=0; i<size; i++ )
        {
            qcs_complex_div(vv->m+i, &sqrt_dv, &dv);
            qcs_set_cell_at_matrix_complex(u, i, k-1, &dv);
        }
    }

    qcs_delete_matrix(vv);
    qcs_delete_matrix(um);
    qcs_delete_matrix(vv_out);

    return u;
}

tf_qcs_matrix *qcs_create_unitary_random_real_matrix_by_qr(int size)
{
    return NULL;
}

tf_qcs_matrix *qcs_create_unitary_random_matrix_by_qr(int size)
{
    return NULL;
}


// % Generate a random matrix with elements with normal distributions
//      s=randn(d^N,d^N)+i*randn(d^N,d^N);
// % Create a positive matrix
//      r=s'*s;
// % Normalize it
//      r=r/trace(r);

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_density_random_real_matrix(int size)
{
    tf_qcs_matrix *tmp1 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp_out = qcs_create_matrix(size, size);
    tf_qcs_complex tr, v;

    int i,j;

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            qcs_set_cell_at_matrix_direct(tmp1, i, j, qcs_randn(), 0);
        }
    }

    tf_qcs_matrix *tmp2 = qcs_clone_matrix(tmp1);

    qcs_transpose_matrix(tmp1);

    qcs_mul_matrix(tmp1, tmp2, tmp_out);

    tr = qcs_trace_matrix(tmp_out);

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            v=*qcs_get_cell_at_matrix_complex(tmp_out, i,j);
            v.re=v.re / tr.re;
            qcs_set_cell_at_matrix_direct(tmp_out, i, j, v.re, 0);
        }
    }

    qcs_delete_matrix(tmp1);
    return tmp_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_density_random_matrix(int size)
{
    tf_qcs_matrix *tmp1 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp_out = qcs_create_matrix(size, size);
    tf_qcs_complex tr, v, o;

    int i,j;

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            qcs_set_cell_at_matrix_direct(tmp1, i, j, qcs_randn(), qcs_randn());
        }
    }

    tf_qcs_matrix *tmp2 = qcs_clone_matrix(tmp1);

    qcs_transpose_matrix(tmp1);

    qcs_mul_matrix(tmp1, tmp2, tmp_out);

    tr = qcs_trace_matrix(tmp_out);

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            v=*qcs_get_cell_at_matrix_complex(tmp_out, i,j);
            o.re=0; o.im=0;
            qcs_complex_div(&v,  &tr, &o);
            qcs_set_cell_at_matrix_direct(tmp_out, i, j, o.re, o.im);
        }
    }

    qcs_delete_matrix(tmp1);
    return tmp_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_hermitian_random_real_matrix(int size)
{
    tf_qcs_matrix *tmp1 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp2 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp_out = qcs_create_matrix(size, size);
    tf_qcs_complex tr, v;

    qcs_fill_random_real_values( tmp1 );
    tmp2=qcs_clone_matrix( tmp1 );
    qcs_transpose_matrix( tmp2 );

    qcs_add_matrix(tmp1, tmp2, tmp_out);

    qcs_delete_matrix(tmp1);
    qcs_delete_matrix(tmp2);

    return tmp_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_hermitian_random_matrix(int size)
{
    tf_qcs_matrix *tmp1 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp2 = qcs_create_matrix(size, size);
    tf_qcs_matrix *tmp_out = qcs_create_matrix(size, size);
    tf_qcs_complex tr, v;

    qcs_fill_random_complex_values( tmp1 );
    tmp2=qcs_clone_matrix( tmp1 );
    qcs_transpose_matrix( tmp2 );

    qcs_add_matrix(tmp1, tmp2, tmp_out);

    qcs_delete_matrix( tmp1 );
    qcs_delete_matrix( tmp2 );

    return tmp_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_arange_operation(tf_qcs_real_number _start, tf_qcs_real_number _end, tf_qcs_real_number _step)
{
    tf_qcs_matrix *m;
    int n_count = ceil ((_end - _start) / _step) ;
    int idx = 0;

#if 0
    printf("qcs_create_matrix_arange_operation:\n");
    printf("\t n_count = %d\n", n_count);
#endif

    m = qcs_create_matrix(1, n_count);

    tf_qcs_real_number v = _start;
    while( v < _end )
    {
        qcs_set_cell_at_matrix_direct( m, 0, idx, v, 0);
        v = v + _step;
        idx++;
    }

    return m;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_linspace_operation_with_endpoint(tf_qcs_real_number _start, tf_qcs_real_number _end, int _n_count)
{
    tf_qcs_matrix *m;
    tf_qcs_real_number _step = (_end - _start) / (_n_count-1)*1.0;
    int idx=0;

    m = qcs_create_matrix(1, _n_count);

    tf_qcs_real_number v = 0.0f;
    while( idx < _n_count )
    {
        qcs_set_cell_at_matrix_direct( m, 0, idx, (v * _step) + _start, 0);
        v = v + 1.0;
        idx++;
    }

    return m;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_linspace_operation_without_endpoint(tf_qcs_real_number _start, tf_qcs_real_number _end, int _n_count)
{
    tf_qcs_matrix *m;
    tf_qcs_real_number _step = (_end - _start) / (_n_count)*1.0;
    int idx=0;

    m = qcs_create_matrix(1, _n_count);

    tf_qcs_real_number v = 0.0f;
    while( idx < _n_count )
    {
        qcs_set_cell_at_matrix_direct( m, 0, idx, (v * _step) + _start, 0);
        v = v + 1.0;
        idx++;
    }

    return m;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_arange_operation_with_float_args(float _start, float _end, float _step)
{
    return qcs_create_matrix_arange_operation( (tf_qcs_real_number)_start, (tf_qcs_real_number)_end, (tf_qcs_real_number)_step);
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_linspace_operation_with_endpoint_with_float_args(float _start, float _end, int _n_count)
{
    return qcs_create_matrix_linspace_operation_with_endpoint( (tf_qcs_real_number)_start, (tf_qcs_real_number)_end, _n_count);
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_create_matrix_linspace_operation_without_endpoint_with_float_args(float _start, float _end, int _n_count)
{
    return qcs_create_matrix_linspace_operation_without_endpoint( (tf_qcs_real_number)_start, (tf_qcs_real_number)_end, _n_count);
}

tf_qcs_matrix *qcs_create_create_operator(int N)
{
    int i,j;
    tf_qcs_matrix *m;

    m = qcs_create_matrix(N, N);

    i=0;
    for(j=1;j<m->cols;j++)
    {
            (m->m+i+(j*m->cols))->re = sqrt((double)j);
            (m->m+i+(j*m->cols))->im = 0;
            i++;
    }

    return m;
}

tf_qcs_matrix *qcs_create_destroy_operator(int N)
{
    int i,j;
    tf_qcs_matrix *m;

    m = qcs_create_matrix(N, N);

    i=0;
    for(j=1;j<m->cols;j++)
    {
            (m->m+j+(i*m->cols))->re = sqrt((double)j);
            (m->m+j+(i*m->cols))->im = 0;
            i++;
    }

    return m;
}


DYNAMIC_LIB_DECORATION void qcs_delete_matrix(tf_qcs_matrix *a_out)
{
    if(a_out->m != NULL) free((void*)a_out->m);
    a_out->m=NULL;

    (a_out)->rows=0;
    (a_out)->cols=0;

    if(a_out != NULL) free((void*)a_out);

    a_out=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_copy_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_out)
{
    int i,j;

    b_out->q = a_in->q;
    b_out->freedom_level = a_in->freedom_level;

    for(i=0;i<a_in->rows;i++)
    {
        for(j=0;j<a_in->cols;j++)
        {
            (b_out->m+j+(i*b_out->cols))->re = (a_in->m+j+(i*a_in->cols))->re;
            (b_out->m+j+(i*b_out->cols))->im = (a_in->m+j+(i*a_in->cols))->im;
        }
    }
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_clone_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_matrix *m;

    m = qcs_create_matrix(a_in->rows, a_in->cols);

    qcs_copy_matrix(a_in, m);

    return m;
}

DYNAMIC_LIB_DECORATION void qcs_eye_matrix(tf_qcs_matrix *a_in)
{
     int i,j;
     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             if(i!=j)
             {
                (a_in->m+j+(i*a_in->cols))->re=0;
                (a_in->m+j+(i*a_in->cols))->im=0;
             }
             else
             {
                (a_in->m+j+(i*a_in->cols))->re=1;
                (a_in->m+j+(i*a_in->cols))->im=0;
             }
         }
     }
}

DYNAMIC_LIB_DECORATION void qcs_eye_matrix_with_param(tf_qcs_matrix *a_in, tf_qcs_real_number p)
{
     int i,j;
     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             if(i!=j)
             {
                (a_in->m+j+(i*a_in->cols))->re=0;
                (a_in->m+j+(i*a_in->cols))->im=0;
             }
             else
             {
                (a_in->m+j+(i*a_in->cols))->re=1 * p;
                (a_in->m+j+(i*a_in->cols))->im=0;
             }
         }
     }
}

DYNAMIC_LIB_DECORATION void qcs_zero_matrix(tf_qcs_matrix *a_in)
{
     int i,j;

     /*for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             (a_in->m+j+(i*a_in->cols))->re=0;
             (a_in->m+j+(i*a_in->cols))->im=0;
         }
     }*/

	memset(a_in->m, 0, sizeof(tf_qcs_complex)*a_in->rows*a_in->cols);
}

DYNAMIC_LIB_DECORATION void qcs_fill_matrix(tf_qcs_matrix *a_in)
{
     int i,j,c=1;

     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
//             (a_in->m+j+(i*a_in->cols))->re=10+(1+j+(i)*10);
//             (a_in->m+j+(i*a_in->cols))->im=0;

//             (a_in->m+j+(i*a_in->cols))->re=i+1;
//             (a_in->m+j+(i*a_in->cols))->im=j+1;

             (a_in->m+j+(i*a_in->cols))->re=c+10;
             (a_in->m+j+(i*a_in->cols))->im=0;

             c++;
         }
     }
}

DYNAMIC_LIB_DECORATION void qcs_fill_random_real_values(tf_qcs_matrix *a_in)
{
    int i,j;

    for(i=0;i<a_in->rows;i++)
    {
        for(j=0;j<a_in->cols;j++)
        {
            (a_in->m+j+(i*a_in->cols))->re=qcs_randu();
            (a_in->m+j+(i*a_in->cols))->im=0;
        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_fill_random_complex_values(tf_qcs_matrix *a_in)
{
    int i,j;

    for(i=0;i<a_in->rows;i++)
    {
        for(j=0;j<a_in->cols;j++)
        {
            (a_in->m+j+(i*a_in->cols))->re=qcs_randn();
            (a_in->m+j+(i*a_in->cols))->im=qcs_randn();
        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_make_projector_from_vector(tf_qcs_matrix *vec_in, tf_qcs_matrix *prj_out)
{
    int size, i, j;
    tf_qcs_complex a, b, b_conj, out;

    size = vec_in->rows;

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            a = *qcs_get_cell_at_matrix_complex(vec_in, i, 0);

            b = *qcs_get_cell_at_matrix_complex(vec_in, j, 0);
            qcs_conj_complex(&b, &b_conj);

            qcs_complex_mul( &a, &b_conj, &out );

            qcs_set_cell_at_matrix_complex(prj_out, i, j, &out);

        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_normalise_real_matrix(tf_qcs_matrix *a_in)
{
}

DYNAMIC_LIB_DECORATION void qcs_roundn_matrix(tf_qcs_matrix *a_in, int n)
{
     int i,j;

     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             qcs_complex_roundn( a_in->m+j+(i*a_in->cols) , n);
         }
     }
}

DYNAMIC_LIB_DECORATION void qcs_chop_matrix(tf_qcs_matrix *a_in)
{
     int i,j;

     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
              if ( fabsf((a_in->m+j+(i*a_in->cols))->re) < QCS_EPS ) (a_in->m+j+(i*a_in->cols))->re = 0;
              if ( fabsf((a_in->m+j+(i*a_in->cols))->im) < QCS_EPS ) (a_in->m+j+(i*a_in->cols))->im = 0;
         }
     }
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_get_column_from_matrix(tf_qcs_matrix *a_in, int c)
{
    int i;
    tf_qcs_matrix *m = qcs_create_matrix(a_in->rows, 1);

    for(i=0; i<a_in->rows; i++)
    {
        qcs_set_cell_at_matrix_complex(m, i, 0, qcs_get_cell_at_matrix_complex(a_in, i, c));
    }

    return m;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qcs_get_row_from_matrix(tf_qcs_matrix *a_in, int r)
{
    int i;
    tf_qcs_matrix *m = qcs_create_matrix(1, a_in->cols);

    for(i=0; i<a_in->cols; i++)
    {
        qcs_set_cell_at_matrix_complex(m, 0, i, qcs_get_cell_at_matrix_complex(a_in, r, i));
    }

    return m;
}

DYNAMIC_LIB_DECORATION void qcs_get_column_from_matrix_into_vector(tf_qcs_matrix *a_in, int c, tf_qcs_matrix *vec_out)
{
    int i;

    for(i=0; i<a_in->rows; i++)
    {
        qcs_set_cell_at_matrix_complex(vec_out, i, 0, qcs_get_cell_at_matrix_complex(a_in, i, c));
    }
}


DYNAMIC_LIB_DECORATION void qcs_set_cell_at_matrix_complex(tf_qcs_matrix *a_in, int r, int c, tf_qcs_complex *b_in)
{
     (a_in->m + (r*a_in->cols) + c)->re = b_in->re;
     (a_in->m + (r*a_in->cols) + c)->im = b_in->im;
}

DYNAMIC_LIB_DECORATION void qcs_set_cell_at_matrix_direct(tf_qcs_matrix *a_in, int r, int c, tf_qcs_real_number re, tf_qcs_real_number im)
{
     (a_in->m + (r*a_in->cols) + c)->re = re;
     (a_in->m + (r*a_in->cols) + c)->im = im;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex* qcs_get_cell_at_matrix_complex(tf_qcs_matrix *a_in, int r, int c)
{
     return (a_in->m + (r*a_in->cols) + c);
}

DYNAMIC_LIB_DECORATION void qcs_add_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out)
{
    int i, j;

    c_out->q = a_in->q;
    c_out->freedom_level = a_in->freedom_level;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_add((a_in->m + (i*a_in->cols) + j),
            (b_in->m + (i*a_in->cols) + j),
            (c_out->m + (i*a_in->cols) + j));
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_sub_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out)
{
     int i, j;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_sub((a_in->m + (i*a_in->cols) + j),
            (b_in->m + (i*a_in->cols) + j),
            (c_out->m + (i*a_in->cols) + j));
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_mul_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out)
{
     int i, j, k;

     tf_qcs_complex num, num_tmp1, num_tmp2;

     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<b_in->cols;j++)
         {
             for(k=0;k<a_in->cols;k++)
             {
               num.re=(c_out->m + (i*c_out->cols) + j)->re;
               num.im=(c_out->m + (i*c_out->cols) + j)->im;

               qcs_complex_mul((a_in->m + (i*a_in->cols) + k),
                           (b_in->m + (k*b_in->cols) + j), &num_tmp2);

               num_tmp1.re=0;
               num_tmp1.im=0;

               qcs_complex_add(&num_tmp2, &num, &num_tmp1);

               num.re=num_tmp1.re;
               num.im=num_tmp1.im;

               (c_out->m + (i*c_out->cols) + j)->re=num.re;
               (c_out->m + (i*c_out->cols) + j)->im=num.im;
             }
         }
     }
}

// A = A + cB
void qcs_add_isM_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_in)
{
    int i, j;
    tf_qcs_complex tmp1, tmp2;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            tmp1.re=0; tmp1.im=0;
            qcs_complex_mul(b_in, (c_in->m + (i*a_in->cols) + j), &tmp1 );

            tmp2.re=0; tmp2.im=0;
            qcs_complex_add((a_in->m + (i*a_in->cols) + j), &tmp1, &tmp2);

            (a_in->m + (i*a_in->cols) + j)->re = tmp2.re;
            (a_in->m + (i*a_in->cols) + j)->im = tmp2.im;
        }
     }
}

// A = A - cB
void qcs_sub_isM_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_in)
{
    int i, j;
    tf_qcs_complex tmp1, tmp2;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            tmp1.re=0; tmp1.im=0;
            qcs_complex_mul(b_in, (c_in->m + (i*a_in->cols) + j), &tmp1 );

            tmp2.re=0; tmp2.im=0;
            qcs_complex_sub((a_in->m + (i*a_in->cols) + j), &tmp1, &tmp2);

            (a_in->m + (i*a_in->cols) + j)->re = tmp2.re;
            (a_in->m + (i*a_in->cols) + j)->im = tmp2.im;
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_add_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out)
{
     int i, j;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_add((a_in->m + (i*a_in->cols) + j),
                        b_in,
                        (c_out->m + (i*a_in->cols) + j));
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_sub_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out)
{
     int i, j;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_sub((a_in->m + (i*a_in->cols) + j),
                        b_in,
                        (c_out->m + (i*a_in->cols) + j));
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_mul_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out)
{
     int i, j;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_mul((a_in->m + (i*a_in->cols) + j),
                        b_in,
                        (c_out->m + (i*a_in->cols) + j));
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_div_scalar_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_in, tf_qcs_matrix *c_out)
{
     int i, j;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            qcs_complex_div((a_in->m + (i*a_in->cols) + j),
                        b_in,
                        (c_out->m + (i*a_in->cols) + j));
        }
     }
}

#if 0
extern void cgetrf_(int* M, int *N, Complex* A, int* LDA, int* IPIV, int* INFO);

int cgetrf(int M, int N, Complex* A, int LDA, int* IPIV)
{
    int INFO;

    cgetrf_ (&M, &N, A, &LDA, IPIV, &INFO);

    return INFO;
}

extern void cgetri_(int *N, Complex *A, int *LDA, int *IPIV, Complex *WORK, int *LWORK, int *INFO) ;

int cgetri(int N, Complex *A, int LDA, int *IPIV, Complex *WORK, int LWORK)
{
    int INFO;

    cgetri_( &N, A, &LDA, IPIV, WORK, &LWORK, &INFO);

    return INFO;
}
#endif

DYNAMIC_LIB_DECORATION void qcs_inv_matrix(tf_qcs_matrix *a_in)
{
    Complex *A, *WORK;

    int N, info, i, j, k, LDA, LWORK, *IPIV;

    N = a_in->rows;
    LDA=N;
    LWORK=N*N;

    IPIV = malloc(sizeof(int)*(N+1));
    memset(IPIV, 0, sizeof(int)*(N+1));

    WORK = malloc(sizeof(Complex)*LWORK);
    memset(WORK, 0, sizeof(Complex)*LWORK);

    A = a_in->m;

    //info = cgetrf(N, N, A, N, IPIV );
    //printf("cgetrf: info %d\n", info);

    //info = cgetri(N, A, LDA, IPIV, WORK, LWORK );
    //printf("cgetri: info %d\n", info);

    free(WORK);
    free(IPIV);
}

DYNAMIC_LIB_DECORATION void qcs_inv_2x2_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_matrix *invm;
    tf_qcs_complex a,b,c,d, ad, bc, adbc, one, r;

    one.re = 1.0;
    one.im = 0.0;

    a = *qcs_get_cell_at_matrix_complex(a_in, 0, 0);
    b = *qcs_get_cell_at_matrix_complex(a_in, 0, 1);
    c = *qcs_get_cell_at_matrix_complex(a_in, 1, 0);
    d = *qcs_get_cell_at_matrix_complex(a_in, 1, 1);

    qcs_complex_mul(&a, &d, &ad);
    qcs_complex_mul(&b, &c, &bc);
    qcs_complex_sub(&ad, &bc, &adbc);

    qcs_complex_div(&one, &adbc, &r);

    invm = qcs_create_matrix(2, 2);

    b.re=-1.0f * b.re;
    b.im=-1.0f * b.im;

    c.re=-1.0f * c.re;
    c.im=-1.0f * c.im;

    qcs_set_cell_at_matrix_complex(invm, 0, 0, &d);
    qcs_set_cell_at_matrix_complex(invm, 0, 1, &b);
    qcs_set_cell_at_matrix_complex(invm, 1, 0, &c);
    qcs_set_cell_at_matrix_complex(invm, 1, 1, &a);

    qcs_scalar_mul_matrix(invm, &r);

    qcs_set_cell_at_matrix_complex(a_in, 0, 0, qcs_get_cell_at_matrix_complex(invm, 0, 0));
    qcs_set_cell_at_matrix_complex(a_in, 0, 1, qcs_get_cell_at_matrix_complex(invm, 0, 1));
    qcs_set_cell_at_matrix_complex(a_in, 1, 0, qcs_get_cell_at_matrix_complex(invm, 1, 0));
    qcs_set_cell_at_matrix_complex(a_in, 1, 1, qcs_get_cell_at_matrix_complex(invm, 1, 1));


    qcs_delete_matrix( invm );
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_calc_dot_of_vector(tf_qcs_matrix *vec_in)
{
    int i;
    tf_qcs_complex f, a, out, r;

    r.re=0;
    r.im=0;

    for(i=0; i<vec_in->rows; i++)
    {
        qcs_conj_complex(vec_in->m+i, &f);
        qcs_complex_mul(&f, vec_in->m+i, &out);
        qcs_complex_add(&r, &out, &a);
        r=a;
    }

    return r;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_calc_dot_of_two_vector(tf_qcs_matrix *vec_in_1, tf_qcs_matrix *vec_in_2)
{
    int i;
    tf_qcs_complex f, a, out, r;

    r.re=0;
    r.im=0;

    for(i=0; i<vec_in_1->rows; i++)
    {
        qcs_conj_complex(vec_in_1->m+i, &f);
        qcs_complex_mul(&f, vec_in_2->m+i, &out);
        qcs_complex_add(&r, &out, &a);
        r=a;
    }

    return r;
}

DYNAMIC_LIB_DECORATION int qcs_non_zero_elements_of_matrix(tf_qcs_matrix *a_in)
{
     int i, j, n=0;

     for(i=0;i<a_in->rows;i++)
     {
        for(j=0;j<a_in->cols;j++)
        {
            if ( (a_in->m + (i*a_in->cols) + j)->re!=0 || (a_in->m + (i*a_in->cols) + j)->im!=0 )
                n++;
        }
     }
     return n;
}

DYNAMIC_LIB_DECORATION void qcs_det_matrix(tf_qcs_matrix *a_in, tf_qcs_complex *b_out)
{
}

DYNAMIC_LIB_DECORATION void qcs_transpose_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_complex a, b;
    int i,j;

    if(a_in->rows == a_in->cols)
    {
        for(i=0;i<a_in->rows;i++)
        {
            for(j=i;j<a_in->cols;j++)
            {
                a.re = qcs_get_cell_at_matrix_complex(a_in,i,j)->re;
                a.im = qcs_get_cell_at_matrix_complex(a_in,i,j)->im;

                qcs_conj_complex(&a, &a);
                b.re = qcs_get_cell_at_matrix_complex(a_in,j,i)->re;
                b.im = qcs_get_cell_at_matrix_complex(a_in,j,i)->im;

                qcs_conj_complex(&b, &b);
                qcs_set_cell_at_matrix_complex(a_in, i, j, &b);
                qcs_set_cell_at_matrix_complex(a_in, j, i, &a);
            }
        }
    }
    else
    {
        tf_qcs_matrix *n;

        n=qcs_create_matrix(a_in->cols, a_in->rows);

        for(i=0;i<a_in->rows;i++)
        {
            for(j=0;j<a_in->cols;j++)
            {
                a.re = qcs_get_cell_at_matrix_complex(a_in,i,j)->re;
                a.im = qcs_get_cell_at_matrix_complex(a_in,i,j)->im;

                qcs_conj_complex(&a, &a);

                qcs_set_cell_at_matrix_complex(n, j, i, &a);
            }
        }
        a_in->rows=n->rows;
        a_in->cols=n->cols;

        free((void*)a_in->m);

        a_in->m = n->m;
        n->m = NULL;

        qcs_delete_matrix(n);

    }
}

DYNAMIC_LIB_DECORATION void qcs_scalar_mul_matrix(tf_qcs_matrix *a_mat, tf_qcs_complex *b_in)
{
     int i, j;
        tf_qcs_complex t;

     for(i=0;i<a_mat->rows;i++)
     {
        for(j=0;j<a_mat->cols;j++)
        {
            qcs_complex_mul((a_mat->m + (i*a_mat->cols) + j),
                        b_in,
                        &t);
            *(a_mat->m + (i*a_mat->cols) + j) = t;
        }
     }
}

DYNAMIC_LIB_DECORATION void qcs_basic_sum_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out)
{
     int i, j, ii, jj, tcols, trows;

     tcols=a_in->cols + b_in->cols;
     trows=a_in->rows + b_in->rows;

/* fill matrix c with zeros */

     for(i=0;i<c_out->rows;i++)
     {
         for(j=0;j<c_out->cols;j++)
         {
             (c_out->m+j+(i*tcols))->re=0;
             (c_out->m+j+(i*tcols))->im=0;
         }
     }

/* copy matrix a to matrix c*/

     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             (c_out->m+j+(i*tcols))->re=(a_in->m+j+(i*a_in->cols))->re;
             (c_out->m+j+(i*tcols))->im=(a_in->m+j+(i*a_in->cols))->im;
         }
     }

/* copy matrix b to matrix c*/

     ii=0;
     jj=0;

     for(i=a_in->rows;i<trows;i++)
     {
         for(j=a_in->cols;j<tcols;j++)
         {
             (c_out->m+j+(i*tcols))->re=(b_in->m+jj+(ii*b_in->cols))->re;
             (c_out->m+j+(i*tcols))->im=(b_in->m+jj+(ii*b_in->cols))->im;
             jj++;
         }
         jj=0;
         ii++;
     }
}

tf_qcs_matrix *qcs_basic_sum(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in)
{
    int tcols, trows;

    tcols=a_in->cols + b_in->cols;
    trows=a_in->rows + b_in->rows;

    tf_qcs_matrix *out_m;

    out_m = qcs_create_matrix( a_in->rows * b_in->rows, a_in->cols * a_in->cols );

    qcs_basic_sum_matrix(a_in, b_in, out_m);

    return out_m;

}

DYNAMIC_LIB_DECORATION void qcs_tensor_matrix(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in, tf_qcs_matrix *c_out)
{
     int x,y,i,j,ii,jj;
     tf_qcs_complex num;

     assert(a_in!=NULL);
     assert(b_in!=NULL);
     assert(c_out!=NULL);

/* fill matrix c with zeros */

     for(i=0;i<c_out->rows;i++)
     {
         for(j=0;j<c_out->cols;j++)
         {
             (c_out->m+j+(i*c_out->cols))->re=0;
             (c_out->m+j+(i*c_out->cols))->im=0;
         }
     }

     ii=0;
     jj=0;
     for(x=0;x<a_in->rows;x++)
     {
         for(y=0;y<a_in->cols;y++)
         {
             ii=(x*b_in->rows);
             jj=(y*b_in->cols);
             for(i=0;i<b_in->rows;i++)
             {
                 for(j=0;j<b_in->cols;j++)
                 {
                     num.re=0;
                     num.im=0;

                     qcs_complex_mul(
                            (a_in->m+y+(x*a_in->cols)),
                            (b_in->m+j+(i*b_in->cols)),
                            &num);

                     (c_out->m+jj+(ii*c_out->cols))->re=num.re;
                     (c_out->m+jj+(ii*c_out->cols))->im=num.im;
                     jj++;
                 }
                ii++;
                jj-=b_in->cols;
             }
         }
     }
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix * qcs_tensor_product(tf_qcs_matrix *a_in, tf_qcs_matrix *b_in)
{
    tf_qcs_matrix *out_m;

    out_m = qcs_create_matrix( a_in->rows * b_in->rows, a_in->cols * a_in->cols );

    qcs_tensor_matrix( a_in, b_in, out_m );

    return out_m;
}

/* --------------------------------------------------------------------- */

DYNAMIC_LIB_DECORATION void qcs_partial_transpose_matrix(tf_qcs_matrix *a_in, int n)
{
    int is, js, i, j, sizeN, bitmask=0, shift, ii, delta=0;
    tf_qcs_complex ctmp;

    sizeN=floor(log2(a_in->rows)/log2(2)+0.5);

    delta=(pow(2,sizeN-(n+1)));

    shift=sizeN-n-1;
    bitmask=bitmask ^ (1 << (shift));

    //printf("n=%d sizeN=%d bitmask=%d shift=%d delta=%d\n", n, sizeN, bitmask, shift, delta);

    for(i=0;i<a_in->rows;i++)
    {
         for(j=0;j<a_in->cols;j++)
         {
              if( ((i & bitmask) >> shift)==0 && ((j & bitmask) >> shift)==1 )
              {
                   is=i+delta;
                   js=j-delta;

                   //printf("swap [%d, %d] - [%d, %d]\n", i,j,is,js);

                   ctmp.re=(a_in->m+j+(i*a_in->cols))->re;
                   ctmp.im=(a_in->m+j+(i*a_in->cols))->im;

                   (a_in->m+j+(i*a_in->cols))->re=(a_in->m+js+(is*a_in->cols))->re;
                   (a_in->m+j+(i*a_in->cols))->im=(a_in->m+js+(is*a_in->cols))->im;

                   (a_in->m+js+(is*a_in->cols))->re=ctmp.re;
                   (a_in->m+js+(is*a_in->cols))->im=ctmp.im;
              } // if( ((i & bitmask) >> shift)==0 && ((j & bitmask) >> shift)==1 )
         } // for(j=0;j<a_in->cols;j++)
    } // for(i=0;i<a_in->rows;i++)
}

// realize D * D' operations
tf_qcs_matrix *qcs_calculate_d_dot_dt_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_matrix *at_in, *a_out;

    at_in = qcs_clone_matrix(a_in);
    qcs_transpose_matrix( at_in );

    a_out=qcs_create_matrix(a_in->rows, at_in->cols);

    qcs_mul_matrix(a_in, at_in, a_out);

    qcs_delete_matrix( at_in);

    return a_out;
}

void qcs_matrix_realignment(tf_qcs_matrix *a_in, tf_qcs_matrix *out_in)
{
    int ox=0, oy=0, x, y, i, j, d, dsqr;
    tf_qcs_complex t;

   if(a_in->q == 2) // for q=2 and freedom > 1
    {
        d=a_in->freedom_level;
        dsqr=(d*d);

        x=0;
        y=0;
        for(i=0;i<d;i++)
        {
            for(j=0;j<d;j++)
            {
                for(y=0;y<dsqr;y++)
                {
                    ox=(i * d) + (y / d);
                    oy=(j * d) + (y % d);
                    //printf("{[%d,%d] <- [%d,%d]} ", ox, oy, y+1, x+1);
                    //t =
                    qcs_set_cell_at_matrix_complex(
                        out_in,
                        ox, oy, qcs_get_cell_at_matrix_complex(a_in, y, x));
                }
                //printf("\n");
                x++;
            }
            //printf("\n");
        }
    }
    else
    {
       if((a_in->q % 2) == 0)
        {
            //printf("other case (q %% 2 == 0)\n");
            d = (int)sqrtf(a_in->cols);
            dsqr = a_in->cols;
            //printf("d=%d\ndsqr=%d\n",d,dsqr);
            x=0;
            for(i=0;i<d;i++)
            {
                for(j=0;j<d;j++)
                {
                    //printf("[%3d,%3d] ", (i*d)+1, (j*d)+1 );
                    for(y=0;y<dsqr;y++)
                    {
                        ox=(i * d) + (y / d);
                        oy=(j * d) + (y % d);
                        //printf("{[%d,%d] <- [%d,%d]} ", ox+1, oy+1, y+1, x+1);
                        //t = *qcs_get_cell_at_matrix_complex(a_in, y, x);
                        //printf("[%d,%d] [%2.2f,%2.2f] ", y+1, x+1, t.re, t.im);

                        qcs_set_cell_at_matrix_complex(
                            out_in,
                            ox, oy, qcs_get_cell_at_matrix_complex(a_in, y, x));
                    }
                    //printf("\n");
                    x++;
                }
                //printf("\n");
            }
        } // if((a_in->q % 2) == 0)
    }
}

void qcs_matrix_realignment2(tf_qcs_matrix *a_in, tf_qcs_matrix *out_in)
{
    int ox=0, oy=0, x, y, i, j, d, dsqr;
    tf_qcs_complex t;

        if((a_in->q % 2) == 0)
        {
            printf("other case (q %% 2 == 0)\n");
            d = (int)sqrtf(a_in->cols);
            dsqr = a_in->cols;
            printf("d=%d\ndsqr=%d\n",d,dsqr);
            x=0;
            for(i=0;i<d;i++)
            {
                for(j=0;j<d;j++)
                {
                    //printf("[%3d,%3d] ", (i*d)+1, (j*d)+1 );
                    for(y=0;y<dsqr;y++)
                    {
                        ox=(i * d) + (y / d);
                        oy=(j * d) + (y % d);
                        //printf("{[%d,%d] <- [%d,%d]} ", ox+1, oy+1, y+1, x+1);
                        //t = *qcs_get_cell_at_matrix_complex(a_in, y, x);
                        //printf("[%d,%d] [%2.2f,%2.2f] ", y+1, x+1, t.re, t.im);

                        qcs_set_cell_at_matrix_complex(
                            out_in,
                            ox, oy, qcs_get_cell_at_matrix_complex(a_in, y, x));
                    }
                    //printf("\n");
                    x++;
                }
                //printf("\n");
            }
        } // if((a_in->q % 2) == 0)
}

DYNAMIC_LIB_DECORATION void qcs_swap_block_matrix(tf_qcs_matrix *a_in, int r1, int c1, int r2, int c2, int size_r, int size_c)
{
    int r,c;
    tf_qcs_complex t1,t2;

    for(r=0;r<size_r;r++)
    {
        for(c=0;c<size_c;c++)
        {
            t1=*qcs_get_cell_at_matrix_complex(a_in, r1+r, c1+c);
            t2=*qcs_get_cell_at_matrix_complex(a_in, r2+r, c2+c);
            qcs_set_cell_at_matrix_complex(a_in, r2+r, c2+c, &t1);
            qcs_set_cell_at_matrix_complex(a_in, r1+r, c1+c, &t2);
        }
    }
}

DYNAMIC_LIB_DECORATION void qcs_partial_transpose_matrix_qudit(tf_qcs_matrix *a_in, int i)
{

    int h, block, d, q, row, col, s_row, s_col, row_m, col_m;
    int max_row, max_col;

    d = a_in->freedom_level;
    q = a_in->q;

    h = (int)powf(d, q-i);
    block = (int)powf(d, q-(i-1));

   // printf("h=%d block=%d\n", h, block);

 /*   for(row=0;row<a_in->rows/h;row++)
    {
        printf("r %d c: ", row);
        for(col=0;col<a_in->cols/h;col++)
        {
            printf("%d", col);
            row_m = row % d;
            col_m = col % d;
            if((row_m) > (col_m))
            {
                //printf("s[%d,%d][%d,%d]", row_m, col_m, col_m, row_m);
                printf("s");
            }
            else
            {
                printf(" ");
            }
            printf(" ");
        }
        printf("\n");
    }*/

    max_row = (int)powf(d,q);
    max_col = (int)powf(d,q);

    for(row=0;row<max_row;row+=block)
    {
        for(col=0;col<max_col;col+=block)
        {
            for(s_row=0;s_row<block;s_row+=h)
            {
                for(s_col=0;s_col<block;s_col+=h)
                {
                    if(s_row > s_col)
                    {
                        //printf("s[%d,%d][%d,%d]", row+s_row, col+s_col, row+s_col, col+s_row);
                        qcs_swap_block_matrix(a_in, row+s_row, col+s_col, row+s_col, col+s_row, h, h);

                    }
                    else
                    {
                        //printf(" ");
                    }
                }
            }
        }
        //printf("\n");
    }
}

/*

    n qubit to traced out -- (0, 1, 2, 3, ...,q-1)

*/
DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_partial_trace_matrix_1_qubit(tf_qcs_matrix *a_in, int n)
{
     tf_qcs_matrix *c_out=NULL;
     tf_qcs_complex ctmp;
     int ival, jval, i,j,sizeN, traceoutN, bitmask=0, shift, ii;

     sizeN=floor(log2(a_in->rows)/log2(2)+0.5);
     traceoutN=(int)pow(2,sizeN-1);

     shift=sizeN-n-1;
     bitmask=bitmask ^ (1 << (shift));

     c_out=qcs_create_matrix(traceoutN, traceoutN);
     c_out->q = -1;
     c_out->freedom_level = 2;

     //printf("sizeN=%d zero_bitmask=%d shift=%d\n", sizeN, bitmask, shift);

     ii=0;
     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             ival=( i & bitmask ) >> shift;
             jval=( j & bitmask ) >> shift;
             if( ival==0 && jval==0)
             {
                 //printf("zero a%d%d\n", i, j);

                 (c_out->m+ii)->re=(a_in->m+j+(i*a_in->cols))->re;
                 (c_out->m+ii)->im=(a_in->m+j+(i*a_in->cols))->im;

                 ii++;

             }
         } // for(j=0;j<a_in->cols;j++)
     } // for(i=0;i<a_in->rows;i++)

     ii=0;
     for(i=0;i<a_in->rows;i++)
     {
         for(j=0;j<a_in->cols;j++)
         {
             ival=( i & bitmask ) >> shift;
             jval=( j & bitmask ) >> shift;
             if(ival==1 && jval==1)
             {
                 //printf("one a%d%d\n", i, j);

                 ctmp.re=(c_out->m+ii)->re;
                 ctmp.im=(c_out->m+ii)->im;

                 qcs_complex_add(&ctmp, (a_in->m+j+(i*a_in->cols)), (c_out->m+ii));

                 ii++;
             }
         } // for(j=0;j<a_in->cols;j++)
     } // for(i=0;i<a_in->rows;i++)


     return c_out;
}

/*

*/
DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_partial_trace_matrix_n_qubit(tf_qcs_matrix *a_in, int _from, int _to)
{
     tf_qcs_matrix *c_out=NULL;
     tf_qcs_complex ctmp;
     int ival, jval, i, j, sizeN, traceoutN, bitmask=0, shift, ii, n, q, maxq, bits;

     n=_to - _from +1;
     sizeN=floor(log2(a_in->rows)/log2(2)+0.5);
     traceoutN=(int)pow(2,sizeN-n);
     maxq=(int)pow(2, n);


     c_out=qcs_create_matrix(traceoutN, traceoutN);

     //printf("sizeN=%d traceoutN=%d maxq=%d bits=%d\n", sizeN, traceoutN, maxq, n);

     shift=((sizeN-1)-_to);

     for(q=(sizeN-1 - _from) ; q >= (sizeN-1 - _to) ; q--) bitmask=bitmask ^ (1 << q);

     //printf("bitmask=%d shift=%d\n", bitmask, shift);

     for(q=0;q<maxq;q++)
     {
         ii=0;

         for(i=0;i<a_in->rows;i++)
         {
             for(j=0;j<a_in->cols;j++)
             {
                  ival=(i & bitmask) >> shift;
                  jval=(j & bitmask) >> shift;
                  if(ival==q && jval==q)
                  {
                       //printf("[%d] a%d%d ", q, i, j);

                       ctmp.re=(c_out->m+ii)->re;
                       ctmp.im=(c_out->m+ii)->im;

                       qcs_complex_add(&ctmp, (a_in->m+j+(i*a_in->cols)), (c_out->m+ii));

                       ii++;
                  }
             } // for(j=0;j<a_in->cols;j++)
         } // for(i=0;i<a_in->rows;i++)
         //printf("\n");
     } // for(q=0;q<maxq;q++)

     return c_out;
}

// d -- fredom level
// q -- count of qudits
// i -- number of qudit traced out (1,2,3,...,q)
DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_partial_trace_matrix_1_qudit(tf_qcs_matrix *a_in, int i)
{
    tf_qcs_matrix *o;
    tf_qcs_complex tmp1, tmp2, tmp3;
    int d, q;
    int *indx, h, nbit, max_rows, number, ii,x,y;
    char txt[128];

    d = a_in->freedom_level;
    q = a_in->q;

    h = (int)powf(d,q-1);
    o=qcs_create_matrix(h, h);
    o->q = q-1;
    o->freedom_level = d;

    max_rows = (int)powf(d,q);

    //printf("h=%d\n", h);

    indx=(int*)malloc(sizeof(int)*d*h);

    ii=0;
    for(number=0;number<d;number++)
    {
        for(nbit=0;nbit<max_rows;nbit++)
        {
            qcs_dec2base_d(nbit, q, d, &txt[0]);

            if(qcs_numbase_d2int(txt[i-1]) == number)
            {
               //printf("%d [%s]", nbit, &txt[0]); printf(" \"%d\" at %d", number, nbit);printf("\n");
               indx[ii] = nbit;
               ii++;
            }
        }
    }

    for(x=0;x<h;x++)
    {
        for(y=0;y<h;y++)
        {
            //printf("set [%d,%d] with ",x,y);
            for(number=0;number<d;number++)
            {
                //printf("[%d,%d]", indx[x+number*h], indx[y+number*h]);
                tmp1=*qcs_get_cell_at_matrix_complex(a_in, indx[x+number*h], indx[y+number*h]);
                tmp2=*qcs_get_cell_at_matrix_complex(o, x, y);
                qcs_complex_add(&tmp1, &tmp2, &tmp3);
                qcs_set_cell_at_matrix_complex(o,x,y,&tmp3);
            }
            //printf(" ");
        }
        //printf(" \n");
    }

    free(indx);

    return o;
}

/* --------------------------------------------------------------------- */


extern void cheev_ (char *JOBZ, char *UPLO, int *N, tf_qcs_complex *A, int *LDA, tf_qcs_real_number *W, tf_qcs_complex *WORK, int *LWORK, tf_qcs_real_number *RWORK, int *INFO);

int cheev ( char JOBZ, char UPLO, int N, tf_qcs_complex *A, int LDA, tf_qcs_real_number *W, tf_qcs_complex *WORK, int LWORK, tf_qcs_real_number *RWORK)

{
  int INFO;

  cheev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &INFO);

  return INFO;
}

extern void cgesvd_( char *JOBU, char *JOBVT, int *M, int *N, tf_qcs_complex *A, int *LDA, tf_qcs_real_number *S, tf_qcs_complex *U, int *LDU, tf_qcs_complex *VT, int *LDVT, tf_qcs_complex *WORK, int *LWORK, tf_qcs_real_number *RWORK, int *INFO);

int cgesvd ( char JOBU, char JOBVT, int M, int N, tf_qcs_complex *A, int LDA, tf_qcs_real_number *S, tf_qcs_complex *U, int LDU, tf_qcs_complex *VT, int LDVT, tf_qcs_complex *WORK, int LWORK, tf_qcs_real_number *RWORK)

{
  int INFO;

  cgesvd_( &JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, &INFO);

  return INFO;
}

DYNAMIC_LIB_DECORATION void qcs_spectral_decompose_of_matrix(tf_qcs_matrix *a_mat, tf_qcs_matrix *eigenvalues, tf_qcs_matrix *eigenvectors)
{
    Complex *A, *WORK;
    tf_qcs_real_number *W, *RWORK;

    int N, info, i, j, k, LDA, LWORK;

    qcs_copy_matrix(a_mat, eigenvectors);

    N=eigenvectors->rows;
    LDA=N;
    LWORK=2*N-1;
    W = malloc(N*sizeof(tf_qcs_real_number));
    WORK = malloc(LWORK*sizeof(tf_qcs_complex));
    RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

    A=eigenvectors->m;

    info=0;
    info=cheev('V', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

    qcs_transpose_matrix(eigenvectors);
    //qcs_print_matrix( eigenvectors );
    //printf("\n");

    //eigenvalues=qcs_create_matrix(N,1);
    eigenvalues->q=-1;
    eigenvalues->freedom_level=-1;

    for(i=0;i<N;i++)
    {
        qcs_set_cell_at_matrix_direct(eigenvalues, i, 0, *(W+i), 0);
    }

    free(RWORK);
    free(WORK);
    free(W);
}

// sprawdzic czy b1 i b2 potrzebuja dodatkowej transpozycji
DYNAMIC_LIB_DECORATION void qcs_svd_decompose_of_matrix(tf_qcs_matrix *state, tf_qcs_matrix *out_coeff, tf_qcs_matrix *out_base1, tf_qcs_matrix *out_base2)
{
    int idx,i, j;
    int info = 0, lda, ldu, ldvt, lwork, m, n;
    tf_qcs_real_number *RWORK, *S;
    Complex *A, *U, *VT, *WORK;

    m = state->rows;
    n = state->cols;

    lda=m;
    A = state->m;

    S = malloc(qcs_min(m,n)*sizeof(tf_qcs_real_number));

    ldu = m;
    U = malloc(ldu*m*sizeof(tf_qcs_complex));

    ldvt = n;
    VT = malloc(ldvt * n * sizeof(tf_qcs_complex));

    lwork = 2 * qcs_min(m,n) + qcs_max(m,n);
    WORK = malloc(lwork * sizeof(tf_qcs_complex));

    RWORK = malloc(5*qcs_min(m,n)*sizeof(tf_qcs_real_number)); memset(RWORK, 0, 5*qcs_min(m,n)*sizeof(tf_qcs_real_number));

    //info = cgesvd( 'A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, WORK, lwork, RWORK);

    //printf("Info = %d\nSchmidt coff:\n", info);
    //for(i=0;i<qcs_min(m,n);i++) printf("-> %f\n", *(S+i));

    out_coeff->rows=1;
    out_coeff->cols=qcs_min(m,n);
    out_coeff->m=malloc(sizeof(tf_qcs_matrix)*qcs_min(m,n));

    for(i=0;i<qcs_min(m,n);i++)
    {
        qcs_set_cell_at_matrix_direct( out_coeff, 0, i, *(S+i), 0);
    }

    if(out_base1!=NULL)
    {
        out_base1->rows=ldu;
        out_base1->cols=m;
        out_base1->m=U;
        //qcs_print_matrix(&phantom);
        //printf("\n");
    }
    else
        free(U);

    if(out_base2!=NULL)
    {
        out_base2->rows=ldvt;
        out_base2->cols=n;
        out_base2->m=VT;
        //qcs_print_matrix(&phantom);
    }
    else
        free(VT);
    free(RWORK);
    free(WORK);
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_square_root_of_operator_matrix(tf_qcs_matrix *a_mat)
{
    tf_qcs_matrix *tmp_mat=NULL, *tmp1=NULL, *diag=NULL;

    Complex *A, *WORK;
    tf_qcs_real_number *W, *RWORK;
#if defined(__GNUC__)
    _Complex sqrt_cc, in_data;
#endif

    int N, info, i, j, k, LDA, LWORK;

    N=a_mat->rows;
    LDA=N;
    LWORK=2*N-1;
    W = malloc(N*sizeof(tf_qcs_real_number));
    WORK = malloc(LWORK*sizeof(tf_qcs_complex));
    RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

    A=a_mat->m;

    info=0;
    //info=cheev('V', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

    diag=qcs_create_matrix(a_mat->rows, a_mat->cols);

    for(i=0;i<N;i++)
    {
        if (*(W+i) > 0) qcs_set_cell_at_matrix_direct(diag, i, i, sqrt(*(W+i)), 0);
        if (*(W+i) == 0) qcs_set_cell_at_matrix_direct(diag, i, i, 0, 0);
#if defined(__GNUC__)
        if (*(W+i) < 0)
        {
           __real__ in_data=*(W+i);
           __imag__ in_data=0;
           sqrt_cc=csqrtf(in_data);
           qcs_set_cell_at_matrix_direct(diag, i, i, crealf(sqrt_cc), cimagf(sqrt_cc));
        }
#endif
    }

    free(RWORK);
    free(WORK);
    free(W);

    qcs_transpose_matrix(a_mat);

    tmp_mat=qcs_create_matrix(a_mat->rows, a_mat->cols);
    tmp1=qcs_create_matrix(a_mat->rows, a_mat->cols);

    qcs_mul_matrix(a_mat, diag, tmp1);
    qcs_transpose_matrix(a_mat);
    qcs_mul_matrix(tmp1, a_mat, tmp_mat);

    qcs_delete_matrix(tmp1);
    qcs_delete_matrix(diag);

    return tmp_mat;
}

DYNAMIC_LIB_DECORATION void qcs_square_root_of_operator_matrix_self(tf_qcs_matrix *a_mat)
{
    tf_qcs_matrix *sqrt_m;

    sqrt_m=qcs_square_root_of_operator_matrix(a_mat);

    if(a_mat->m!=NULL) free(a_mat->m);
    a_mat->m=sqrt_m->m;
    sqrt_m->m=NULL;

    qcs_delete_matrix(sqrt_m);
}

DYNAMIC_LIB_DECORATION void qcs_calculate_norm1_of_column_vector(tf_qcs_matrix *a_mat, tf_qcs_real_number *norm_val_out)
{
    tf_qcs_real_number tmpval=(tf_qcs_real_number)0.0, 
                       nval=(tf_qcs_real_number)0.0;
    int i;

    assert(a_mat!=NULL);
    assert(a_mat->cols == 1);

    for(i=0;i<a_mat->rows;i++)
    {
        qcs_mod_complex( qcs_get_cell_at_matrix_complex(a_mat, i, 0), &tmpval );
        nval = nval + tmpval;
    }

    *norm_val_out = nval;
}

DYNAMIC_LIB_DECORATION void qcs_calculate_norm2_of_column_vector(tf_qcs_matrix *a_mat, tf_qcs_real_number *norm_val_out)
{
    tf_qcs_real_number tmpval=(tf_qcs_real_number)0.0, 
                       nval=(tf_qcs_real_number)0.0;
    int i;

    assert(a_mat!=NULL);
    assert(a_mat->cols == 1);

    for(i=0;i<a_mat->rows;i++)
    {
        qcs_mod_complex( qcs_get_cell_at_matrix_complex(a_mat, i, 0), &tmpval );
        nval = nval + (tmpval*tmpval);
    }

    *norm_val_out = (tf_qcs_real_number)sqrt((double)nval);
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_norm_of_matrix(tf_qcs_matrix *a_mat)
{
    tf_qcs_matrix *tmp, *out;

    tmp = qcs_clone_matrix(a_mat);
    out = qcs_create_matrix(a_mat->rows, a_mat->cols);

    qcs_transpose_matrix(tmp);
    qcs_mul_matrix(tmp, a_mat, out);
    qcs_delete_matrix(tmp);

    qcs_square_root_of_operator_matrix_self(out);

    return out;
}

DYNAMIC_LIB_DECORATION void qcs_norm_of_matrix_self(tf_qcs_matrix *a_mat)
{
    tf_qcs_matrix *tmp, *out;

    tmp = qcs_clone_matrix(a_mat);
    out = qcs_create_matrix(a_mat->rows, a_mat->cols);

    qcs_transpose_matrix(tmp);
    qcs_mul_matrix(tmp, a_mat, out);

    qcs_copy_matrix(out, a_mat);

    qcs_delete_matrix(tmp);
    qcs_delete_matrix(out);

    qcs_square_root_of_operator_matrix_self(a_mat);
}

DYNAMIC_LIB_DECORATION void qcs_exp_of_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_matrix *a_out, *E, *D, *X, *cX, *Eshadow, *Dshadow, *Xshadow;


    tf_qcs_real_number f,v,c;
    tf_qcs_complex ctmp;
    int e,s,k,q,p;

    a_out = qcs_clone_matrix( a_in );

    v = qcs_infinity_norm_of_matrix( a_in ).re;
    f = (float)frexp( v, &e );
    s=qcs_max(0,e+1);

//    printf("v=%f\n", v);
//    printf("f=%f, e=%d, s=%d\n", f, e, s);

    ctmp.re = powf(2.0f,(float)s);
    ctmp.im=0;
    qcs_div_scalar_matrix(a_in, &ctmp, a_out );

    Xshadow = qcs_create_matrix( a_in->rows, a_in->cols);
    X = qcs_clone_matrix( a_out );
    c = 0.5f;

    Eshadow = qcs_create_matrix( a_in->rows, a_in->cols);
    E = qcs_create_matrix( a_in->rows, a_in->cols);
    qcs_eye_matrix( E );
    ctmp.re=c;
    ctmp.im=0;
    qcs_add_isM_matrix(E, &ctmp, a_out);

    Dshadow = qcs_create_matrix( a_in->rows, a_in->cols);
    D = qcs_create_matrix( a_in->rows, a_in->cols);
    qcs_eye_matrix( D );
    ctmp.re=c;
    ctmp.im=0;
    qcs_sub_isM_matrix(D, &ctmp, a_out);

    q = 6;
    p = 1;

/*    printf("X\n");
    qcs_print_matrix( X );

    printf("E\n");
    qcs_print_matrix( E );

    printf("D\n");
    qcs_print_matrix( D );*/

    cX = qcs_create_matrix( a_in->rows, a_in->cols);

    for(k=2; k<=q; k++)
    {
        c = c * (q-k+1) / (k*(2*q-k+1));
        //printf("k = %d, c = %f\n", k, c);

        // X = A * X
        qcs_zero_matrix(Xshadow);
        qcs_mul_matrix(a_out,X,Xshadow);
        qcs_copy_matrix(Xshadow, X);

        // cX = c * X
        ctmp.re=c;
        ctmp.im=0;
        qcs_zero_matrix(cX);
        qcs_mul_scalar_matrix(X, &ctmp, cX);

        // E = E + cX
        qcs_zero_matrix(Eshadow);
        qcs_add_matrix(E, cX, Eshadow);
        qcs_copy_matrix(Eshadow, E);

        if(p) {
            // D = D + cX
            qcs_zero_matrix(Dshadow);
            qcs_add_matrix(D, cX, Dshadow);
            qcs_copy_matrix(Dshadow, D);
        }
        else {
            // D = D - cX
            qcs_zero_matrix(Dshadow);
            qcs_sub_matrix(D, cX, Dshadow);
            qcs_copy_matrix(Dshadow, D);
        }

        p=!p;

    } // for(k=2; k<=q; k++)

/*    printf("D\n");
    qcs_print_matrix( D );

    printf("E\n");
    qcs_print_matrix( E );*/

    // E = D\E or E = inv(D) * E
    qcs_zero_matrix(Eshadow);
    qcs_inv_matrix( D );
    qcs_mul_matrix(D,E,Eshadow);
    qcs_copy_matrix(Eshadow, E);

    for(k=1;k<=s;k++)
    {
        qcs_zero_matrix(Eshadow);
        qcs_mul_matrix(E,E,Eshadow);
        qcs_copy_matrix(Eshadow, E);
    }

    qcs_copy_matrix(E, a_in);

    qcs_delete_matrix( cX );

    qcs_delete_matrix( D );
    qcs_delete_matrix( Dshadow );

    qcs_delete_matrix( E );
    qcs_delete_matrix( Eshadow );

    qcs_delete_matrix( X );
    qcs_delete_matrix( Xshadow );

    qcs_delete_matrix( a_out );
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_infinity_norm_of_matrix(tf_qcs_matrix *a)
{
    tf_qcs_complex t;
    tf_qcs_real_number rnorm, r_tmp;
    int r,c;

    t.re=0; t.im=0;

    for( r = 0 ; r < a->rows ; r++)
    {
        // calulate sum or row "r"
        rnorm=0;
        for( c = 0 ; c < a->cols; c++)
        {
            r_tmp=0;
            qcs_mod_complex(a->m + (r*a->cols) + c, &r_tmp);
            rnorm = rnorm + r_tmp;
        }

        if (rnorm > t.re) t.re=rnorm;

    }

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_entropy_of_matrix(tf_qcs_matrix *a)
{
	Complex *A, *WORK;
	tf_qcs_real_number *W, *RWORK, entropy;

	int N, info, i, j, LDA, LWORK;

    N=a->rows;
  	LDA=N;
  	LWORK=2*N-1;
  	W = malloc(N*sizeof(tf_qcs_real_number));
  	WORK = malloc(LWORK*sizeof(tf_qcs_complex));
  	RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

	A=a->m;

    info=0;
    //info=cheev('N', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

    entropy=0;
    for(i=0;i<N;i++)
	{
		if(W[i]!=0) entropy+=W[i]*(log(W[i])/log(2));
	}

	free(RWORK);
	free(WORK);
	free(W);

	return -entropy;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_linear_entropy_of_matrix(tf_qcs_matrix *a_in)
{
	tf_qcs_complex a;

	a=qcs_trace_square_matrix(a_in);

	return (4.0/3.0)*(1 - a.re);
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_negativity_of_matrix(tf_qcs_matrix *a)
{
	Complex *A, *WORK;
	tf_qcs_real_number *W, *RWORK, negativity;

	int N, info, i, j, LDA, LWORK;

    N=a->rows;
  	LDA=N;
  	LWORK=2*N-1;
  	W = malloc(N*sizeof(tf_qcs_real_number));
  	WORK = malloc(LWORK*sizeof(tf_qcs_complex));
  	RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

	A=a->m;

    info=0;
    //info=cheev('N', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

    negativity=0;
    for(i=0;i<N;i++)
	{
		if(W[i]<0) negativity+=W[i];
	}

	free(RWORK);
	free(WORK);
	free(W);

	return -negativity;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp;
    tf_qcs_matrix *m;

    tmp.re=0;
    tmp.im=0;

    m=qcs_create_matrix(sigma->rows, sigma->cols);

    qcs_square_root_of_operator_matrix_self(rho);

    qcs_mul_matrix(rho, sigma, m);
    qcs_zero_matrix(sigma);
    qcs_mul_matrix(m, rho, sigma);

    qcs_square_root_of_operator_matrix_self(sigma);
    tmp=qcs_trace_matrix(sigma);

    qcs_delete_matrix(m);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_square_of_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp;

    tmp.re=0;
    tmp.im=0;

    tmp = qcs_fidelity(sigma, rho);

    tmp.re=tmp.re * tmp.re;

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_super_fidelity(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp, tmp1, tmp2, tmp3;
    tf_qcs_matrix *m;

    tmp.re=0;  tmp.im=0;
    tmp1.re=0; tmp1.im=0;
    tmp2.re=0; tmp2.im=0;
    tmp3.re=0; tmp3.im=0;

    m=qcs_create_matrix(sigma->rows, sigma->cols);
    qcs_zero_matrix(m);
    qcs_mul_matrix(sigma, rho, m);

    tmp1 = qcs_trace_matrix( m );

    tmp2 = qcs_trace_square_matrix_fast( sigma );
    tmp3 = qcs_trace_square_matrix_fast( rho );

    tmp.re=tmp1.re + sqrtf( (1-tmp2.re)*(1-tmp3.re) );

    qcs_delete_matrix(m);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_trace_distance(tf_qcs_matrix *rho, tf_qcs_matrix *sigma)
{
    tf_qcs_complex final_out, tmp, one_per_two;
    tf_qcs_matrix *m_tmp1, *m_tmp2, *m_tmp3;

    tmp.re=0;
    tmp.im=0;

    final_out.re=0;
    final_out.im=0;

    one_per_two.re=0.5;
    one_per_two.im=0;

    m_tmp1=qcs_create_matrix(sigma->rows, sigma->cols);
    qcs_zero_matrix( m_tmp1);

    qcs_sub_matrix(rho, sigma, m_tmp1);

    m_tmp1 = qcs_norm_of_matrix( m_tmp1 );

    tmp=qcs_trace_matrix(m_tmp1);

    qcs_complex_mul(&one_per_two, &tmp, &final_out);

    return final_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_hilbert_schmidt_distance(tf_qcs_matrix *rho, tf_qcs_matrix *sigma)
{
    tf_qcs_complex tmp;
    tf_qcs_matrix *m_tmp1, *m_tmp2, *m_tmp3;

    tmp.re=0;
    tmp.im=0;

    m_tmp1=qcs_create_matrix(sigma->rows, sigma->cols);
    qcs_zero_matrix( m_tmp1);
    m_tmp3=qcs_create_matrix(sigma->rows, sigma->cols);
    qcs_zero_matrix( m_tmp3);

    qcs_sub_matrix(rho, sigma, m_tmp1);
    m_tmp2 = qcs_clone_matrix( m_tmp1 );
    qcs_transpose_matrix( m_tmp2 );

    qcs_mul_matrix(m_tmp2, m_tmp1, m_tmp3);

    tmp=qcs_trace_matrix(m_tmp3);

    tmp.re=sqrtf(tmp.re);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_bures_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp;

    tmp.re=0;
    tmp.im=0;

    tmp=qcs_fidelity(sigma, rho);
    tmp.re=sqrt(2-2*sqrt(tmp.re));

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_angle_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp;

    tmp.re=0;
    tmp.im=0;

    tmp=qcs_fidelity(sigma, rho);
    tmp.re=acos(sqrt(tmp.re));

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_sine_metric(tf_qcs_matrix *sigma, tf_qcs_matrix *rho)
{
    tf_qcs_complex tmp;

    tmp.re=0;
    tmp.im=0;

    tmp=qcs_fidelity(sigma, rho);
    tmp.re=sqrt(1 - tmp.re);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex* qcs_eigenvalues_of_matrix(tf_qcs_matrix *a_in)
{
    return NULL;
}

/* --------------------------------------------------------------------- */

DYNAMIC_LIB_DECORATION void qcs_add_noise_to_matrix(tf_qcs_matrix *a_in, tf_qcs_real_number n)
{
     // rho_noisy=(1-p)*rho+p*eye(2^N)/2^N;
     // N=log2(sx)/log2(d);
     // N=floor(N+0.5);
     int i,j, N;
     tf_qcs_complex tt, out;

     N=floor(log2(a_in->rows)/log2(2)+0.5);

     tt.re=1-n;
     tt.im=0;

    for (i=0; i<a_in->rows; i++) {
      for (j=0; j<a_in->cols; j++) {
            qcs_complex_mul((a_in->m+i*a_in->cols+j), &tt, &out);
            (a_in->m+i*a_in->cols+j)->re=out.re;
            (a_in->m+i*a_in->cols+j)->im=out.im;
        }
    }

     tt.re=n*(1.0/powf(2.0,(tf_qcs_real_number)N));
     tt.im=0;

    for(i=0; i<a_in->rows; i++) {
             qcs_complex_add((a_in->m+i*a_in->cols+i), &tt, &out);
            (a_in->m+i*a_in->cols+i)->re=out.re;
            (a_in->m+i*a_in->cols+i)->im=out.im;
    }
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_trace_matrix(tf_qcs_matrix *a_in)
{
      tf_qcs_complex f,t;

      f.re=0;
      f.im=0;

      int i;

      for(i=0;i<a_in->rows;i++)
      {
           qcs_complex_add(&f, qcs_get_cell_at_matrix_complex(a_in, i, i), &t);

           f.re=t.re;
           f.im=t.im;
      }

      return f;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_trace_square_matrix(tf_qcs_matrix *a_in) /* purity */
{
    tf_qcs_complex f;

    tf_qcs_matrix *t = qcs_create_matrix(a_in->rows, a_in->cols);

    qcs_mul_matrix(a_in, a_in, t);

    f = qcs_trace_matrix( t );

    qcs_delete_matrix( t );

    return f;
}

DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_trace_square_matrix_fast(tf_qcs_matrix *a_in) /* purity */
{
     tf_qcs_complex two_tmp_const, a_tmp, b_tmp, c_tmp, d_tmp, e_tmp, out_tmp;
     int x,y;

     two_tmp_const.re=2;
     two_tmp_const.im=0;

     out_tmp.re=0;
     out_tmp.im=0;

     // first calc the square of a purity a11^2 + b22^2 + c22^2+ ....
     for(x=0;x<a_in->rows;x++)
     {
         a_tmp.re=qcs_get_cell_at_matrix_complex(a_in, x, x)->re;
         a_tmp.im=qcs_get_cell_at_matrix_complex(a_in, x, x)->im;

         b_tmp.re=a_tmp.re;
         b_tmp.im=a_tmp.im;
         qcs_complex_mul(&a_tmp, &b_tmp, &c_tmp);

         d_tmp.re=out_tmp.re;
         d_tmp.im=out_tmp.im;
         qcs_complex_add(&d_tmp, &c_tmp, &out_tmp);
     }

     // second calc the rest of a purity
     for(x=0;x<a_in->rows;x++)
     {
         for(y=x+1;y<a_in->rows;y++)
         {
                 a_tmp.re=qcs_get_cell_at_matrix_complex(a_in, x, y)->re;
                 a_tmp.im=qcs_get_cell_at_matrix_complex(a_in, x, y)->im;

                 b_tmp.re=qcs_get_cell_at_matrix_complex(a_in, y, x)->re;
                 b_tmp.im=qcs_get_cell_at_matrix_complex(a_in, y, x)->im;

                 qcs_complex_mul(&a_tmp, &b_tmp, &c_tmp);
                 qcs_complex_mul(&two_tmp_const, &c_tmp, &d_tmp);

                 e_tmp.re=out_tmp.re;
                 e_tmp.im=out_tmp.im;
                 qcs_complex_add(&e_tmp, &d_tmp, &out_tmp);
         }
     }

     return out_tmp;
}


DYNAMIC_LIB_DECORATION tf_qcs_complex qcs_calculate_sum_of_matrix(tf_qcs_matrix *a_in)
{
    tf_qcs_complex a, b, c;
    int i,j;

    a.re=0;
    a.im=0;
    b.re=0;
    b.im=0;
    c.re=0;
    c.im=0;

    for(i=0;i<a_in->rows;i++)
    {
        for(j=0;j<a_in->cols;j++)
        {
            b=*qcs_get_cell_at_matrix_complex(a_in, i, j);
            qcs_complex_add(&c,&b,&a);
            c=a;
        }
    }
    return c;
}

/* --------------------------------------------------------------------- */

DYNAMIC_LIB_DECORATION void qcs_print_eigenvalue_of_matrix(tf_qcs_matrix *a)
{
	Complex *A, *WORK;
	tf_qcs_real_number *W, *RWORK;

	int N, info, i, j, LDA, LWORK;

    N=a->rows;
  	LDA=N;
  	LWORK=2*N-1;
  	W = malloc(N*sizeof(tf_qcs_real_number));
  	WORK = malloc(LWORK*sizeof(tf_qcs_complex));
  	RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

	A=a->m;

    info=0;
    //info=cheev('N', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

    for(i=0;i<N;i++)
	{
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("lambda %d = %f\n", i, W[i]);
#else
     printf("lambda %d = %f\n", i, W[i]);
#endif
	}

	free(RWORK);
	free(WORK);
	free(W);
}

DYNAMIC_LIB_DECORATION void qcs_make_eigenvectors_of_matrix(tf_qcs_matrix *a)
{
	Complex *A, *WORK;
	tf_qcs_real_number *W, *RWORK;

	int N, info, i, j, LDA, LWORK;

    N=a->rows;
  	LDA=N;
  	LWORK=2*N-1;
  	W = malloc(N*sizeof(tf_qcs_real_number));
  	WORK = malloc(LWORK*sizeof(tf_qcs_complex));
  	RWORK = malloc((3*N-2)*sizeof(tf_qcs_real_number));

	A=a->m;

    info=0;
    //info=cheev('V', 'U', N, A, LDA, W, WORK, LWORK, RWORK);

	free(RWORK);
	free(WORK);
	free(W);
}

/* --------------------------------------------------------------------- */

DYNAMIC_LIB_DECORATION void qcs_print_matrix(tf_qcs_matrix *a_in)
{
    int i, j;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++)
    {
        for (j=0; j<a_in->cols; j++)
        {
            //qcs_print_complex(a_in->m+i*a_in->cols+j);
            tf_qcs_real_number re=(a_in->m+(i*a_in->cols)+j)->re;
            tf_qcs_real_number im=(a_in->m+(i*a_in->cols)+j)->im;
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("(%3.3f+%3.3fi) ", re, im);
#else
            printf("(%3.3f+%3.3fi) ", re, im);
#endif
        }
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("\n");
#else
     printf("\n");
#endif
    }
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_only_real_part(tf_qcs_matrix *a_in)
{
    int i, j;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++)
    {
        for (j=0; j<a_in->cols; j++)
        {
            //qcs_print_complex(a_in->m+i*a_in->cols+j);
            tf_qcs_real_number re=(a_in->m+i*a_in->cols+j)->re;
            //tf_qcs_real_number im=(a_in->m+i*a_in->cols+j)->im;
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("%2.2f ", re);
#else
            printf("%2.2f ", re);
#endif
        }
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("\n");
#else
     printf("\n");
#endif
    }
}


DYNAMIC_LIB_DECORATION void qcs_print_matrix_dot(tf_qcs_matrix *a_in, int bits, int base)
{
    int i, j, x;
    char tmp[512];

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for(x=0;x<bits+7;x++)
#ifdef PYTHON_SCRIPT
        PySys_WriteStdout(" ");
#else
        printf(" ");
#endif


    for (i=0; i<a_in->rows; i++)
    {
        qcs_dec2base_d(i, bits, base, &tmp[0]);
#ifdef PYTHON_SCRIPT
        PySys_WriteStdout("%s ", &tmp[0]);
#else
        printf("%s ", &tmp[0]);
#endif
    }
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("\n");
#else
     printf("\n");
#endif
    for (i=0; i<a_in->rows; i++)
    {
        qcs_dec2base_d(i, bits, base, &tmp[0]);
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("(%3d) %s ", i ,&tmp[0]);
#else
        printf("(%3d) %s ", i ,&tmp[0]);
#endif
        for (j=0; j<a_in->cols; j++)
        {
            //qcs_print_complex_dot(a_in->m+i*a_in->cols+j);
            if((a_in->m+i*a_in->cols+j)->re == 1.0)
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout("1");
#else
                printf("1");
#endif
            else
#ifdef PYTHON_SCRIPT
                PySys_WriteStdout(".");
#else
                printf(".");
#endif
            for(x=0;x<bits;x++)
#ifdef PYTHON_SCRIPT
        PySys_WriteStdout(" ");
#else
        printf(" ");
#endif
        }
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("\n");
#else
     printf("\n");
#endif
    }
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_in_matlab_format(tf_qcs_matrix *a_in)
{
    int i,j;
    tf_qcs_complex r;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

//#ifdef PYTHON_SCRIPT
//     PySys_WriteStdout("[ \n");
//#else
     _PRINT("[ \n");
//#endif
    for (i=0; i<a_in->rows; i++) {
      for (j=0; j<a_in->cols; j++) {
            //qcs_print_complex(a_in->m+i*a_in->cols+j);
            r = *qcs_get_cell_at_matrix_complex(a_in, i, j);
//#ifdef PYTHON_SCRIPT
//            PySys_WriteStdout("(%2.7f+%2.7fi) ", r.re, r.im);
//#else
            _PRINT("(%2.7f+%2.7fi) ", r.re, r.im);
//#endif

        }
//#ifdef PYTHON_SCRIPT
//     PySys_WriteStdout("; ... \n");
//#else
     _PRINT("; ... \n");
//#endif
    }
//#ifdef PYTHON_SCRIPT
//     PySys_WriteStdout("] \n");
//#else
     _PRINT("] \n");
//#endif
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_in_matlab_format_with_prefix(tf_qcs_matrix *a_in, char *prefix_str)
{
    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    _PRINT("%s ", prefix_str);
    qcs_print_matrix_in_matlab_format( a_in );
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_in_mathematica_format(tf_qcs_matrix *a_in)
{
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("please write me \n");
#else
     printf("please write me \n");
#endif
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_to_file(tf_qcs_matrix *a_in, FILE *b_out)
{
    int i,j;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++) {
        for (j=0; j<a_in->cols; j++) {
            qcs_print_complex_to_file(a_in->m+i*a_in->cols+j, b_out);
        }
        fprintf(b_out,"\n");
    }
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_to_file_sqr(tf_qcs_matrix *a_in, FILE *b_out)
{
    int i,j;
    tf_qcs_real_number l;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++) {
        for (j=0; j<a_in->cols; j++) {
            qcs_mod_complex(a_in->m+i*a_in->cols+j, &l);
            fprintf(b_out, "%2.2f ",l);
        }
        fprintf(b_out,"\n");
    }
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_to_file_sqr_01(tf_qcs_matrix *a_in, FILE *b_out)
{
    int i,j;
    tf_qcs_real_number l;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++) {
        for (j=0; j<a_in->cols; j++) {
            qcs_mod_complex(a_in->m+i*a_in->cols+j, &l);
            if(l==0) fprintf(b_out, "0");
            if(l==1) fprintf(b_out, "1");
        }
        fprintf(b_out,"\n");
    }
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix_to_file_sqr_mathematica(tf_qcs_matrix *a_in, FILE *b_out)
{
    int i,j;
    tf_qcs_real_number l;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    fprintf(b_out, "{");
    for ( i=0; i<a_in->rows; i++ )
    {
        fprintf(b_out, "{");
        for ( j=0; j<a_in->cols; j++ )
        {
            qcs_mod_complex( a_in->m + i * a_in->cols + j, &l );
            // fprintf(b_out, "%f2.2",l);
            if( l==0 ) fprintf(b_out, "0,");
            if( l==1 ) fprintf(b_out, "1,");
        }
        if(i != a_in->rows-1)
              fprintf(b_out,"\b},\n");
        else
              fprintf(b_out,"\b}\n");
    }
    fprintf(b_out, "}\n");
}

// TODO (marek#1#): dokonczyc implementacje qcs_print_matrix_to_file_sqr_matlab
DYNAMIC_LIB_DECORATION void qcs_print_matrix_to_file_sqr_matlab(tf_qcs_matrix *a_in, FILE *b_out)
{
    int i,j;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++) {
        for (j=0; j<a_in->cols; j++) {
           // qcs_print_complex_to_file(a_in->m+i*a_in->cols+j, b_out);
           fprintf(b_out, "%2.17f + %2.17fi\t", (a_in->m+i*a_in->cols+j)->re, (a_in->m+i*a_in->cols+j)->im);
        }
        fprintf(b_out,";\n");
    }
}

DYNAMIC_LIB_DECORATION void qcs_save_matrix_as_flat_text_file_for_gnuplot(tf_qcs_matrix *a_in, char *fname_out)
{
    int i,j;
    tf_qcs_real_number mod_value;
    FILE *f;

    f=fopen(fname_out,"w");
    fprintf(f, "{");
    for (i=0; i<a_in->rows; i++)
    {
        for (j=0; j<a_in->cols; j++)
        {
           qcs_mod_complex( (a_in->m+i*a_in->cols+j), &mod_value);
           fprintf(f, "{%d, %d, %2.3f},\n", i, j, mod_value);
        }
    }
    fprintf(f, "}");
    fclose( f );
}

DYNAMIC_LIB_DECORATION tf_qcs_schmidt_decomposition* qcs_create_schmidt_decomposition_empty()
{
    tf_qcs_schmidt_decomposition *tmp;

    tmp=malloc(sizeof(tf_qcs_schmidt_decomposition));

    tmp->base1 = NULL;
    tmp->base2 = NULL;

    tmp->schmidt_coeff = NULL;

    return tmp;
}

tf_qcs_schmidt_decomposition* qcs_create_schmidt_decomposition(tf_qcs_matrix *b1, tf_qcs_matrix *b2, tf_qcs_matrix *sc)
{
    tf_qcs_schmidt_decomposition *tmp;

    tmp=malloc(sizeof(tf_qcs_schmidt_decomposition));

    tmp->base1 = b1;
    tmp->base2 = b2;

    tmp->schmidt_coeff = sc;

    return tmp;
}


DYNAMIC_LIB_DECORATION void qcs_delete_schmidt_decomposition( tf_qcs_schmidt_decomposition *sd)
{
    if(sd->base1 != NULL)
    {
        qcs_delete_matrix(sd->base1);
        sd->base1 = NULL;
    }

    if(sd->base2 != NULL)
    {
        qcs_delete_matrix(sd->base2);
        sd->base2 = NULL;
    }

    if(sd->schmidt_coeff != NULL)
    {
        qcs_delete_matrix(sd->schmidt_coeff);
        sd->schmidt_coeff = NULL;

    }
    free(sd);
    sd=NULL;
}

int qcs_compare_schmidt_decompositions(tf_qcs_schmidt_decomposition *a, tf_qcs_schmidt_decomposition *b, int test_type)
{
    int i, j, a_non_zeros, b_non_zeros;
    tf_qcs_complex *t, a_sum, b_sum, tmp;
    tf_qcs_real_number l1, l2;

    a_non_zeros = 0;
    b_non_zeros = 0;

    a_sum.re=0;
    a_sum.im=0;

    b_sum.re=0;
    b_sum.im=0;

    for(i=0;i<a->schmidt_coeff->rows;i++)
    {
        for(j=0;j<a->schmidt_coeff->cols;j++)
        {
            t=qcs_get_cell_at_matrix_complex(a->schmidt_coeff, i, j);
            if(t->re!=0 || t->im!=0)
            {
                a_non_zeros++;
                qcs_complex_add(t, &a_sum, &tmp);
                a_sum = tmp;
            }
        }
    }

    for(i=0;i<b->schmidt_coeff->rows;i++)
    {
        for(j=0;j<b->schmidt_coeff->cols;j++)
        {
            t=qcs_get_cell_at_matrix_complex(b->schmidt_coeff, i, j);
            if(t->re!=0 || t->im!=0)
            {
                b_non_zeros++;
                qcs_complex_add(t, &b_sum, &tmp);
                b_sum = tmp;
            }
        }
    }

    //printf("qcs_compare_schmidt_decompositions debug: %d %d\n", a_non_zeros, b_non_zeros);

    if(a_non_zeros > b_non_zeros)
        return 1;
    if(b_non_zeros > a_non_zeros)
        return 2;

    if( a_non_zeros == b_non_zeros && test_type == 0 )
    {
        if(a_sum.re > b_sum.re) return 100;
        if(a_sum.re < b_sum.re) return 200;
    }

    if( a_non_zeros == b_non_zeros && test_type == 1 )
    {
        l1=qcs_sqr(qcs_get_cell_at_matrix_complex(a->schmidt_coeff, 0, 0)->re);
        l2=qcs_sqr(qcs_get_cell_at_matrix_complex(b->schmidt_coeff, 0, 0)->re);
        //printf("qcs_compare_schmidt_decompositions debug: %f %f %d %d\n", l1, l2, 0, a_non_zeros);
        i=1;
        while ( (l1 <= l2) && (i < a_non_zeros) )
        {
            t=qcs_get_cell_at_matrix_complex(a->schmidt_coeff, 0, i);
            l1=l1+qcs_sqr(t->re);

            t=qcs_get_cell_at_matrix_complex(b->schmidt_coeff, 0, i);
            l2=l2+qcs_sqr(t->re);

            //printf("qcs_compare_schmidt_decompositions debug: %f %f %d %d\n", l1, l2, i, a_non_zeros);
            i++;
        }

        if(l1 <= l2)
            return 10;
        else
            return -1;
    }

    return -1;
}

DYNAMIC_LIB_DECORATION tf_qcs_spectral_decomposition* qcs_create_spectral_decomposition()
{
    tf_qcs_spectral_decomposition *tmp;

    tmp=malloc(sizeof(tf_qcs_spectral_decomposition));

    tmp->eigenvectors = NULL;
    tmp->eigenvalues = NULL;

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_spectral_decomposition( tf_qcs_spectral_decomposition *sd)
{
    if(sd->eigenvectors != NULL)
    {
        qcs_delete_matrix(sd->eigenvectors);
        sd->eigenvectors = NULL;
    }

    if(sd->eigenvalues != NULL)
    {
        qcs_delete_matrix(sd->eigenvalues);
        sd->eigenvalues = NULL;

    }

    free(sd);

    sd=NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_svd_decomposition* qcs_create_svd_decomposition()
{
    tf_qcs_svd_decomposition *tmp;

    tmp=malloc(sizeof(tf_qcs_svd_decomposition));

    tmp->U = NULL;
    tmp->S = NULL;
    tmp->V = NULL;

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_svd_decomposition( tf_qcs_svd_decomposition *sd)
{
    if(sd->U != NULL)
    {
        qcs_delete_matrix( sd->U );
        sd->U = NULL;
    }

    if(sd->S != NULL)
    {
        qcs_delete_matrix( sd->S );
        sd->S = NULL;
    }

    if(sd->V != NULL)
    {
        qcs_delete_matrix( sd->V );
        sd->V = NULL;
    }


    free(sd);

    sd=NULL;
}


DYNAMIC_LIB_DECORATION int qcs_matrix_ispure(tf_qcs_matrix *a)
{
    int ispure = -1 ;

    tf_qcs_complex trace, trace_sqr;

    trace = qcs_trace_matrix(a);

    trace_sqr = qcs_trace_square_matrix_fast(a);

    if(trace.re == trace_sqr.re)
        ispure = 1;
    else
        ispure = 0;

    return ispure;
}

/*************************************************************************************/

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_mixing_two_density_matrix(tf_qcs_real_number p1, tf_qcs_matrix *m1, tf_qcs_real_number p2, tf_qcs_matrix* m2)
{
    tf_qcs_matrix* c_out = NULL;

    int i, j;

    c_out = qcs_create_matrix(m1->rows, m1->cols);

    c_out->q = m1->q;
    c_out->freedom_level = m1->freedom_level;

     for(i=0;i<m1->rows;i++)
     {
        for(j=0;j<m1->cols;j++)
        {
            //qcs_complex_add((a_in->m + (i*a_in->cols) + j),(b_in->m + (i*a_in->cols) + j),

            (c_out->m + (i*m1->cols) + j)->re= p1*((m1->m + (i*m1->cols) + j)->re) + p2*(m2->m + (i*m2->cols) + j)->re;
            (c_out->m + (i*m1->cols) + j)->im= p1*((m1->m + (i*m1->cols) + j)->im) + p2*(m2->m + (i*m2->cols) + j)->im;
        }
     }


    return c_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_horodecky_9x9_state(tf_qcs_real_number a)
{
    tf_qcs_matrix* tmp;
    tf_qcs_real_number scale;

    tmp = qcs_create_matrix(9,9);

    tmp->q=2;
    tmp->freedom_level=3;

    scale = 1.0 / (8*a + 1.0);

    qcs_set_cell_at_matrix_direct(tmp, 0, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 6, 6, scale*((1+a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 6, 8, scale*(sqrt(1-a*a)/2.0), 0);

    qcs_set_cell_at_matrix_direct(tmp, 7, 7, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 6, scale*(sqrt(1-a*a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8, scale*((1+a)/2.0), 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_horodecky_9x9_state_with_param(tf_qcs_real_number a, tf_qcs_real_number p)
{
    tf_qcs_matrix* tmp;
    tf_qcs_real_number scale;

    tmp = qcs_create_matrix(9,9);

    tmp->q=2;
    tmp->freedom_level=3;

    scale = 1.0 / (8*a + 1.0);

    qcs_set_cell_at_matrix_direct(tmp, 0, 0, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, p*scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, p*scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, p*scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5, p*scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 6, 6, p*scale*((1+a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 6, 8, p*scale*(sqrt(1-a*a)/2.0), 0);

    qcs_set_cell_at_matrix_direct(tmp, 7, 7, p*scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, p*scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 6, p*scale*(sqrt(1-a*a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8, p*scale*((1+a)/2.0), 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_horodecki_4x4_state(tf_qcs_real_number p, tf_qcs_real_number a, tf_qcs_real_number b)
{
    tf_qcs_matrix* tmp;

    tmp = qcs_create_matrix(4,4);

    tmp->q=2;
    tmp->freedom_level=2;


    qcs_set_cell_at_matrix_direct(tmp, 0, 0, p*a*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 3, p*a*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, (1-p)*a*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 1, 2, (1-p)*a*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 2, 1, (1-p)*a*b, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, (1-p)*b*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 3, 0, p*a*b, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, p*b*b, 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_ha_9x9_state(tf_qcs_real_number g)
{
    tf_qcs_matrix* tmp;
    tf_qcs_real_number n, a, b;

    tmp = qcs_create_matrix(9, 9);
    tmp->q=2;
    tmp->freedom_level=3;

    a=(1.0/3.0)*(g*g + 2);
    b=(1.0/3.0)*((1.0/(g*g)) + 2);
    n=1.0 / (7 + g*g + ((1.0)/(g*g)));

    qcs_set_cell_at_matrix_direct(tmp, 0, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, n*1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, n*b, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, n*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, n*1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5, n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 6, 6, n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 7, 7, n*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8, n*1, 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_ha_9x9_state_with_param(tf_qcs_real_number g, tf_qcs_real_number p)
{
    tf_qcs_matrix* tmp;
    tf_qcs_real_number n, a, b;

    tmp = qcs_create_matrix(9, 9);
    tmp->q=2;
    tmp->freedom_level=3;

    a=(1.0/3.0)*(g*g + 2);
    b=(1.0/3.0)*((1.0/(g*g)) + 2);
    n=1.0 / (7 + g*g + ((1.0)/(g*g)));

    qcs_set_cell_at_matrix_direct(tmp, 0, 0, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, p*n*1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, p*n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, p*n*b, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, p*n*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, p*n*1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5, p*n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 6, 6, p*n*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 7, 7, p*n*b, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, p*n*1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8, p*n*1, 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_w0_9x9_matrix()
{
    tf_qcs_matrix* tmp;

    tmp = qcs_create_matrix(9, 9);

    tmp->q=2;
    tmp->freedom_level=3;

    qcs_set_cell_at_matrix_direct(tmp, 0, 0,  1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, -1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, -1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1,  1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, -1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4,  1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, -1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5,  1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 6, 6,  1, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, -1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, -1, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8,  1, 0);

    return tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_maximally_mixed_state(int x)
{
    tf_qcs_matrix* tmp;
    int i;

    tmp = qcs_create_matrix(x,x);

    for(i=0;i<x;i++)
    {
        qcs_set_cell_at_matrix_direct(tmp, i, i, (1.0 / (tf_qcs_real_number)x), 0);
    }

    return tmp;
}


DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_maximally_mixed_state_with_param(int x, tf_qcs_real_number p)
{
    tf_qcs_matrix* tmp;
    int i;

    tmp = qcs_create_matrix(x,x);

    for(i=0;i<x;i++)
    {
        qcs_set_cell_at_matrix_direct(tmp, i, i, p * (1.0/(tf_qcs_real_number)x), 0);
    }

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_update_horodecky_9x9_state(tf_qcs_matrix *tmp, tf_qcs_real_number a)
{
    tf_qcs_real_number scale;

    scale = 1.0 / (8*a + 1.0);

    qcs_set_cell_at_matrix_direct(tmp, 0, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 0, 8, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 1, 1, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 2, 2, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 3, 3, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 4, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 4, 8, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 5, 5, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 6, 6, scale*((1+a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 6, 8, scale*(sqrt(1-a*a)/2.0), 0);

    qcs_set_cell_at_matrix_direct(tmp, 7, 7, scale*a, 0);

    qcs_set_cell_at_matrix_direct(tmp, 8, 0, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 4, scale*a, 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 6, scale*(sqrt(1-a*a)/2.0), 0);
    qcs_set_cell_at_matrix_direct(tmp, 8, 8, scale*((1+a)/2.0), 0);
}

DYNAMIC_LIB_DECORATION void qcs_update_horodecki_4x4_state(tf_qcs_matrix *m, tf_qcs_real_number p, tf_qcs_real_number a, tf_qcs_real_number b)
{
    qcs_set_cell_at_matrix_direct(m, 0, 0, p*a*a, 0);
    qcs_set_cell_at_matrix_direct(m, 0, 3, p*a*b, 0);

    qcs_set_cell_at_matrix_direct(m, 1, 1, (1-p)*a*a, 0);
    qcs_set_cell_at_matrix_direct(m, 1, 2, (1-p)*a*b, 0);

    qcs_set_cell_at_matrix_direct(m, 2, 1, (1-p)*a*b, 0);
    qcs_set_cell_at_matrix_direct(m, 2, 2, (1-p)*b*b, 0);

    qcs_set_cell_at_matrix_direct(m, 3, 0, p*a*b, 0);
    qcs_set_cell_at_matrix_direct(m, 3, 3, p*b*b, 0);
}

DYNAMIC_LIB_DECORATION void qcs_update_ha_9x9_state(tf_qcs_matrix *m, tf_qcs_real_number g)
{
    tf_qcs_real_number n, a, b;

    a=(1.0/3.0)*(g*g + 2);
    b=(1.0/3.0)*((1.0/(g*g)) + 2);
    n=1.0 / (7 + g*g + ((1.0)/(g*g)));

    qcs_set_cell_at_matrix_direct(m, 0, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 0, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 0, 8, n*1, 0);

    qcs_set_cell_at_matrix_direct(m, 1, 1, n*a, 0);
    qcs_set_cell_at_matrix_direct(m, 2, 2, n*b, 0);
    qcs_set_cell_at_matrix_direct(m, 3, 3, n*b, 0);

    qcs_set_cell_at_matrix_direct(m, 4, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 4, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 4, 8, n*1, 0);

    qcs_set_cell_at_matrix_direct(m, 5, 5, n*a, 0);
    qcs_set_cell_at_matrix_direct(m, 6, 6, n*a, 0);
    qcs_set_cell_at_matrix_direct(m, 7, 7, n*b, 0);

    qcs_set_cell_at_matrix_direct(m, 8, 0, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 8, 4, n*1, 0);
    qcs_set_cell_at_matrix_direct(m, 8, 8, n*1, 0);
}

/* entangled detection */

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_witnesses_application_to_matrix(tf_qcs_matrix *w, tf_qcs_matrix *m)
{
    tf_qcs_matrix *yy;
    tf_qcs_complex t;

    yy = qcs_create_matrix(w->rows, w->cols);

    qcs_mul_matrix(w, m, yy);

    t = qcs_trace_matrix( yy );

    qcs_delete_matrix(yy);

    return t.re;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_ppt_criterion(tf_qcs_matrix *m)
{
    tf_qcs_matrix *yy;
    tf_qcs_real_number t;

    yy = qcs_create_matrix(m->rows, m->cols);

    qcs_copy_matrix(m, yy);

    qcs_partial_transpose_matrix(yy, 1);
    t = qcs_negativity_of_matrix( yy );

    qcs_delete_matrix(yy);

    return t;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_ccnr_criterion(tf_qcs_matrix *m)
{
    tf_qcs_matrix *yy, *s;
    tf_qcs_complex t;

    yy = qcs_create_matrix(m->rows, m->cols);

    qcs_matrix_realignment(m, yy);

    s = qcs_create_empty_matrix();
    qcs_svd_decompose_of_matrix(yy, s, NULL, NULL);

    t = qcs_calculate_sum_of_matrix(s);

    qcs_delete_matrix(s);
    qcs_delete_matrix(yy);

    return t.re;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_depolarizing_channel_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;
    tf_qcs_real_number pp;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;


    if( m == 0 )
    {
        pp = sqrtf( 1.0 - (3.0*p)/4.0 );

        qcs_set_cell_at_matrix_direct(mat, 0, 0, pp, 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, pp, 0);
    }

    if( m == 1 )
    {
        pp=sqrtf(p/4.0);
        qcs_set_cell_at_matrix_direct(mat, 0, 1, pp, 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 0, pp, 0);
    }

    if( m == 2 )
    {
        pp=sqrtf(p/4.0);
        qcs_set_cell_at_matrix_direct(mat, 0, 1, 0, -pp);
        qcs_set_cell_at_matrix_direct(mat, 1, 0, 0,  pp);
    }

    if( m == 3 )
    {
        pp=sqrtf(p/4.0);
        qcs_set_cell_at_matrix_direct(mat, 0, 0, pp, 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, -pp, 0);
    }

    return mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_depolarizing_channel_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_amplitude_damping_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, 1, 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(1-p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 1, sqrtf(p), 0);
    }

    return mat;
}


DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_amplitude_damping_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_phase_damping_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, 1, 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(1-p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(p), 0);
    }

    return mat;
}

tf_qcs_matrix* qcs_create_phase_damping_v2_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(1+p)/sqrtf(2.0), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(1+p)/sqrtf(2.0), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(1-p)/sqrtf(2.0), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, -sqrtf(1-p)/sqrtf(2.0), 0);
    }

    return mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_phase_damping_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    return NULL;
}


DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_phase_flip_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(1-p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, -sqrtf(1-p), 0);
    }

    return mat;
}

#if 0
DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_phase_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    tf_qcs_complex pp;
    tf_qcs_matrix *z, *out_mat;

    out_mat = qcs_create_matrix(d,d);

    pp.re = sqrtf( qcs_newton_symbol(d-1.0,m) * powf((1.0-p)/2.0, m) * powf((1.0+p)/2.0, d-1.0-m) );
    pp.im=0;

    z = get_qudit_power_of_u_gate(get_qudit_pauli_z_d_gate( d ), m);
    qcs_mul_scalar_matrix(z, &pp, out_mat);

    qcs_delete_matrix(z);

    return out_mat;
}
#endif

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_bit_flip_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 1, sqrtf(1-p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 0, sqrtf(1-p), 0);
    }

    return mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_bit_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_bit_phase_flip_operator_for_qubit(int m, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 1, 0, -sqrtf(1-p));
        qcs_set_cell_at_matrix_direct(mat, 1, 0, 0,  sqrtf(1-p));
    }

    return mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_bit_phase_flip_operator_for_qudit(int m, int d, tf_qcs_real_number p)
{
    tf_qcs_matrix *mat;

    mat=qcs_create_matrix(2,2);
    mat->q=1;
    mat->freedom_level=2;

    if( m == 0 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 0, sqrtf(1-p), 0);
        qcs_set_cell_at_matrix_direct(mat, 1, 1, sqrtf(1-p), 0);
    }

    if( m == 1 )
    {
        qcs_set_cell_at_matrix_direct(mat, 0, 1, 0, -sqrtf(p));
        qcs_set_cell_at_matrix_direct(mat, 1, 0, 0,  sqrtf(p));
    }

    return mat;
}


/* fully correlated phase flip for qubit */

#if 0
DYNAMIC_LIB_DECORATION tf_qcs_matrix* qcs_create_fully_correlated_phase_flip_operator_for_qubit_register(int m, int q, tf_qcs_real_number p)
{
    int i;
    char gates[128];
    tf_qcs_matrix *mat;
    tf_qcs_complex c;


    if( m == 0 )
    {
        mat = qcs_create_matrix( (int)powf(2, q), (int)powf(2, q) );
        mat->q=q;
        mat->freedom_level=2;

        for(i=0;i<(int)powf(2, q);i++)
            qcs_set_cell_at_matrix_direct(mat, i, i, sqrtf(1-p), 0);
    }

    if( m == 1 )
    {
        for(i=0;i<q;i++) gates[i]='z';
        gates[i]=0;

        mat = qcs_build_matrix_from_tensor( gates );
        mat->q=q;
        mat->freedom_level=2;

        for(i=0;i<(int)powf(2, q);i++)
        {
            c = *qcs_get_cell_at_matrix_complex(mat, i, i);
            c.re = c.re * sqrtf(p) ;
            qcs_set_cell_at_matrix_direct(mat, i, i, c.re, 0);
        }
    }

    return mat;
}
#endif

/* four dimensional matrix of any size which contains complex numbers */

DYNAMIC_LIB_DECORATION tf_qcs_matrix4x4 *qcs_create_matrix4x4(int x1, int x2, int x3, int x4)
{
    tf_qcs_matrix4x4 *tmp=(tf_qcs_matrix4x4 *)malloc(sizeof(tf_qcs_matrix4x4));

    tmp->x1=x1;
    tmp->x2=x2;
    tmp->x3=x3;
    tmp->x4=x4;

    tmp->w4 = ( tmp->x4 ) ;
    tmp->w3 = ( tmp->x3 * tmp->w4) ;
    tmp->w2 = ( tmp->x2 * tmp->w3) ;
    tmp->w1 = ( tmp->x1 * tmp->w2) ;

    tmp->m=malloc(x1*x2*x3*x4*sizeof(tf_qcs_complex));

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_matrix4x4(tf_qcs_matrix4x4 *a_in)
{
    if(a_in->m != NULL)
    {
        free((void*)a_in->m);
        a_in->m = NULL;
    }

    free((void*)a_in);
    a_in = NULL;
}

DYNAMIC_LIB_DECORATION void qcs_set_cell_at_matrix4x4(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, tf_qcs_complex *m)
{
    int delta = x1*a_in->w2+x2*a_in->w3+x3*a_in->w4+x4;

    /*printf("%d delta [%d,%d,%d,%d] w1=%d w2=%d w3=%d w4=%d\n",
               delta,
               x1, x2, x3, x4,
               a_in->w1, a_in->w2, a_in->w3, a_in->w4);
    */
    *(a_in->m + delta) = *m;

}

DYNAMIC_LIB_DECORATION void qcs_set_cell_at_matrix4x4_direct(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, tf_qcs_real_number re, tf_qcs_real_number im)
{
    tf_qcs_complex m;

    m.re=re; m.im=im;

    qcs_set_cell_at_matrix4x4(a_in, x1,x2,x3,x4, &m);
}

DYNAMIC_LIB_DECORATION tf_qcs_complex* qcs_get_cell_at_matrix4x4(tf_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4)
{
    int delta = x1*a_in->w2 + x2*a_in->w3 + x3*a_in->w4 + x4;

    return (a_in->m + delta);
}

DYNAMIC_LIB_DECORATION void qcs_zero_matrix4x4(tf_qcs_matrix4x4 *a_in)
{
    memset(a_in->m, 0, a_in->x1 * a_in->x2 * a_in->x3 * a_in->x4 * sizeof(tf_qcs_complex));
}

DYNAMIC_LIB_DECORATION void qcs_print_matrix4x4(tf_qcs_matrix4x4 *a_in)
{
    int i,j,k,l;
    tf_qcs_complex a;

    for(i=0;i<a_in->x1;i++)
    {
        printf(":x1: %d\n", i);
        for(j=0;j<a_in->x2;j++)
        {
            printf("\t:x2: %d\n", j);
            for(k=0;k<a_in->x3;k++)
            {
                printf("\t\t:x3: %d\n", k);
                printf("\t\t\t");
                for(l=0;l<a_in->x4;l++)
                {
                    a = *qcs_get_cell_at_matrix4x4(a_in, i, j, k, l);
                    printf("(%2.2f,%2.2f) ", a.re, a.im);
                }
                printf("\n");
            }
        }
    }
}

/* matrix of any size which contains integer numbers */

DYNAMIC_LIB_DECORATION ti_qcs_matrix *qcs_create_int_matrix(int rows, int cols)
{
    ti_qcs_matrix *t;

    t=(ti_qcs_matrix*)malloc(sizeof(ti_qcs_matrix));

    t->rows=rows;
    t->cols=cols;

    t->m=(int*)malloc(sizeof(int)*rows*cols);

    memset(t->m, 0, sizeof(int)*rows*cols);

    return t;
}

DYNAMIC_LIB_DECORATION void qcs_delete_int_matrix(ti_qcs_matrix *a_out)
{
    if(a_out->m != NULL) free((void*)a_out->m);
    a_out->m=NULL;

    (a_out)->rows=0;
    (a_out)->cols=0;

    if(a_out != NULL) free((void*)a_out);

    a_out=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_set_cell_at_int_matrix(ti_qcs_matrix *a_in, int r, int c, int v)
{
     *(a_in->m + (r*a_in->cols) + c) = v;
}

DYNAMIC_LIB_DECORATION int qcs_get_cell_at_int_matrix(ti_qcs_matrix *a_in, int r, int c)
{
    return *(a_in->m + (r*a_in->cols) + c);
}

DYNAMIC_LIB_DECORATION void qcs_print_int_matrix(ti_qcs_matrix *a_in)
{
    int v, i, j;

    assert(a_in != NULL);
    assert(a_in->rows > 0 && a_in->cols > 0);

    for (i=0; i<a_in->rows; i++)
    {
        for (j=0; j<a_in->cols; j++)
        {
            //qcs_print_complex(a_in->m+i*a_in->cols+j);
            v=*(a_in->m+i*a_in->cols+j);
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("%d ", v);
#else
            printf("%d ", v);
#endif
        }
#ifdef PYTHON_SCRIPT
     PySys_WriteStdout("\n");
#else
     printf("\n");
#endif
    }
}

DYNAMIC_LIB_DECORATION ti_qcs_matrix4x4 *qcs_create_int_matrix4x4(int x1, int x2, int x3, int x4)
{
    ti_qcs_matrix4x4 *tmp=(ti_qcs_matrix4x4 *)malloc(sizeof(ti_qcs_matrix4x4));

    tmp->x1=x1;
    tmp->x2=x2;
    tmp->x3=x3;
    tmp->x4=x4;

    tmp->w4 = ( tmp->x4 ) ;
    tmp->w3 = ( tmp->x3 * tmp->w4) ;
    tmp->w2 = ( tmp->x2 * tmp->w3) ;
    tmp->w1 = ( tmp->x1 * tmp->w2) ;

    tmp->m=malloc(x1*x2*x3*x4*sizeof(tf_qcs_complex));

    return tmp;
}

DYNAMIC_LIB_DECORATION void qcs_delete_int_matrix4x4(ti_qcs_matrix4x4 *a_in)
{
    if(a_in->m != NULL)
    {
        free((void*)a_in->m);
        a_in->m = NULL;
    }

    free((void*)a_in);
    a_in=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_set_cell_at_int_matrix4x4(ti_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4, int v)
{
    int delta = x1*a_in->w2+x2*a_in->w3+x3*a_in->w4+x4;

    /**printf("%d delta [%d,%d,%d,%d] w1=%d w2=%d w3=%d w4=%d\n",
               delta,
               x1, x2, x3, x4,
               a_in->w1, a_in->w2, a_in->w3, a_in->w4);
    */
    *(a_in->m + delta) = v;
}

DYNAMIC_LIB_DECORATION int qcs_get_cell_at_int_matrix4x4(ti_qcs_matrix4x4 *a_in, int x1, int x2, int x3, int x4)
{
    int delta = x1*a_in->w2+x2*a_in->w3+x3*a_in->w4+x4;

    return *(a_in->m + delta);
}

DYNAMIC_LIB_DECORATION void qcs_zero_int_matrix4x4(tf_qcs_matrix4x4 *a_in)
{
}

DYNAMIC_LIB_DECORATION void qcs_print_int_matrix4x4(ti_qcs_matrix4x4 *a_in)
{
    int i,j,k,l;
    int a;

    for(i=0;i<a_in->x1;i++)
    {
        printf(":x1: %d\n", i);
        for(j=0;j<a_in->x2;j++)
        {
            printf("\t:x2: %d\n", j);
            for(k=0;k<a_in->x3;k++)
            {
                printf("\t\t:x3: %d\n", k);
                printf("\t\t\t");
                for(l=0;l<a_in->x4;l++)
                {
                    a = qcs_get_cell_at_int_matrix4x4(a_in, i, j, k, l);
                    printf("%d ", a);
                }
                printf("\n");
            }
        }
    }
}
