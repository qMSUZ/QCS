/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
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
#include <stdlib.h>
#include <string.h>

#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif

#include "qcs.h"
#include "qcs_matrix.h"
#include "qcs_qubit.h"
#include "qcs_rand.h"
#include "qcs_misc.h"

DYNAMIC_LIB_DECORATION tf_qcs_qubit* qcs_new_qubit()
{
    tf_qcs_qubit* t;

    t=(tf_qcs_qubit*)malloc(sizeof(tf_qcs_qubit));

    t->alpha.re=0;
    t->alpha.im=0;

    t->beta.re=0;
    t->beta.im=0;

    return t;
}

DYNAMIC_LIB_DECORATION void qcs_delete_qubit(tf_qcs_qubit *q)
{
    (q)->alpha.re=0;
    (q)->alpha.im=0;

    (q)->beta.re=0;
    (q)->beta.im=0;

    free((void*)q);
    q=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_make_projector_in_base(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_complex tmp;

    qcs_complex_mul(&base->v0_a , &base->v0_a, &tmp );
    base->gate_s0->m[0]=tmp;

    qcs_complex_mul(&base->v0_a , &base->v0_b, &tmp );
    base->gate_s0->m[2]=tmp;

    qcs_complex_mul(&base->v0_a , &base->v0_b, &tmp );
    base->gate_s0->m[1]=tmp;

    qcs_complex_mul(&base->v0_b , &base->v0_b, &tmp );
    base->gate_s0->m[3]=tmp;


    qcs_complex_mul(&base->v1_a , &base->v1_a, &tmp );
    base->gate_s1->m[0]=tmp;

    qcs_complex_mul(&base->v1_a , &base->v1_b, &tmp );
    base->gate_s1->m[2]=tmp;

    qcs_complex_mul(&base->v1_a , &base->v1_b, &tmp );
    base->gate_s1->m[1]=tmp;

    qcs_complex_mul(&base->v1_b , &base->v1_b, &tmp );
    base->gate_s1->m[3]=tmp;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_std_base()
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_std_base(b);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_hadamard_base()
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_hadamard_base(b);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_pauli_x_base()
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_pauli_x_base(b);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_pauli_y_base()
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_pauli_y_base(b);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_pauli_z_base()
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_pauli_z_base(b);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_m_base(tf_qcs_real_number theta)
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_m_base(b, theta);

    return b;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_base_desc *qcs_new_b_base(tf_qcs_real_number theta)
{
    tf_qcs_qubit_base_desc *b;

    b=(tf_qcs_qubit_base_desc *)malloc(sizeof(tf_qcs_qubit_base_desc));
    qcs_make_b_base(b, theta);

    return b;
}

DYNAMIC_LIB_DECORATION void qcs_set_base(tf_qcs_qubit_base_desc *base, tf_qcs_complex v0_a, tf_qcs_complex v0_b, tf_qcs_complex v1_a, tf_qcs_complex v1_b)
{
    tf_qcs_complex tmp;

    base->v0_a = v0_a;
    base->v0_b = v0_b;

    base->v1_a = v1_a;
    base->v1_b = v1_b;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    base->gate_s0->m[0]=base->v0_a;
    base->gate_s0->m[2]=base->v0_b;

    base->gate_s1->m[1]=base->v1_a;
    base->gate_s1->m[3]=base->v1_b;


    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_zero_base(tf_qcs_qubit_base_desc *base)
{
    base->v0_a.re = 0;
    base->v0_a.im = 0;

    base->v0_b.re = 0;
    base->v0_b.im = 0;

    base->v1_a.re = 0;
    base->v1_a.im = 0;

    base->v1_b.re = 0;
    base->v1_b.im = 0;

    base->gate_s0 = NULL;
    base->gate_s1 = NULL;
}

DYNAMIC_LIB_DECORATION void qcs_make_std_base(tf_qcs_qubit_base_desc *base)
{
    base->v0_a.re=1; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=0; base->v1_a.im=0;
    base->v1_b.re=1; base->v1_b.im=0;

    base->gate_s0=qcs_create_matrix(2,2);
    base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_hadamard_base(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_complex tmp;

    base->v0_a.re=0.707107; base->v0_a.im=0;
    base->v0_b.re=0.707107; base->v0_b.im=0;

    base->v1_a.re=0.707107;  base->v1_a.im=0;
    base->v1_b.re=-0.707107; base->v1_b.im=0;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_pauli_x_base(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_complex tmp;

    base->v0_a.re=1; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=0; base->v1_a.im=0;
    base->v1_b.re=1; base->v1_b.im=0;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_pauli_y_base(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_complex tmp;

    base->v0_a.re=0; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=1;

    base->v1_a.re=0; base->v1_a.im=-1;
    base->v1_b.re=0; base->v1_b.im=0;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s0 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_pauli_z_base(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_complex tmp;

    base->v0_a.re=1; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=0; base->v1_a.im=0;
    base->v1_b.re=-1; base->v1_b.im=0;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_m_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta)
{
    tf_qcs_complex tmp, tmp_e;

    base->v0_a.re=1; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=0; base->v1_a.im=0;
    base->v1_b.re=1; base->v1_b.im=0;

    tmp.re=0;
    tmp.im=theta;

    qcs_exp_complex(&tmp, &tmp_e);

    qcs_complex_mul(&tmp_e, &base->v1_b, &tmp);
    base->v1_b.re=tmp.re;
    base->v1_b.im=tmp.im;

    if ( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if ( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_b_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta)
{
    tf_qcs_complex tmp, tmp2, tmp_e, tmp_sqrt_e, tmp_one, tmp_zero;

    base->v0_a.re=1/sqrt(2); base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=1/sqrt(2); base->v1_a.im=0;
    base->v1_b.re=0; base->v1_b.im=0;

    tmp.re=theta;
    tmp.im=0;

    tmp_one.re = 0;
    tmp_one.im = 1;

    qcs_complex_mul(&tmp_one, &tmp, &tmp2);
    qcs_exp_complex(&tmp2, &tmp_e);
//    printf("debug -> %f+i%f\n", tmp_e.re, tmp_e.im);

    tmp_one.re=sqrt(2.0);
    tmp_one.im=0;

    qcs_complex_div(&tmp_e, &tmp_one, &tmp_sqrt_e);
//    printf("debug -> %f+i%f\n", tmp_sqrt_e.re, tmp_sqrt_e.im);


    base->v0_b.re = tmp_sqrt_e.re; base->v0_b.im = tmp_sqrt_e.im;
    base->v1_b.re = -tmp_sqrt_e.re; base->v1_b.im = -tmp_sqrt_e.im;

    if( base->gate_s0 == NULL ) base->gate_s0=qcs_create_matrix(2,2);
    if( base->gate_s1 == NULL ) base->gate_s1=qcs_create_matrix(2,2);

    qcs_make_projector_in_base(base);
}

DYNAMIC_LIB_DECORATION void qcs_make_b_base_float(tf_qcs_qubit_base_desc *base, float theta)
{
    qcs_make_b_base( base, theta);
}
// to delete
/*void qcs_update_m_base(tf_qcs_qubit_base_desc *base, tf_qcs_real_number theta)
{
    tf_qcs_complex tmp, tmp_e;

    base->v0_a.re=1; base->v0_a.im=0;
    base->v0_b.re=0; base->v0_b.im=0;

    base->v1_a.re=0; base->v1_a.im=0;
    base->v1_b.re=1; base->v1_b.im=0;

    tmp.re=0;
    tmp.im=theta;

    qcs_exp_complex(&tmp, &tmp_e);

    qcs_complex_mul(&tmp_e, &base->v1_b, &tmp);
    base->v1_b.re=tmp.re;
    base->v1_b.im=tmp.im;

    qcs_make_projector_in_base(base);
}
*/

DYNAMIC_LIB_DECORATION void qcs_base_print(tf_qcs_qubit_base_desc *base)
{
    tf_qcs_real_number fmod;

    printf("a=\n");

    qcs_mod_complex(&base->v0_a, &fmod);
    printf("%f %fi (|v0_a|=%f)\n", base->v0_a.re, base->v0_a.im, fmod*fmod);

    qcs_mod_complex(&base->v0_b, &fmod);
    printf("%f %fi (|v1_b|=%f)\n", base->v0_b.re, base->v0_b.im, fmod*fmod);

    printf("b=\n");

    qcs_mod_complex(&base->v1_a, &fmod);
    printf("%f %fi (|v1_a|=%f)\n", base->v1_a.re, base->v1_a.im, fmod*fmod);
    qcs_mod_complex(&base->v1_b, &fmod);
    printf("%f %fi (|v1_b|=%f)\n", base->v1_b.re, base->v1_b.im, fmod*fmod);

    printf("gate s0\n");
    qcs_print_matrix(base->gate_s0);
    printf("gate s1\n");
    qcs_print_matrix(base->gate_s1);
}

DYNAMIC_LIB_DECORATION void qcs_delete_base(tf_qcs_qubit_base_desc *base)
{
    if(base->gate_s0!=NULL)
    {
        qcs_delete_matrix(base->gate_s0);
    }

    if(base->gate_s1!=NULL)
    {
        qcs_delete_matrix(base->gate_s1);
    }

    free((void*)base);

    base=NULL;
}

DYNAMIC_LIB_DECORATION void qcs_set_random_real_state_qubit(tf_qcs_qubit *q)
{
    int i;
    double sum = 0, tmp_r;

    tmp_r = qcs_randu(); sum = sum + ( tmp_r * tmp_r );
    q->alpha.re=(tf_qcs_real_number)tmp_r;
    q->alpha.im=0;

    tmp_r = qcs_randu(); sum = sum + ( tmp_r * tmp_r );
    q->beta.re=(tf_qcs_real_number)tmp_r;
    q->beta.im=0;

    q->alpha.re= q->alpha.re / sqrtf(sum);
    q->beta.re= q->beta.re / sqrtf(sum);
}

DYNAMIC_LIB_DECORATION void qcs_set_superposition_state(tf_qcs_qubit *q)
{
    tf_qcs_real_number theta, phi;
    tf_qcs_complex a, b, img_one, half_phi, e, e_out, t1, t2;

    theta=(tf_qcs_real_number)qcs_randu() * QCS_PI;
    phi=(tf_qcs_real_number)qcs_randu() * 2 * QCS_PI;

    a.re=cos(theta/2);
    a.im=0;

    b.re=sin(theta/2);
    b.im=0;

    e.re=0;
    e.im=-1*(phi/2);

    qcs_exp_complex(&e, &e_out);

    qcs_complex_mul(&a, &e_out, &t1);
    qcs_complex_mul(&b, &e_out, &t2);

    q->alpha.re=t1.re;
    q->alpha.im=t1.im;

    q->beta.re=t2.re;
    q->beta.im=t2.im;
}

DYNAMIC_LIB_DECORATION void qcs_set_superposition_state_per(tf_qcs_qubit *q, tf_qcs_real_number per)
{
    qcs_set_superposition_state_theta_psi(q, (int)(180 * per), 0);
}

DYNAMIC_LIB_DECORATION void qcs_set_superposition_state_theta_psi(tf_qcs_qubit *q, int _theta, int _psi)
{
  tf_qcs_complex a, b, img_one, half_phi, e, e_out, t1, t2;

  tf_qcs_real_number theta, phi;

  theta=_theta * QCS_PI/180;
  phi=_psi * QCS_PI/180;

  a.re=cos(theta/2);
  a.im=0;

  b.re=sin(theta/2);
  b.im=0;

  e.re=0;
  e.im=-1*(phi/2);

  qcs_exp_complex(&e, &e_out);

  qcs_complex_mul(&a, &e_out, &t1);
  qcs_complex_mul(&b, &e_out, &t2);

 (q)->alpha.re=t1.re;
 (q)->alpha.im=t1.im;

 (q)->beta.re=t2.re;
 (q)->beta.im=t2.im;
}

DYNAMIC_LIB_DECORATION void qcs_set_state0_qubit(tf_qcs_qubit *q)
{
    q->alpha.re=1;
    q->alpha.im=0;

    q->beta.re=0;
    q->beta.im=0;
}

DYNAMIC_LIB_DECORATION void qcs_set_state1_qubit(tf_qcs_qubit *q)
{
    q->alpha.re=0;
    q->alpha.im=0;

    q->beta.re=1;
    q->beta.im=0;
}

DYNAMIC_LIB_DECORATION int qcs_qubit_measure(tf_qcs_qubit *q)
{
    tf_qcs_real_number r_num, m_num;

    r_num=rand()/(tf_qcs_real_number)RAND_MAX;
    qcs_mod_complex(&q->alpha, &m_num);

    if (r_num > m_num) {
      return 1;
    }
    else {
      return 0;
    }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qubit(tf_qcs_qubit *q)
{
     if(q!=NULL)
     {
#ifdef PYTHON_SCRIPT
     		PySys_WriteStdout("(%-f + %-fi) |0> + (%-f + %-fi) |1>\n",
			     q->alpha.re, q->alpha.im,
			     q->beta.re,  q->beta.im);
#else
     		printf("(%-f + %-fi) |0> + (%-f + %-fi) |1>\n",
			     q->alpha.re, q->alpha.im,
			     q->beta.re,  q->beta.im);
#endif
     }
     else
     {
#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("qubit is null\n");
#else
         printf("qubit is null\n");
#endif
     }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qubit_sqr(tf_qcs_qubit *q)
{
     tf_qcs_real_number a, b;

     if(q!=NULL)
     {
         qcs_mod_complex(&q->alpha,&a);
         qcs_mod_complex(&q->beta, &b);

#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("%2.5f |0> + %2.5f |1> (%2.5f)\n", a*a, b*b, (a*a)+(b*b));
#else
         printf("%2.5f |0> + %2.5f |1> (%2.5f)\n", a*a, b*b, (a*a)+(b*b));
#endif
     }
     else
     {
#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("qubit is null\n");
#else
         printf("qubit is null\n");
#endif
     }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qubit_vec(tf_qcs_qubit *q)
{
     if(q!=NULL)
     {
#ifdef PYTHON_SCRIPT
        PySys_WriteStdout("/                    \\\n");
        PySys_WriteStdout("| %2.6f + %2.6fi |\n", q->alpha.re, q->alpha.im);
        PySys_WriteStdout("|                    |\n");
        PySys_WriteStdout("| %2.6f + %2.6fi |\n", q->beta.re,  q->beta.im);
        PySys_WriteStdout("\\                    /\n");
#else
        printf("/                    \\\n");
   		printf("| %2.6f + %2.6fi |\n", q->alpha.re, q->alpha.im);
        printf("|                    |\n");
        printf("| %2.6f + %2.6fi |\n", q->beta.re,  q->beta.im);
        printf("\\                    /\n");
#endif
     }
     else
     {
#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("qubit is null\n");
#else
         printf("qubit is null\n");
#endif
     }
}

DYNAMIC_LIB_DECORATION void qcs_printf_qubit_vec_sqr(tf_qcs_qubit *q)
{
     tf_qcs_real_number a, b;

     if(q!=NULL)
     {
         qcs_mod_complex(&q->alpha,&a);
         qcs_mod_complex(&q->beta, &b);

#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("/          \\\n");
         PySys_WriteStdout("| %2.6f |\n",a);
         PySys_WriteStdout("|          |\n");
         PySys_WriteStdout("| %2.6f |\n", b);
         PySys_WriteStdout("\\           /\n");
#else
         printf("/          \\\n");
         printf("| %2.6f |\n",a);
         printf("|          |\n");
         printf("| %2.6f |\n", b);
         printf("\\           /\n");
#endif
     }
     else
     {
#ifdef PYTHON_SCRIPT
         PySys_WriteStdout("qubit is null\n");
#else
         printf("qubit is null\n");
#endif
     }
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit_array* qcs_new_qubit_array(int n)
{
    tf_qcs_qubit_array* t;
    int i;

    t = (tf_qcs_qubit_array*)malloc( sizeof(tf_qcs_qubit_array) );
    t->arr = (tf_qcs_qubit*)malloc( n*sizeof(tf_qcs_qubit) );
    t->n = n;

    for(i=0;i<n;i++)
    {
        (t->arr + i)->alpha.re=0;
        (t->arr + i)->alpha.im=0;
        (t->arr + i)->beta.re=0;
        (t->arr + i)->beta.im=0;
    }

    return t;
}

DYNAMIC_LIB_DECORATION void qcs_set_qubit_array_n(tf_qcs_qubit_array *qa, int i, tf_qcs_qubit *q)
{
//    if (qa==NULL) printf("NULL\n");
//    printf("qcs_set_qubit_array_n %d at %d:\n", qa->n, i);
//    qcs_printf_qubit( q );

    (qa->arr + i)->alpha.re = q->alpha.re;
    (qa->arr + i)->alpha.im = q->alpha.im;
    (qa->arr + i)->beta.re  = q->beta.re;
    (qa->arr + i)->beta.im  = q->beta.im;
}

DYNAMIC_LIB_DECORATION tf_qcs_qubit *qcs_get_qubit_array_n(tf_qcs_qubit_array *qa, int i)
{
            return (qa->arr + i);
}

DYNAMIC_LIB_DECORATION void qcs_delete_qubit_array(tf_qcs_qubit_array *qa)
{
     free(qa->arr);
     qa->n=0;
     free(qa);
     qa=NULL;
}
