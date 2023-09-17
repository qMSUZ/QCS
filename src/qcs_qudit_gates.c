/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Copyright (C) 2007 -- 2008 by Przemysław Ratajczak                   *
 *                                                                         *
 *   Part of the Quantum Computing Simulator:                              *
  *   https://github.com/qMSUZ/QCS                                         *
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
#include <string.h>

#include "qcs.h"
#include "qcs_matrix.h"
#include "qcs_gates.h"
#include "qcs_qudit.h"
#include "qcs_qudit_gates.h"

static tf_qcs_matrix *qudit_zero_reset_gate=NULL;
static tf_qcs_matrix *qudit_id_gate=NULL;
static tf_qcs_matrix *qudit_not_gate=NULL;
static tf_qcs_matrix *qudit_pauli_x_gate=NULL;
static tf_qcs_matrix *qudit_pauli_y_gate=NULL;
static tf_qcs_matrix *qudit_pauli_z_gate=NULL;
static tf_qcs_matrix *qudit_hadamard_gate=NULL;
static tf_qcs_matrix *qudit_hadamardv2_gate=NULL;
static tf_qcs_matrix *qudit_square_root_gate=NULL;
static tf_qcs_matrix *qudit_x_rot90_gate=NULL;
static tf_qcs_matrix *qudit_y_rot90_gate=NULL;
static tf_qcs_matrix *qudit_z_rot90_gate=NULL;
static tf_qcs_matrix *qudit_minus_x_rot90_gate=NULL;
static tf_qcs_matrix *qudit_minus_y_rot90_gate=NULL;
static tf_qcs_matrix *qudit_minus_z_rot90_gate=NULL;
static tf_qcs_matrix *qudit_rotate_alpha_gate=NULL;
static tf_qcs_matrix *qudit_rotate_theta_gate=NULL;
static tf_qcs_matrix *qudit_t_gate=NULL;
static tf_qcs_matrix *qudit_s_gate=NULL;
static tf_qcs_matrix *qudit_phase_gate=NULL;
static tf_qcs_matrix *qudit_phase_f_gate=NULL;
static tf_qcs_matrix *qudit_cnot_gate=NULL;
static tf_qcs_matrix *qudit_cnot_swaped_gate=NULL; //control & target qubits swaped
static tf_qcs_matrix *qudit_toffoli_gate=NULL;
static tf_qcs_matrix *qudit_swap_gate=NULL;

/* QUDIT GATES */

DYNAMIC_LIB_DECORATION char *create_qudit_gates_cache_key(enum _qcs_gate_type gate_type, int freedom)
{
    static char key[128];

    memset(&key[0], 0, 128);
    sprintf(&key[0], "%s_%d", get_gate_name(gate_type), freedom);
    return &key[0];
}

DYNAMIC_LIB_DECORATION pt_qcs_qudit_gate_cache_element make_qudit_gates_cache_element( enum _qcs_gate_type gate_type, int freedom, tf_qcs_matrix *m)
{
    pt_qcs_qudit_gate_cache_element e = NULL;

    e = (pt_qcs_qudit_gate_cache_element) malloc (sizeof (t_qcs_qudit_gate_cache_element));

    e->gate_type = gate_type;
    e->freedom = freedom;
    e->m = m;

    return e;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_zero_reset_gate(int freedom_level)
{
    return qudit_zero_reset_gate;
}

/*tf_qcs_matrix *get_qudit_id_gate(int freedom_level)
{
    return qudit_id_gate;
}*/


DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_id_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e = NULL;

    tf_qcs_matrix *local_qudit_id_gate;


    key = create_qudit_gates_cache_key( QCS_ID_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e =  qcs_hash_table_get ( qudit_gates_cache, key);
    if ( e != NULL )
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element exists\n" );
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element does not exist!\ngeneration a new matrix\n" );
#endif

        local_qudit_id_gate=qcs_create_matrix(freedom_level, freedom_level);

        for(i=0; i<freedom_level*freedom_level; i++)
        {
            if (i%(freedom_level+1)==0)  //not optimal for large  //%=fmod
            {
                (local_qudit_id_gate->m+i)->re=1; (local_qudit_id_gate->m+i)->im=0;
            }
            else
            {
                (local_qudit_id_gate->m+i)->re=0; (local_qudit_id_gate->m+i)->im=0;
            }
        }

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element( QCS_ID_GATE, freedom_level, local_qudit_id_gate );
        qcs_hash_table_insert ( qudit_gates_cache, key, (void*)e );
    }
    return local_qudit_id_gate;
}


DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_not_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e =  NULL;

    tf_qcs_matrix *local_qudit_not_gate;


    key = create_qudit_gates_cache_key( QCS_NOT_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif

    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if ( e != NULL )
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element exists\n" );
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element does not exist!\ngeneration a new matrix\n" );
#endif

        local_qudit_not_gate=qcs_create_matrix(freedom_level,freedom_level);

        for(i=0; i<freedom_level*freedom_level; i++)
        {
            if (i%(freedom_level+1)==freedom_level)  //not optimal for large
            {
                (local_qudit_not_gate->m+i)->re=1; (local_qudit_not_gate->m+i)->im=0;
            }
            else
            {
                (local_qudit_not_gate->m+i)->re=0; (local_qudit_not_gate->m+i)->im=0;
            }
        }
        (local_qudit_not_gate->m+freedom_level-1)->re=1; (local_qudit_not_gate->m+freedom_level-1)->im=0;

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element( QCS_NOT_GATE, freedom_level, local_qudit_not_gate );
        qcs_hash_table_insert ( qudit_gates_cache, key, (void*)e );
    }

    return local_qudit_not_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_x_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e =  NULL;

    tf_qcs_matrix *local_qudit_pauli_x_gate;


    key = create_qudit_gates_cache_key( QCS_PAULIX_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif

    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if ( e != NULL )
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element exists\n" );
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element does not exist!\ngeneration a new matrix\n" );
#endif
        local_qudit_pauli_x_gate=qcs_create_matrix(freedom_level,freedom_level);

        for(i=0; i<freedom_level*freedom_level; i++)
        {
            if (i%(freedom_level+1)==freedom_level)  //not optimal for large
            {
                (local_qudit_pauli_x_gate->m+i)->re=1; (local_qudit_pauli_x_gate->m+i)->im=0;
            }
            else
            {
                (local_qudit_pauli_x_gate->m+i)->re=0; (local_qudit_pauli_x_gate->m+i)->im=0;
            }
        }
        (local_qudit_pauli_x_gate->m+freedom_level-1)->re=1; (local_qudit_pauli_x_gate->m+freedom_level-1)->im=0;

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element( QCS_PAULIX_GATE, freedom_level, local_qudit_pauli_x_gate );
        qcs_hash_table_insert ( qudit_gates_cache, key, (void*)e );
    }

    return local_qudit_pauli_x_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_y_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e =  NULL;

    tf_qcs_matrix *local_qudit_pauli_y_gate;


    key = create_qudit_gates_cache_key( QCS_PAULIY_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if ( e != NULL )
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element exists\n" );
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element does not exist!\ngeneration a new matrix\n" );
#endif
        local_qudit_pauli_y_gate=qcs_create_matrix(freedom_level, freedom_level);

        temp=-1;

        for(i=freedom_level*freedom_level-2; i>0; i--)
        {
            if (i%(freedom_level-1)==0)  //not optimal for large
            {
                temp++;
                (local_qudit_pauli_y_gate->m+i)->re=sin(2*QCS_PI*temp/freedom_level); (local_qudit_pauli_y_gate->m+i)->im=cos(2*QCS_PI*temp/freedom_level);//cos(2*QCS_PI*temp/freedom_level); (qudit_pauli_z_gate->m+i)->im=sin(2*QCS_PI*temp/freedom_level);
            }
            else
            {
                (local_qudit_pauli_y_gate->m+i)->re=0; (local_qudit_pauli_y_gate->m+i)->im=0;
            }
        }
        (local_qudit_pauli_y_gate->m+0)->re=0; (local_qudit_pauli_y_gate->m+0)->im=0;
        (local_qudit_pauli_y_gate->m+freedom_level*freedom_level-1)->re=0; (local_qudit_pauli_y_gate->m+freedom_level*freedom_level-1)->im=0;
        (local_qudit_pauli_y_gate->m+freedom_level*(freedom_level-1))->re=0; (local_qudit_pauli_y_gate->m+freedom_level*(freedom_level-1))->im=1;


        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element( QCS_PAULIY_GATE, freedom_level, local_qudit_pauli_y_gate );
        qcs_hash_table_insert ( qudit_gates_cache, key, (void*)e );
    }

    return local_qudit_pauli_y_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_z_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e =  NULL;

    tf_qcs_matrix *local_qudit_pauli_z_gate;


    key = create_qudit_gates_cache_key( QCS_PAULIZ_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if ( e != NULL )
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element exists\n" );
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf( "element does not exist!\ngeneration a new matrix\n" );
#endif
        local_qudit_pauli_z_gate=qcs_create_matrix(freedom_level,freedom_level);

        temp=0;
        (local_qudit_pauli_z_gate->m+0)->re=1;
        (local_qudit_pauli_z_gate->m+0)->im=0;
        for(i=1; i<freedom_level*freedom_level; i++)
        {
            if (i%(freedom_level+1)==0)  //not optimal for large
            {
                temp++;
                (local_qudit_pauli_z_gate->m+i)->re=cos(2*QCS_PI*temp/freedom_level);
                (local_qudit_pauli_z_gate->m+i)->im=sin(2*QCS_PI*temp/freedom_level);
            }
            else
            {
                (local_qudit_pauli_z_gate->m+i)->re=0;
                (local_qudit_pauli_z_gate->m+i)->im=0;
            }
        }

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element( QCS_PAULIZ_GATE, freedom_level, local_qudit_pauli_z_gate );
        qcs_hash_table_insert ( qudit_gates_cache, key, (void*)e );
    }

    return local_qudit_pauli_z_gate;
}

/*
tf_qcs_matrix *get_qudit_hadamard_gate(int freedom_level)
{
    return qudit_hadamard_gate;
}
*/

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_hadamard_gate(int freedom_level)
{
    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e = NULL;

    tf_qcs_matrix *local_qudit_hadamard_gate;


    key = create_qudit_gates_cache_key( QCS_HADAMARD_GATE, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if( e != NULL)
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element exists\n");
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element does not exist!\ngeneration a new matrix\n");
#endif
        local_qudit_hadamard_gate=qcs_create_matrix(freedom_level,freedom_level);

        multiplication_factor = 1/sqrt(freedom_level);
        for(i=0; i<freedom_level*freedom_level; i++)
        {
            if (i < freedom_level)
            {
                (local_qudit_hadamard_gate->m+i)->re = multiplication_factor; (local_qudit_hadamard_gate->m+i)->im=0;
            }
            else
            {
                if (i % freedom_level==0)
                {
                    (local_qudit_hadamard_gate->m+i)->re = multiplication_factor; (local_qudit_hadamard_gate->m+i)->im=0;
                }
                else
                {
                    (local_qudit_hadamard_gate->m+i)->re = cos(2*QCS_PI*(i%freedom_level)*(freedom_level-i/freedom_level)/freedom_level)*multiplication_factor;
                    (local_qudit_hadamard_gate->m+i)->im = sin(2*QCS_PI*(i%freedom_level)*(freedom_level-i/freedom_level)/freedom_level)*multiplication_factor;
                }
            }
        }
        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element(QCS_HADAMARD_GATE, freedom_level, local_qudit_hadamard_gate);
        qcs_hash_table_insert (qudit_gates_cache, key, (void*)e);
    }

    return local_qudit_hadamard_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_square_root_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_x_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_y_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_z_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_minus_x_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_minus_y_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_minus_z_rot90_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_rotate_alpha_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_rotate_theta_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_t_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_s_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_phase_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_phase_f_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_cnot_gate(int freedom_level)
{
    return qudit_cnot_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_cnot_swaped_gate(int freedom_level)
{
    return qudit_cnot_swaped_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_toffoli_gate(int freedom_level)
{
    return NULL;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_swap_gate(int freedom_level)
{
    return qudit_swap_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_hadamard_d_gate(int freedom_level)
{
    int temp=0;
    int k=0;
    int l=0;
    tf_qcs_real_number multiplication_factor=0;
    char *key;
    void *e = NULL;

    tf_qcs_matrix *local_qudit_hadamard_gate;
    tf_qcs_complex v1,v2,v3;

    key = create_qudit_gates_cache_key( QCS_HADAMARD_GATE_d, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if( e != NULL)
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element exists\n");
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element does not exist!\ngeneration a new matrix\n");
#endif
        local_qudit_hadamard_gate=qcs_create_matrix(freedom_level, freedom_level);
        multiplication_factor = 1/sqrt(freedom_level);

        for(k=1;k<=freedom_level;k++)
        {
            for(l=1;l<=freedom_level;l++)
            {
                v1.re = 0;
                v1.im = 2*(k-1)*(l-1) * QCS_PI;

                v2.re = 1.0/(tf_qcs_real_number)freedom_level;
                v2.im = 0;

                qcs_complex_mul(&v1,&v2,&v3);
                qcs_exp_complex(&v3, &v2);

                v1.re=multiplication_factor;
                v1.im=0;

                qcs_complex_mul(&v1, &v2, &v3);
                qcs_set_cell_at_matrix_complex(local_qudit_hadamard_gate, k-1, l-1, &v3);
            }
        }

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element(QCS_HADAMARD_GATE_d, freedom_level, local_qudit_hadamard_gate);
        qcs_hash_table_insert (qudit_gates_cache, key, (void*)e);
    }

    return local_qudit_hadamard_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_x_d_gate(int freedom_level)
{
    tf_qcs_matrix *local_qudit_pauli_x_d_gate;
    int new_idx, k=0;
    int l=0;
    tf_qcs_complex v1;
    char *key;
    void *e = NULL;

    v1.re=1;
    v1.im=0;

    key = create_qudit_gates_cache_key( QCS_PAULIX_GATE_d, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if( e != NULL)
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element exists\n");
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element does not exist!\ngeneration a new matrix\n");
#endif
        local_qudit_pauli_x_d_gate=qcs_create_matrix(freedom_level, freedom_level);
        for(k=1;k<=freedom_level;k++)
        {
            new_idx=k-1-1;
            if(new_idx < 0) new_idx = freedom_level - 1;
            qcs_set_cell_at_matrix_complex(local_qudit_pauli_x_d_gate, (k-1) , new_idx, &v1);
        }

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element(QCS_PAULIX_GATE_d, freedom_level, local_qudit_pauli_x_d_gate);
        qcs_hash_table_insert (qudit_gates_cache, key, (void*)e);
    }

    return local_qudit_pauli_x_d_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_z_d_gate(int freedom_level)
{
    tf_qcs_matrix *local_qudit_pauli_z_d_gate;
    int k=0;
    int l=0;
    tf_qcs_complex v1, v2, v3;
    char *key;
    void *e = NULL;
    tf_qcs_real_number multiplication_factor=0;

    v1.re=1;
    v1.im=0;

    key = create_qudit_gates_cache_key( QCS_PAULIZ_GATE_d, freedom_level);
#ifdef __qcs_core_library_debug_mode__
    printf("key = %s\n", key);
#endif
    e = qcs_hash_table_get ( qudit_gates_cache, key);
    if( e != NULL)
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element exists\n");
#endif
        return ((pt_qcs_qudit_gate_cache_element)e)->m;
    }
    else
    {
#ifdef __qcs_core_library_debug_mode__
        printf("element does not exist!\ngeneration a new matrix\n");
#endif
        local_qudit_pauli_z_d_gate=qcs_create_matrix(freedom_level, freedom_level);
        //multiplication_factor = 1.0/(tf_qcs_real_number)(freedom_level);
        for(k=1;k<=freedom_level;k++)
        {
                v1.re = 0;
                v1.im = 2*(k-1) * QCS_PI;

                v2.re = 1.0/(tf_qcs_real_number)freedom_level;
                v2.im = 0;

                qcs_complex_mul(&v1, &v2, &v3);
                qcs_exp_complex(&v3, &v2);

                //v1.re=multiplication_factor;
                //v1.im=0;

                //qcs_complex_mul(&v1, &v2, &v3);

            qcs_set_cell_at_matrix_complex(local_qudit_pauli_z_d_gate, (k-1) , (k-1), &v2);
        }

        pt_qcs_qudit_gate_cache_element e = make_qudit_gates_cache_element(QCS_PAULIZ_GATE_d, freedom_level, local_qudit_pauli_z_d_gate);
        qcs_hash_table_insert (qudit_gates_cache, key, (void*)e);
    }

    return local_qudit_pauli_z_d_gate;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_power_of_x_gate(int freedom_level, int power)
{
    if(power == 0)
    {
         tf_qcs_matrix *out_math = qcs_clone_matrix(get_qudit_pauli_x_gate(freedom_level));

         qcs_eye_matrix( out_math );

         return out_math;
    }

    int i;

    tf_qcs_matrix *x_gate = get_qudit_pauli_x_gate(freedom_level);

    tf_qcs_matrix *tmp = qcs_clone_matrix(x_gate);
    tf_qcs_matrix *out_math = qcs_create_matrix(freedom_level, freedom_level);

    qcs_copy_matrix(x_gate, out_math);


    for(i=1;i<power;i++)
    {
        qcs_zero_matrix(out_math);

        qcs_mul_matrix(tmp, x_gate, out_math);

        qcs_copy_matrix(out_math, tmp);
    }

    qcs_delete_matrix(tmp);

    return out_math;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_power_of_z_gate(int freedom_level, int power)
{
    if(power == 0)
    {
         tf_qcs_matrix *out_math = qcs_clone_matrix(get_qudit_pauli_z_d_gate(freedom_level));

         qcs_eye_matrix( out_math );

         return out_math;
    }

    int i;

    tf_qcs_matrix *z_gate = get_qudit_pauli_z_d_gate( freedom_level );

    tf_qcs_matrix *tmp=qcs_clone_matrix( z_gate );
    tf_qcs_matrix *out_math=qcs_create_matrix(freedom_level, freedom_level);


    qcs_copy_matrix(z_gate, out_math);


    for(i=1;i<power;i++)
    {
        qcs_zero_matrix(out_math);

        qcs_mul_matrix(tmp, z_gate, out_math);

        qcs_copy_matrix(out_math, tmp);
    }

    qcs_delete_matrix(tmp);

    return out_math;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_power_of_u_gate(tf_qcs_matrix *u, int power)
{

    if(power == 0)
    {
         tf_qcs_matrix *out_math=qcs_clone_matrix(u);
         qcs_eye_matrix( out_math);

         return out_math;
    }

    int i;
    tf_qcs_matrix *u_gate = u;

    tf_qcs_matrix *tmp=qcs_clone_matrix(u);
    tf_qcs_matrix *out_math=qcs_clone_matrix(u);

    qcs_zero_matrix(out_math);
    qcs_copy_matrix(u_gate, out_math);


    for(i=1;i<power;i++)
    {
        qcs_zero_matrix(out_math);

        qcs_mul_matrix(tmp, u_gate, out_math);

        qcs_copy_matrix(out_math, tmp);
    }

    qcs_delete_matrix(tmp);

    return out_math;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_qudit_pauli_operator_gate(int freedom_level, int r, int s)
{
    tf_qcs_matrix *x, *z, *out_mat;

    x = get_qudit_power_of_x_gate( freedom_level, r);
    z = get_qudit_power_of_z_gate( freedom_level, s);

    out_mat = qcs_create_matrix( freedom_level, freedom_level);

    qcs_mul_matrix(x, z, out_mat);

    qcs_delete_matrix( z );
    qcs_delete_matrix( x );

    return out_mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *make_arbitrary_matrix_for_one_qudit(tf_qcs_matrix *gate, int n, int fd, int t)
{
    tf_qcs_matrix *u=NULL;

    int dx, dy, x, y, i, j, minimatrix=0, step=0;
    char dits[128];

    t++; // poprawka dotyczaca numeru quditu

    u=qcs_create_matrix(pow(fd,n), pow(fd,n));

    minimatrix=pow(fd,n-1);
    step=pow(fd,n-t)-1;

    printf("minimatrices=%d, s=%d\n", minimatrix, step);

    j=1;
    dx=dy=0;
    for (i=0;i<pow(fd,n);i++)
    {
        qcs_dec2base_d(i, n, fd, &dits[0]);
        dits[n]=0;

        if (dits[t-1]=='0' && j<=minimatrix)
        {
            for(x=0;x<fd;x++)
            {
                if(x==0)
                    dx=0;
                else
                    dx=x+step*x;

                for(y=0;y<fd;y++)
                {
                    if(y==0)
                        dy=0;
                    else
                        dy=y+step*y;

                    //printf("%d %d from %d %d\n",i+dx,i+dy,x,y);
                    qcs_set_cell_at_matrix_complex(u, i+dx, i+dy, qcs_get_cell_at_matrix_complex(gate, x, y));
                }
            }
            printf("\n");
            //qcs_set_cell_at_matrix_complex(u, i, i,        (gate->m+0)); qcs_set_cell_at_matrix_complex(u, i, i+1+step  (gate->m+1));
            //qcs_set_cell_at_matrix_complex(u, i+1+step, i, (gate->m+2)); qcs_set_cell_at_matrix_complex(u, i+1+step, i+1+step, (gate->m+3));
            j++;
        }
    }

    return u;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_general_e_matrix(int fd, int k, int j)
{
    int v, u, value;
    tf_qcs_matrix *e_mat;

    e_mat = qcs_create_matrix(fd, fd);

    for(v=1 ; v<=fd ; v++)
    {
        for(u=1 ; u<=fd ; u++)
        {
            value = QKDF(v,j) * QKDF(u,k); // QKDF -- kronecker delta function
            qcs_set_cell_at_matrix_direct(e_mat, v-1, u-1, value, 0);
        }
    }

    return e_mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_general_theta_matrix(int fd, int k, int j)
{
    tf_qcs_matrix *e_mat, *e1, *e2;

    e_mat = qcs_create_matrix(fd, fd);

    e1 = get_general_e_matrix( fd, k, j);
    e2 = get_general_e_matrix( fd, j, k);

    qcs_add_matrix(e1, e2, e_mat);

    return e_mat;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_general_beta_matrix(int fd, int k, int j) /* 1 <= k < j <= d */
{
    tf_qcs_matrix *e_mat, *e1, *e2, *e_out;
    tf_qcs_complex minus_i;

    minus_i.re =  0;
    minus_i.im = -1;

    e_mat = qcs_create_matrix(fd, fd);
    e_out = qcs_create_matrix(fd, fd);

    e1 = get_general_e_matrix( fd, j, k);
    e2 = get_general_e_matrix( fd, k, j);

    qcs_sub_matrix(e1, e2, e_mat);
    qcs_mul_scalar_matrix(e_mat, &minus_i, e_out);

    qcs_delete_matrix( e_mat );

    return e_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *get_general_eta_matrix(int fd, int r)
{
    int j;
    tf_qcs_real_number sqrt_r;
    tf_qcs_matrix *e1, *e2, *e_tmp1, *e_tmp2, *e_tmp3, *e_out;
    tf_qcs_complex rcmx;

    e_tmp1 = qcs_create_matrix( fd, fd);
    e_tmp2 = qcs_create_matrix( fd, fd);
    e_tmp3 = qcs_create_matrix( fd, fd);

    e_out = qcs_create_matrix( fd, fd);

    rcmx.re=r; rcmx.im=0;
    for(j=1;j<=r;j++)
    {
        e1 = get_general_e_matrix(fd, j, j);

        qcs_add_matrix(e1, e_tmp1, e_tmp2);
        qcs_copy_matrix(e_tmp2, e_tmp1);
    }

    e2 = get_general_e_matrix(fd, r+1, r+1);
    qcs_mul_scalar_matrix(e2, &rcmx, e_tmp2);

    qcs_sub_matrix(e_tmp1, e_tmp2, e_tmp3);

    rcmx.re = sqrtf(2.0 / (r*(r+1)) );
    rcmx.im = 0;

    qcs_mul_scalar_matrix(e_tmp3, &rcmx, e_out);

    qcs_delete_matrix( e_tmp1 );
    qcs_delete_matrix( e_tmp2 );
    qcs_delete_matrix( e_tmp3 );

    return e_out;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *qudit_d_level_teleportation_correction_matrix(int n, int m, int d)
{
    int k, kk, i;
    tf_qcs_matrix *tmp_matrix = NULL;
    tf_qcs_complex v1, v2;

    tmp_matrix = qcs_create_matrix(d, d);

    i=m;
    for(k=0;k<d;k++)
    {
        v1.re=0;
        v1.im=(2.0 * QCS_PI * n * (tf_qcs_real_number)k)/(tf_qcs_real_number)d;
        qcs_exp_complex(&v1, &v2);

        //kk=((k+m) % d);
        kk=qcs_mod(i,d);
        //printf("|%d><%d| ", k, kk);
        qcs_set_cell_at_matrix_complex(tmp_matrix, k, kk, &v2);
        i--;
    }

    return tmp_matrix;
}

/******/

DYNAMIC_LIB_DECORATION void qcs_build_qudit_gates_matrix(freedom_level)
{
    qudit_zero_reset_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_id_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_not_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_pauli_x_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_pauli_y_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_pauli_z_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_hadamard_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_square_root_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_x_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_y_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_z_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_minus_x_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_minus_y_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_minus_z_rot90_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_rotate_alpha_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_rotate_theta_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_t_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_s_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_phase_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_phase_f_gate=qcs_create_matrix(freedom_level,freedom_level);
    qudit_cnot_gate=qcs_create_matrix(freedom_level*freedom_level,freedom_level*freedom_level);
    qudit_cnot_swaped_gate=qcs_create_matrix(freedom_level*freedom_level,freedom_level*freedom_level); //control & target qubits swaped
    qudit_toffoli_gate=qcs_create_matrix(freedom_level*freedom_level*freedom_level,freedom_level*freedom_level*freedom_level);
    qudit_swap_gate=qcs_create_matrix(freedom_level*freedom_level,freedom_level*freedom_level);

    int temp=0;
    int i=0;
    int j=0;
    tf_qcs_real_number multiplication_factor=0;

/* Qudit Zero matrix */
    for(i=0; i<freedom_level*freedom_level; i++)
    {
        (qudit_zero_reset_gate->m+i)->re=0; (qudit_zero_reset_gate->m+i)->im=0;
    }

/* Qudit Id matrix */

    for(i=0; i<freedom_level*freedom_level; i++)
    {
        if (i%(freedom_level+1)==0)  //not optimal for large  //%=fmod
        {
            (qudit_id_gate->m+i)->re=1; (qudit_id_gate->m+i)->im=0;
        }
        else
        {
            (qudit_id_gate->m+i)->re=0; (qudit_id_gate->m+i)->im=0;
        }
    }

/* Qudit NOT matrix */

    for(i=0; i<freedom_level*freedom_level; i++)
    {
        if (i%(freedom_level+1)==freedom_level)  //not optimal for large
        {
            (qudit_not_gate->m+i)->re=1; (qudit_not_gate->m+i)->im=0;
        }
        else
        {
            (qudit_not_gate->m+i)->re=0; (qudit_not_gate->m+i)->im=0;
        }
    }
    (qudit_not_gate->m+freedom_level-1)->re=1; (qudit_not_gate->m+freedom_level-1)->im=0;

/* Qudit Pauli X matrix */

    for(i=0; i<freedom_level*freedom_level; i++)
    {
        if (i%(freedom_level+1)==freedom_level)  //not optimal for large
        {
            (qudit_pauli_x_gate->m+i)->re=1; (qudit_pauli_x_gate->m+i)->im=0;
        }
        else
        {
            (qudit_pauli_x_gate->m+i)->re=0; (qudit_pauli_x_gate->m+i)->im=0;
        }
    }
    (qudit_pauli_x_gate->m+freedom_level-1)->re=1; (qudit_pauli_x_gate->m+freedom_level-1)->im=0;

/* Qudit Pauli Y matrix */

/* Qudit Pauli Z matrix */

    temp=0;
    (qudit_pauli_z_gate->m+0)->re=1; (qudit_pauli_z_gate->m+0)->im=0;
    for(i=1; i<freedom_level*freedom_level; i++)
    {
        if (i%(freedom_level+1)==0)  //not optimal for large
        {
            temp++;
            (qudit_pauli_z_gate->m+i)->re=cos(2*QCS_PI*temp/freedom_level); (qudit_pauli_z_gate->m+i)->im=sin(2*QCS_PI*temp/freedom_level);
        }
        else
        {
            (qudit_pauli_z_gate->m+i)->re=0; (qudit_pauli_z_gate->m+i)->im=0;
        }
    }

/* Qudit Hadamard matrix */

    multiplication_factor = 1/sqrt(freedom_level);
    for(i=0; i<freedom_level*freedom_level; i++)
    {
        if (i<freedom_level)
        {
            (qudit_hadamard_gate->m+i)->re=multiplication_factor; (qudit_hadamard_gate->m+i)->im=0;
        }
        else
        {
            if (i%freedom_level==0)
            {
                (qudit_hadamard_gate->m+i)->re=multiplication_factor; (qudit_hadamard_gate->m+i)->im=0;
            }
            else
            {
                (qudit_hadamard_gate->m+i)->re=cos(2*QCS_PI*(i%freedom_level)*(freedom_level-i/freedom_level)/freedom_level)*multiplication_factor;
                (qudit_hadamard_gate->m+i)->im=sin(2*QCS_PI*(i%freedom_level)*(freedom_level-i/freedom_level)/freedom_level)*multiplication_factor;
            }
        }
    }
    //elimination of numeric errors (when value is close to 0 (~10^(-16)), but it should be equal to 0)
    for(i=0; i<freedom_level*freedom_level; i++)
    {
        if (((qudit_hadamard_gate->m+i)->re > -0.000000000000001) & ((qudit_hadamard_gate->m+i)->re < 0.000000000000001)) (qudit_hadamard_gate->m+i)->re = 0;
        if (((qudit_hadamard_gate->m+i)->im > -0.000000000000001) & ((qudit_hadamard_gate->m+i)->im < 0.000000000000001)) (qudit_hadamard_gate->m+i)->im = 0;
    }



/* Qudit Square Root matrix */

/* Qudit X Rot90 matrix */

/* Qudit Y Rot90 matrix */

/* Qudit Z Rot90 matrix */

/* Qudit Minus X Rot90 matrix */

/* Qudit Minus Y Rot90 matrix */

/* Qudit Minus Z Rot90 matrix */

/* Qudit Rotate Alpha matrix */

/* Qudit Rotate Theta matrix */

/* Qudit T gate matrix */

/* Qudit S gate matrix */

/* Qudit Phase gate matrix */

/* Qudit Phase F gate matrix */

/* Qudit CNOT gate matrix */

    for(i=0; i<freedom_level*freedom_level*freedom_level*freedom_level; i++)
    {
        (qudit_cnot_gate->m+i)->re=0; (qudit_cnot_gate->m+i)->im=0;
    }
    for(i=0; i<freedom_level; i++)
    {
        for(j=0; j<freedom_level; j++)
        {
            //x=i*freedom_level+j y=i*freedom_level+(i+j)%freedom_level
            (qudit_cnot_gate->m+(i*freedom_level+(i+j)%freedom_level)*freedom_level*freedom_level+(i*freedom_level+j))->re=1;
            (qudit_cnot_gate->m+(i*freedom_level+(i+j)%freedom_level)*freedom_level*freedom_level+(i*freedom_level+j))->im=0;
        }
    }

/* Qudit CNOT SWAPED gate matrix */ //control & target qubits swaped

    for(i=0; i<freedom_level*freedom_level*freedom_level*freedom_level; i++)
    {
        (qudit_cnot_swaped_gate->m+i)->re=0; (qudit_cnot_swaped_gate->m+i)->im=0;
    }
    for(i=0; i<freedom_level*freedom_level; i++)
    {
       if (i%freedom_level==0)
       {
            (qudit_cnot_swaped_gate->m+i*freedom_level*freedom_level+i)->re=1;
            (qudit_cnot_swaped_gate->m+i*freedom_level*freedom_level+i)->im=0;

       }
       else
       {
            (qudit_cnot_swaped_gate->m+((freedom_level+1)*i)%(freedom_level*freedom_level)*freedom_level*freedom_level+i)->re=1;
            (qudit_cnot_swaped_gate->m+((freedom_level+1)*i)%(freedom_level*freedom_level)*freedom_level*freedom_level+i)->im=0;
       }
    }

/* Qudit Toffoli gate matrix */

/* Qudit SWAP gate matrix */
    int ind=0;
    int ind2=0;
    for(i=0; i<freedom_level*freedom_level*freedom_level*freedom_level; i++)
    {
        (qudit_swap_gate->m+i)->re=0; (qudit_swap_gate->m+i)->im=0;
    }
    for(i=0; i<freedom_level*freedom_level; i++)
    {
       if (i%freedom_level==0)
       {
            (qudit_swap_gate->m+ind*freedom_level*freedom_level+ind*freedom_level)->re=1;
            (qudit_swap_gate->m+ind*freedom_level*freedom_level+ind*freedom_level)->im=0;
            ind++;
            ind2=0;
       }
       else
       {
            ind2++;
            (qudit_swap_gate->m+(ind+freedom_level*ind2-1)*freedom_level*freedom_level+i)->re=1;
            (qudit_swap_gate->m+(ind+freedom_level*ind2-1)*freedom_level*freedom_level+i)->im=0;
       }
    }
}

/******/

DYNAMIC_LIB_DECORATION void qcs_destroy_qudit_gates_matrix()
{
     qcs_delete_matrix(qudit_zero_reset_gate);
     qcs_delete_matrix(qudit_id_gate);
     qcs_delete_matrix(qudit_not_gate);
     qcs_delete_matrix(qudit_pauli_x_gate);
     qcs_delete_matrix(qudit_pauli_y_gate);
     qcs_delete_matrix(qudit_pauli_z_gate);
     qcs_delete_matrix(qudit_hadamard_gate);
     qcs_delete_matrix(qudit_square_root_gate);
     qcs_delete_matrix(qudit_x_rot90_gate);
     qcs_delete_matrix(qudit_y_rot90_gate);
     qcs_delete_matrix(qudit_z_rot90_gate);
     qcs_delete_matrix(qudit_minus_x_rot90_gate);
     qcs_delete_matrix(qudit_minus_y_rot90_gate);
     qcs_delete_matrix(qudit_minus_z_rot90_gate);
     qcs_delete_matrix(qudit_rotate_alpha_gate);
     qcs_delete_matrix(qudit_rotate_theta_gate);
     qcs_delete_matrix(qudit_t_gate);
     qcs_delete_matrix(qudit_s_gate);
     qcs_delete_matrix(qudit_phase_gate);
     qcs_delete_matrix(qudit_phase_f_gate);
     qcs_delete_matrix(qudit_cnot_gate);
     qcs_delete_matrix(qudit_cnot_swaped_gate);
     qcs_delete_matrix(qudit_toffoli_gate);
}

///!!!!! QUDIT PART

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qudit_syntesis_one_control_one_target(int n, int d, int c1, int t)
{
    // n - number of qudits in circuit (lines in circuit)
    // d - dimension of qudits (freedom_level)
    // c - number of control qudit [1..n]
    // t - number of target qudit [1..n]

    c1=c1+1; //in [0..n-1]
    t=t+1;


     tf_qcs_matrix *cnot = NULL;
     cnot = qcs_create_matrix(pow(d, n), pow(d, n));

     long i = 0;
     int tabb[n]; //for representation in d-number system
     long ind = 0;
     long liczba = 0; //max (qudit dimension)^(number of qudits) = pow(d, n)
     double ww = 0;
     double liczba2 = 0;
     int a=0;
     int control_bit_value = 0;
     double result = 0;

     //zeros in matrix
     for(i = 0; i < pow(d, n) * pow(d, n); i++)
     {
        (cnot->m+i)->re=0; (cnot->m+i)->im=0;
     }

     for(i = 0; i < pow(d, n); i++) //rows
     {
        //finding d-number representation (matlab == dec2base(i,d))
        liczba = i;
        liczba2 = liczba;
        for (ind=0; ind < n; ind++)
        {
            if (pow(d, n-ind-1) > liczba2)
            {
                tabb[ind]=0;
            }
            else
            {
                ww = 1;
                for (a=0; a < n-ind-1; a++) //=pow()
                    ww=ww*d;
                tabb[ind] = floor(liczba / ww);
                liczba = liczba % (long)ww;
            }
        }
        tabb[ind] = liczba;

        //finding new value of target qudit
        control_bit_value = tabb[c1 - 1];
        tabb[t - 1] = (tabb[t - 1] + control_bit_value) % d;

        //finding decimal representation of d-number value (matlab == base2dec(' ',d))
        result = 0;
        for (ind = n-1; ind >= 0; ind--)
            result = result + pow(d, n-ind-1) * tabb[ind];

        //changing one zero in cnot matrix (result row, i column)
        ind = result*pow(d, n);
        (cnot->m+i+ind)->re = 1;
     }
     return cnot;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qudit_syntesis_two_control_one_target(int n, int d, int c1, int c2, int t)
{
    // n - number of qudits in circuit (lines in circuit)
    // d - dimension of qudits (freedom_level)
    // c - number of control qudit [1..n]
    // t - number of target qudit [1..n]

     tf_qcs_matrix *cnot = NULL;
     cnot = qcs_create_matrix(pow(d, n), pow(d, n));

     long i = 0;
     int tabb[n]; //for representation in d-number system
     long ind = 0;
     long liczba = 0; //max (qudit dimension)^(number of qudits) = pow(d, n)
     double ww = 0;
     double liczba2 = 0;
     int a=0;
     int control_bit_value = 0;
     double result = 0;

     //zeros in matrix
     for(i = 0; i < pow(d, n) * pow(d, n); i++)
     {
        (cnot->m+i)->re=0; (cnot->m+i)->im=0;
     }

     for(i = 0; i < pow(d, n); i++) //rows
     {
        //finding d-number representation (matlab == dec2base(i,d))
        liczba = i;
        liczba2 = liczba;
        for (ind=0; ind < n; ind++)
        {
            if (pow(d, n-ind-1) > liczba2)
            {
                tabb[ind]=0;
            }else
            {
                ww = 1;
                for (a=0; a < n-ind-1; a++) //=pow()
                    ww=ww*d;
                tabb[ind] = floor(liczba / ww);
                liczba = liczba % (long)ww;
            }
        }
        tabb[ind] = liczba;

        //finding new value of target qudit
        control_bit_value = tabb[c1 - 1] + tabb[c2 - 1];
        tabb[t - 1] = (tabb[t - 1] + control_bit_value) % d;

        //finding decimal representation of d-number value (matlab == base2dec(' ',d))
        result = 0;
        for (ind = n-1; ind >= 0; ind--)
            result = result + pow(d, n-ind-1) * tabb[ind];

        //changing one zero in cnot matrix (result row, i column)
        ind = result*pow(d, n);
        (cnot->m+i+ind)->re = 1;
     }
     return cnot;
}

DYNAMIC_LIB_DECORATION tf_qcs_matrix *cnot_qudit_syntesis_three_control_one_target(int n, int d, int c1, int c2, int c3, int t)
{
    // n - number of qudits in circuit (lines in circuit)
    // d - dimension of qudits (freedom_level)
    // c - number of control qudit [1..n]
    // t - number of target qudit [1..n]

     tf_qcs_matrix *cnot = NULL;
     cnot = qcs_create_matrix(pow(d, n), pow(d, n));

     long i = 0;
     int tabb[n]; //for representation in d-number system
     long ind = 0;
     long liczba = 0; //max (qudit dimension)^(number of qudits) = pow(d, n)
     double ww = 0;
     double liczba2 = 0;
     int a=0;
     int control_bit_value = 0;
     double result = 0;

     //zeros in matrix
     for(i = 0; i < pow(d, n) * pow(d, n); i++)
     {
        (cnot->m+i)->re=0; (cnot->m+i)->im=0;
     }

     for(i = 0; i < pow(d, n); i++) //rows
     {
        //finding d-number representation (matlab == dec2base(i,d))
        liczba = i;
        liczba2 = liczba;
        for (ind=0; ind < n; ind++)
        {
            if (pow(d, n-ind-1) > liczba2)
            {
                tabb[ind]=0;
            }else
            {
                ww = 1;
                for (a=0; a < n-ind-1; a++) //=pow()
                    ww=ww*d;
                tabb[ind] = floor(liczba / ww);
                liczba = liczba % (long)ww;
            }
        }
        tabb[ind] = liczba;

        //finding new value of target qudit
        control_bit_value = tabb[c1 - 1] + tabb[c2 - 1] + tabb[c3 - 1];
        tabb[t - 1] = (tabb[t - 1] + control_bit_value) % d;

        //finding decimal representation of d-number value (matlab == base2dec(' ',d))
        result = 0;
        for (ind = n-1; ind >= 0; ind--)
            result = result + pow(d, n-ind-1) * tabb[ind];

        //changing one zero in cnot matrix (result row, i column)
        ind = result*pow(d, n);
        (cnot->m+i+ind)->re = 1;
     }
     return cnot;
}
