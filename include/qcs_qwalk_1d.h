/***************************************************************************
 *   Copyright (C) 2005 -- 2011 by Marek Sawerwain                         *
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
 
 #ifndef __qcs_qwalk_1d_h__
#define __qcs_qwalk_1d_h__

#include "qcs_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    QWalk1D_Line_Lattice = 1000,
    QWalk1D_Segment_Lattice = 1001,
    QWalk1D_Cycle_Lattice = 1002
};

typedef struct {
    int mt_version;
    int type;
    int steps, total_size, extra_size, max_index;
    int coin_size;
    int left, right;
    int iteration, experiments;
    tf_qcs_real_number rbl_prob;
    tf_qcs_matrix *coin;
    ti_qcs_matrix *RBL; // pointer for random broken links matrix
    tf_qcs_matrix *org_A;
    tf_qcs_matrix *copy_A;
    tf_qcs_matrix *probability;
} QWalk1D;

typedef QWalk1D tf_qwalk_1d;
typedef tf_qwalk_1d* pf_qwalk_1d;

tf_qwalk_1d* qcs_create_qwalk_1d( tf_qcs_matrix *coin, int type, int steps );
tf_qwalk_1d* qcs_create_qwalk_1d_mt( tf_qcs_matrix *coin, int type, int steps );

void qcs_qwalk_1d_reset( tf_qwalk_1d *q );
void qcs_qwalk_1d_set_coin( tf_qwalk_1d *q, tf_qcs_matrix *coin);
void qcs_destroy_qwalk_1d( tf_qwalk_1d *q );

void qcs_qwalk_1d_initialise_random_broken_links( tf_qwalk_1d *q );
void qcs_qwalk_1d_set_hadamard_state( tf_qwalk_1d *q );

void qcs_qwalk_1d_random_broken_link( tf_qwalk_1d *q );

void qcs_qwalk_1d_measure_state( tf_qwalk_1d *q );
void qcs_qwalk_1d_random_measure_state( tf_qwalk_1d *q );

void qcs_qwalk_1d_one_iteration( tf_qwalk_1d *q );

void qcs_qwalk_1d_set_probability( tf_qwalk_1d *q );
void qcs_qwalk_1d_update_probability( tf_qwalk_1d *q );
tf_qcs_real_number qcs_qwalk_1d_get_actual_probability( tf_qwalk_1d *q, int n);
tf_qcs_real_number qcs_qwalk_1d_get_average_probability( tf_qwalk_1d *q, int n);

void qcs_qwalk_1d_save_actual_probatility( tf_qwalk_1d *q, char *fname );
void qcs_qwalk_1d_save_average_probatility( tf_qwalk_1d *q, char *fname );

tf_qcs_complex* qcs_qwalk_1d_get_org_a_element( tf_qwalk_1d *q, int x, int y);
void qcs_qwalk_1d_set_org_a_element( tf_qwalk_1d *q, int x, int y, tf_qcs_real_number re, tf_qcs_real_number im);

#ifdef __cplusplus
};
#endif

#endif // __qcs_qwalk_1d_h__
