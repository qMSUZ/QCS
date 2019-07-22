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

#ifndef __qcs_gates_h__
#define __qcs_gates_h__

#include "qcs_hash_table.h"
#include "qcs_complex.h"
#include "qcs_matrix.h"

extern t_qcs_hash_table *qudit_gates_cache;

typedef struct {
        int gate_type;
        int freedom;
        tf_qcs_matrix *m;
} t_qcs_qudit_gate_cache_element;

typedef t_qcs_qudit_gate_cache_element *pt_qcs_qudit_gate_cache_element;

enum _qcs_gate_type {
	QCS_ZERO_GATE = 0,
	QCS_ID_GATE,
	QCS_NOT_GATE,
	QCS_PAULIX_GATE,
	QCS_PAULIY_GATE,
	QCS_PAULIZ_GATE,
	QCS_HADAMARD_GATE,
	QCS_SQUARE_ROOT_GATE,
	QCS_X_ROT90_GATE,
	QCS_Y_ROT90_GATE,
	QCS_Z_ROT90_GATE,
	QCS_MX_ROT90_GATE,
	QCS_MY_ROT90_GATE,
	QCS_MZ_ROT90_GATE,
	QCS_ROT_ALPHA_GATE,
	QCS_ROT_THETA_GATE,
	QCS_T_GATE,
	QCS_S_GATE,
	QCS_PHASE_GATE,
	QCS_PHASE_F_GATE,
    QCS_HADAMARD_GATE_d,
    QCS_PAULIX_GATE_d,
    QCS_PAULIZ_GATE_d
};

tf_qcs_matrix *get_gate(enum _qcs_gate_type t);
char *get_gate_name(enum _qcs_gate_type t);



#endif /* __qcs_gates_h__ */
