/***************************************************************************
 *   Copyright (C) 2006 -- 2010 by Marek Sawerwain                         *
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

#include "qcs.h"

#include "qcs_hash_table.h"
#include "qcs_gates.h"


t_qcs_hash_table *qudit_gates_cache;

DYNAMIC_LIB_DECORATION char *get_gate_name(enum _qcs_gate_type t)
{
    switch(t)
    {
        case QCS_ZERO_GATE:
            return "zero";
        break;

        case QCS_ID_GATE:
            return "id";
        break;

        case QCS_NOT_GATE:
            return "not";
        break;

	    case QCS_PAULIX_GATE:
            return "pauli_x";
        break;

	    case QCS_PAULIY_GATE:
            return "pauli_y";
        break;

	    case QCS_PAULIZ_GATE:
            return "pauli_z";
        break;

	    case QCS_HADAMARD_GATE:
            return "hadamard";
        break;

	    case QCS_X_ROT90_GATE:
            return "x_rot90";
        break;

	    case QCS_Y_ROT90_GATE:
            return "y_rot90";
        break;

	    case QCS_Z_ROT90_GATE:
            return "z_rot90";
        break;

	    case QCS_MX_ROT90_GATE:
            return "minus_x_rot90";
        break;

	    case QCS_MY_ROT90_GATE:
            return "minus_y_rot90";
        break;

	    case QCS_MZ_ROT90_GATE:
            return "minus_z_rot90";
        break;

	    case QCS_ROT_ALPHA_GATE:
            return "rotate_alpha";
        break;

	    case QCS_ROT_THETA_GATE:
            return "rotate_theta";
        break;

	    case QCS_T_GATE:
            return "t";
        break;

	    case QCS_S_GATE:
            return "s";
        break;

	    case QCS_PHASE_GATE:
            return "phase";
        break;

	    case QCS_PHASE_F_GATE:
            return "phase_f";
        break;

	    case QCS_HADAMARD_GATE_d:
            return "hadamard_d";
        break;

	    case QCS_PAULIX_GATE_d:
            return "pauli_x_d";
        break;

	    case QCS_PAULIZ_GATE_d:
            return "pauli_z_d";
        break;

        default:
            return NULL;
    }
}

