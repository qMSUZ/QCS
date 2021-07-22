/***************************************************************************
 *   Copyright (C) 2005 -- 2010 by Marek Sawerwain <qcs@gmail.com>         *
 *                                                                         *
 *   Part of the Quantum Computing Simulator:                              *
 *   http://code.google.com/p/qcs                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
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

#ifndef __qcs_info_h__
#define __qcs_info_h__

#ifdef __cplusplus
extern "C" {
#endif


#define ADT_SHORT 0
#define ADT_FULL 1

extern int _amplitude_disp_type;
extern int no_state_vector_max_qubits;
extern int use_state_vector_max_qubits;
extern int use_graph_state_desc_max_qubits;
extern int pqc_mode_max_qubits;
extern int chp_mode_max_qubits;

char* version();
char* compile_system();
char* compilator_name();

void amplitude_display_format(int _amplitude_disp_type);
void info();
void qcs_core_library_initialization();
void qcs_core_library_deinitialization();

#ifdef __cplusplus
}
#endif

#endif /* __qcs_info_h__ */
