/***************************************************************************
 *   Copyright (C) 2005 -- 2011 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
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

#ifndef __qcs_hash_table_h__
#define __qcs_hash_table_h__

#include <stdio.h>
#include <stdlib.h>

typedef size_t qcs_hash_table_size;

typedef struct s_qcs_hash_table_node {
	char *key;
	void *data;
	struct s_qcs_hash_table_node *next;
} t_qcs_hash_table_node ;


typedef struct {
	qcs_hash_table_size size;
	t_qcs_hash_table_node **nodes;
	qcs_hash_table_size (*qcs_hash_table_func)(const char *);
} t_qcs_hash_table;

t_qcs_hash_table *qcs_hash_table_create(qcs_hash_table_size, qcs_hash_table_size (*hash_func)(const char *));
void qcs_hash_table_destroy(t_qcs_hash_table *hash_tbl);
int qcs_hash_table_insert(t_qcs_hash_table *hash_tbl, const char *key, void *data);
int qcs_hash_table_remove(t_qcs_hash_table *hash_tbl, const char *key);
void *qcs_hash_table_get(t_qcs_hash_table *hash_tbl, const char *key);

void qcs_hash_table_basic_display(t_qcs_hash_table *hash_tbl, FILE *f);

#endif /* __qcs_hash_table_h__ */
