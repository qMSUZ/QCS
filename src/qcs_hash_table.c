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

#include<string.h>
#include<stdio.h>

#include "qcs.h"
#include "qcs_hash_table.h"

static char *mystrdup(const char *s)
{
	char *b;

	if(!(b=malloc(strlen(s)+1))) return NULL;

	strcpy(b, s);

	return b;
}

static qcs_hash_table_size def_hashfunc(const char *key)
{
	qcs_hash_table_size hash=0;

	while( *key )
        hash += (unsigned char)*key++;

	return hash;
}

DYNAMIC_LIB_DECORATION t_qcs_hash_table *qcs_hash_table_create(qcs_hash_table_size size, qcs_hash_table_size (*hash_func)(const char *))
{
	t_qcs_hash_table *hashtbl;

	if(!(hashtbl=malloc(sizeof(t_qcs_hash_table)))) return NULL;

	if(!(hashtbl->nodes=calloc(size, sizeof(struct t_qcs_hash_table_node *))))
	{
		free(hashtbl);
		return NULL;
	}

	hashtbl->size = size;

	if(hash_func)
        hashtbl->qcs_hash_table_func = hash_func ;
	else
        hashtbl->qcs_hash_table_func = def_hashfunc ;

	return hashtbl;
}


DYNAMIC_LIB_DECORATION void qcs_hash_table_destroy(t_qcs_hash_table *hash_tbl)
{
	qcs_hash_table_size n;
	t_qcs_hash_table_node *node, *oldnode;

	for( n = 0 ; n < hash_tbl->size ; n++ )
	{
		node = hash_tbl->nodes[n];

		while(node)
		{
			free(node->key);
			oldnode=node;
			node=node->next;
			free(oldnode);
		}
	}

	free(hash_tbl->nodes);
	free(hash_tbl);
}

DYNAMIC_LIB_DECORATION int qcs_hash_table_insert(t_qcs_hash_table *hash_tbl, const char *key, void *data)
{
	t_qcs_hash_table_node *node;

	qcs_hash_table_size hash=hash_tbl->qcs_hash_table_func(key) % hash_tbl->size;


//	fprintf(stderr, "hashtbl_insert() key=%s, hash=%d, data=%s\n", key, hash, (char*)data);

	node = hash_tbl->nodes[hash];
	while(node)
	{
		if(!strcmp(node->key, key))
		{
			node->data=data;
			return 0;
		}
		node = node->next;
	}


	if(!(node=malloc(sizeof(t_qcs_hash_table_node)))) return -1;
	if(!(node->key=mystrdup(key)))
	{
		free(node);
		return -1;
	}
	node->data=data;
	node->next=hash_tbl->nodes[hash];
	hash_tbl->nodes[hash]=node;


	return 0;
}

DYNAMIC_LIB_DECORATION int qcs_hash_table_remove(t_qcs_hash_table *hash_tbl, const char *key)
{
	t_qcs_hash_table_node *node, *prev_node=NULL;
	qcs_hash_table_size hash=hash_tbl->qcs_hash_table_func(key) % hash_tbl->size;

	node=hash_tbl->nodes[hash];
	while(node)
	{
		if(!strcmp(node->key, key))
		{
			free(node->key);

			if(prev_node)
                prev_node->next=node->next;
			else
                hash_tbl->nodes[hash]=node->next;

			free(node);

			return 0;
		}
		prev_node=node;
		node=node->next;
	}

	return -1;
}

DYNAMIC_LIB_DECORATION void *qcs_hash_table_get(t_qcs_hash_table *hash_tbl, const char *key)
{
	t_qcs_hash_table_node *node;
	qcs_hash_table_size hash = hash_tbl->qcs_hash_table_func(key) % hash_tbl->size;

//	fprintf(stderr, "hashtbl_get() key=%s, hash=%d\n", key, hash);

	node=hash_tbl->nodes[hash];
	while(node)
	{
		if(!strcmp(node->key, key)) return node->data;
		node  = node->next;
	}

	return NULL;
}

DYNAMIC_LIB_DECORATION void qcs_hash_table_basic_display(t_qcs_hash_table *hash_tbl, FILE *f)
{
	qcs_hash_table_size n;
	t_qcs_hash_table_node *node, *oldnode;

	for( n = 0 ; n < hash_tbl->size ; n++ )
	{
	    fprintf(f, "entry: %ld [\n", n);
		node = hash_tbl->nodes[n];
		while(node)
		{
			fprintf(f, " %s", node->key);
			node=node->next;
		}
		fprintf(f, "\n]\n");
	}

}
