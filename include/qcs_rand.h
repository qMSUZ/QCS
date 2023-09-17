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

#ifndef __qcs_rand_h__
#define __qcs_rand_h__

#ifdef __cplusplus
extern "C" {
#endif

#define MT_N 624

/* === Mersenne Twister === */
void qcs_mt_init_by_int (unsigned int s);
void qcs_mt_init_by_array (unsigned int init_key[], int key_length);
void qcs_mt_init_by_entropy (void);
void qcs_mt_set_state (unsigned int save[]);
void qcs_mt_get_state (unsigned int save[]);

/* === Scalar Generators === */
double qcs_randu (void);
double qcs_randn (void);
double qcs_rande (void);

/* === Array Generators === */
void qcs_fill_randu (int n, double *p);
void qcs_fill_randn (int n, double *p);
void qcs_fill_rande (int n, double *p);


/* another random number system from: http://local.wasp.uwa.edu.au/~pbourke/other/random/ by Paul Burke*/

void   RandomInitialise(int,int);
double RandomUniform(void);
double RandomGaussian(double,double);
int    RandomInt(int,int);
double RandomDouble(double,double);

#ifdef __cplusplus
};
#endif

#endif /* __qcs_rand_h__ */
