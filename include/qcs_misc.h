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

#ifndef __qcs_all_h__
#define __qcs_all_h__

#include <time.h>
#include <stdlib.h>
//#include <mm_malloc.h>


//#define tf_qcs_real_number double
//typedef float tf_qcs_real_number;
typedef double tf_qcs_real_number;

#define QCS_ACCURACY -7
#define QCS_EPS 0.000001
#define QCS_PI  3.14159265358979323846


#define QCS_TRUE 1
#define QCS_FALSE 0

#ifdef DEBUG
#define BADFLAG  0xBADB
#define GOODFLAG 0x6666
#endif

#ifdef BUILD_DYNAMIC_DLL
#define DYNAMIC_LIB_DECORATION __declspec(dllexport)
#endif

#ifdef USE_DYNAMIC_DLL
#define DYNAMIC_LIB_DECORATION __declspec(dllimport)
#endif

#ifdef USE_STATIC_LIB
#define DYNAMIC_LIB_DECORATION
#endif

#if !defined(BUILD_DYNAMIC_DLL) && !defined(USE_DYNAMIC_DLL) && !defined(USE_STATIC_LIB)
#define DYNAMIC_LIB_DECORATION
#endif

#ifdef __cplusplus
extern "C" {
#endif

int qcs_max(int a, int b);
int qcs_min(int a, int b);

int qcs_max_float(float a, float b);
int qcs_min_float(float a, float b);

int qcs_max_double(double a, double b);
int qcs_min_double(double a, double b);

tf_qcs_real_number qcs_max_tf_qcs_real_number(tf_qcs_real_number a, tf_qcs_real_number b);
tf_qcs_real_number qcs_min_tf_qcs_real_number(tf_qcs_real_number a, tf_qcs_real_number b);

tf_qcs_real_number qcs_sqr(tf_qcs_real_number a);

void qcs_dec2bin(int n, int nbits, char *t);
void qcs_dec2base_d(int n, int nbits, int base, char *t);
int qcs_numbase_d2int(char t);

int qcs_bin2dec(char *t);
int qcs_base_d2dec(char *t, int base);

void qcs_str_insert_at(int pos, char bit, char *str);

int qcs_mod(int value, int modulo);

int qcs_newton_symbol(int n, int k);

int qcs_power_mod(int a, int k, int n);
int qcs_fast_power_mod(int b, int x, int m);

#define QKDF qcs_kronecker_delta_function
int qcs_kronecker_delta_function(int k, int j);

double qcs_diffclock(clock_t clock1, clock_t clock2);

void reset_and_start_timer();
double get_elapsed_mcycles();
double get_elapsed_cycles_raw();

int CPUSupportsSSE2();
int CPUSupportsSSE4();
int CPUSupportsAVX();

#ifdef __cplusplus
}
#endif

#endif /* __qcs_all_h__ */
