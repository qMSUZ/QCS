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

#include <time.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "qcs.h"
#include "qcs_misc.h"

#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif


static void local_strreverse(char* begin, char* end)
{

    char aux;

    while(end>begin)
	aux=*end, *end--=*begin, *begin++=aux;

}

static void local_itoa(int value, char* str, int base)
{

    static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char* wstr=str;
    int sign;

    // Validate base
    if (base < 2 || base > 35)
    {
        *wstr='\0';
        return;
    }

    // Take care of sign
    if ((sign=value) < 0) value = -value;

    // Conversion. Number is reversed.
    do *wstr++ = num[value%base]; while(value/=base);

    if(sign<0) *wstr++='-';
    *wstr='\0';

    // Reverse string
    local_strreverse(str,wstr-1);
}

DYNAMIC_LIB_DECORATION int qcs_max(int a, int b)
{
    return a > b ? a : b ;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_max_tf_qcs_real_number(tf_qcs_real_number a, tf_qcs_real_number b)
{
    return a > b ? a : b ;
}

DYNAMIC_LIB_DECORATION int qcs_min(int a, int b)
{
    return a > b ? b : a ;
}

DYNAMIC_LIB_DECORATION int qcs_max_float(float a, float b)
{
    return a > b ? a : b ;
}

DYNAMIC_LIB_DECORATION int qcs_min_float(float a, float b)
{
    return a > b ? b : a ;
}

DYNAMIC_LIB_DECORATION int qcs_max_double(double a, double b)
{
    return a > b ? a : b ;
}

DYNAMIC_LIB_DECORATION int qcs_min_double(double a, double b)
{
    return a > b ? b : a ;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_min_tf_qcs_real_number(tf_qcs_real_number a, tf_qcs_real_number b)
{
    return a > b ? b : a ;
}

DYNAMIC_LIB_DECORATION tf_qcs_real_number qcs_sqr(tf_qcs_real_number a)
{
    return a*a;
}

DYNAMIC_LIB_DECORATION void qcs_dec2bin(int n, int nbits, char *t)
{
    int ind=0;
    tf_qcs_real_number l=(tf_qcs_real_number)n;
    tf_qcs_real_number base;

    while(nbits>0)
    {
         nbits--;

         base=1 << nbits; //base=pow(2, nbits);
         //printf("%d %f %f\n", nbits, base, l-base);
         if ( (l-base) >= 0 )
         {
               t[ind]='1';
               l=l-base;
         }
         else
         {
              t[ind]='0';
         }
         ind++;
    }
}

DYNAMIC_LIB_DECORATION void qcs_dec2base_d(int n, int nbits, int base, char *t)
{
    int i;
    char tmp[512];

    memset(&tmp[0], 0, sizeof(tmp));
    local_itoa(n, &tmp[0], base);
    memset(t, '0', nbits);
    strcpy(t+nbits-strlen(tmp), &tmp[0]);
}

DYNAMIC_LIB_DECORATION int qcs_numbase_d2int(char t)
{
    if(t>='0' && t<='9') return t-'0';
    if(t>='a' && t<='z') return t-'a'+10;
    if(t>='A' && t<='Z') return t-'A'+10;

}

DYNAMIC_LIB_DECORATION int qcs_bin2dec(char *t)
{
    int n, i, len, v=1;

    n = 0;
    len = strlen(t);
    for(i = len - 1 ; i >= 0 ; i--)
    {
        if(t[i] == '1')
        {
            n = n + ( v << ( len - i - 1 ));
            //printf("qcs_bin2dec: %d %d\n", v << (len-i-1),len-i-1);
        }
    }
    return n;
}

DYNAMIC_LIB_DECORATION int qcs_base_d2dec(char *t, int base)
{
    int i, v = 0, len;

    len = strlen(t);

    for(i = len - 1; i >= 0  ; i--)
    {
        //printf("-> %d -> pow=%d c=[%c|%d]\n", i, (int)powf(base, i), t[len - i - 1], qcs_numbase_d2int(t[len - i - 1]));
        v = v + (int)powf(base, i) * qcs_numbase_d2int(t[len - i - 1]);
    }

    return v;
}

DYNAMIC_LIB_DECORATION void qcs_str_insert_at(int pos, char bit, char *str)
{
    int i,len;

    len=strlen(str);
    if(pos < len)
    {
        for(i=len+1; i>=pos; i--)
        {
            str[i+1]=str[i];
        }
        str[pos]=bit;
        str[len+1]=0;
    }
    else
    {
        str[pos]=bit;
        str[pos+1]=0;
    }
}

DYNAMIC_LIB_DECORATION int qcs_mod(int n, int p)
{
        if ( n >= 0 )
        {
            if ((n-(n/p)*p) == 0)
                return 0;
            else
                return  n-(n/p)*p;
        }
        else
        {
            if ((n-(n/p)*p) == 0)
                return 0;
            else
                return  p - ( (-n)-((-n)/p)*p);
        }
}

DYNAMIC_LIB_DECORATION int qcs_newton_symbol(int n, int k)
{
    if(n == k) return 1;
    if(k == 0) return 1;

    tf_qcs_real_number t=1,i;

    for( i=1 ; i<=k ; i++ )
    {
        t = t * ( ( n - i + 1 ) / i );
    }

    return (int)t;
}

DYNAMIC_LIB_DECORATION int qcs_power_mod(int a, int k, int n)
{
    int v;

    v=(int)powf(a,k);

    return qcs_mod(v, n);
}

DYNAMIC_LIB_DECORATION int qcs_fast_power_mod(int b, int x, int m)
{
	int erg = 1;
	while ( x > 0 )
	{
		if ( x & 1 ) { erg = ( erg * b ) % m; }
		b = ( b * b ) % m;
		x = x >> 1;
	}
	return erg;
}

int qcs_kronecker_delta_function(int k, int j)
{
    if( k==j )
        return 1;
    else
        return 0;
}

DYNAMIC_LIB_DECORATION double qcs_diffclock(clock_t clock1, clock_t clock2)
{
	double diffticks=clock1-clock2;

        double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
	//double diffms=(diffticks)/(CLOCKS_PER_SEC);

	return diffms;
}

static inline void __cpuid(int info[4], int infoType)
{
    __asm__  __volatile__ ("cpuid" : "=a" (info[0]), "=b" (info[1]), "=c" (info[2]), "=d" (info[3])
                      : "0" (infoType));
}

static inline uint64_t __rdtsc()
 {
        uint32_t low, high;

        __asm__ __volatile__ (
            "xorl %%eax,%%eax \n    cpuid"
            ::: "%rax", "%rbx", "%rcx", "%rdx" );
        __asm__ __volatile__ (
                              "rdtsc" : "=a" (low), "=d" (high));

        return (uint64_t)high << 32 | low;
}

static uint64_t start, end;

DYNAMIC_LIB_DECORATION inline void reset_and_start_timer()
{
    start = __rdtsc();
}

DYNAMIC_LIB_DECORATION inline double get_elapsed_mcycles()
{
    end = __rdtsc();

    return (end-start) / (1000.0 * 1000.0);
}

DYNAMIC_LIB_DECORATION inline double get_elapsed_cycles_raw()
{
        end = __rdtsc();

    return (end-start);
}

DYNAMIC_LIB_DECORATION inline int CPUSupportsSSE2()
{
    int info[4];

    __cpuid(info, 1);
    return (info[3] & (1 << 26));
}

DYNAMIC_LIB_DECORATION inline int CPUSupportsSSE4()
{
    int info[4];

    __cpuid(info, 1);
    return (info[2] & (1 << 19));
}

DYNAMIC_LIB_DECORATION inline int CPUSupportsAVX()
{
    int info[4];

    __cpuid(info, 1);
    return (info[2] & (1 << 28));
}
