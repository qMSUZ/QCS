#include "qcs.h"
#include "qcs_misc.h"

/*
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.
   David Bateman added normal and exponential distributions following
   Marsaglia and Tang's Ziggurat algorithm.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   Copyright (C) 2004, David Bateman
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp

   * 2006-04-01 David Bateman
   * * convert for use in octave, declaring static functions only used
   *   here and adding oct_ to functions visible externally
   * * inverse sense of ALLBITS
   * 2004-01-19 Paul Kienzle
   * * comment out main
   * add init_by_entropy, get_state, set_state
   * * converted to allow compiling by C++ compiler
   *
   * 2004-01-25 David Bateman
   * * Add Marsaglia and Tsang Ziggurat code
   *
   * 2004-07-13 Paul Kienzle
   * * make into an independent library with some docs.
   * * introduce new main and test code.
   *
   * 2004-07-28 Paul Kienzle & David Bateman
   * * add -DALLBITS flag for 32 vs. 53 bits of randomness in mantissa
   * * make the naming scheme more uniform
   * * add -DHAVE_X86 for faster support of 53 bit mantissa on x86 arch.
   *
   * 2005-02-23 Paul Kienzle
   * * fix -DHAVE_X86_32 flag and add -DUSE_X86_32=0|1 for explicit control
   */

/*
   === Build instructions ===

   Compile with -DHAVE_GETTIMEOFDAY if the gettimeofday function is
   available.  This is not necessary if your architecture has
   /dev/urandom defined.

   Compile with -DALLBITS to disable 53-bit random numbers. This is about
   50% slower than using 32-bit random numbers.

   Uses implicit -Di386 or explicit -DHAVE_X86_32 to determine if CPU=x86.
   You can force X86 behaviour with -DUSE_X86_32=1, or suppress it with
   -DUSE_X86_32=0. You should also consider -march=i686 or similar for
   extra performance. Check whether -DUSE_X86_32=0 is faster on 64-bit
   x86 architectures.

   If you want to replace the Mersenne Twister with another
   generator then redefine randi32 appropriately.

   === Usage instructions ===
   Before using any of the generators, initialize the state with one of
   oct_init_by_int, oct_init_by_array or oct_init_by_entropy.

   All generators share the same state vector.

   === Mersenne Twister ===
   void qcs_init_by_int(unsigned int s)           32-bit initial state
   void qcs_init_by_array(unsigned int k[],int m) m*32-bit initial state
   void qcs_init_by_entropy(void)             random initial state
   void qcs_get_state(unsigned int save[MT_N+1])  saves state in array
   void qcs_set_state(unsigned int save[MT_N+1])  restores state from array
   static unsigned int randmt(void)               returns 32-bit unsigned int

   === inline generators ===
   static unsigned int randi32(void)   returns 32-bit unsigned int
   static unsigned long long randi53(void)   returns 53-bit unsigned int
   static unsigned long long randi54(void)   returns 54-bit unsigned int
   static unsigned long long randi64(void)   returns 64-bit unsigned int
   static double randu32(void)     returns 32-bit uniform in (0,1)
   static double randu53(void)     returns 53-bit uniform in (0,1)

   double qcs_randu(void)       returns M-bit uniform in (0,1)
   double qcs_randn(void)       returns M-bit standard normal
   double qcs_rande(void)       returns N-bit standard exponential

   === Array generators ===
   void oct_fill_randi32(octave_idx_type, unsigned int [])
   void oct_fill_randi64(octave_idx_type, unsigned long long [])
   void oct_fill_randu(octave_idx_type, double [])
   void oct_fill_randn(octave_idx_type, double [])
   void oct_fill_rande(octave_idx_type, double [])

*/

#if defined (HAVE_CONFIG_H)
#include <config.h>
#endif

#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "qcs_rand.h"

#if !defined(USE_X86_32)
# if defined(i386) || defined(HAVE_X86_32)
#  define USE_X86_32 1
# else
#  define USE_X86_32 0
# endif
#endif

/* ===== Mersenne Twister 32-bit generator ===== */

#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

static unsigned int *next;
static unsigned int state[MT_N]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static int initt = 1;

/* initializes state[MT_N] with a seed */
DYNAMIC_LIB_DECORATION void qcs_mt_init_by_int (unsigned int s)
{
    int j;
    state[0] = s & 0xffffffffUL;
    for (j = 1; j < MT_N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    left = 1;
    initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
DYNAMIC_LIB_DECORATION void qcs_mt_init_by_array (unsigned int *init_key, int key_length)
{
  int i, j, k;
  qcs_mt_init_by_int (19650218UL);
  i = 1;
  j = 0;
  k = (MT_N > key_length ? MT_N : key_length);
  for (; k; k--)
    {
      state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
	+ init_key[j] + j; /* non linear */
      state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      j++;
      if (i >= MT_N)
	{
	  state[0] = state[MT_N-1];
	  i = 1;
	}
      if (j >= key_length)
	j = 0;
    }
  for (k = MT_N - 1; k; k--)
    {
      state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
	- i; /* non linear */
      state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i >= MT_N)
	{
	  state[0] = state[MT_N-1];
	  i = 1;
	}
    }

  state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
  left = 1;
  initf = 1;
}

DYNAMIC_LIB_DECORATION void qcs_mt_init_by_entropy (void)
{
    unsigned int entropy[MT_N];
    int n = 0;

    /* Look for entropy in /dev/urandom */
    FILE* urandom =fopen("/dev/urandom", "rb");
    if (urandom)
      {
	while (n < MT_N)
	  {
	    unsigned char word[4];
	    if (fread(word, 4, 1, urandom) != 1)
	      break;
	    entropy[n++] = word[0]+(word[1]<<8)+(word[2]<<16)+(word[3]<<24);
	  }
	fclose(urandom);
      }

    /* If there isn't enough entropy, gather some from various sources */
    if (n < MT_N)
      entropy[n++] = time(NULL); /* Current time in seconds */
    if (n < MT_N)
      entropy[n++] = clock();    /* CPU time used (usec) */
#ifdef HAVE_GETTIMEOFDAY
    if (n < MT_N)
      {
	struct timeval tv;
	if (gettimeofday(&tv, NULL) != -1)
	  entropy[n++] = tv.tv_usec;   /* Fractional part of current time */
      }
#endif
    /* Send all the entropy into the initial state vector */
    qcs_mt_init_by_array(entropy,n);
}

DYNAMIC_LIB_DECORATION void qcs_mt_set_state (unsigned int *save)
{
  int i;
  for (i = 0; i < MT_N; i++)
    state[i] = save[i];
  left = save[MT_N];
  next = state + (MT_N - left + 1);
}

DYNAMIC_LIB_DECORATION void qcs_mt_get_state (unsigned int *save)
{
  int i;
  for (i = 0; i < MT_N; i++)
    save[i] = state[i];
  save[MT_N] = left;
}

static void next_state (void)
{
  unsigned int *p = state;
  int j;

  /* if init_by_int() has not been called, */
  /* a default initial seed is used         */
  /* if (initf==0) init_by_int(5489UL); */
  /* Or better yet, a random seed! */
  if (initf == 0)
    qcs_mt_init_by_entropy();

  left = MT_N;
  next = state;

  for (j = MT_N - MT_M + 1; --j; p++)
    *p = p[MT_M] ^ TWIST(p[0], p[1]);

  for (j = MT_M; --j; p++)
    *p = p[MT_M-MT_N] ^ TWIST(p[0], p[1]);

  *p = p[MT_M-MT_N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
static unsigned int randmt (void)
{
  register unsigned int y;

  if (--left == 0)
    next_state();
  y = *next++;

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  return (y ^ (y >> 18));
}

/* ===== Uniform generators ===== */

/* Select which 32 bit generator to use */
#define randi32 randmt

static unsigned long long randi53 (void)
{
  const unsigned int lo = randi32();
  const unsigned int hi = randi32()&0x1FFFFF;
#if HAVE_X86_32
  unsigned long long u;
  unsigned int *p = (unsigned int *)&u;
  p[0] = lo;
  p[1] = hi;
  return u;
#else
  return (((unsigned long long)hi<<32)|lo);
#endif
}

static unsigned long long randi54 (void)
{
  const unsigned int lo = randi32();
  const unsigned int hi = randi32()&0x3FFFFF;
#if HAVE_X86_32
  unsigned long long u;
  unsigned int *p = (unsigned int *)&u;
  p[0] = lo;
  p[1] = hi;
  return u;
#else
  return (((unsigned long long)hi<<32)|lo);
#endif
}

#if 0
// FIXME -- this doesn't seem to be used anywhere; should it be removed?
static unsigned long long randi64 (void)
{
  const unsigned int lo = randi32();
  const unsigned int hi = randi32();
#if HAVE_X86_32
  unsigned long long u;
  unsigned int *p = (unsigned int *)&u;
  p[0] = lo;
  p[1] = hi;
  return u;
#else
  return (((unsigned long long)hi<<32)|lo);
#endif
}
#endif

#ifdef ALLBITS /* generates a random number on (0,1)-real-interval */
static double randu32 (void)
{
  return ((double)randi32() + 0.5) * (1.0/4294967296.0);
  /* divided by 2^32 */
}
#else
/* generates a random number on (0,1) with 53-bit resolution */
static double randu53 (void)
{
  const unsigned int a=randi32()>>5;
  const unsigned int b=randi32()>>6;
  return (a*67108864.0+b+0.4) * (1.0/9007199254740992.0);
}
#endif

/* Determine mantissa for uniform doubles */
DYNAMIC_LIB_DECORATION double qcs_randu (void)
{
#ifdef ALLBITS
  return randu32 ();
#else
  return randu53 ();
#endif
}

/* ===== Ziggurat normal and exponential generators ===== */
#ifdef ALLBITS
# define ZIGINT unsigned int
# define EMANTISSA 4294967296.0 /* 32 bit mantissa */
# define ERANDI randi32() /* 32 bits for mantissa */
# define NMANTISSA 2147483648.0 /* 31 bit mantissa */
# define NRANDI randi32() /* 31 bits for mantissa + 1 bit sign */
# define RANDU randu32()
#else
# define ZIGINT unsigned long long
# define EMANTISSA 9007199254740992.0  /* 53 bit mantissa */
# define ERANDI randi53() /* 53 bits for mantissa */
# define NMANTISSA EMANTISSA
# define NRANDI randi54() /* 53 bits for mantissa + 1 bit sign */
# define RANDU randu53()
#endif

#define ZIGGURAT_TABLE_SIZE 256

#define ZIGGURAT_NOR_R 3.6541528853610088
#define ZIGGURAT_NOR_INV_R 0.27366123732975828
#define NOR_SECTION_AREA 0.00492867323399

#define ZIGGURAT_EXP_R 7.69711747013104972
#define ZIGGURAT_EXP_INV_R 0.129918765548341586
#define EXP_SECTION_AREA 0.0039496598225815571993

static ZIGINT ki[ZIGGURAT_TABLE_SIZE];
static double wi[ZIGGURAT_TABLE_SIZE], fi[ZIGGURAT_TABLE_SIZE];
static ZIGINT ke[ZIGGURAT_TABLE_SIZE];
static double we[ZIGGURAT_TABLE_SIZE], fe[ZIGGURAT_TABLE_SIZE];

/*
This code is based on the paper Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software. Code was
presented in this paper for a Ziggurat of 127 levels and using a 32 bit
integer random number generator. This version of the code, uses the
Mersenne Twister as the integer generator and uses 256 levels in the
Ziggurat. This has several advantages.

  1) As Marsaglia and Tsang themselves states, the more levels the few
     times the expensive tail algorithm must be called
  2) The cycle time of the generator is determined by the integer
     generator, thus the use of a Mersenne Twister for the core random
     generator makes this cycle extremely long.
  3) The license on the original code was unclear, thus rewriting the code
     from the article means we are free of copyright issues.
  4) Compile flag for full 53-bit random mantissa.

It should be stated that the authors made my life easier, by the fact that
the algorithm developed in the text of the article is for a 256 level
ziggurat, even if the code itself isn't...

One modification to the algorithm developed in the article, is that it is
assumed that 0 <= x < Inf, and "unsigned long"s are used, thus resulting in
terms like 2^32 in the code. As the normal distribution is defined between
-Inf < x < Inf, we effectively only have 31 bit integers plus a sign. Thus
in Marsaglia and Tsang, terms like 2^32 become 2^31. We use NMANTISSA for
this term.  The exponential distribution is one sided so we use the
full 32 bits.  We use EMANTISSA for this term.

It appears that I'm slightly slower than the code in the article, this
is partially due to a better generator of random integers than they
use. But might also be that the case of rapid return was optimized by
inlining the relevant code with a #define. As the basic Mersenne
Twister is only 25% faster than this code I suspect that the main
reason is just the use of the Mersenne Twister and not the inlining,
so I'm not going to try and optimize further.
*/

static void create_ziggurat_tables (void)
{
  int i;
  double x, x1;

  /* Ziggurat tables for the normal distribution */
  x1 = ZIGGURAT_NOR_R;
  wi[255] = x1 / NMANTISSA;
  fi[255] = exp (-0.5 * x1 * x1);

  /* Index zero is special for tail strip, where Marsaglia and Tsang
   * defines this as
   * k_0 = 2^31 * r * f(r) / v, w_0 = 0.5^31 * v / f(r), f_0 = 1,
   * where v is the area of each strip of the ziggurat.
   */
  ki[0] = (ZIGINT) (x1 * fi[255] / NOR_SECTION_AREA * NMANTISSA);
  wi[0] = NOR_SECTION_AREA / fi[255] / NMANTISSA;
  fi[0] = 1.;

  for (i = 254; i > 0; i--)
    {
      /* New x is given by x = f^{-1}(v/x_{i+1} + f(x_{i+1})), thus
       * need inverse operator of y = exp(-0.5*x*x) -> x = sqrt(-2*ln(y))
       */
      x = sqrt(-2. * log(NOR_SECTION_AREA / x1 + fi[i+1]));
      ki[i+1] = (ZIGINT)(x / x1 * NMANTISSA);
      wi[i] = x / NMANTISSA;
      fi[i] = exp (-0.5 * x * x);
      x1 = x;
    }

  ki[1] = 0;

  /* Zigurrat tables for the exponential distribution */
  x1 = ZIGGURAT_EXP_R;
  we[255] = x1 / EMANTISSA;
  fe[255] = exp (-x1);

  /* Index zero is special for tail strip, where Marsaglia and Tsang
   * defines this as
   * k_0 = 2^32 * r * f(r) / v, w_0 = 0.5^32 * v / f(r), f_0 = 1,
   * where v is the area of each strip of the ziggurat.
   */
  ke[0] = (ZIGINT) (x1 * fe[255] / EXP_SECTION_AREA * EMANTISSA);
  we[0] = EXP_SECTION_AREA / fe[255] / EMANTISSA;
  fe[0] = 1.;

  for (i = 254; i > 0; i--)
    {
      /* New x is given by x = f^{-1}(v/x_{i+1} + f(x_{i+1})), thus
       * need inverse operator of y = exp(-x) -> x = -ln(y)
       */
      x = - log(EXP_SECTION_AREA / x1 + fe[i+1]);
      ke[i+1] = (ZIGINT)(x / x1 * EMANTISSA);
      we[i] = x / EMANTISSA;
      fe[i] = exp (-x);
      x1 = x;
    }
  ke[1] = 0;

  initt = 0;
}

/*
 * Here is the guts of the algorithm. As Marsaglia and Tsang state the
 * algorithm in their paper
 *
 * 1) Calculate a random signed integer j and let i be the index
 *     provided by the rightmost 8-bits of j
 * 2) Set x = j * w_i. If j < k_i return x
 * 3) If i = 0, then return x from the tail
 * 4) If [f(x_{i-1}) - f(x_i)] * U < f(x) - f(x_i), return x
 * 5) goto step 1
 *
 * Where f is the functional form of the distribution, which for a normal
 * distribution is exp(-0.5*x*x)
 */

DYNAMIC_LIB_DECORATION double qcs_randn (void)
{
  if (initt)
    create_ziggurat_tables();

  while (1)
    {
      /* The following code is specialized for 32-bit mantissa.
       * Compared to the arbitrary mantissa code, there is a performance
       * gain for 32-bits:  PPC: 2%, MIPS: 8%, x86: 40%
       * There is a bigger performance gain compared to using a full
       * 53-bit mantissa:  PPC: 60%, MIPS: 65%, x86: 240%
       * Of course, different compilers and operating systems may
       * have something to do with this.
       */
#if !defined(ALLBITS)
# if HAVE_X86_32
      /* 53-bit mantissa, 1-bit sign, x86 32-bit architecture */
      double x;
      int si,idx;
      register unsigned int lo, hi;
      long long rabs;
      unsigned int *p = (unsigned int *)&rabs;
      lo = randi32();
      idx = lo&0xFF;
      hi = randi32();
      si = hi&UMASK;
      p[0] = lo;
      p[1] = hi&0x1FFFFF;
      x = ( si ? -rabs : rabs ) * wi[idx];
# else /* !HAVE_X86_32 */
      /* arbitrary mantissa (selected by NRANDI, with 1 bit for sign) */
      const unsigned long long r = NRANDI;
      const long long rabs=r>>1;
      const int idx = (int)(rabs&0xFF);
      const double x = ( r&1 ? -rabs : rabs) * wi[idx];
# endif /* !HAVE_X86_32 */
      if (rabs < (long long)ki[idx])
#else /* ALLBITS */
      /* 32-bit mantissa */
      const unsigned int r = randi32();
      const unsigned int rabs = r&LMASK;
      const int idx = (int)(r&0xFF);
      const double x = ((int32_t)r) * wi[idx];
      if (rabs < ki[idx])
#endif /* ALLBITS */
	return x;        /* 99.3% of the time we return here 1st try */
      else if (idx == 0)
	{
	  /* As stated in Marsaglia and Tsang
	   *
	   * For the normal tail, the method of Marsaglia[5] provides:
	   * generate x = -ln(U_1)/r, y = -ln(U_2), until y+y > x*x,
	   * then return r+x. Except that r+x is always in the positive
	   * tail!!!! Any thing random might be used to determine the
	   * sign, but as we already have r we might as well use it
	   *
	   * [PAK] but not the bottom 8 bits, since they are all 0 here!
	   */
	  double xx, yy;
	  do
	    {
	      xx = - ZIGGURAT_NOR_INV_R * log (RANDU);
	      yy = - log (RANDU);
	    }
	  while ( yy+yy <= xx*xx);
	  return (rabs&0x100 ? -ZIGGURAT_NOR_R-xx : ZIGGURAT_NOR_R+xx);
	}
      else if ((fi[idx-1] - fi[idx]) * RANDU + fi[idx] < exp(-0.5*x*x))
	return x;
    }
}

DYNAMIC_LIB_DECORATION double qcs_rande (void)
{
  if (initt)
    create_ziggurat_tables();

  while (1)
    {
      ZIGINT ri = ERANDI;
      const int idx = (int)(ri & 0xFF);
      const double x = ri * we[idx];
      if (ri < ke[idx])
	return x;		// 98.9% of the time we return here 1st try
      else if (idx == 0)
	{
	  /* As stated in Marsaglia and Tsang
	   *
	   * For the exponential tail, the method of Marsaglia[5] provides:
           * x = r - ln(U);
	   */
	  return ZIGGURAT_EXP_R - log(RANDU);
	}
      else if ((fe[idx-1] - fe[idx]) * RANDU + fe[idx] < exp(-x))
	return x;
    }
}

/* Array generators */
DYNAMIC_LIB_DECORATION void qcs_fill_randu (int n, double *p)
{
  int i;
  for (i = 0; i < n; i++)
    p[i] = qcs_randu();
}

DYNAMIC_LIB_DECORATION void qcs_fill_randn (int n, double *p)
{
  int i;
  for (i = 0; i < n; i++)
    p[i] = qcs_randn();
}

DYNAMIC_LIB_DECORATION void qcs_fill_rande (int n, double *p)
{
  int i;
  for (i = 0; i < n; i++)
    p[i] = qcs_rande();
}

/*
   This Random Number Generator is based on the algorithm in a FORTRAN
   version published by George Marsaglia and Arif Zaman, Florida State
   University; ref.: see original comments below.
   At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
   Science, we have written sources in further languages (C, Modula-2
   Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
   results compared with the original FORTRAN version.
   April 1989
   Karl-L. Noell <NOELL@DWIFH1.BITNET>
      and  Helmut  Weber <WEBER@DWIFH1.BITNET>

   This random number generator originally appeared in "Toward a Universal
   Random Number Generator" by George Marsaglia and Arif Zaman.
   Florida State University Report: FSU-SCRI-87-50 (1987)
   It was later modified by F. James and published in "A Review of Pseudo-
   random Number Generators"
   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
   (However, a newly discovered technique can yield
   a period of 10^600. But that is still in the development stage.)
   It passes ALL of the tests for random number generators and has a period
   of 2^144, is completely portable (gives bit identical results on all
   machines with at least 24-bit mantissas in the floating point
   representation).
   The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).

   Use IJ = 1802 & KL = 9373 to test the random number generator. The
   subroutine RANMAR should be used to generate 20000 random numbers.
   Then display the next six random numbers generated multiplied by 4096*4096
   If the random number generator is working properly, the random numbers
   should be:
           6533892.0  14220222.0  7275067.0
           6172232.0  8354498.0   10633180.0
*/

/* Globals */
static double u[97],c,cd,cm;
static int i97,j97;
static int test = QCS_FALSE;

/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
   each subsequence having a length of approximately 10^30.
*/
DYNAMIC_LIB_DECORATION void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
		ij = 1802;
		kl = 9373;
   }

   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);

   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }

   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = QCS_TRUE;
}

/*
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
DYNAMIC_LIB_DECORATION double RandomUniform(void)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test)
   	RandomInitialise(1802,9373);

   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;

   return(uni);
}

/*
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number
  with a given mean and standard devaiation.  Calls are made to a
  function subprogram which must return independent random
  numbers uniform in the interval (0,1).
  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  and J.F. Monahan augmented with quadratic bounding curves.
*/
DYNAMIC_LIB_DECORATION double RandomGaussian(double mean,double stddev)
{
   double  q,u,v,x,y;

	/*
		Generate P = (u,v) uniform in rect. enclosing acceptance region
      Make sure that any random numbers <= 0 are rejected, since
      gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
	*/
   do {
      u = RandomUniform();
      v = RandomUniform();
   	if (u <= 0.0 || v <= 0.0) {
       	u = 1.0;
       	v = 1.0;
   	}
      v = 1.7156 * (v - 0.5);

      /*  Evaluate the quadratic form */
      x = u - 0.449871;
   	y = fabs(v) + 0.386595;
      q = x * x + y * (0.19600 * y - 0.25472 * x);

      /* Accept P if inside inner ellipse */
      if (q < 0.27597)
			break;

      /*  Reject P if outside outer ellipse, or outside acceptance region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}

/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
DYNAMIC_LIB_DECORATION int RandomInt(int lower,int upper)
{
   return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random double  within a range, lower -> upper
*/

DYNAMIC_LIB_DECORATION double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}

