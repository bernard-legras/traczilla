/***********************************************************************
This file defines typedefs and symbols for interfacing C IEEE 754
support functions to Fortran code.

The model representation of a t-digit floating-point number is

	x = (-1)**s * 0.d1 d2 d3 ... dt * beta**e

where the digits dk satisfy

	0 <= dk < beta

The fractional part, which we call the significand, is defined to lie
in the range

	1/beta <= significand < 1

For IEEE floating-point, with its hidden bit and denormalized numbers,
we adjust parameters to conform to our model.  Denormalized numbers
are normalized by expanding their exponent range.

IEEE floating point arithmetic has these formats, where s is the sign
bit, e is an exponent bit, and f is a fraction bit.

Single precision:
	seee eeee efff ffff ffff ffff ffff ffff

	significand = 1.fff ffff ffff ffff ffff ffff (1 + 23 bits)

	exponent bias = 127

Double precision:
	seee eeee eeee ffff ffff ffff ffff ffff
	ffff ffff ffff ffff ffff ffff ffff ffff

	significand = 1.ffff ffff ffff ffff ffff ffff ffff ffff
			ffff ffff ffff ffff ffff (1 + 52 bits)

	exponent bias = 1023

Here are some sample IEEE bit patterns:

========================================================================
			    LITTLE ENDIAN

Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	        Infinity	0x7f800000
	       -Infinity	0xff800000
	             NaN	0x7f80ffff

Double precision
	               0	0x00000000 00000000
	               1	0x00000000 3ff00000
	              -1	0x00000000 bff00000
	               2	0x00000000 40000000
	              -2	0x00000000 c0000000
	     1.11022e-16	0x00000002 3ca00000	eps(1.0)
	    -1.11022e-16	0x00000002 bca00000
	    2.22507e-308	0x00000000 00100000	smallest normal
	   -2.22507e-308	0x00000000 80100000
	    4.94066e-324	0x00000001 00000000	smallest subnormal
	   -4.94066e-324	0x00000001 80000000
	    1.79769e+308	0xffffffff 7fefffff	largest normal
	   -1.79769e+308	0xffffffff ffefffff
	        Infinity	0x00000000 7ff00000
	       -Infinity	0x00000000 fff00000
	             NaN	0xffffffff 7ff7ffff

========================================================================
			     BIG ENDIAN
Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	             Inf	0x7f800000
	            -Inf	0xff800000
	             NaN	0x7fffffff

Double precision
	               0	0x00000000 00000000
	               1	0x3ff00000 00000000
	              -1	0xbff00000 00000000
	               2	0x40000000 00000000
	              -2	0xc0000000 00000000
	     1.11022e-16	0x3ca00000 00000002	eps(1.0)
	    -1.11022e-16	0xbca00000 00000002
	    2.22507e-308	0x00100000 00000000	smallest normal
	   -2.22507e-308	0x80100000 00000000
	    4.94066e-324	0x00000000 00000001	smallest subnormal
	   -4.94066e-324	0x80000000 00000001
	    1.79769e+308	0x7fefffff ffffffff	largest normal
	   -1.79769e+308	0xffefffff ffffffff
	             Inf	0x7ff00000 00000000
	            -Inf	0xfff00000 00000000
	             NaN	0x7fffffff ffffffff

========================================================================

***********************************************************************/

/* type mappings from Fortran to C */
typedef double DOUBLEPRECISION;

typedef float REAL;

#if defined(__alpha)
typedef int INTEGER;	/* need 32-bit integers, not 64-bit ones */
#else
typedef long INTEGER;
#endif

typedef int LOGICAL;	/* use int, not long, to avoid conflicts on HP-UX with */
			/* system header file declarations of isinf(), isnan() */
typedef union
{
    REAL r;
    INTEGER i;
} REAL_PARTS;				/* decomposition of REAL */

typedef union
{
    DOUBLEPRECISION r;
    INTEGER i[2];
} DOUBLEPRECISION_PARTS;		/* decomposition of DOUBLEPRECISION */

/* Fortran LOGICAL values -- compiler dependent!  Most (all?) UNIX */
/* Fortran compilers use 0 for .FALSE. and non-zero for .TRUE., like C. */

#ifdef _FALSE_
#undef _FALSE_
#endif

#ifdef _TRUE_
#undef _TRUE_
#endif

#define _FALSE_				((LOGICAL)0)
#define _TRUE_				((LOGICAL)1)

#define BASE				2 /* number base */

/* stored significand bits in single and double precision */

#define	T_SP				23
#define	T_DP				52

#define BASE_TO_THE_T_SP		((REAL)8388608.0)
#define BASE_TO_THE_T_DP		((DOUBLEPRECISION)4503599627370496.0)

#define _SHIFTED_EXPONENT_MASK_SP	0xff
#define EXPONENT_MASK_SP		0x7f800000L
#define EXPONENT_MASK_DP		0x7ff00000L

#define _SHIFTED_EXPONENT_MASK_DP	0x7ff
#define _SIGNIFICAND_MASK_SP		0x007fffffL
#define _SIGNIFICAND_MASK_DP		0x000fffffL

/* Exponent biases such that significand lies in (1/beta) <= significand < 1.
These are 1 less than the IEEE biases, because its stored significand
lies in 1 <= significand < beta due to the hidden bit.  We define them
with a leading underscore because they are for internal use only. */

#define _BIAS_SP			126
#define _BIAS_DP			1022

/* Indexes into two-word INTEGER array to account for addressing order */

#define i386 				1

#if i386 || sun386 || __i386__ || __sun386__ || msdos || MIPSEL || __alpha
					/* Intel 80xxx or MIPS little endian */
#define DP_LOW				0
#define DP_HIGH				1
#else				/* big endian (MIPS, Motorola, SPARC, ...) */
#define DP_LOW				1
#define DP_HIGH				0
#endif

/* macros to extract (high-order) significand and exponent as integer values */

#define GET_EXPONENT_SP(x) ((((x) >> T_SP) & \
				_SHIFTED_EXPONENT_MASK_SP) - _BIAS_SP)
#define GET_EXPONENT_DP(x) ((((x) >> (T_DP - 32)) & \
				_SHIFTED_EXPONENT_MASK_DP) - _BIAS_DP)

#define SET_EXPONENT_SP(x)	(((x) + _BIAS_SP) << T_SP)
#define SET_EXPONENT_DP(x)	(((x) + _BIAS_DP) << (T_DP - 32))

#define EXPONENT_DENORM_SP		(-_BIAS_SP)
#define EXPONENT_DENORM_DP		(-_BIAS_DP)

#define EXPONENT_INFNAN_SP		(255 - _BIAS_SP)
#define EXPONENT_INFNAN_DP		(2047 - _BIAS_DP)

#define SIGNIFICAND_SP(x)		(((x) & _SIGNIFICAND_MASK_SP))
#define SIGNIFICAND_DP(x)		(((x) & _SIGNIFICAND_MASK_DP))

#define MAX_NORMAL_SP			0x7f7fffffL
#define MAX_NORMAL_DP			0x7fefffffL
#define MAX_NORMAL_Low_DP		0xffffffffL

#define MIN_NORMAL_SP			0x00800000L
#define MIN_DENORMAL_SP			0x00000001L

#define MIN_NORMAL_DP			0x00100000L
#define MIN_NORMAL_Low_DP		0x00000000L
#define MIN_DENORMAL_DP			0x00000000L
#define MIN_DENORMAL_Low_DP		0x00000001L

#define Inf_SP				0x7f800000L
#define NegInf_SP			0xff800000L
#define NaN_SP				0x7fffffffL /* significand is */
						    /* arbitrary non-zero */

/* High-order words for double-precision Infinity and NaN. */
#define Inf_DP				0x7ff00000L
#define Inf_Low_DP			0x00000000L
#define NegInf_DP			0xfff00000L
#define NegInf_Low_DP			0x00000000L
#define NaN_DP				0x7fffffffL /* significand is */
#define NaN_Low_DP			0xffffffffL /* arbitrary non-zero */

#define ISNEG_SP(x)			((x) & 0x80000000L)
#define ISNEG_DP(x)			((x) & 0x80000000L)

#ifndef ABS
#define ABS(x)				(((x) < 0.0) ? -(x) : (x))
#endif

/* Map external names onto the conventions of the
local Fortran compiler.  */

#if _AIX || __hppa
#define _FTN_NAME_TYPE	1	/* name unchanged */
#endif

#if ardent
#define _FTN_NAME_TYPE	2	/* name -> NAME */
#endif

#ifndef _FTN_NAME_TYPE
#define _FTN_NAME_TYPE	3	/* name -> name_ */
#endif

#if (_FTN_NAME_TYPE == 2)
#define adx	ADX
#define dadx	DADX
#define deps	DEPS
#define deps2	DEPS2
#define dintxp	DINTXP
#define disden	DISDEN
#define disinf	DISINF
#define disnan	DISNAN
#define dsetxp	DSETXP
#define eps	EPS
#define eps2	EPS2
#define intxp	INTXP
#define isden	ISDEN
#define isinf	ISINF
#define isnan	ISNAN
#define setxp	SETXP
#endif /* (_FTN_NAME_TYPE == 2) */

#if (_FTN_NAME_TYPE == 3)
#define adx	adx_
#define dadx	dadx_
#define deps	deps_
#define deps2	deps2_
#define dintxp	dintxp_
#define disden	disden_
#define disinf	disinf_
#define disnan	disnan_
#define dsetxp	dsetxp_
#define eps	eps_
#define eps2	eps2_
#define intxp	intxp_
#define isden	isden_
#define isinf	isinf_
#define isnan	isnan_
#define setxp	setxp_
#endif /* (_FTN_NAME_TYPE == 3) */

#if defined(__STDC__) || defined(__cplusplus)
#define ARGS(plist)	plist
#define STDC		1
#define VOID_ARG	void
#else
#define const
#define ARGS(plist)	()
#define STDC		0
#define VOID_ARG
#endif

#if defined(__cplusplus)
extern "C" {
#endif

DOUBLEPRECISION	dadx ARGS((DOUBLEPRECISION *x__, INTEGER *n__));
DOUBLEPRECISION	deps ARGS((DOUBLEPRECISION *x__));
INTEGER		dintxp ARGS((DOUBLEPRECISION *x__));
DOUBLEPRECISION	dsetxp ARGS((DOUBLEPRECISION *x__, INTEGER *n__));

LOGICAL		disden ARGS((DOUBLEPRECISION *x__));
LOGICAL		disinf ARGS((DOUBLEPRECISION *x__));
LOGICAL		disnan ARGS((DOUBLEPRECISION *x__));

REAL		adx ARGS((REAL *x__, INTEGER *n__));
REAL		eps ARGS((REAL *x__));
INTEGER		intxp ARGS((REAL *x__));
REAL		setxp ARGS((REAL *x__, INTEGER *n__));

LOGICAL		isden ARGS((REAL *x__));
LOGICAL		isinf ARGS((REAL *x__));
LOGICAL		isnan ARGS((REAL *x__));

#if defined(__cplusplus)
};
#endif

#if STDC
#include <stdlib.h>
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
