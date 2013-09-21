#include "ieeeftn.h"

#if STDC
LOGICAL
isinf(REAL *x)
#else /* NOT STDC */
LOGICAL
isinf(x)
REAL *x;
#endif /* STDC */
{
    REAL_PARTS w;

    w.r = *x;
    return (((GET_EXPONENT_SP(w.i) == EXPONENT_INFNAN_SP) &&
	     (SIGNIFICAND_SP(w.i) == 0))
	    ? _TRUE_ : _FALSE_);
}
#include "ieeeftn.h"

#if STDC
LOGICAL isnan(REAL *x)
#else /* NOT STDC */
LOGICAL isnan(x)
REAL *x;
#endif /* STDC */
{
#if defined(__alpha)
    if (isden(x))
	return (_FALSE_);
    else
	return ((*x != *x) ? _TRUE_ : _FALSE_);
#else
    return ((*x != *x) ? _TRUE_ : _FALSE_);
#endif
}
/*
   quicksort.c
*/

/*
Quick Sort (Sorting array A[size])

    While Low is less than High
    {
        Choose Pivot as the element at position A[Low]
        While A[High] is greater than Pivot, decrement High; else move A[High] to A[Low]
        While A[Low] is less than Pivot, increment Low; else move A[Low] to A[High]
    }
    Move Pivot into A[High], see Pivot position as High.
    If Low is less than Pivot point, recursively call Quick Sort with Low =
        Low, High = Pivot point - 1
    If High is greater than Pivot point, recursively call Quick Sort with Low =
        Pivot point + 1, High = High.
*/

#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
        
void quiksort(float *sort, int low, int high){
   float temp;
   float pivot;
   int m;
   int i;
   if(low < high){
      SWAP(sort[low], sort[(high+low)/2]);
      pivot = sort[low];
      m = low;
      for (i = low + 1; i <= high; i++){
			if (sort[i] <  pivot) {
				m++;
                SWAP(sort[m], sort[i]);
			}
      }
        SWAP(sort[low], sort[m]);
		quiksort(sort, low, m - 1);
		quiksort(sort, m + 1, high);
   }
}
void quicksort_(void *sort, int *lowf, int *highf){
   int low;
   int high;
   low = *lowf;
   high = *highf;
   quiksort(sort,low,high);
}
