#ifndef MULTIGS_UTILS_H
#define MULTIGS_UTILS_H

#if defined(_NDEBUG)
#define DASSERT(expr)
#else
#include <assert.h>
#define DASSERT(expr) assert(expr);
#endif

/* Uniform randomly choose m elements in [0, n-1]
 * Store them into out.
 * Return 0 in this time
 */
int urandint_m(int n, int m, int *out);

/* Uniform randomly choose one element in [0, n-1]
 * Return the element
 */
int urandint_1(int n);

/* Weighted Randomly choose one element in [0, n-1].
 * w is n weights.
 * Return the element.
 */
int wrandint_1(int n, const double *w);



#endif
