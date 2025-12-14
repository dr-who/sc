/* newton.h - Newton-Raphson root finder
 * C89 portable for DOS, Linux
 * 
 * Solves f(x) = 0 using Newton-Raphson iteration
 */

#ifndef NEWTON_H
#define NEWTON_H

#include "config.h"

#ifdef HAVE_NEWTON

#include "apf.h"
#include "apfc.h"

/* Maximum iterations */
#ifndef NEWTON_MAX_ITER
#define NEWTON_MAX_ITER 100
#endif

/* Result codes */
#define NEWTON_OK 0
#define NEWTON_NO_CONVERGE 1
#define NEWTON_ZERO_DERIV 2

/* Solve f(x) = 0 starting from x0
 * func_idx: index of user function (0-25 for a-z)
 * x0: initial guess
 * result: solution
 * Returns: NEWTON_OK on success, error code otherwise
 */
int newton_solve(int func_idx, const apf *x0, apf *result);

/* Process solve command */
void cmd_newton(const char *args);

#endif /* HAVE_NEWTON */
#endif /* NEWTON_H */
