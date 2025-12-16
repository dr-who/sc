/* optim.h - Optimization functions for scalc
 * Implements quadratic programming for portfolio optimization
 * C89 compliant for Watcom C / DOS
 */
#ifndef OPTIM_H
#define OPTIM_H

#include "config.h"

#ifdef HAVE_OPTIM

#include "apf.h"
#include "matrix.h"

/* Quadratic Programming Solver
 * Minimize: (1/2) * x' * H * x + f' * x
 * Subject to: Aeq * x = beq  (equality constraints)
 *             lb <= x <= ub   (bound constraints)
 *
 * Parameters:
 *   x      - output: optimal solution vector (n x 1)
 *   fval   - output: optimal objective value
 *   H      - Hessian matrix (n x n, symmetric positive semi-definite)
 *   f      - linear term (n x 1), can be NULL for zero
 *   Aeq    - equality constraint matrix (m x n), can be NULL
 *   beq    - equality constraint RHS (m x 1), can be NULL
 *   lb     - lower bounds (n x 1), can be NULL for -inf
 *   ub     - upper bounds (n x 1), can be NULL for +inf
 *   n      - number of variables
 *
 * Returns: 1 on success, 0 on failure
 */
int quadprog_solve(matrix_t *x, apf *fval,
                   const matrix_t *H, const matrix_t *f,
                   const matrix_t *Aeq, const matrix_t *beq,
                   const matrix_t *lb, const matrix_t *ub,
                   int n);

/* Simplified portfolio optimization
 * Minimize: x' * Sigma * x  (portfolio variance)
 * Subject to: sum(x) = 1    (fully invested)
 *             x >= 0        (no short selling)
 *
 * Parameters:
 *   weights - output: optimal weights (n x 1)
 *   variance - output: optimal portfolio variance
 *   Sigma   - covariance matrix (n x n)
 *   n       - number of assets
 *
 * Returns: 1 on success, 0 on failure
 */
int minvar_portfolio(matrix_t *weights, apf *variance,
                     const matrix_t *Sigma, int n);

#endif /* HAVE_OPTIM */
#endif /* OPTIM_H */
