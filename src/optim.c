/* optim.c - Optimization functions for scalc
 * Implements quadratic programming for portfolio optimization
 * C89 compliant for Watcom C / DOS
 */

#include "config.h"

#ifdef HAVE_OPTIM

#include <stdio.h>
#include "optim.h"
#include "apf.h"
#include "apfc.h"
#include "matrix.h"

/* Project vector onto probability simplex (sum=1, all >= 0)
 * Uses the algorithm from "Projection onto the probability simplex"
 * by Chen & Ye (2011) - O(n log n) sorting-based method
 * 
 * For small n, we use a simpler O(n^2) method
 */
static void project_simplex(matrix_t *x, int n)
{
    apf sum, one, zero, t, tmp;
    int i, rho;
    apf sorted[MAT_MAX_ELEM];
    
    apf_from_int(&one, 1);
    apf_zero(&zero);
    
    /* First check if already feasible */
    apf_zero(&sum);
    for (i = 0; i < n; i++) {
        apf_add(&sum, &sum, &MAT_AT(x, i, 0).re);
    }
    
    /* Simple clipping + renormalization for non-negative constraint */
    /* Clip to [0, 1] */
    for (i = 0; i < n; i++) {
        if (MAT_AT(x, i, 0).re.sign) {
            apf_zero(&MAT_AT(x, i, 0).re);
        }
        if (apf_cmp(&MAT_AT(x, i, 0).re, &one) > 0) {
            apf_copy(&MAT_AT(x, i, 0).re, &one);
        }
    }
    
    /* Compute sum after clipping */
    apf_zero(&sum);
    for (i = 0; i < n; i++) {
        apf_add(&sum, &sum, &MAT_AT(x, i, 0).re);
    }
    
    /* If sum is zero, use equal weights */
    if (apf_iszero(&sum)) {
        apf equal;
        apf_from_int(&equal, n);
        apf_div(&equal, &one, &equal);
        for (i = 0; i < n; i++) {
            apf_copy(&MAT_AT(x, i, 0).re, &equal);
        }
        return;
    }
    
    /* Normalize to sum to 1 */
    for (i = 0; i < n; i++) {
        apf_div(&tmp, &MAT_AT(x, i, 0).re, &sum);
        apf_copy(&MAT_AT(x, i, 0).re, &tmp);
    }
    
    /* Now project properly using Duchi's algorithm */
    /* Copy and sort in descending order */
    for (i = 0; i < n; i++) {
        apf_copy(&sorted[i], &MAT_AT(x, i, 0).re);
    }
    
    /* Simple insertion sort (fine for small n) */
    {
        int j;
        for (i = 1; i < n; i++) {
            apf_copy(&tmp, &sorted[i]);
            j = i - 1;
            while (j >= 0 && apf_cmp(&sorted[j], &tmp) < 0) {
                apf_copy(&sorted[j+1], &sorted[j]);
                j--;
            }
            apf_copy(&sorted[j+1], &tmp);
        }
    }
    
    /* Find rho */
    apf_zero(&sum);
    rho = 0;
    for (i = 0; i < n; i++) {
        apf rhs;
        apf_add(&sum, &sum, &sorted[i]);
        apf_from_int(&rhs, i + 1);
        apf_mul(&rhs, &rhs, &sorted[i]);
        apf_sub(&tmp, &sum, &one);
        if (apf_cmp(&rhs, &tmp) > 0) {
            rho = i + 1;
        }
    }
    
    if (rho == 0) rho = 1;
    
    /* Compute theta */
    apf_zero(&sum);
    for (i = 0; i < rho; i++) {
        apf_add(&sum, &sum, &sorted[i]);
    }
    apf_sub(&sum, &sum, &one);
    apf_from_int(&tmp, rho);
    apf_div(&t, &sum, &tmp);
    
    /* Project: x_i = max(x_i - t, 0) */
    for (i = 0; i < n; i++) {
        apf_sub(&tmp, &MAT_AT(x, i, 0).re, &t);
        if (tmp.sign) {
            apf_zero(&MAT_AT(x, i, 0).re);
        } else {
            apf_copy(&MAT_AT(x, i, 0).re, &tmp);
        }
        apf_zero(&MAT_AT(x, i, 0).im);
    }
}

/* Compute gradient of (1/2) * x' * H * x + f' * x
 * grad = H * x + f
 */
static void compute_gradient(matrix_t *grad, const matrix_t *H, 
                            const matrix_t *x, const matrix_t *f, int n)
{
    int i, j;
    apf sum, tmp;
    
    for (i = 0; i < n; i++) {
        apf_zero(&sum);
        for (j = 0; j < n; j++) {
            apf_mul(&tmp, &MAT_AT(H, i, j).re, &MAT_AT(x, j, 0).re);
            apf_add(&sum, &sum, &tmp);
        }
        if (f != NULL && f->rows > 0) {
            apf_add(&sum, &sum, &MAT_AT(f, i, 0).re);
        }
        apf_copy(&MAT_AT(grad, i, 0).re, &sum);
        apf_zero(&MAT_AT(grad, i, 0).im);
    }
}

/* Compute objective: (1/2) * x' * H * x + f' * x */
static void compute_objective(apf *fval, const matrix_t *H, 
                             const matrix_t *x, const matrix_t *f, int n)
{
    int i, j;
    apf sum, tmp, half;
    
    apf_from_str(&half, "0.5");
    
    /* Compute x' * H * x */
    apf_zero(&sum);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            apf_mul(&tmp, &MAT_AT(x, i, 0).re, &MAT_AT(H, i, j).re);
            apf_mul(&tmp, &tmp, &MAT_AT(x, j, 0).re);
            apf_add(&sum, &sum, &tmp);
        }
    }
    apf_mul(&sum, &sum, &half);
    
    /* Add f' * x */
    if (f != NULL && f->rows > 0) {
        for (i = 0; i < n; i++) {
            apf_mul(&tmp, &MAT_AT(f, i, 0).re, &MAT_AT(x, i, 0).re);
            apf_add(&sum, &sum, &tmp);
        }
    }
    
    apf_copy(fval, &sum);
}

/* Quadratic Programming via Projected Gradient Descent
 * Minimize: (1/2) * x' * H * x + f' * x
 * Subject to: sum(x) = 1, lb <= x <= ub
 */
int quadprog_solve(matrix_t *x, apf *fval,
                   const matrix_t *H, const matrix_t *f,
                   const matrix_t *Aeq, const matrix_t *beq,
                   const matrix_t *lb, const matrix_t *ub,
                   int n)
{
    matrix_t grad = {0, 0, NULL};
    apf step, one, tmp, prev_fval, diff;
    int iter, i, max_iter;
    int has_bounds;
    
    MAT_INIT_EMPTY(grad);
    
    if (n <= 0 || n > MAT_MAX_ROWS) {
        printf("Error: invalid problem size\n");
        return 0;
    }
    
    /* Check for bounds */
    has_bounds = (lb != NULL && lb->rows > 0) || (ub != NULL && ub->rows > 0);
    
    /* Initialize x with equal weights */
    mat_zero(x, n, 1);
    if (!x->data) return 0;
    {
        apf equal;
        apf_from_int(&equal, n);
        apf_from_int(&one, 1);
        apf_div(&equal, &one, &equal);
        for (i = 0; i < n; i++) {
            apf_copy(&MAT_AT(x, i, 0).re, &equal);
            apf_zero(&MAT_AT(x, i, 0).im);
        }
    }
    
    /* Initialize gradient matrix */
    mat_zero(&grad, n, 1);
    if (!grad.data) return 0;
    
    /* Step size - use 1 / (2 * max eigenvalue estimate) */
    /* For simplicity, use a small fixed step */
    apf_from_str(&step, "0.1");
    
    max_iter = 1000;
    apf_from_int(&prev_fval, 0);
    
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute gradient */
        compute_gradient(&grad, H, x, f, n);
        
        /* Gradient step: x = x - step * grad */
        for (i = 0; i < n; i++) {
            apf_mul(&tmp, &step, &MAT_AT(&grad, i, 0).re);
            apf_sub(&MAT_AT(x, i, 0).re, &MAT_AT(x, i, 0).re, &tmp);
        }
        
        /* Apply bounds if specified */
        if (has_bounds) {
            for (i = 0; i < n; i++) {
                if (lb != NULL && lb->rows > i) {
                    if (apf_cmp(&MAT_AT(x, i, 0).re, &MAT_AT(lb, i, 0).re) < 0) {
                        apf_copy(&MAT_AT(x, i, 0).re, &MAT_AT(lb, i, 0).re);
                    }
                }
                if (ub != NULL && ub->rows > i) {
                    if (apf_cmp(&MAT_AT(x, i, 0).re, &MAT_AT(ub, i, 0).re) > 0) {
                        apf_copy(&MAT_AT(x, i, 0).re, &MAT_AT(ub, i, 0).re);
                    }
                }
            }
        }
        
        /* Project onto simplex (sum = 1, non-negative) if we have equality constraint */
        if (Aeq != NULL && beq != NULL && Aeq->rows > 0) {
            project_simplex(x, n);
        }
        
        /* Check convergence */
        compute_objective(fval, H, x, f, n);
        apf_sub(&diff, fval, &prev_fval);
        apf_abs(&diff, &diff);
        
        if (iter > 0 && diff.exp < fval->exp - AP_BITS + 10) {
            break;  /* Converged */
        }
        
        apf_copy(&prev_fval, fval);
        
        /* Adaptive step size - decrease if oscillating */
        if (iter % 100 == 99) {
            apf half;
            apf_from_str(&half, "0.5");
            apf_mul(&step, &step, &half);
        }
    }
    
    /* Final objective */
    compute_objective(fval, H, x, f, n);
    
    return 1;
}

/* Simplified minimum variance portfolio
 * Minimize: x' * Sigma * x
 * Subject to: sum(x) = 1, x >= 0
 */
int minvar_portfolio(matrix_t *weights, apf *variance,
                     const matrix_t *Sigma, int n)
{
    matrix_t Aeq = {0,0,NULL}, beq = {0,0,NULL}, lb = {0,0,NULL}, ub = {0,0,NULL};
    int i;
    
    /* Set up constraints */
    /* Aeq = ones(1, n), beq = 1 */
    mat_zero(&Aeq, 1, n);
    mat_zero(&beq, 1, 1);
    if (!Aeq.data || !beq.data) return 0;
    
    apf_from_int(&MAT_AT(&beq, 0, 0).re, 1);
    apf_zero(&MAT_AT(&beq, 0, 0).im);
    for (i = 0; i < n; i++) {
        apf_from_int(&MAT_AT(&Aeq, 0, i).re, 1);
        apf_zero(&MAT_AT(&Aeq, 0, i).im);
    }
    
    /* lb = zeros(n, 1), ub = ones(n, 1) */
    mat_zero(&lb, n, 1);
    mat_zero(&ub, n, 1);
    if (!lb.data || !ub.data) return 0;
    
    for (i = 0; i < n; i++) {
        apf_zero(&MAT_AT(&lb, i, 0).re);
        apf_zero(&MAT_AT(&lb, i, 0).im);
        apf_from_int(&MAT_AT(&ub, i, 0).re, 1);
        apf_zero(&MAT_AT(&ub, i, 0).im);
    }
    
    return quadprog_solve(weights, variance, Sigma, NULL, &Aeq, &beq, &lb, &ub, n);
}

#endif /* HAVE_OPTIM */
