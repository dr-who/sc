/* solver.h - Equation solvers header
 * C89 portable for DOS, Linux
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "config.h"

#ifdef HAVE_SOLVER

#include "apfc.h"

/* Solve quadratic equation ax² + bx + c = 0
 * Returns: 0 = no solution, 1 = one root (double), 2 = two roots
 */
int solve_quadratic(const apfc *a, const apfc *b, const apfc *c,
                    apfc *x1, apfc *x2);

/* Command handlers */
void cmd_quadratic(const char *args);
void cmd_solve(const char *args);

/* Format a complex result for display */
void format_complex_result(const apfc *val);

#endif /* HAVE_SOLVER */

#endif /* SOLVER_H */
