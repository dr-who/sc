/*
 * decomp_funcs.h - Table-driven decomposition function dispatch
 * 
 * Handles functions that can return multiple matrices:
 * - Single output: svd(A) -> singular values
 * - Multi output: [U,S,V] = svd(A)
 * 
 * C89 compliant for Watcom C / DOS
 */
#ifndef DECOMP_FUNCS_H
#define DECOMP_FUNCS_H

#include "matrix.h"
#include "apfc.h"

/* Maximum number of output matrices from a decomposition */
#define DECOMP_MAX_OUTPUTS 4

/* Decomposition function info */
typedef struct {
    const char *name;
    int min_outputs;      /* Minimum outputs (for single-output call) */
    int max_outputs;      /* Maximum outputs */
    int min_inputs;       /* Minimum input arguments */
    int max_inputs;       /* Maximum input arguments */
    /* Function pointer set at dispatch time based on output count */
} DecompFuncEntry;

/* Result structure for decompositions */
typedef struct {
    int num_outputs;
    matrix_t outputs[DECOMP_MAX_OUTPUTS];
} DecompResult;

/* Initialize decomposition function table */
void decomp_funcs_init(void);

/* Look up a decomposition function by name */
DecompFuncEntry *decomp_func_lookup(const char *name);

/* Execute a decomposition with specified number of outputs
 * Returns 1 on success, 0 on error
 * nout = number of outputs requested (1 for single, 2+ for multi)
 * result->num_outputs set to actual outputs produced
 */
int decomp_func_exec(const char *name, int nout, 
                     DecompResult *result, 
                     const matrix_t *input, 
                     const matrix_t *input2);  /* input2 can be NULL */

#endif /* DECOMP_FUNCS_H */
