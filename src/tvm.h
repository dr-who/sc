/* tvm.h - Time Value of Money functions
 * C89 portable for DOS, Linux
 * 
 * HP-12C style TVM solver:
 *   N   - number of periods
 *   I%  - interest rate per period (percent)
 *   PV  - present value
 *   PMT - payment per period
 *   FV  - future value
 * 
 * Formula: PV + PMT * ((1-(1+i)^-n)/i) + FV * (1+i)^-n = 0
 */

#ifndef TVM_H
#define TVM_H

#include "config.h"

#ifdef HAVE_TVM

#include "apf.h"

/* TVM registers */
typedef struct {
    apf n;      /* Number of periods */
    apf i;      /* Interest rate per period (as decimal, not %) */
    apf pv;     /* Present value */
    apf pmt;    /* Payment */
    apf fv;     /* Future value */
    int beg;    /* Begin mode (1) vs End mode (0) */
} tvm_regs_t;

extern tvm_regs_t tvm_regs;

/* Clear TVM registers */
void tvm_clear(void);

/* Set TVM values */
void tvm_set_n(const apf *val);
void tvm_set_i(const apf *val);   /* Takes percent, stores decimal */
void tvm_set_pv(const apf *val);
void tvm_set_pmt(const apf *val);
void tvm_set_fv(const apf *val);
void tvm_set_beg(int beg);

/* Calculate unknown value */
int tvm_calc_n(apf *result);
int tvm_calc_i(apf *result);      /* Returns percent */
int tvm_calc_pv(apf *result);
int tvm_calc_pmt(apf *result);
int tvm_calc_fv(apf *result);

/* Show TVM status */
void tvm_show(void);

/* Process TVM command */
void cmd_tvm(const char *args);

#endif /* HAVE_TVM */
#endif /* TVM_H */
