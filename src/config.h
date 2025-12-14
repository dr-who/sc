/* config.h - Build configuration for scalc
 * 
 * PRESETS:
 *   (default)     - Full build, all features (Linux)
 *   SCALC_MEDIUM  - Most features, no matrices (DOS)
 *   SCALC_TINY    - Basic scientific calculator
 *   SCALC_MINIMAL - BODMAS + sqrt + exp/log + pi/e
 *   SCALC_VIC20   - Use bigint library instead (separate build)
 */

#ifndef CONFIG_H
#define CONFIG_H

/* ========== Precision ========== */
#ifndef AP_LIMBS
#define AP_LIMBS 8    /* 128 bits = 8 x 16-bit limbs */
#endif

/* ========== Platform Presets ========== */

#if defined(SCALC_MINIMAL)
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #ifndef MAX_INPUT
    #define MAX_INPUT 80
  #endif

#elif defined(SCALC_TINY)
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_POW 1
  #define HAVE_FACTORIAL 1
  #define HAVE_VARIABLES 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_HEXBIN 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1

#elif defined(SCALC_MEDIUM)
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_HYPER 1
  #define HAVE_POW 1
  #define HAVE_COMPLEX 1
  #define HAVE_FACTORIAL 1
  #define HAVE_ROMAN 1
  #define HAVE_HEXBIN 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_VARIABLES 1
  #define HAVE_LOOPS 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_SOLVER 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1
  #define HAVE_BITWISE 1
  #define HAVE_STATS 1
  #define HAVE_ORBITAL 1

#else
  /* Full build (default) - everything ON */
  #define SCALC_FULL 1
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_HYPER 1
  #define HAVE_POW 1
  #define HAVE_COMPLEX 1
  #define HAVE_MATRIX 1
  #define HAVE_PLOT 1
  #define HAVE_FACTORIAL 1
  #define HAVE_ROMAN 1
  #define HAVE_HEXBIN 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_VARIABLES 1
  #define HAVE_LOOPS 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_SOLVER 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1
  #define HAVE_BITWISE 1
  #define HAVE_STATS 1
  #define HAVE_NEWTON 1
  #define HAVE_TVM 1
  #define HAVE_ORBITAL 1
#endif

/* ========== Dependencies ========== */

#if defined(HAVE_TRIG) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_TRIG) && !defined(HAVE_CONST)
  #define HAVE_CONST 1
#endif
#if defined(HAVE_HYPER) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_POW) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_COMPLEX) && !defined(HAVE_TRIG)
  #define HAVE_TRIG 1
#endif
#if defined(HAVE_MATRIX) && !defined(HAVE_COMPLEX)
  #define HAVE_COMPLEX 1
#endif
#if defined(HAVE_PLOT) && !defined(HAVE_USER_FUNCS)
  #define HAVE_USER_FUNCS 1
#endif
/* Solver needs sqrt for discriminant, complex for complex roots */
#if defined(HAVE_SOLVER) && !defined(HAVE_SQRT)
  #define HAVE_SQRT 1
#endif
#if defined(HAVE_SOLVER) && !defined(HAVE_COMPLEX)
  #define HAVE_COMPLEX 1
#endif

/* ========== Buffer Sizes ========== */

#ifndef MAX_INPUT
  #define MAX_INPUT 256
#endif

#ifdef HAVE_USER_FUNCS
  #ifndef MAX_FUNCTIONS
    #define MAX_FUNCTIONS 26
  #endif
  #ifndef MAX_FUNC_BODY
    #define MAX_FUNC_BODY 64
  #endif
#else
  #define MAX_FUNCTIONS 0
  #define MAX_FUNC_BODY 0
#endif

#ifdef HAVE_VARIABLES
  #ifndef MAX_SCALAR_VARS
    #define MAX_SCALAR_VARS 26
  #endif
#else
  #define MAX_SCALAR_VARS 0
#endif

#ifdef HAVE_MATRIX
  #ifndef MAX_MATRIX_VARS
    #define MAX_MATRIX_VARS 8
  #endif
#else
  #define MAX_MATRIX_VARS 1  /* Minimum 1 for C89 compliance */
#endif

#endif /* CONFIG_H */
