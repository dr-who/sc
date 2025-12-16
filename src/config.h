/* config.h - Build configuration for scalc
 * 
 * PRESETS (define ONE before including):
 *   (default)     - Full build, all features (Linux/modern)
 *   SCALC_MEDIUM  - Most features, no matrices (DOS 16-bit)
 *   SCALC_TINY    - Basic scientific calculator (8-bit friendly)
 *   SCALC_MINIMAL - BODMAS + sqrt + exp/log + pi/e (embedded)
 *   SCALC_VIC20   - 8-bit VIC-20 build (5KB RAM, CORDIC)
 *
 * PHASES (from research):
 *   Phase 1: Core foundation - arithmetic, trig, exp/log, RPN
 *   Phase 2: Scientific essentials - gamma, erf, bessel, elliptic
 *   Phase 3: Professional - units, integration, ODE, distributions
 *   Phase 4: Advanced - Lambert W, eigenvalues, QR, symbolic
 */

#ifndef CONFIG_H
#define CONFIG_H

/* ========== Platform Detection ========== */
#if defined(__WATCOMC__)
  #define PLATFORM_DOS 1
  #define PLATFORM_16BIT 1
#elif defined(__CC65__) || defined(__C64__) || defined(__VIC20__)
  #define PLATFORM_8BIT 1
  #define PLATFORM_VIC20 1
#elif defined(_WIN32) || defined(_WIN64)
  #define PLATFORM_WINDOWS 1
#else
  #define PLATFORM_UNIX 1
#endif

/* ========== Precision ========== */
#if defined(PLATFORM_8BIT) || defined(SCALC_VIC20)
  /* 8-bit: 64 bits = 4 x 16-bit limbs (saves RAM) */
  #ifndef AP_LIMBS
    #define AP_LIMBS 4
  #endif
#else
  /* 16-bit+: 128 bits = 8 x 16-bit limbs */
  #ifndef AP_LIMBS
    #define AP_LIMBS 8
  #endif
#endif

/* ========== Algorithm Selection ========== */
/* CORDIC: shift-add only, no multiply - best for 8-bit */
/* TAYLOR: faster with hardware multiply - best for 16-bit+ */
#if defined(PLATFORM_8BIT)
  #define USE_CORDIC 1
#else
  #define USE_TAYLOR 1
#endif

/* ========== Platform Presets ========== */

#if defined(SCALC_VIC20)
  /* VIC-20: 5KB RAM, 6502 @ 1MHz - absolute minimum */
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1        /* CORDIC-based */
  #define HAVE_RPN 1         /* 4-level stack with LASTX */
  #define HAVE_DISPLAY_MODES 1
  #ifndef MAX_INPUT
    #define MAX_INPUT 40
  #endif
  #ifndef RPN_STACK_SIZE
    #define RPN_STACK_SIZE 4  /* Classic HP 4-level */
  #endif

#elif defined(SCALC_MINIMAL)
  /* Embedded: minimal footprint */
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_DISPLAY_MODES 1
  #ifndef MAX_INPUT
    #define MAX_INPUT 80
  #endif

#elif defined(SCALC_TINY)
  /* Phase 1: Core foundation */
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_HYPER 1
  #define HAVE_POW 1
  #define HAVE_FACTORIAL 1
  #define HAVE_VARIABLES 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_HEXBIN 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1
  #define HAVE_DISPLAY_MODES 1

#elif defined(SCALC_MEDIUM)
  /* DOS 16-bit target - all features with reduced sizes (64KB data limit) */
  /* Core arithmetic */
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_HYPER 1
  #define HAVE_POW 1
  #define HAVE_FACTORIAL 1
  #define HAVE_VARIABLES 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_HEXBIN 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1
  #define HAVE_DISPLAY_MODES 1
  #define HAVE_COMPLEX 1
  #define HAVE_BITWISE 1
  /* #define HAVE_ROMAN 1 - not implemented */
  #define HAVE_LOOPS 1
  #define HAVE_SOLVER 1
  /* Phase 2: Scientific functions */
  #define HAVE_GAMMA 1
  #define HAVE_ERF 1
  /* These are not yet implemented:
  #define HAVE_BESSEL 1
  #define HAVE_ELLIPTIC 1
  #define HAVE_LAMBERTW 1
  */
  #define HAVE_STATS 1
  #define HAVE_DISTRIBUTIONS 1
  #define HAVE_MATRIX 1          /* 2x2 matrices, 2 variables */
  /* Phase 3: Professional */
  #define HAVE_TVM 1
  #define HAVE_NEWTON 1
  #define HAVE_ORBITAL 1
  /* Buffer limits for 64KB DGROUP */
  #ifndef MAX_INPUT
    #define MAX_INPUT 128
  #endif
  #ifndef MAX_HISTORY
    #define MAX_HISTORY 10
  #endif
  #ifndef MAX_VARIABLES
    #define MAX_VARIABLES 26
  #endif
  #ifndef MAX_FUNCTIONS
    #define MAX_FUNCTIONS 10
  #endif
  #ifndef MAX_MATRIX_VARS
    #define MAX_MATRIX_VARS 4    /* 4 matrices with __far memory */
  #endif

#else
  /* Full build (default) - Linux/modern - ALL PHASES */
  #define SCALC_FULL 1
  
  /* Phase 1: Core foundation */
  #define HAVE_SQRT 1
  #define HAVE_EXP 1
  #define HAVE_CONST 1
  #define HAVE_TRIG 1
  #define HAVE_HYPER 1
  #define HAVE_POW 1
  #define HAVE_FACTORIAL 1
  #define HAVE_VARIABLES 1
  #define HAVE_USER_FUNCS 1
  #define HAVE_HEXBIN 1
  #define HAVE_SCI_NOTATION 1
  #define HAVE_RPN 1
  #define HAVE_COMB 1
  #define HAVE_GCD 1
  #define HAVE_DISPLAY_MODES 1
  #define HAVE_COMPLEX 1
  /* #define HAVE_ROMAN 1 - not implemented */
  #define HAVE_LOOPS 1
  #define HAVE_SOLVER 1
  #define HAVE_BITWISE 1
  
  /* Phase 2: Scientific essentials */
  #define HAVE_GAMMA 1       /* gamma, lgamma, factorial extension */
  #define HAVE_ERF 1         /* erf, erfc */
  #define HAVE_DISTRIBUTIONS 1  /* normal, t, chi2, F, binomial, poisson */
  #define HAVE_BESSEL 1      /* Bessel functions J, Y, I, K */
  /* #define HAVE_ELLIPTIC 1 - not implemented */
  #define HAVE_MATRIX 1
  #define HAVE_STATS 1
  
  /* Phase 3: Professional features */
  #define HAVE_UNITS 1       /* Physical units with dimensional analysis */
  #define HAVE_POLYROOTS 1   /* Polynomial root finding */
  #define HAVE_INTEGRATE 1   /* Gauss-Legendre numerical integration */
  #define HAVE_ODE 1         /* ODE solver (RK4) */
  #define HAVE_PLOT 1
  #define HAVE_NEWTON 1
  #define HAVE_TVM 1
  #define HAVE_ORBITAL 1
  
  /* Phase 4: Advanced */
  /* #define HAVE_LAMBERTW 1 - not implemented */
  #define HAVE_EIGENVALUES 1 /* Eigenvalue computation */
  #define HAVE_QR 1          /* QR decomposition */
  #define HAVE_SYMBOLIC 1    /* Limited symbolic differentiation */
  #define HAVE_NAMED_VARS 1  /* Multi-character variable names (MATLAB compat) */
  #define HAVE_OPTIM 1       /* Optimization: quadprog, fmincon */
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
#if defined(HAVE_SOLVER) && !defined(HAVE_SQRT)
  #define HAVE_SQRT 1
#endif
#if defined(HAVE_SOLVER) && !defined(HAVE_COMPLEX)
  #define HAVE_COMPLEX 1
#endif

/* Phase 2-4 dependencies */
#if defined(HAVE_GAMMA) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_ERF) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_DISTRIBUTIONS) && !defined(HAVE_ERF)
  #define HAVE_ERF 1
#endif
#if defined(HAVE_DISTRIBUTIONS) && !defined(HAVE_GAMMA)
  #define HAVE_GAMMA 1
#endif
#if defined(HAVE_BESSEL) && !defined(HAVE_TRIG)
  #define HAVE_TRIG 1
#endif
#if defined(HAVE_ELLIPTIC) && !defined(HAVE_SQRT)
  #define HAVE_SQRT 1
#endif
#if defined(HAVE_LAMBERTW) && !defined(HAVE_EXP)
  #define HAVE_EXP 1
#endif
#if defined(HAVE_EIGENVALUES) && !defined(HAVE_MATRIX)
  #define HAVE_MATRIX 1
#endif
#if defined(HAVE_QR) && !defined(HAVE_MATRIX)
  #define HAVE_MATRIX 1
#endif
#if defined(HAVE_INTEGRATE) && !defined(HAVE_USER_FUNCS)
  #define HAVE_USER_FUNCS 1
#endif
#if defined(HAVE_ODE) && !defined(HAVE_USER_FUNCS)
  #define HAVE_USER_FUNCS 1
#endif

/* ========== Buffer Sizes ========== */

#ifndef MAX_INPUT
  #if defined(PLATFORM_8BIT)
    #define MAX_INPUT 40
  #elif defined(PLATFORM_DOS)
    #define MAX_INPUT 128
  #else
    #define MAX_INPUT 256
  #endif
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
    #if defined(PLATFORM_8BIT)
      #define MAX_MATRIX_VARS 2
    #elif defined(PLATFORM_DOS)
      #define MAX_MATRIX_VARS 4
    #else
      #define MAX_MATRIX_VARS 8
    #endif
  #endif
#else
  #define MAX_MATRIX_VARS 1  /* Minimum 1 for C89 compliance */
#endif

/* ========== RPN Stack Configuration ========== */
#ifdef HAVE_RPN
  #ifndef RPN_STACK_SIZE
    #if defined(PLATFORM_8BIT)
      #define RPN_STACK_SIZE 4   /* Classic HP 4-level */
    #else
      #define RPN_STACK_SIZE 16  /* Extended stack */
    #endif
  #endif
#endif

/* ========== Display Mode Configuration ========== */
#ifdef HAVE_DISPLAY_MODES
  #define DISPLAY_MODE_ALL 0   /* Show all significant digits */
  #define DISPLAY_MODE_FIX 1   /* Fixed decimal places */
  #define DISPLAY_MODE_SCI 2   /* Scientific notation */
  #define DISPLAY_MODE_ENG 3   /* Engineering notation (exp % 3) */
#endif

/* ========== Stats Data Limits ========== */
#ifdef HAVE_STATS
  #ifndef STATS_MAX_DATA
    #if defined(PLATFORM_8BIT)
      #define STATS_MAX_DATA 50
    #elif defined(PLATFORM_DOS)
      #define STATS_MAX_DATA 500
    #else
      #define STATS_MAX_DATA 1000
    #endif
  #endif
#endif

/* ========== Unit System Configuration ========== */
#ifdef HAVE_UNITS
  #ifndef MAX_UNITS
    #define MAX_UNITS 150
  #endif
  #ifndef MAX_UNIT_EXPR
    #define MAX_UNIT_EXPR 32
  #endif
#endif

/* ========== Integration Configuration ========== */
#ifdef HAVE_INTEGRATE
  /* Gauss-Legendre nodes (5, 10, or 15 point) */
  #ifndef GAUSS_LEGENDRE_POINTS
    #if defined(PLATFORM_8BIT)
      #define GAUSS_LEGENDRE_POINTS 5
    #else
      #define GAUSS_LEGENDRE_POINTS 15
    #endif
  #endif
#endif

#endif /* CONFIG_H */
