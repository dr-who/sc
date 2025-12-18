/*
 * func_registry.h - Unified function registry with documentation and dispatch
 * 
 * This provides a single source of truth for all built-in functions.
 * The parser looks up functions here, and the help system uses the same data.
 *
 * C89 compliant
 */
#ifndef FUNC_REGISTRY_H
#define FUNC_REGISTRY_H

#include <stddef.h>

/* ========== Handler IDs ========== */
/* Each function has a unique handler ID that the parser uses in a switch */

typedef enum {
    /* Invalid/not found */
    FH_NONE = 0,
    
    /* Constants */
    FH_CONST_PI,
    FH_CONST_E,
    FH_CONST_I,
    FH_CONST_INF,
    FH_CONST_NAN,
    FH_CONST_TRUE,
    FH_CONST_FALSE,
    FH_CONST_EPS,
    FH_CONST_REALMAX,
    FH_CONST_REALMIN,
    FH_CONST_PHI,
    FH_CONST_TAU,
    FH_CONST_EPSP,       /* Native epsilon for current precision */
    FH_CONST_PRECISION,  /* Current precision in bits */
    FH_SILVER,           /* Silver ratio */
    
    /* Trigonometric - basic */
    FH_SIN, FH_COS, FH_TAN,
    FH_ASIN, FH_ACOS, FH_ATAN, FH_ATAN2,
    FH_SINH, FH_COSH, FH_TANH,
    FH_ASINH, FH_ACOSH, FH_ATANH,
    
    /* Trigonometric - reciprocal */
    FH_SEC, FH_CSC, FH_COT,
    FH_ASEC, FH_ACSC, FH_ACOT,
    FH_SECH, FH_CSCH, FH_COTH,
    FH_ASECH, FH_ACSCH, FH_ACOTH,
    
    /* Trigonometric - degree variants */
    FH_SIND, FH_COSD, FH_TAND,
    FH_ASIND, FH_ACOSD, FH_ATAND, FH_ATAN2D,
    FH_SECD, FH_CSCD, FH_COTD,
    FH_ASECD, FH_ACSCD, FH_ACOTD,
    
    /* Trigonometric - special */
    FH_SINC, FH_SINPI, FH_COSPI,
    
    /* Exponential and powers */
    FH_EXP, FH_EXP2, FH_EXPM1,
    FH_POW, FH_POW2,
    FH_SQRT, FH_CBRT, FH_NTHROOT,
    FH_HYPOT,
    FH_REALSQRT, FH_REALPOW,
    
    /* Logarithmic */
    FH_LOG, FH_LN, FH_LOG10, FH_LOG2,
    FH_LOG1P, FH_LOGB,
    FH_REALLOG,
    
    /* Complex number operations */
    FH_REAL, FH_IMAG, FH_CONJ,
    FH_ABS, FH_ARG, FH_ANGLE, FH_PHASE,
    FH_COMPLEX,
    FH_CART2POL, FH_POL2CART,
    
    /* Rounding and sign */
    FH_FLOOR, FH_CEIL, FH_ROUND,
    FH_TRUNC, FH_FIX, FH_FRAC,
    FH_SIGN, FH_SIGNUM, FH_NEG, FH_COPYSIGN,
    FH_MOD, FH_REM, FH_FMOD, FH_FDIM,
    FH_CLAMP, FH_CLIP,
    FH_EXP10, FH_POW10, FH_FABS, FH_CPOW,
    FH_ISEVEN, FH_ISODD, FH_MAX2, FH_MIN2,
    FH_MODINV, FH_NORM1, FH_NORMINF, FH_CTRANSPOSE,
    FH_CUBE_NUM, FH_SQUARE_NUM, FH_MEAN_COLS, FH_STD_COLS,
    FH_SUM_COLS, FH_SUM_ROWS, FH_WEEKS, FH_XNOR,
    
    /* Number theory */
    FH_GCD, FH_LCM,
    FH_ISPRIME, FH_NEXTPRIME, FH_PREVPRIME,
    FH_PRIMEPI, FH_NTHPRIME,
    FH_FACTOR, FH_DIVISORS, FH_DIVISORSUM,
    FH_TOTIENT, FH_MOBIUS,
    FH_FACT, FH_FACTORIAL, FH_FACTORIAL2,
    FH_NCR, FH_NPR, FH_COMB, FH_PERM,
    FH_FIBONACCI, FH_LUCAS, FH_CATALAN,
    FH_EVEN, FH_ODD,
    FH_NDIGITS, FH_DIGITSUM, FH_DIGITROOT,
    FH_ISPERFECT, FH_ISABUNDANT, FH_ISDEFICIENT,
    FH_ISSQUAREFREE,
    
    /* Special functions */
    FH_GAMMA, FH_LGAMMA, FH_TGAMMA,
    FH_BETA, FH_DIGAMMA, FH_PSI,
    FH_ERF, FH_ERFC, FH_ERFINV,
    FH_BESSELJ, FH_BESSELY,
    FH_ZETA, FH_HARMONIC,
    FH_SIGMOID, FH_SOFTPLUS,
    FH_HEAVISIDE, FH_STEP, FH_RECT, FH_TRI,
    FH_POCHHAMMER,
    
    /* Statistics - basic */
    FH_SUM, FH_PROD, FH_MEAN, FH_MEDIAN, FH_MODE,
    FH_STD, FH_VAR, FH_MIN, FH_MAX, FH_RANGE,
    FH_GEOMEAN, FH_HARMMEAN,
    FH_SKEWNESS, FH_KURTOSIS,
    FH_MAD, FH_IQR, FH_RMS, FH_SUMSQ, FH_MEANSQ,
    FH_PRCTILE, FH_ZSCORE,
    FH_COV, FH_CORRCOEF,
    
    /* Probability distributions */
    FH_NORMPDF, FH_NORMCDF, FH_NORMINV,
    FH_BINOPDF, FH_BINOCDF, FH_UNIFPDF, FH_UNIFCDF,
    FH_FPDF, FH_FCDF, FH_FINV,
    FH_TCDF, FH_TINV,
    FH_CHI2CDF, FH_CHI2INV,
    FH_POISSCDF, FH_POISSPDF,
    
    /* Random numbers */
    FH_RAND, FH_RANDN, FH_RANDI, FH_RANDPERM,
    
    /* Linear algebra */
    FH_DET, FH_INV, FH_TRACE, FH_TRANS, FH_TRANSPOSE,
    FH_RANK, FH_COND, FH_NORM,
    FH_LINSOLVE, FH_MLDIVIDE,
    FH_CHOL, FH_QR, FH_SVD, FH_LU, FH_EIG,
    FH_NULL, FH_SCHUR, FH_PINV,
    FH_RREF, FH_ORTH, FH_POLY, FH_EXPM, FH_LOGM,
    FH_HESS, FH_BALANCE,
    
    /* Signal processing */
    FH_FFT, FH_IFFT, FH_XCORR, FH_XCOV,
    
    /* Data analysis */
    FH_MESHGRID, FH_HISTCOUNTS, FH_ISOUTLIER,
    
    /* Matrix creation */
    FH_ZEROS, FH_ONES, FH_EYE, FH_DIAG,
    FH_LINSPACE, FH_LOGSPACE,
    FH_MAGIC, FH_PASCAL, FH_HILB, FH_INVHILB,
    FH_VANDER, FH_TOEPLITZ, FH_HANKEL, FH_COMPAN,
    FH_RESHAPE, FH_REPMAT,
    FH_ROT90, FH_FLIP, FH_FLIPLR, FH_FLIPUD,
    FH_TRIU, FH_TRIL, FH_BLKDIAG,
    FH_CAT, FH_HORZCAT, FH_VERTCAT,
    FH_CIRCSHIFT, FH_SORT, FH_SORTROWS, FH_UNIQUE,
    FH_KRON, FH_DOT, FH_CROSS,
    
    /* Matrix query */
    FH_SIZE, FH_LENGTH, FH_ROWS, FH_COLS,
    FH_NUMEL, FH_NDIMS,
    FH_ISEMPTY, FH_ISSCALAR, FH_ISVECTOR,
    FH_ISROW, FH_ISCOLUMN, FH_ISMATRIX, FH_ISSQUARE,
    FH_ISSORTED, FH_ISSYMMETRIC,
    FH_NNZ, FH_FIND,
    
    /* Signal processing */
    FH_CONV, FH_DECONV,
    FH_DIFF, FH_GRADIENT,
    FH_CUMSUM, FH_CUMPROD, FH_CUMMIN, FH_CUMMAX,
    FH_MOVMEAN, FH_MOVSUM, FH_MOVMAX, FH_MOVMIN, FH_MOVSTD,
    FH_INTERP1, FH_TRAPZ,
    
    /* Polynomials */
    FH_POLYVAL, FH_POLYDER, FH_POLYINT,
    FH_ROOTS, FH_ROOTS2, FH_QUADROOTS,
    
    /* Angle conversion */
    FH_DEG2RAD, FH_RAD2DEG,
    FH_WRAPTOPI, FH_WRAPTO2PI,
    FH_WRAP180, FH_WRAP360,
    FH_WRAPTO180, FH_WRAPTO360,
    
    /* Set operations */
    FH_UNION, FH_INTERSECT, FH_SETDIFF, FH_SETXOR,
    FH_ISMEMBER, FH_UNIQUE_SET,
    
    /* Logical */
    FH_ANY, FH_ALL,
    FH_EQ, FH_NE, FH_LT, FH_LE, FH_GT, FH_GE,
    FH_AND, FH_OR, FH_NOT, FH_XOR,
    FH_ISNAN, FH_ISINF, FH_ISFINITE, FH_ISREAL,
    FH_NAND, FH_NOR, FH_IMPLIES,
    FH_APPROXEQ, FH_ISEQUAL,
    FH_ISINTEGER, FH_ISNUMERIC, FH_ISLOGICAL,
    FH_ISNEGATIVE, FH_ISPOSITIVE,
    
    /* Bitwise */
    FH_BITAND, FH_BITOR, FH_BITXOR, FH_BITNOT,
    FH_BITSHIFT, FH_SHL, FH_SHR,
    FH_BITGET, FH_BITSET,
    
    /* Coordinates */
    FH_CART2SPH, FH_SPH2CART,
    FH_CART2CYL, FH_CYL2CART,
    
    /* Utility */
    FH_NEXTPOW2, FH_NORMALIZE, FH_RESCALE,
    FH_DISP, FH_FORMAT,
    FH_LERP,
    FH_TIC, FH_TOC, FH_VER, FH_CLC, FH_CLEAR,
    FH_HELP, FH_DEMO, FH_BENCH, FH_FUNCS, FH_VARS,
    FH_QUIT, FH_EXIT,
    FH_FORMAT_SCI, FH_FORMAT_ENG, FH_FORMAT_HEX, FH_FORMAT_BIN, FH_FORMAT_OCT,
    FH_FORMAT_SHORT, FH_FORMAT_LONG, FH_FORMAT_DIGITS,
    FH_ANGLE_RAD, FH_ANGLE_DEG, FH_ANGLE_GRAD,
    FH_SAVE, FH_LOAD, FH_PLOT, FH_SOLVE, FH_CONST_CMD,
    FH_PREC, FH_BASE, FH_FLOAT, FH_ASCII, FH_SET, FH_UNSET,
    FH_QUAD, FH_INTEGRATE, FH_ORBIT, FH_TLE, FH_RPN,
    FH_ENG1000, FH_ENG24,
    
    /* Data science */
    FH_PCA, FH_PCAREDUCE,
    FH_KMEANS, FH_SILHOUETTE, FH_PDIST,
    FH_MANHATTAN,
    
    /* Text analysis */
    FH_BOW, FH_TFIDF, FH_WORDVEC,
    
    /* Number theory - additional */
    FH_BELL, FH_STIRLING2, FH_DERANGEMENTS,
    FH_OMEGA, FH_BIGOMEGA,
    FH_FALLING, FH_PARTITION,
    FH_GENHARMONIC,
    
    /* Datetime */
    FH_DATENUM, FH_DATETIME, FH_DATESHIFT,
    FH_YEAR, FH_MONTH, FH_DAY,
    FH_HOUR, FH_MINUTE, FH_SECOND,
    FH_DAYS, FH_HOURS, FH_MINUTES, FH_SECONDS,
    FH_STARTOFMONTH, FH_STARTOFYEAR, FH_STARTOFDAY,
    
    /* Tables */
    FH_TIMETABLE, FH_TABLE, FH_HEIGHT, FH_WIDTH,
    FH_GROUPSUMMARY, FH_GROUPCOUNTS,
    FH_SORTROWS_TBL, FH_HEAD, FH_TAIL,
    
    /* Forecasting */
    FH_MOVAVG, FH_EWMA, FH_FORECAST,
    FH_TREND, FH_DETREND, FH_SEASON,
    FH_DIFF_TS, FH_LAG, FH_LEAD,
    FH_AUTOCORR, FH_CROSSCORR,
    FH_ARIMA, FH_EXPSMOOTH,
    
    /* Additional math */
    FH_TRIANGULAR, FH_PENTAGONAL, FH_HEXAGONAL,
    FH_ISTRIANGULAR, FH_ISPERFECTSQUARE, FH_ISPERFECTPOWER, FH_ISPOWER,
    FH_REVERSEDIGITS, FH_ISPALINDROME,
    FH_SUBFACTORIAL, FH_STIRLING1,
    
    /* Data cleaning */
    FH_FILLMISSING, FH_RMMISSING, FH_SMOOTHDATA,
    FH_ISNAN_VEC, FH_ANYNAN, FH_ALLNAN,
    
    /* Financial */
    FH_NPV, FH_IRR, FH_PAYBACK, FH_CAGR, FH_COMPOUND, FH_ROI,
    FH_SHARPE, FH_DRAWDOWN, FH_CUMRET, FH_WINSORIZE,
    
    /* Signal analysis */
    FH_PEAKS, FH_VALLEYS, FH_CROSSOVER,
    FH_STANDARDIZE, FH_MINMAX_SCALE,
    FH_BOUNDS, FH_MINMAX,
    
    /* Curve fitting and integration */
    FH_POLYFIT, FH_CUMTRAPZ, FH_DIFF2, FH_CENTER,
    FH_LINREG, FH_QUANTILE,
    
    /* SaaS Metrics */
    FH_MRR, FH_ARPU, FH_MRRBRIDGE,
    FH_CUSTOMERCOUNT, FH_NEWCUSTOMERS, FH_CHURN,
    FH_REACTIVATED, FH_CHURNRATE,
    FH_NRR, FH_GRR, FH_RETENTION, FH_TENURE,
    
    /* Sentinel */
    FH_MAX_HANDLER
} FuncHandler;

/* ========== Function Flags ========== */

#define FF_NONE         0x0000
#define FF_CONST        0x0001  /* Is a constant (no parens) */
#define FF_SCALAR       0x0002  /* Operates on scalars */
#define FF_MATRIX       0x0004  /* Operates on matrices */
#define FF_VARIADIC     0x0008  /* Variable number of arguments */
#define FF_ANGLE        0x0010  /* Affected by angle mode */
#define FF_DEGREE       0x0020  /* Degree variant (sind, cosd, etc.) */
#define FF_COMPLEX      0x0040  /* Handles complex numbers */
#define FF_REAL_ONLY    0x0080  /* Requires real arguments */
#define FF_NO_PARENS    0x0100  /* Can be called without parens */
#define FF_DEPRECATED   0x0200  /* Deprecated function */

/* ========== Function Registry Entry ========== */

typedef struct {
    const char *name;           /* Function name */
    FuncHandler handler;        /* Handler ID for dispatch */
    int min_args;               /* Minimum arguments (-1 = constant) */
    int max_args;               /* Maximum arguments (-1 = unlimited) */
    unsigned int flags;         /* Function flags */
    const char *category;       /* Category for grouping */
    const char *syntax;         /* Usage syntax */
    const char *description;    /* Full description */
    const char *examples[6];    /* Example expressions */
    const char *see_also;       /* Related functions */
} FuncEntry;

/* ========== API Functions ========== */

/* Initialize function registry (call once at startup) */
void func_registry_init(void);

/* Look up function by name - returns NULL if not found */
const FuncEntry *func_lookup(const char *name);

/* Check if name is a known function */
int func_exists(const char *name);

/* Get handler ID for function (FH_NONE if not found) */
FuncHandler func_get_handler(const char *name);

/* Get total number of registered functions */
int func_count(void);

/* Iterate through all functions */
const FuncEntry *func_get_all(int *count);

/* Show help for a function */
void func_help(const char *name);

/* Run demo for a function */
void func_demo(const char *name);

/* List all functions by category */
void func_list(void);

#endif /* FUNC_REGISTRY_H */
