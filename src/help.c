/*
 * help.c - Function help, demo, and benchmark system
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "help.h"
#include "config.h"

/* Forward declarations */
extern const char *input_ptr;
extern void next_token(void);
extern int parse_expr(void *result);

/* Helper to compare strings case-insensitively */
static int str_eq(const char *a, const char *b) {
    while (*a && *b) {
        char ca = *a, cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
        a++; b++;
    }
    return *a == *b;
}

/* Help entry structure */
typedef struct {
    const char *name;
    const char *category;
    const char *syntax;
    const char *description;
    const char *examples[4];
} HelpEntry;

/* All function help entries */
static const HelpEntry help_entries[] = {
    /* ===== TRIGONOMETRY ===== */
    {"sin", "Trigonometry", "sin(x)", "Sine of x (x in radians)",
     {"sin(0) = 0", "sin(pi/2) = 1", "sin(pi/6) = 0.5", NULL}},
    {"cos", "Trigonometry", "cos(x)", "Cosine of x (x in radians)",
     {"cos(0) = 1", "cos(pi/3) = 0.5", "cos(pi) = -1", NULL}},
    {"tan", "Trigonometry", "tan(x)", "Tangent of x (x in radians)",
     {"tan(0) = 0", "tan(pi/4) = 1", NULL, NULL}},
    {"asin", "Trigonometry", "asin(x)", "Inverse sine, returns radians. Domain: [-1, 1]",
     {"asin(0) = 0", "asin(0.5) = pi/6", "asin(1) = pi/2", NULL}},
    {"acos", "Trigonometry", "acos(x)", "Inverse cosine, returns radians. Domain: [-1, 1]",
     {"acos(1) = 0", "acos(0.5) = pi/3", "acos(0) = pi/2", NULL}},
    {"atan", "Trigonometry", "atan(x)", "Inverse tangent, returns radians",
     {"atan(0) = 0", "atan(1) = pi/4", NULL, NULL}},
    {"atan2", "Trigonometry", "atan2(y, x)", "Two-argument arctangent, returns angle in radians",
     {"atan2(1, 1) = pi/4", "atan2(1, 0) = pi/2", "atan2(-1, -1) = -3*pi/4", NULL}},
    {"sinh", "Trigonometry", "sinh(x)", "Hyperbolic sine",
     {"sinh(0) = 0", "sinh(1) = 1.175", NULL, NULL}},
    {"cosh", "Trigonometry", "cosh(x)", "Hyperbolic cosine",
     {"cosh(0) = 1", "cosh(1) = 1.543", NULL, NULL}},
    {"tanh", "Trigonometry", "tanh(x)", "Hyperbolic tangent",
     {"tanh(0) = 0", "tanh(1) = 0.762", NULL, NULL}},
    {"asinh", "Trigonometry", "asinh(x)", "Inverse hyperbolic sine",
     {"asinh(0) = 0", "asinh(1) = 0.881", NULL, NULL}},
    {"acosh", "Trigonometry", "acosh(x)", "Inverse hyperbolic cosine. Domain: [1, inf)",
     {"acosh(1) = 0", "acosh(2) = 1.317", NULL, NULL}},
    {"atanh", "Trigonometry", "atanh(x)", "Inverse hyperbolic tangent. Domain: (-1, 1)",
     {"atanh(0) = 0", "atanh(0.5) = 0.549", NULL, NULL}},
    {"sec", "Trigonometry", "sec(x)", "Secant: 1/cos(x)",
     {"sec(0) = 1", "sec(pi/3) = 2", NULL, NULL}},
    {"csc", "Trigonometry", "csc(x)", "Cosecant: 1/sin(x)",
     {"csc(pi/2) = 1", "csc(pi/6) = 2", NULL, NULL}},
    {"cot", "Trigonometry", "cot(x)", "Cotangent: 1/tan(x)",
     {"cot(pi/4) = 1", "cot(pi/2) = 0", NULL, NULL}},
    {"asec", "Trigonometry", "asec(x)", "Inverse secant. Domain: |x| >= 1",
     {"asec(1) = 0", "asec(2) = pi/3", NULL, NULL}},
    {"acsc", "Trigonometry", "acsc(x)", "Inverse cosecant. Domain: |x| >= 1",
     {"acsc(1) = pi/2", "acsc(2) = pi/6", NULL, NULL}},
    {"acot", "Trigonometry", "acot(x)", "Inverse cotangent",
     {"acot(1) = pi/4", "acot(0) = pi/2", NULL, NULL}},
    {"sech", "Trigonometry", "sech(x)", "Hyperbolic secant: 1/cosh(x)",
     {"sech(0) = 1", "sech(1) = 0.648", NULL, NULL}},
    {"csch", "Trigonometry", "csch(x)", "Hyperbolic cosecant: 1/sinh(x)",
     {"csch(1) = 0.851", NULL, NULL, NULL}},
    {"coth", "Trigonometry", "coth(x)", "Hyperbolic cotangent: 1/tanh(x)",
     {"coth(1) = 1.313", NULL, NULL, NULL}},
    {"asech", "Trigonometry", "asech(x)", "Inverse hyperbolic secant. Domain: (0, 1]",
     {"asech(1) = 0", "asech(0.5) = 1.317", NULL, NULL}},
    {"acsch", "Trigonometry", "acsch(x)", "Inverse hyperbolic cosecant",
     {"acsch(1) = 0.881", NULL, NULL, NULL}},
    {"acoth", "Trigonometry", "acoth(x)", "Inverse hyperbolic cotangent. Domain: |x| > 1",
     {"acoth(2) = 0.549", NULL, NULL, NULL}},
    {"sind", "Trigonometry", "sind(x)", "Sine of x degrees",
     {"sind(0) = 0", "sind(30) = 0.5", "sind(90) = 1", NULL}},
    {"cosd", "Trigonometry", "cosd(x)", "Cosine of x degrees",
     {"cosd(0) = 1", "cosd(60) = 0.5", "cosd(90) = 0", NULL}},
    {"tand", "Trigonometry", "tand(x)", "Tangent of x degrees",
     {"tand(0) = 0", "tand(45) = 1", NULL, NULL}},
    {"asind", "Trigonometry", "asind(x)", "Inverse sine, returns degrees",
     {"asind(0) = 0", "asind(0.5) = 30", "asind(1) = 90", NULL}},
    {"acosd", "Trigonometry", "acosd(x)", "Inverse cosine, returns degrees",
     {"acosd(1) = 0", "acosd(0.5) = 60", "acosd(0) = 90", NULL}},
    {"atand", "Trigonometry", "atand(x)", "Inverse tangent, returns degrees",
     {"atand(0) = 0", "atand(1) = 45", NULL, NULL}},
    {"atan2d", "Trigonometry", "atan2d(y, x)", "Two-argument arctangent, returns degrees",
     {"atan2d(1, 1) = 45", "atan2d(1, 0) = 90", NULL, NULL}},
    {"sinc", "Trigonometry", "sinc(x)", "Normalized sinc: sin(pi*x)/(pi*x), sinc(0)=1",
     {"sinc(0) = 1", "sinc(1) = 0", "sinc(0.5) = 0.637", NULL}},

    /* ===== EXPONENTIAL & LOGARITHMIC ===== */
    {"exp", "Exponential", "exp(x)", "Exponential function e^x",
     {"exp(0) = 1", "exp(1) = e = 2.718...", "exp(ln(5)) = 5", NULL}},
    {"ln", "Exponential", "ln(x)", "Natural logarithm (base e)",
     {"ln(1) = 0", "ln(e) = 1", "ln(exp(5)) = 5", NULL}},
    {"log", "Exponential", "log(x)", "Natural logarithm (alias for ln)",
     {"log(1) = 0", "log(e) = 1", NULL, NULL}},
    {"log10", "Exponential", "log10(x)", "Base-10 logarithm",
     {"log10(1) = 0", "log10(10) = 1", "log10(1000) = 3", NULL}},
    {"log2", "Exponential", "log2(x)", "Base-2 logarithm",
     {"log2(1) = 0", "log2(8) = 3", "log2(1024) = 10", NULL}},
    {"exp2", "Exponential", "exp2(x)", "Base-2 exponential: 2^x",
     {"exp2(0) = 1", "exp2(10) = 1024", NULL, NULL}},
    {"pow", "Exponential", "x^y or pow(x,y)", "x raised to power y",
     {"2^10 = 1024", "3^4 = 81", "(-1)^0.5 = i", NULL}},
    {"sqrt", "Exponential", "sqrt(x)", "Square root",
     {"sqrt(4) = 2", "sqrt(2) = 1.414...", "sqrt(-1) = i", NULL}},
    {"cbrt", "Exponential", "cbrt(x)", "Cube root",
     {"cbrt(8) = 2", "cbrt(27) = 3", "cbrt(-8) = -2", NULL}},
    {"nthroot", "Exponential", "nthroot(x, n)", "n-th root of x",
     {"nthroot(16, 4) = 2", "nthroot(32, 5) = 2", NULL, NULL}},
    {"log1p", "Exponential", "log1p(x)", "log(1+x), accurate for small x",
     {"log1p(0) = 0", "log1p(1) = ln(2)", NULL, NULL}},
    {"expm1", "Exponential", "expm1(x)", "exp(x)-1, accurate for small x",
     {"expm1(0) = 0", "expm1(1) = e-1", NULL, NULL}},
    {"reallog", "Exponential", "reallog(x)", "Real-valued logarithm (error for x<=0)",
     {"reallog(e) = 1", "reallog(10) = 2.303", NULL, NULL}},
    {"realsqrt", "Exponential", "realsqrt(x)", "Real-valued sqrt (error for x<0)",
     {"realsqrt(4) = 2", "realsqrt(2) = 1.414", NULL, NULL}},
    {"realpow", "Exponential", "realpow(x, y)", "Real-valued power",
     {"realpow(2, 3) = 8", NULL, NULL, NULL}},

    /* ===== COMPLEX NUMBERS ===== */
    {"re", "Complex", "re(z) or real(z)", "Real part of complex number",
     {"re(3+4i) = 3", "re(5) = 5", NULL, NULL}},
    {"real", "Complex", "real(z)", "Real part (alias for re)",
     {"real(3+4i) = 3", NULL, NULL, NULL}},
    {"im", "Complex", "im(z) or imag(z)", "Imaginary part of complex number",
     {"im(3+4i) = 4", "im(5) = 0", NULL, NULL}},
    {"imag", "Complex", "imag(z)", "Imaginary part (alias for im)",
     {"imag(3+4i) = 4", NULL, NULL, NULL}},
    {"conj", "Complex", "conj(z)", "Complex conjugate",
     {"conj(3+4i) = 3-4i", "conj(5) = 5", NULL, NULL}},
    {"abs", "Complex", "abs(x) or abs(z)", "Absolute value or complex magnitude",
     {"abs(-5) = 5", "abs(3+4i) = 5", NULL, NULL}},
    {"arg", "Complex", "arg(z)", "Argument (phase angle) in radians",
     {"arg(1+i) = pi/4", "arg(-1) = pi", NULL, NULL}},
    {"angle", "Complex", "angle(z)", "Phase angle in radians (alias for arg)",
     {"angle(1+i) = pi/4", NULL, NULL, NULL}},
    {"phase", "Complex", "phase(z)", "Phase angle (alias for arg)",
     {"phase(1+i) = pi/4", NULL, NULL, NULL}},
    {"complex", "Complex", "complex(a, b)", "Create complex number a + bi",
     {"complex(3, 4) = 3+4i", NULL, NULL, NULL}},
    {"cart2pol", "Complex", "cart2pol(x, y)", "Cartesian to polar: returns [r, theta]",
     {"cart2pol(3, 4) = [5, 0.927]", NULL, NULL, NULL}},
    {"pol2cart", "Complex", "pol2cart(r, theta)", "Polar to Cartesian: returns [x, y]",
     {"pol2cart(5, 0) = [5, 0]", NULL, NULL, NULL}},

    /* ===== ROUNDING & SIGN ===== */
    {"floor", "Rounding", "floor(x)", "Round toward negative infinity",
     {"floor(3.7) = 3", "floor(-3.7) = -4", NULL, NULL}},
    {"ceil", "Rounding", "ceil(x)", "Round toward positive infinity",
     {"ceil(3.2) = 4", "ceil(-3.2) = -3", NULL, NULL}},
    {"trunc", "Rounding", "trunc(x)", "Round toward zero",
     {"trunc(3.7) = 3", "trunc(-3.7) = -3", NULL, NULL}},
    {"round", "Rounding", "round(x)", "Round to nearest integer",
     {"round(3.5) = 4", "round(3.4) = 3", "round(-3.5) = -4", NULL}},
    {"fix", "Rounding", "fix(x)", "Round toward zero (alias for trunc)",
     {"fix(3.7) = 3", "fix(-3.7) = -3", NULL, NULL}},
    {"frac", "Rounding", "frac(x)", "Fractional part: x - trunc(x)",
     {"frac(3.7) = 0.7", "frac(-3.7) = -0.7", NULL, NULL}},
    {"sign", "Rounding", "sign(x)", "Sign function: -1, 0, or 1",
     {"sign(5) = 1", "sign(-5) = -1", "sign(0) = 0", NULL}},
    {"signum", "Rounding", "signum(x)", "Sign function (alias for sign)",
     {"signum(5) = 1", NULL, NULL, NULL}},
    {"mod", "Rounding", "mod(x, y)", "Modulo (floored division remainder)",
     {"mod(17, 5) = 2", "mod(-17, 5) = 3", NULL, NULL}},
    {"rem", "Rounding", "rem(x, y)", "Remainder (truncated division)",
     {"rem(17, 5) = 2", "rem(-17, 5) = -2", NULL, NULL}},

    /* ===== NUMBER THEORY ===== */
    {"gcd", "Number Theory", "gcd(a, b)", "Greatest common divisor",
     {"gcd(12, 18) = 6", "gcd(48, 18) = 6", NULL, NULL}},
    {"lcm", "Number Theory", "lcm(a, b)", "Least common multiple",
     {"lcm(4, 6) = 12", "lcm(12, 18) = 36", NULL, NULL}},
    {"isprime", "Number Theory", "isprime(n)", "Test if n is prime",
     {"isprime(7) = true", "isprime(8) = false", "isprime(2) = true", NULL}},
    {"divisors", "Number Theory", "divisors(n)", "List all divisors of n",
     {"divisors(12) = [1,2,3,4,6,12]", NULL, NULL, NULL}},
    {"totient", "Number Theory", "totient(n)", "Euler's totient function phi(n)",
     {"totient(10) = 4", "totient(12) = 4", NULL, NULL}},
    {"mobius", "Number Theory", "mobius(n)", "Mobius function mu(n)",
     {"mobius(1) = 1", "mobius(6) = 1", "mobius(4) = 0", NULL}},
    {"even", "Number Theory", "even(n)", "Test if n is even",
     {"even(4) = true", "even(5) = false", NULL, NULL}},
    {"odd", "Number Theory", "odd(n)", "Test if n is odd",
     {"odd(5) = true", "odd(4) = false", NULL, NULL}},
    {"ncr", "Number Theory", "ncr(n, r)", "Binomial coefficient C(n,r) = n!/(r!(n-r)!)",
     {"ncr(5, 2) = 10", "ncr(10, 3) = 120", NULL, NULL}},
    {"npr", "Number Theory", "npr(n, r)", "Permutations P(n,r) = n!/(n-r)!",
     {"npr(5, 2) = 20", "npr(10, 3) = 720", NULL, NULL}},
    {"fact", "Number Theory", "fact(n) or factorial(n)", "Factorial n!",
     {"fact(5) = 120", "fact(0) = 1", "fact(10) = 3628800", NULL}},
    {"factorial", "Number Theory", "factorial(n)", "Factorial n! (alias for fact)",
     {"factorial(5) = 120", NULL, NULL, NULL}},
    {"factorial2", "Number Theory", "factorial2(n)", "Double factorial n!!",
     {"factorial2(10) = 3840", "factorial2(9) = 945", NULL, NULL}},
    {"gamma", "Number Theory", "gamma(x)", "Gamma function: (x-1)! for integers",
     {"gamma(5) = 24", "gamma(0.5) = sqrt(pi)", NULL, NULL}},
    {"lgamma", "Number Theory", "lgamma(x)", "Log-gamma function: ln(gamma(x))",
     {"lgamma(5) = ln(24)", "lgamma(10) = 12.802", NULL, NULL}},
    {"fibonacci", "Number Theory", "fibonacci(n)", "n-th Fibonacci number",
     {"fibonacci(0) = 0", "fibonacci(10) = 55", "fibonacci(20) = 6765", NULL}},
    {"lucas", "Number Theory", "lucas(n)", "n-th Lucas number",
     {"lucas(0) = 2", "lucas(1) = 1", "lucas(10) = 123", NULL}},
    {"catalan", "Number Theory", "catalan(n)", "n-th Catalan number",
     {"catalan(0) = 1", "catalan(5) = 42", "catalan(10) = 16796", NULL}},

    /* ===== SPECIAL FUNCTIONS ===== */
    {"erf", "Special", "erf(x)", "Error function",
     {"erf(0) = 0", "erf(1) = 0.843", "erf(inf) = 1", NULL}},
    {"erfc", "Special", "erfc(x)", "Complementary error function: 1-erf(x)",
     {"erfc(0) = 1", "erfc(1) = 0.157", NULL, NULL}},
    {"sigmoid", "Special", "sigmoid(x)", "Logistic sigmoid: 1/(1+exp(-x))",
     {"sigmoid(0) = 0.5", "sigmoid(10) = 0.9999...", NULL, NULL}},
    {"softplus", "Special", "softplus(x)", "Softplus: ln(1+exp(x))",
     {"softplus(0) = ln(2)", "softplus(10) = 10", NULL, NULL}},
    {"step", "Special", "step(x) or heaviside(x)", "Unit step function",
     {"step(1) = 1", "step(-1) = 0", "step(0) = 0.5", NULL}},
    {"heaviside", "Special", "heaviside(x)", "Heaviside step function",
     {"heaviside(1) = 1", "heaviside(-1) = 0", NULL, NULL}},
    {"rect", "Special", "rect(x)", "Rectangle function: 1 if |x|<0.5",
     {"rect(0) = 1", "rect(0.6) = 0", NULL, NULL}},
    {"tri", "Special", "tri(x)", "Triangle function: max(0, 1-|x|)",
     {"tri(0) = 1", "tri(0.5) = 0.5", "tri(1) = 0", NULL}},
    {"beta", "Special", "beta(a, b)", "Beta function: gamma(a)*gamma(b)/gamma(a+b)",
     {"beta(2, 3) = 0.0833", "beta(1, 1) = 1", NULL, NULL}},
    {"besselj", "Special", "besselj(n, x)", "Bessel function of first kind J_n(x)",
     {"besselj(0, 1) = 0.765", "besselj(1, 1) = 0.440", NULL, NULL}},
    {"bessely", "Special", "bessely(n, x)", "Bessel function of second kind Y_n(x)",
     {"bessely(0, 1) = 0.088", "bessely(1, 1) = -0.781", NULL, NULL}},
    {"digamma", "Special", "digamma(x) or psi(x)", "Digamma function: d/dx ln(Gamma(x))",
     {"digamma(1) = -0.5772 (Euler-Mascheroni)", NULL, NULL, NULL}},
    {"zeta", "Special", "zeta(s)", "Riemann zeta function for s > 1",
     {"zeta(2) = pi^2/6 = 1.6449", "zeta(4) = pi^4/90", NULL, NULL}},
    {"harmonic", "Special", "harmonic(n)", "Harmonic number H_n = 1 + 1/2 + ... + 1/n",
     {"harmonic(1) = 1", "harmonic(10) = 2.9289", NULL, NULL}},
    
    /* ===== DISTRIBUTIONS ===== */
    {"tcdf", "Probability", "tcdf(t, df)", "Student's t cumulative distribution",
     {"tcdf(0, 10) = 0.5", "tcdf(1.812, 10) = 0.95", NULL, NULL}},
    {"chi2cdf", "Probability", "chi2cdf(x, df)", "Chi-squared cumulative distribution",
     {"chi2cdf(0, 5) = 0", "chi2cdf(11.07, 5) = 0.95", NULL, NULL}},

    /* ===== STATISTICS ===== */
    {"sum", "Statistics", "sum(v)", "Sum of vector elements",
     {"sum([1,2,3,4,5]) = 15", NULL, NULL, NULL}},
    {"mean", "Statistics", "mean(v)", "Arithmetic mean",
     {"mean([1,2,3,4,5]) = 3", NULL, NULL, NULL}},
    {"median", "Statistics", "median(v)", "Median value",
     {"median([1,2,3,4,5]) = 3", "median([1,2,3,4]) = 2.5", NULL, NULL}},
    {"std", "Statistics", "std(v)", "Sample standard deviation",
     {"std([1,2,3,4,5]) = 1.581", NULL, NULL, NULL}},
    {"var", "Statistics", "var(v)", "Sample variance",
     {"var([1,2,3,4,5]) = 2.5", NULL, NULL, NULL}},
    {"min", "Statistics", "min(v) or min(a,b,...)", "Minimum value",
     {"min([3,1,4,1,5]) = 1", "min(3, 7) -> use min([3,7])", NULL, NULL}},
    {"max", "Statistics", "max(v) or max(a,b,...)", "Maximum value",
     {"max([3,1,4,1,5]) = 5", "max(3, 7) -> use max([3,7])", NULL, NULL}},
    {"range", "Statistics", "range(v)", "Range: max - min",
     {"range([1,2,3,4,5]) = 4", NULL, NULL, NULL}},
    {"geomean", "Statistics", "geomean(v)", "Geometric mean",
     {"geomean([1,2,4,8]) = 2.828", NULL, NULL, NULL}},
    {"harmmean", "Statistics", "harmmean(v)", "Harmonic mean",
     {"harmmean([1,2,4]) = 1.714", NULL, NULL, NULL}},
    {"skewness", "Statistics", "skewness(v)", "Sample skewness",
     {"skewness([1,2,3,4,5]) = 0", NULL, NULL, NULL}},
    {"kurtosis", "Statistics", "kurtosis(v)", "Sample kurtosis",
     {"kurtosis([1,2,3,4,5])", NULL, NULL, NULL}},
    {"mad", "Statistics", "mad(v)", "Mean absolute deviation",
     {"mad([1,2,3,4,5]) = 1.2", NULL, NULL, NULL}},
    {"iqr", "Statistics", "iqr(v)", "Interquartile range: Q3-Q1",
     {"iqr([1,2,3,4,5,6,7,8,9,10])", NULL, NULL, NULL}},
    {"rms", "Statistics", "rms(v)", "Root mean square",
     {"rms([3,4]) = 3.536", NULL, NULL, NULL}},
    {"sumsq", "Statistics", "sumsq(v)", "Sum of squares",
     {"sumsq([1,2,3]) = 14", NULL, NULL, NULL}},
    {"meansq", "Statistics", "meansq(v)", "Mean of squares",
     {"meansq([1,2,3]) = 4.667", NULL, NULL, NULL}},
    {"prctile", "Statistics", "prctile(v, p)", "p-th percentile",
     {"prctile([1:10], 50) = median", NULL, NULL, NULL}},
    {"zscore", "Statistics", "zscore(v)", "Standardized z-scores",
     {"zscore([1,2,3,4,5])", NULL, NULL, NULL}},

    /* ===== PROBABILITY DISTRIBUTIONS ===== */
    {"normpdf", "Probability", "normpdf(x) or normpdf(x,mu,sigma)", "Normal PDF",
     {"normpdf(0) = 0.399", "normpdf(0, 0, 1) = 0.399", NULL, NULL}},
    {"normcdf", "Probability", "normcdf(x) or normcdf(x,mu,sigma)", "Normal CDF",
     {"normcdf(0) = 0.5", "normcdf(1.96) = 0.975", NULL, NULL}},
    {"norminv", "Probability", "norminv(p)", "Inverse normal CDF (quantile)",
     {"norminv(0.5) = 0", "norminv(0.975) = 1.96", NULL, NULL}},
    {"rand", "Probability", "rand or rand(n) or rand(m,n)", "Uniform random [0,1)",
     {"rand = random scalar", "rand(3) = 3x3 random matrix", NULL, NULL}},
    {"randn", "Probability", "randn or randn(n) or randn(m,n)", "Standard normal random",
     {"randn = random scalar N(0,1)", "randn(3) = 3x3 random matrix", NULL, NULL}},
    {"randi", "Probability", "randi(max) or randi(max,n)", "Random integer [1,max]",
     {"randi(6) = dice roll", "randi(10, 3) = 3x3 ints", NULL, NULL}},
    {"randperm", "Probability", "randperm(n)", "Random permutation of 1:n",
     {"randperm(5) = shuffled [1,2,3,4,5]", NULL, NULL, NULL}},

    /* ===== LINEAR ALGEBRA ===== */
    {"det", "Linear Algebra", "det(A)", "Matrix determinant",
     {"det([1,2;3,4]) = -2", "det(eye(3)) = 1", NULL, NULL}},
    {"inv", "Linear Algebra", "inv(A)", "Matrix inverse",
     {"inv([1,2;3,4])", "A * inv(A) = eye", NULL, NULL}},
    {"trace", "Linear Algebra", "trace(A)", "Sum of diagonal elements",
     {"trace([1,2;3,4]) = 5", "trace(eye(n)) = n", NULL, NULL}},
    {"trans", "Linear Algebra", "trans(A) or A'", "Matrix transpose",
     {"trans([1,2;3,4]) = [1,3;2,4]", "[1,2,3]' = column", NULL, NULL}},
    {"eig", "Linear Algebra", "eig(A)", "Eigenvalues of matrix",
     {"eig([1,2;3,4])", NULL, NULL, NULL}},
    {"rank", "Linear Algebra", "rank(A)", "Matrix rank",
     {"rank([1,2;3,4]) = 2", "rank([1,2;2,4]) = 1", NULL, NULL}},
    {"cond", "Linear Algebra", "cond(A)", "Condition number",
     {"cond([1,2;3,4])", NULL, NULL, NULL}},
    {"norm", "Linear Algebra", "norm(v) or norm(A)", "Vector/matrix norm",
     {"norm([3,4]) = 5", "norm([1,2;3,4])", NULL, NULL}},
    {"linsolve", "Linear Algebra", "linsolve(A, b)", "Solve Ax = b",
     {"linsolve([1,2;3,4], [5;6])", NULL, NULL, NULL}},
    {"mldivide", "Linear Algebra", "A \\\\ b or mldivide(A, b)", "Matrix left division",
     {"[1,2;3,4] \\\\ [5;6]", NULL, NULL, NULL}},
    {"chol", "Linear Algebra", "chol(A) or [L] = chol(A)", "Cholesky decomposition (A = L*L')",
     {"chol([4,2;2,5])", "For positive definite matrices", NULL, NULL}},
    {"qr", "Linear Algebra", "[Q,R] = qr(A)", "QR factorization (A = Q*R)",
     {"[Q,R] = qr([1,2;3,4])", "Q*R = A", NULL, NULL}},
    {"svd", "Linear Algebra", "[U,S,V] = svd(A)", "Singular value decomposition (A = U*S*V')",
     {"[U,S,V] = svd([1,2;3,4])", "U*S*V' = A", NULL, NULL}},
    {"lu", "Linear Algebra", "[L,U] = lu(A) or [L,U,P] = lu(A)", "LU factorization",
     {"[L,U] = lu([1,2;3,4])", "L*U = P*A", NULL, NULL}},
    {"eig", "Linear Algebra", "eig(A) or [V,D] = eig(A)", "Eigenvalues/eigenvectors",
     {"[V,D] = eig([1,2;2,1])", "A*V = V*D", NULL, NULL}},
    {"null", "Linear Algebra", "null(A)", "Null space basis vectors",
     {"null([1,2,3;2,4,6])", "Returns basis for Ax=0", NULL, NULL}},
    {"schur", "Linear Algebra", "[Q,T] = schur(A)", "Schur decomposition (A = Q*T*Q')",
     {"[Q,T] = schur([1,2;3,4])", "T upper triangular", NULL, NULL}},

    /* ===== MATRIX CREATION ===== */
    {"zeros", "Matrix", "zeros(n) or zeros(m,n)", "Matrix of zeros",
     {"zeros(3) = 3x3 zeros", "zeros(2,3) = 2x3 zeros", NULL, NULL}},
    {"ones", "Matrix", "ones(n) or ones(m,n)", "Matrix of ones",
     {"ones(3) = 3x3 ones", "ones(2,3) = 2x3 ones", NULL, NULL}},
    {"eye", "Matrix", "eye(n) or eye(m,n)", "Identity matrix",
     {"eye(3) = 3x3 identity", NULL, NULL, NULL}},
    {"diag", "Matrix", "diag(v) or diag(A)", "Create diagonal or extract diagonal",
     {"diag([1,2,3]) = 3x3 diagonal", "diag(A) = diagonal of A", NULL, NULL}},
    {"linspace", "Matrix", "linspace(a, b, n)", "n linearly spaced points from a to b",
     {"linspace(0, 10, 5) = [0,2.5,5,7.5,10]", NULL, NULL, NULL}},
    {"logspace", "Matrix", "logspace(a, b, n)", "n log-spaced points: 10^a to 10^b",
     {"logspace(0, 2, 3) = [1, 10, 100]", NULL, NULL, NULL}},
    {"magic", "Matrix", "magic(n)", "n x n magic square",
     {"magic(3) = 3x3 magic square", NULL, NULL, NULL}},
    {"pascal", "Matrix", "pascal(n)", "n x n Pascal matrix",
     {"pascal(4)", NULL, NULL, NULL}},
    {"hilb", "Matrix", "hilb(n)", "n x n Hilbert matrix",
     {"hilb(3) = ill-conditioned matrix", NULL, NULL, NULL}},
    {"vander", "Matrix", "vander(v)", "Vandermonde matrix",
     {"vander([1,2,3])", NULL, NULL, NULL}},
    {"toeplitz", "Matrix", "toeplitz(v)", "Symmetric Toeplitz matrix",
     {"toeplitz([1,2,3])", NULL, NULL, NULL}},
    {"hankel", "Matrix", "hankel(v)", "Hankel matrix",
     {"hankel([1,2,3])", NULL, NULL, NULL}},
    {"compan", "Matrix", "compan(p)", "Companion matrix of polynomial",
     {"compan([1,-6,11,-6])", NULL, NULL, NULL}},

    /* ===== MATRIX MANIPULATION ===== */
    {"reshape", "Matrix", "reshape(A, m, n)", "Reshape matrix to m x n",
     {"reshape([1:6], 2, 3)", NULL, NULL, NULL}},
    {"repmat", "Matrix", "repmat(A, m, n)", "Tile A into m x n copies",
     {"repmat([1,2], 2, 3)", NULL, NULL, NULL}},
    {"rot90", "Matrix", "rot90(A) or rot90(A, k)", "Rotate matrix 90 degrees k times",
     {"rot90([1,2;3,4])", NULL, NULL, NULL}},
    {"flip", "Matrix", "flip(v)", "Flip vector/matrix",
     {"flip([1,2,3]) = [3,2,1]", NULL, NULL, NULL}},
    {"fliplr", "Matrix", "fliplr(A)", "Flip matrix left-right",
     {"fliplr([1,2,3]) = [3,2,1]", NULL, NULL, NULL}},
    {"flipud", "Matrix", "flipud(A)", "Flip matrix up-down",
     {"flipud([1;2;3]) = [3;2;1]", NULL, NULL, NULL}},
    {"triu", "Matrix", "triu(A)", "Upper triangular part",
     {"triu([1,2;3,4]) = [1,2;0,4]", NULL, NULL, NULL}},
    {"tril", "Matrix", "tril(A)", "Lower triangular part",
     {"tril([1,2;3,4]) = [1,0;3,4]", NULL, NULL, NULL}},
    {"blkdiag", "Matrix", "blkdiag(A, B, ...)", "Block diagonal matrix",
     {"blkdiag([1,2;3,4], [5,6;7,8])", NULL, NULL, NULL}},
    {"cat", "Matrix", "cat(dim, A, B)", "Concatenate along dimension",
     {"cat(1, [1,2], [3,4])", NULL, NULL, NULL}},
    {"horzcat", "Matrix", "horzcat(A, B) or [A, B]", "Horizontal concatenation",
     {"horzcat([1;2], [3;4])", NULL, NULL, NULL}},
    {"vertcat", "Matrix", "vertcat(A, B) or [A; B]", "Vertical concatenation",
     {"vertcat([1,2], [3,4])", NULL, NULL, NULL}},
    {"circshift", "Matrix", "circshift(v, k)", "Circular shift by k positions",
     {"circshift([1,2,3,4], 1) = [4,1,2,3]", NULL, NULL, NULL}},
    {"sort", "Matrix", "sort(v)", "Sort elements ascending",
     {"sort([3,1,4,1,5]) = [1,1,3,4,5]", NULL, NULL, NULL}},
    {"sortrows", "Matrix", "sortrows(A)", "Sort rows of matrix",
     {"sortrows([3,1;1,4;2,2])", NULL, NULL, NULL}},
    {"unique", "Matrix", "unique(v)", "Unique elements sorted",
     {"unique([1,2,2,3,3,3]) = [1,2,3]", NULL, NULL, NULL}},

    /* ===== MATRIX QUERIES ===== */
    {"size", "Matrix Query", "size(A)", "Dimensions [rows, cols]",
     {"size([1,2,3;4,5,6]) = [2, 3]", NULL, NULL, NULL}},
    {"length", "Matrix Query", "length(v)", "Length of vector",
     {"length([1,2,3,4,5]) = 5", NULL, NULL, NULL}},
    {"rows", "Matrix Query", "rows(A)", "Number of rows",
     {"rows([1,2,3;4,5,6]) = 2", NULL, NULL, NULL}},
    {"cols", "Matrix Query", "cols(A)", "Number of columns",
     {"cols([1,2,3;4,5,6]) = 3", NULL, NULL, NULL}},
    {"numel", "Matrix Query", "numel(A)", "Number of elements",
     {"numel([1,2,3;4,5,6]) = 6", NULL, NULL, NULL}},
    {"ndims", "Matrix Query", "ndims(A)", "Number of dimensions",
     {"ndims([1,2,3;4,5,6]) = 2", NULL, NULL, NULL}},
    {"isempty", "Matrix Query", "isempty(A)", "Test if empty",
     {"isempty([]) = true", "isempty([1]) = false", NULL, NULL}},
    {"isscalar", "Matrix Query", "isscalar(x)", "Test if scalar (1x1)",
     {"isscalar(5) = true", "isscalar([1,2]) = false", NULL, NULL}},
    {"isvector", "Matrix Query", "isvector(A)", "Test if row or column vector",
     {"isvector([1,2,3]) = true", NULL, NULL, NULL}},
    {"ismatrix", "Matrix Query", "ismatrix(A)", "Test if 2D matrix",
     {"ismatrix([1,2;3,4]) = true", NULL, NULL, NULL}},
    {"issquare", "Matrix Query", "issquare(A)", "Test if square matrix",
     {"issquare([1,2;3,4]) = true", NULL, NULL, NULL}},
    {"isrow", "Matrix Query", "isrow(A)", "Test if row vector",
     {"isrow([1,2,3]) = true", NULL, NULL, NULL}},
    {"iscolumn", "Matrix Query", "iscolumn(A)", "Test if column vector",
     {"iscolumn([1;2;3]) = true", NULL, NULL, NULL}},
    {"issorted", "Matrix Query", "issorted(v)", "Test if sorted ascending",
     {"issorted([1,2,3]) = true", NULL, NULL, NULL}},
    {"issymmetric", "Matrix Query", "issymmetric(A)", "Test if symmetric",
     {"issymmetric([1,2;2,1]) = true", NULL, NULL, NULL}},

    /* ===== SIGNAL PROCESSING ===== */
    {"conv", "Signal", "conv(u, v)", "Convolution of vectors",
     {"conv([1,2,3], [1,1]) = [1,3,5,3]", NULL, NULL, NULL}},
    {"deconv", "Signal", "deconv(u, v)", "Deconvolution",
     {"deconv([1,3,5,3], [1,1])", NULL, NULL, NULL}},
    {"gradient", "Signal", "gradient(v)", "Numerical gradient",
     {"gradient([1,4,9,16,25])", NULL, NULL, NULL}},
    {"movmean", "Signal", "movmean(v, k)", "Moving average with window k",
     {"movmean([1,2,3,4,5], 3)", NULL, NULL, NULL}},
    {"movsum", "Signal", "movsum(v, k)", "Moving sum with window k",
     {"movsum([1,2,3,4,5], 3)", NULL, NULL, NULL}},
    {"diff", "Signal", "diff(v)", "Differences between consecutive elements",
     {"diff([1,4,9,16]) = [3,5,7]", NULL, NULL, NULL}},
    {"cumsum", "Signal", "cumsum(v)", "Cumulative sum",
     {"cumsum([1,2,3,4,5]) = [1,3,6,10,15]", NULL, NULL, NULL}},
    {"cumprod", "Signal", "cumprod(v)", "Cumulative product",
     {"cumprod([1,2,3,4,5]) = [1,2,6,24,120]", NULL, NULL, NULL}},
    {"interp1", "Signal", "interp1(x, y, xi)", "Linear interpolation",
     {"interp1([1,2,3], [1,4,9], 1.5) = 2.5", NULL, NULL, NULL}},
    {"trapz", "Signal", "trapz(y) or trapz(x, y)", "Trapezoidal integration",
     {"trapz([1,2,3,4,5])", NULL, NULL, NULL}},
    {"corrcoef", "Signal", "corrcoef(x, y)", "Correlation coefficient",
     {"corrcoef([1,2,3], [2,4,6]) = 1", NULL, NULL, NULL}},

    /* ===== POLYNOMIALS ===== */
    {"polyval", "Polynomial", "polyval(p, x)", "Evaluate polynomial p at x",
     {"polyval([1,2,3], 2) = 11", "p = [1,2,3] means x^2+2x+3", NULL, NULL}},
    {"polyder", "Polynomial", "polyder(p)", "Polynomial derivative",
     {"polyder([1,2,3]) = [2,2]", NULL, NULL, NULL}},
    {"polyint", "Polynomial", "polyint(p)", "Polynomial integral",
     {"polyint([2,2]) = [1,2,0]", NULL, NULL, NULL}},

    /* ===== ANGLE FUNCTIONS ===== */
    {"deg2rad", "Angle", "deg2rad(x)", "Convert degrees to radians",
     {"deg2rad(180) = pi", "deg2rad(90) = pi/2", NULL, NULL}},
    {"rad2deg", "Angle", "rad2deg(x)", "Convert radians to degrees",
     {"rad2deg(pi) = 180", "rad2deg(pi/2) = 90", NULL, NULL}},
    {"hypot", "Angle", "hypot(x, y)", "Hypotenuse: sqrt(x^2 + y^2)",
     {"hypot(3, 4) = 5", NULL, NULL, NULL}},
    {"wrapToPi", "Angle", "wrapToPi(x)", "Wrap angle to [-pi, pi]",
     {"wrapToPi(4) = 4-2*pi", NULL, NULL, NULL}},
    {"wrap180", "Angle", "wrap180(x)", "Wrap angle to [-180, 180]",
     {"wrap180(270) = -90", NULL, NULL, NULL}},
    {"wrap360", "Angle", "wrap360(x)", "Wrap angle to [0, 360)",
     {"wrap360(370) = 10", "wrap360(-10) = 350", NULL, NULL}},

    /* ===== SET OPERATIONS ===== */
    {"union", "Set", "union(A, B)", "Set union",
     {"union([1,2,3], [2,3,4]) = [1,2,3,4]", NULL, NULL, NULL}},
    {"intersect", "Set", "intersect(A, B)", "Set intersection",
     {"intersect([1,2,3], [2,3,4]) = [2,3]", NULL, NULL, NULL}},
    {"setdiff", "Set", "setdiff(A, B)", "Set difference: A - B",
     {"setdiff([1,2,3], [2,3,4]) = [1]", NULL, NULL, NULL}},
    {"setxor", "Set", "setxor(A, B)", "Symmetric difference",
     {"setxor([1,2,3], [2,3,4]) = [1,4]", NULL, NULL, NULL}},

    /* ===== LOGICAL ===== */
    {"any", "Logical", "any(v)", "True if any element is nonzero",
     {"any([0,0,1]) = true", "any([0,0,0]) = false", NULL, NULL}},
    {"all", "Logical", "all(v)", "True if all elements are nonzero",
     {"all([1,1,1]) = true", "all([1,0,1]) = false", NULL, NULL}},
    {"eq", "Logical", "eq(a, b)", "Element-wise equality",
     {"eq(5, 5) = 1", "5 == 5 = true", NULL, NULL}},
    {"ne", "Logical", "ne(a, b)", "Element-wise not-equal",
     {"ne(5, 6) = 1", "5 ~= 6 = true", NULL, NULL}},
    {"lt", "Logical", "lt(a, b)", "Less than",
     {"lt(5, 6) = 1", "5 < 6 = true", NULL, NULL}},
    {"le", "Logical", "le(a, b)", "Less than or equal",
     {"le(5, 5) = 1", "5 <= 5 = true", NULL, NULL}},
    {"gt", "Logical", "gt(a, b)", "Greater than",
     {"gt(6, 5) = 1", "6 > 5 = true", NULL, NULL}},
    {"ge", "Logical", "ge(a, b)", "Greater than or equal",
     {"ge(5, 5) = 1", "5 >= 5 = true", NULL, NULL}},
    {"isnan", "Logical", "isnan(x)", "Test if NaN",
     {"isnan(NaN) = true", "isnan(1) = false", NULL, NULL}},
    {"isinf", "Logical", "isinf(x)", "Test if infinite",
     {"isinf(Inf) = true", "isinf(1) = false", NULL, NULL}},
    {"isfinite", "Logical", "isfinite(x)", "Test if finite",
     {"isfinite(1) = true", "isfinite(Inf) = false", NULL, NULL}},
    {"isreal", "Logical", "isreal(x)", "Test if real (no imaginary part)",
     {"isreal(5) = true", "isreal(5+3i) = false", NULL, NULL}},

    /* ===== BITWISE ===== */
    {"bitand", "Bitwise", "bitand(a, b)", "Bitwise AND",
     {"bitand(15, 7) = 7", "bitand(255, 15) = 15", NULL, NULL}},
    {"bitor", "Bitwise", "bitor(a, b)", "Bitwise OR",
     {"bitor(8, 4) = 12", "bitor(240, 15) = 255", NULL, NULL}},
    {"bitxor", "Bitwise", "bitxor(a, b)", "Bitwise XOR",
     {"bitxor(255, 170) = 85", NULL, NULL, NULL}},
    {"bitshift", "Bitwise", "bitshift(a, k)", "Shift left (k>0) or right (k<0)",
     {"bitshift(1, 4) = 16", NULL, NULL, NULL}},
    {"shl", "Bitwise", "shl(a, k)", "Shift left by k bits",
     {"shl(1, 4) = 16", NULL, NULL, NULL}},
    {"shr", "Bitwise", "shr(a, k)", "Shift right by k bits",
     {"shr(16, 2) = 4", NULL, NULL, NULL}},
    {"bnot", "Bitwise", "bnot(a) or not(a)", "Bitwise NOT",
     {"bnot(0) = -1", NULL, NULL, NULL}},

    /* ===== COORDINATE TRANSFORMS ===== */
    {"cart2sph", "Coordinates", "cart2sph(x, y, z)", "Cartesian to spherical",
     {"cart2sph(1, 2, 3)", NULL, NULL, NULL}},
    {"sph2cart", "Coordinates", "sph2cart(az, el, r)", "Spherical to Cartesian",
     {"sph2cart(0, 0, 1)", NULL, NULL, NULL}},

    /* ===== UTILITY ===== */
    {"clamp", "Utility", "clamp(x, lo, hi)", "Clamp x to range [lo, hi]",
     {"clamp(5, 0, 10) = 5", "clamp(-5, 0, 10) = 0", "clamp(15, 0, 10) = 10", NULL}},
    {"rescale", "Utility", "rescale(v)", "Scale to [0, 1]",
     {"rescale([1,2,3,4,5])", NULL, NULL, NULL}},
    {"normalize", "Utility", "normalize(v)", "Normalize to unit length",
     {"normalize([3,4]) = [0.6, 0.8]", NULL, NULL, NULL}},
    {"find", "Utility", "find(v)", "Indices of nonzero elements",
     {"find([0,1,0,1,0]) = [2,4]", NULL, NULL, NULL}},
    {"nextpow2", "Utility", "nextpow2(n)", "Exponent of next power of 2",
     {"nextpow2(5) = 3", "nextpow2(9) = 4", NULL, NULL}},
    {"eps", "Utility", "eps", "Machine epsilon (smallest x where 1+x > 1)",
     {"eps = 2.22e-16 (double)", NULL, NULL, NULL}},
    {"realmax", "Utility", "realmax", "Largest finite floating-point number",
     {"realmax = 1.79e+308 (double)", NULL, NULL, NULL}},
    {"realmin", "Utility", "realmin", "Smallest positive normalized number",
     {"realmin = 2.22e-308 (double)", NULL, NULL, NULL}},

    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL}}
};

/* Find help entry by name */
static const HelpEntry *find_help(const char *name) {
    int i;
    for (i = 0; help_entries[i].name != NULL; i++) {
        if (str_eq(help_entries[i].name, name)) {
            return &help_entries[i];
        }
    }
    return NULL;
}

/* Show help for a specific function */
void show_function_help(const char *name) {
    const HelpEntry *h = find_help(name);
    int i;
    
    if (h == NULL) {
        printf("No help available for '%s'\n", name);
        printf("Type 'functions' to see all available functions.\n");
        return;
    }
    
    printf("\n%s - %s\n", h->name, h->category);
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("Syntax:  %s\n", h->syntax);
    printf("Description: %s\n", h->description);
    printf("\nExamples:\n");
    for (i = 0; i < 4 && h->examples[i] != NULL; i++) {
        printf("  %s\n", h->examples[i]);
    }
    printf("\n");
}

/* Check if function exists */
int function_exists(const char *name) {
    return find_help(name) != NULL;
}

/* Show demo for a specific function */
void show_function_demo(const char *name) {
    const HelpEntry *h = find_help(name);
    int i;
    
    if (h == NULL) {
        printf("No demo available for '%s'\n", name);
        return;
    }
    
    printf("\n=== Demo: %s ===\n", h->name);
    printf("Description: %s\n\n", h->description);
    
    for (i = 0; i < 4 && h->examples[i] != NULL; i++) {
        const char *ex = h->examples[i];
        const char *eq = strchr(ex, '=');
        
        if (eq != NULL) {
            /* Extract expression part before '=' */
            char expr[256];
            int len = (int)(eq - ex);
            if (len > 255) len = 255;
            strncpy(expr, ex, len);
            expr[len] = '\0';
            
            /* Trim whitespace */
            while (len > 0 && (expr[len-1] == ' ' || expr[len-1] == '\t')) {
                expr[--len] = '\0';
            }
            
            printf(">>> %s\n", expr);
            /* Actually execute if it looks like a simple expression */
            /* For now, just show expected result */
            printf("    %s\n\n", eq + 1);
        } else {
            printf(">>> %s\n\n", ex);
        }
    }
}

/* Count help entries */
static int count_help_entries(void) {
    int count = 0;
    while (help_entries[count].name != NULL) count++;
    return count;
}

/* Run benchmark on all functions */
void run_function_bench(void) {
    printf("\nFunction benchmark not yet implemented.\n");
    printf("Total documented functions: %d\n\n", count_help_entries());
}
