#!/bin/bash
# Comprehensive test suite for scalc - tests all ~350 functions
SC="${1:-../bin/sc}"
PASS=0
FAIL=0

# Test expression equals expected (handles matrix output with spaces/newlines)
test_expr() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | tr '\n' ' ' | sed 's/  */ /g;s/^ *//;s/ *$//')
    expected=$(echo "$expected" | tr '\n' ' ' | sed 's/  */ /g;s/^ *//;s/ *$//')
    if [ "$result" = "$expected" ]; then
        ((PASS++)); return 0
    else
        echo "  FAIL: $expr = '$result' (expected '$expected')"
        ((FAIL++)); return 1
    fi
}

# Test approximate equality
test_approx() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | head -1 | sed 's/^ *//;s/ *$//')
    if echo "$result" | grep -qi "error"; then
        echo "  FAIL: $expr (error)"; ((FAIL++)); return 1
    fi
    local ok=$(echo "$result $expected" | awk '{d=$1-$2; if(d<0)d=-d; t=$2*0.0001; if(t<0)t=-t; if(t<0.0001)t=0.0001; print (d<=t)?1:0}')
    if [ "$ok" = "1" ]; then
        ((PASS++)); return 0
    else
        echo "  FAIL: $expr = '$result' (expected ~$expected)"
        ((FAIL++)); return 1
    fi
}

# Test no error
test_ok() {
    local expr="$1"
    local result
    result=$(echo "$expr" | $SC 2>&1)
    if echo "$result" | grep -qi "error"; then
        echo "  FAIL: $expr (error)"; ((FAIL++)); return 1
    else
        ((PASS++)); return 0
    fi
}

# Test boolean (accepts 1, 0, true, false)
test_bool() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | head -1 | sed 's/^ *//;s/ *$//')
    if [ "$expected" = "true" ]; then
        if [ "$result" = "1" ] || [ "$result" = "true" ]; then ((PASS++)); return 0; fi
    else
        if [ "$result" = "0" ] || [ "$result" = "false" ]; then ((PASS++)); return 0; fi
    fi
    echo "  FAIL: $expr = '$result' (expected $expected)"
    ((FAIL++)); return 1
}

echo "=== SCALC COMPREHENSIVE TEST SUITE ==="
echo

# CONSTANTS
test_approx "pi" "3.14159265358979"
test_approx "e" "2.71828182845905"
test_expr "i*i" "-1"
test_expr "Inf" "Inf"

# ARITHMETIC
test_expr "2+3" "5"
test_expr "10-7" "3"
test_expr "4*5" "20"
test_expr "15/3" "5"
test_expr "2^10" "1024"
test_expr "5!" "120"
test_expr "mod(17,5)" "2"
test_expr "rem(17,5)" "2"
test_expr "mod(-7,3)" "2"
test_expr "rem(-7,3)" "-1"

# TRIGONOMETRY
test_expr "sin(0)" "0"
test_approx "sin(pi/2)" "1"
test_expr "cos(0)" "1"
test_approx "cos(pi)" "-1"
test_expr "tan(0)" "0"
test_approx "tan(pi/4)" "1"
test_expr "asin(0)" "0"
test_expr "acos(1)" "0"
test_expr "atan(0)" "0"
test_approx "sinh(1)" "1.1752011936438"
test_expr "cosh(0)" "1"
test_expr "tanh(0)" "0"
test_expr "asinh(0)" "0"
test_expr "acosh(1)" "0"
test_expr "atanh(0)" "0"
test_expr "sec(0)" "1"
test_approx "csc(pi/2)" "1"
test_approx "cot(pi/4)" "1"
test_expr "sech(0)" "1"
test_approx "csch(1)" "0.8509181282"
test_approx "coth(1)" "1.3130352855"
test_approx "asech(0.5)" "1.3169578969"
test_approx "acsch(1)" "0.8813735870"
test_approx "acoth(2)" "0.5493061443"
test_expr "sind(30)" "0.5"
test_expr "cosd(60)" "0.5"
test_expr "tand(45)" "1"
test_expr "asind(0.5)" "30"
test_expr "acosd(0.5)" "60"
test_expr "atand(1)" "45"
test_expr "atan2d(1,1)" "45"
test_expr "sinc(0)" "1"

# EXPONENTIAL/LOG
test_expr "exp(0)" "1"
test_approx "exp(1)" "2.71828182845905"
test_expr "ln(1)" "0"
test_expr "log10(100)" "2"
test_expr "log2(8)" "3"
test_expr "exp2(3)" "8"
test_expr "pow2(4)" "16"
test_expr "sqrt(9)" "3"
test_expr "cbrt(27)" "3"
test_expr "nthroot(16,4)" "2"
test_approx "log1p(0.001)" "0.0009995"
test_approx "expm1(0.001)" "0.0010005"
test_approx "reallog(10)" "2.302585"
test_expr "realsqrt(4)" "2"
test_expr "realpow(2,3)" "8"

# COMPLEX
test_expr "re(3+4i)" "3"
test_expr "im(3+4i)" "4"
test_ok "conj(3+4i)"
test_expr "abs(3+4i)" "5"
test_approx "arg(1+i)" "0.78539816"
test_ok "complex(3,4)"

# ROUNDING
test_expr "floor(3.7)" "3"
test_expr "floor(-3.7)" "-4"
test_expr "ceil(3.2)" "4"
test_expr "trunc(3.7)" "3"
test_expr "round(3.5)" "4"
test_expr "fix(-3.7)" "-3"
test_expr "frac(3.7)" "0.7"
test_expr "sign(5)" "1"
test_expr "sign(-5)" "-1"
test_expr "abs(-5)" "5"

# NUMBER THEORY
test_expr "gcd(12,18)" "6"
test_expr "lcm(4,6)" "12"
test_bool "isprime(7)" "true"
test_bool "isprime(8)" "false"
test_ok "prime(5)"
test_bool "even(4)" "true"
test_bool "odd(5)" "true"
test_expr "fact(5)" "120"
test_expr "ncr(5,2)" "10"
test_expr "npr(5,2)" "20"
test_approx "gamma(5)" "24"
test_expr "totient(10)" "4"
test_expr "fibonacci(10)" "55"
test_expr "lucas(10)" "123"
test_expr "catalan(5)" "42"
test_expr "factorial2(10)" "3840"

# SPECIAL
test_approx "erf(1)" "0.84270079"
test_approx "erfc(1)" "0.15729921"
test_expr "sigmoid(0)" "0.5"
test_approx "softplus(0)" "0.69314718"
test_expr "step(1)" "1"
test_expr "step(-1)" "0"
test_expr "heaviside(1)" "1"
test_expr "rect(0)" "1"
test_expr "tri(0)" "1"

# STATISTICS
test_expr "sum([1,2,3,4,5])" "15"
test_expr "mean([1,2,3,4,5])" "3"
test_expr "median([1,2,3,4,5])" "3"
test_expr "min([3,1,4])" "1"
test_expr "max([3,1,4])" "4"
test_expr "prod([1,2,3,4])" "24"
test_expr "range([1,5,3,9])" "8"
test_approx "std([1,2,3,4,5])" "1.5811388"
test_approx "var([1,2,3,4,5])" "2.5"
test_expr "sumsq([1,2,3])" "14"
test_approx "geomean([1,2,4,8])" "2.8284271"
test_approx "harmmean([1,2,4])" "1.7142857"
test_expr "skewness([1,2,3,4,5])" "0"
test_expr "mad([1,2,3,4,5])" "1.2"
test_expr "iqr([1,2,3,4,5,6,7,8])" "3.5"
test_expr "prctile([1,2,3,4,5],50)" "3"

# DISTRIBUTIONS
test_approx "normpdf(0)" "0.39894228"
test_expr "normcdf(0)" "0.5"
test_approx "norminv(0.975)" "1.9604"
test_approx "norminv(0.5)" "0" "0.001"
test_ok "rand"
test_ok "randn"
test_ok "randi(10,3)"
test_ok "randperm(5)"

# LINEAR ALGEBRA
test_expr "det([1,2;3,4])" "-2"
test_expr "trace([1,2;3,4])" "5"
test_ok "inv([1,2;3,4])"
test_ok "trans([1,2;3,4])"
test_ok "rank([1,2;3,4])"
test_approx "norm([3,4])" "5"
test_approx "cond([1,0;0,2])" "2"

# MATRIX CREATION
test_expr "zeros(2,2)" "0 0 0 0"
test_expr "ones(2,2)" "1 1 1 1"
test_expr "eye(2)" "1 0 0 1"
test_expr "linspace(0,10,3)" "0 5 10"
test_ok "magic(3)"
test_ok "magic(4)"
test_ok "pascal(4)"
test_ok "hilb(3)"
test_ok "vander([1,2,3])"
test_ok "toeplitz([1,2,3])"
test_ok "hankel([1,2,3])"

# MATRIX MANIPULATION
test_ok "reshape([1,2,3,4,5,6],2,3)"
test_ok "repmat([1,2],2,2)"
test_ok "rot90([1,2;3,4])"
test_expr "flip([1,2,3])" "3 2 1"
test_ok "triu([1,2,3;4,5,6;7,8,9])"
test_ok "tril([1,2,3;4,5,6;7,8,9])"
test_expr "sort([3,1,4,1,5])" "1 1 3 4 5"

# MATRIX QUERIES
test_ok "size([1,2,3])"
test_expr "length([1,2,3])" "3"
test_expr "rows([1,2;3,4;5,6])" "3"
test_expr "cols([1,2;3,4;5,6])" "2"
test_expr "numel([1,2;3,4])" "4"
test_bool "isempty([])" "true"
test_bool "isscalar(5)" "true"
test_bool "isvector([1,2,3])" "true"
test_bool "ismatrix([1,2;3,4])" "true"
test_bool "issquare([1,2;3,4])" "true"
test_bool "issorted([1,2,3])" "true"
test_bool "issymmetric([1,2;2,1])" "true"

# SIGNAL PROCESSING
test_expr "conv([1,2],[1,1])" "1 3 2"
test_expr "deconv([1,3,2],[1,1])" "1 2"
test_ok "gradient([1,2,4,7,11])"
test_approx "trapz([0,1,2,3])" "4.5"
test_ok "movmean([1,2,3,4,5],3)"
test_ok "diff([1,2,4,7,11])"
test_expr "cumsum([1,2,3,4])" "1 3 6 10"
test_expr "cumprod([1,2,3,4])" "1 2 6 24"
test_expr "interp1([0,1,2],[0,1,4],1.5)" "2.5"
test_approx "corrcoef([1,2,3],[1,2,3])" "1"

# POLYNOMIALS
test_expr "polyval([1,2,3],2)" "11"
test_expr "polyder([3,2,1])" "6 2"
test_expr "polyint([2,1])" "1 1 0"
test_expr "quadratic(1,-5,6)" "3 2"

# ANGLE
test_approx "deg2rad(180)" "3.14159265"
test_expr "rad2deg(pi)" "180"
test_expr "hypot(3,4)" "5"
test_expr "wrap360(370)" "10"
test_expr "wrap180(270)" "-90"

# SET OPERATIONS
test_expr "unique([1,2,2,3,3,3])" "1 2 3"
test_expr "union([1,2,3],[3,4,5])" "1 2 3 4 5"
test_expr "intersect([1,2,3,4],[3,4,5,6])" "3 4"
test_expr "setdiff([1,2,3,4],[3,4,5])" "1 2"
test_expr "setxor([1,2,3],[2,3,4])" "1 4"

# LOGICAL
test_bool "1 && 1" "true"
test_bool "1 && 0" "false"
test_bool "0 || 1" "true"
test_bool "not(0)" "true"
test_bool "any([0,0,1])" "true"
test_bool "all([1,1,1])" "true"
test_bool "eq(5,5)" "true"
test_bool "ne(5,6)" "true"
test_bool "lt(3,5)" "true"
test_bool "isnan(NaN)" "true"
test_bool "isinf(Inf)" "true"
test_bool "isfinite(5)" "true"
test_bool "isreal(5)" "true"

# BITWISE
test_expr "bitand(12,10)" "8"
test_expr "bitor(12,10)" "14"
test_expr "bitxor(12,10)" "6"
test_expr "bitshift(8,2)" "32"
test_expr "shl(1,4)" "16"
test_expr "shr(16,2)" "4"

# COORDINATE TRANSFORMS
test_ok "cart2pol(3,4)"
test_ok "pol2cart(5,0)"
test_ok "cart2sph(1,1,1)"
test_ok "sph2cart(0.785,0.615,1.73)"

# UTILITY
test_expr "clamp(5,0,3)" "3"
test_ok "rescale([1,2,3,4,5])"
test_ok "normalize([3,4])"
test_expr "find([0,1,0,1,1])" "2 4 5"
test_expr "nextpow2(5)" "3"
test_ok "eps"

# MATRIX OPERATIONS
test_ok "cat(1,[1,2],[3,4])"
test_expr "horzcat([1,2],[3,4])" "1 2 3 4"
test_ok "vertcat([1,2],[3,4])"
test_expr "[1,2,3]+[4,5,6]" "5 7 9"
test_expr "[1,2,3].*[4,5,6]" "4 10 18"
test_expr "[1,2,3].^2" "1 4 9"

# COLON
test_expr "1:5" "1 2 3 4 5"
test_expr "0:2:10" "0 2 4 6 8 10"
test_expr "5:-1:1" "5 4 3 2 1"

# LINEAR SOLVE
test_ok "linsolve([1,2;3,4],[5;11])"
test_ok "mldivide([1,2;3,4],[5;11])"

# DECOMPOSITIONS
test_ok "eig([1,2;2,1])"

# SPECIAL
test_expr "kron([1,2],[1,1])" "1 1 2 2"
test_expr "sub2ind([3,3],2,2)" "5"
test_bool "approxeq(1.0001,1.0002,0.001)" "true"
test_ok "printbin(255)"
test_ok "printhex(255)"

# COMMENTS
test_expr "5+3 % comment" "8"

echo
echo "=== RESULTS: $PASS passed, $FAIL failed ==="
exit $FAIL
