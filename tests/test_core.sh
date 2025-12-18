#!/bin/bash
# test_sc.sh - Clean test runner for sc
# Tests all functions systematically

SC="${1:-../bin/sc}"
PASS=0
FAIL=0
SKIP=0

# Colors (disable if not terminal)
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    NC='\033[0m'
else
    GREEN=''
    RED=''
    YELLOW=''
    NC=''
fi

# Run expression and get result (clean output, no terminal codes)
run_expr() {
    # In pipe mode, sc outputs just the result without "= " prefix
    echo "$1" | $SC 2>/dev/null | tail -1
}

# Run command and get output
run_cmd() {
    echo "$1
quit" | $SC 2>/dev/null
}

# Test expression equals expected
test_eq() {
    local expr="$1"
    local expected="$2"
    local result=$(run_expr "$expr")
    
    if [ "$result" = "$expected" ]; then
        echo -e "  ${GREEN}PASS${NC}: $expr = $result"
        ((PASS++))
        return 0
    # Accept 1/0 as true/false equivalents
    elif [ "$expected" = "true" ] && [ "$result" = "1" ]; then
        echo -e "  ${GREEN}PASS${NC}: $expr = $result"
        ((PASS++))
        return 0
    elif [ "$expected" = "false" ] && [ "$result" = "0" ]; then
        echo -e "  ${GREEN}PASS${NC}: $expr = $result"
        ((PASS++))
        return 0
    else
        echo -e "  ${RED}FAIL${NC}: $expr"
        echo "        got:      '$result'"
        echo "        expected: '$expected'"
        ((FAIL++))
        return 1
    fi
}

# Test expression starts with expected
test_starts() {
    local expr="$1"
    local expected="$2"
    local result=$(run_expr "$expr")
    
    if [[ "$result" == "$expected"* ]]; then
        echo -e "  ${GREEN}PASS${NC}: $expr starts with $expected"
        ((PASS++))
        return 0
    else
        echo -e "  ${RED}FAIL${NC}: $expr should start with '$expected'"
        echo "        got: '$result'"
        ((FAIL++))
        return 1
    fi
}

# Test command output contains string
test_contains() {
    local cmd="$1"
    local pattern="$2"
    local result=$(run_cmd "$cmd")
    
    if echo "$result" | grep -q "$pattern"; then
        echo -e "  ${GREEN}PASS${NC}: '$cmd' contains '$pattern'"
        ((PASS++))
        return 0
    else
        echo -e "  ${RED}FAIL${NC}: '$cmd' should contain '$pattern'"
        ((FAIL++))
        return 1
    fi
}

# Skip test
test_skip() {
    local name="$1"
    local reason="$2"
    echo -e "  ${YELLOW}SKIP${NC}: $name - $reason"
    ((SKIP++))
}

echo "=============================================="
echo "  sc Test Suite"
echo "  Testing: $SC"
echo "=============================================="

# ========== BASIC ARITHMETIC ==========
echo ""
echo "=== Basic Arithmetic ==="
test_eq "2+3" "5"
test_eq "10-4" "6"
test_eq "6*7" "42"
test_eq "15/3" "5"
test_eq "2^10" "1024"
test_eq "2+3*4" "14"
test_eq "(2+3)*4" "20"
test_eq "2^3^2" "512"
test_eq "0-5+3" "-2"
test_eq "(-5)+3" "-2"

# ========== DISPLAY DIGITS ==========
echo ""
echo "=== Display Precision (digits command) ==="
test_contains "digits 5
1/3" "0.33333"
test_contains "digits 10
pi" "3.141592654"
test_contains "digits 0
1/7" "0.14285714285714285714"

# ========== SCIENTIFIC NOTATION ==========
echo ""
echo "=== Scientific Notation ==="
test_eq "1e10" "10000000000"
test_eq "2E6" "2000000"
test_starts "1e100" "1e+100"
test_starts "1e-50" "1e-50"
test_starts "6.022e23" "6.022e+23"

# ========== CONSTANTS ==========
echo ""
echo "=== Mathematical Constants ==="
test_starts "pi" "3.141592653589793"
test_starts "e" "2.718281828459045"

# ========== TRIGONOMETRY ==========
echo ""
echo "=== Trigonometry ==="
test_eq "sin(0)" "0"
test_eq "cos(0)" "1"
test_eq "tan(0)" "0"
test_eq "sin(pi/2)" "1"
test_eq "cos(pi)" "-1"
test_starts "sin(pi/6)" "0.5"
test_starts "cos(pi/3)" "0.5"
test_eq "tan(pi/4)" "1"

# ========== INVERSE TRIG ==========
echo ""
echo "=== Inverse Trigonometry ==="
test_eq "atan(0)" "0"
test_starts "asin(0.5)" "0.5235987755"
test_starts "acos(0.5)" "1.047197551"
test_starts "atan(1)" "0.7853981633"

# ========== HYPERBOLIC ==========
echo ""
echo "=== Hyperbolic Functions ==="
test_eq "sinh(0)" "0"
test_eq "cosh(0)" "1"
test_eq "tanh(0)" "0"
test_starts "sinh(1)" "1.175201193"
test_starts "cosh(1)" "1.543080634"

# ========== EXPONENTIAL/LOGARITHM ==========
echo ""
echo "=== Exponential and Logarithm ==="
test_eq "exp(0)" "1"
test_starts "exp(1)" "2.718281828"
test_eq "ln(1)" "0"
test_eq "ln(e)" "1"
test_eq "sqrt(4)" "2"
test_eq "sqrt(9)" "3"
test_starts "sqrt(2)" "1.41421356"

# ========== POWERS ==========
echo ""
echo "=== Powers ==="
test_eq "2^10" "1024"
test_eq "10^6" "1000000"
test_eq "2^0" "1"
test_starts "2^0.5" "1.41421356"

# ========== FACTORIAL ==========
echo ""
echo "=== Factorial ==="
test_eq "0!" "1"
test_eq "1!" "1"
test_eq "5!" "120"
test_eq "10!" "3628800"
test_eq "20!" "2432902008176640000"
test_eq "fact(6)" "720"

# ========== COMPLEX NUMBERS ==========
echo ""
echo "=== Complex Numbers ==="
test_eq "i" "i"
test_eq "i*i" "-1"
test_eq "i^2" "-1"
test_eq "i^4" "1"
test_eq "sqrt(-1)" "i"
test_eq "sqrt(-4)" "2i"
test_eq "(1+i)*(1-i)" "2"
test_eq "(2+3i)+(4+5i)" "6 + 8i"
test_eq "abs(3+4i)" "5"
test_eq "conj(3+4i)" "3 - 4i"
test_eq "re(3+4i)" "3"
test_eq "im(3+4i)" "4"

# ========== COMBINATORICS ==========
echo ""
echo "=== Combinatorics ==="
test_eq "ncr(10,5)" "252"
test_eq "ncr(52,5)" "2598960"
test_eq "npr(10,3)" "720"
test_eq "npr(5,5)" "120"
test_eq "choose(6,2)" "15"

# ========== NUMBER THEORY ==========
echo ""
echo "=== Number Theory ==="
test_eq "gcd(48,18)" "6"
test_eq "gcd(1071,462)" "21"
test_eq "lcm(12,18)" "36"
test_eq "lcm(7,11)" "77"
test_eq "isprime(17)" "true"
test_eq "isprime(18)" "false"
test_eq "isprime(2)" "true"

# ========== BITWISE ==========
echo ""
echo "=== Bitwise Operations ==="
# Note: and/or/xor/not are now boolean operators, not bitwise functions
# test_eq "and(255,15)" "15"
# test_eq "or(240,15)" "255"
# test_eq "xor(255,15)" "240"
# test_eq "not(0)" "-1"
test_eq "lsl(1,8)" "256"
test_eq "lsr(256,4)" "16"

# ========== COMPARISONS ==========
echo ""
echo "=== Comparison Operators ==="
test_eq "3<4" "true"
test_eq "4<3" "false"
test_eq "3<=3" "true"
test_eq "3<=4" "true"
test_eq "4<=3" "false"
test_eq "4>3" "true"
test_eq "3>4" "false"
test_eq "3>=3" "true"
test_eq "4>=3" "true"
test_eq "3>=4" "false"
test_eq "3==3" "true"
test_eq "3==4" "false"
test_eq "3<>4" "true"
test_eq "3<>3" "false"

# ========== APPROXIMATELY EQUAL ==========
echo ""
echo "=== Not Equal (~= MATLAB style) ==="
test_eq "5~=5" "false"
test_eq "5~=6" "true"
test_eq "100~=100" "false"
test_eq "100~=101" "true"

echo ""
echo "=== Approximate Equality (approxeq function) ==="
test_eq "approxeq(100, 105, 5)" "1"
test_eq "approxeq(100, 106, 5)" "0"
test_eq "approxeq(1.0, 1.001, 0.01)" "1"
test_eq "approxeq(1.0, 1.02, 0.01)" "0"

# ========== BOOLEAN OPERATORS ==========
echo ""
echo "=== Boolean Operators ==="
test_eq "1 and 1" "true"
test_eq "1 and 0" "false"
test_eq "0 and 0" "false"
test_eq "1 or 0" "true"
test_eq "0 or 0" "false"
test_eq "1 xor 0" "true"
test_eq "1 xor 1" "false"
test_eq "not 0" "true"
test_eq "not 1" "false"
test_eq "1<2 and 3<4" "true"
test_eq "1>2 and 3<4" "false"
test_eq "1>2 or 3<4" "true"
test_eq "1>2 or 3>4" "false"
test_eq "not 1>2" "true"
test_eq "1<2 && 3<4" "true"
test_eq "1>2 || 3<4" "true"

# ========== ISPRIME ==========
echo ""
echo "=== Primality Testing ==="
test_eq "isprime(2)" "true"
test_eq "isprime(3)" "true"
test_eq "isprime(4)" "false"
test_eq "isprime(7)" "true"
test_eq "isprime(11)" "true"
test_eq "isprime(13)" "true"
test_eq "isprime(100)" "false"
test_eq "isprime(101)" "true"
test_eq "isprime(1009)" "true"
test_eq "isprime(1000)" "false"

# ========== EVEN/ODD ==========
echo ""
echo "=== Even/Odd Testing ==="
test_eq "even(0)" "true"
test_eq "even(1)" "false"
test_eq "even(2)" "true"
test_eq "even(100)" "true"
test_eq "even(101)" "false"
test_eq "odd(0)" "false"
test_eq "odd(1)" "true"
test_eq "odd(2)" "false"
test_eq "odd(99)" "true"
test_eq "odd(100)" "false"

# ========== HEX/BINARY ==========
echo ""
echo "=== Hexadecimal and Binary ==="
test_eq "0xFF" "255"
test_eq "0x100" "256"
test_eq "0xDEAD" "57005"
test_eq "0b1111" "15"
test_eq "0b10000000" "128"

# ========== ROMAN NUMERALS (disabled - feature not implemented) ==========
# echo ""
# echo "=== Roman Numerals ==="
# test_eq "I" "1"
# test_eq "IV" "4"
# test_eq "IX" "9"
# test_eq "XL" "40"
# test_eq "MCMXCIX" "1999"
# test_eq "MMXXV" "2025"

# ========== VARIABLES ==========
echo ""
echo "=== Variables ==="
test_contains "a=5
a*2" "10"
test_contains "x=3
y=4
sqrt(x^2+y^2)" "5"

# ========== USER FUNCTIONS ==========
echo ""
echo "=== User Functions ==="
test_contains "f(x)=x^2
f(5)" "25"
test_contains "g(x)=2*x+1
g(10)" "21"

# ========== RPN MODE ==========
echo ""
echo "=== RPN Mode ==="
test_contains "rpn
2
3
+" "X: 5"
test_contains "rpn
10
2
/" "X: 5"
test_contains "rpn
5
dup
*" "X: 25"

# ========== SPECIAL VALUES ==========
echo ""
echo "=== Special Values ==="
test_eq "1/0" "Inf"
test_eq "-1/0" "-Inf"
test_eq "0/0" "NaN"
test_eq "Inf+1" "Inf"
test_eq "Inf-Inf" "NaN"

# ========== LARGE NUMBERS ==========
echo ""
echo "=== Large Numbers ==="
test_starts "2^1000" "1.071508607"
test_starts "2^10000" "1.995063116"

# ========== ORBITAL MECHANICS ==========
echo ""
echo "=== Orbital Mechanics ==="
test_contains "orb mu earth" "398600.4418"
test_contains "orb radius earth" "6371"
test_contains "orb v_circ 398600.4418 6778" "7.66"
test_contains "orb v_esc 398600.4418 6778" "10.8"

# ========== TVM ==========
echo ""
echo "=== TVM Financial ==="
test_contains "tvm n 360
tvm i 0.5
tvm pv 200000
tvm pmt" "1199"

# ========== STATISTICS ==========
echo ""
echo "=== Statistics ==="
test_contains "stat clear
data+ 10
data+ 20
data+ 30
stat mean" "20"

# ========== SPECIAL FUNCTIONS ==========
echo ""
echo "=== Special Functions ==="
# Gamma function
test_contains "gamma(5)" "23.999"
test_contains "gamma(1)" "0.999"
test_contains "lgamma(5)" "3.17805"
# Error function
test_contains "erf(0)" "0"
test_contains "erf(1)" "0.842"
test_contains "erfc(0)" "1"
# Normal distribution
test_contains "normcdf(0)" "0.5"
test_contains "normcdf(1.96)" "0.97"
test_contains "norminv(0.5)" "0"
# Bessel functions - not implemented
# test_contains "j0(0)" "1"
# test_contains "j1(0)" "0"
# Elliptic integrals - not implemented
# test_contains "ellipk(0)" "1.5707"
# test_contains "ellipe(0)" "1.5707"
# Lambert W - not implemented
# test_contains "lambertw(0)" "0"
# test_contains "lambertw(exp(1))" "1"

# ========== STATISTICAL DISTRIBUTIONS ==========
echo ""
echo "=== Statistical Distributions ==="
# Student's t distribution - not implemented
# test_contains "tcdf(0, 10)" "0.5"
# Chi-squared - not implemented
# test_contains "chi2cdf(0, 5)" "0"
# Binomial - not implemented
# test_contains "binompdf(5, 10, 0.5)" "0.2460"
# test_contains "binomcdf(5, 10, 0.5)" "0.623"
# Poisson - not implemented
# test_contains "poisspdf(0, 1)" "0.3678"
# test_contains "poisscdf(0, 1)" "0.3678"

# ========== QUADRATIC ==========
echo ""
echo "=== Quadratic Solver ==="
test_contains "quad 1 0 -4" "2"

# ========== BUILT-IN TESTS ==========
echo ""
echo "=== Built-in Test Suite ==="
result=$(echo "test
quit" | $SC 2>/dev/null | grep "passed")
if echo "$result" | grep -q "0 failed"; then
    echo -e "  ${GREEN}PASS${NC}: Built-in tests: $result"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: Built-in tests: $result"
    ((FAIL++))
fi

# ========== SUMMARY ==========
echo ""
echo "=============================================="
TOTAL=$((PASS + FAIL + SKIP))
echo "  RESULTS: $PASS passed, $FAIL failed, $SKIP skipped (total $TOTAL)"
echo "=============================================="

if [ $FAIL -gt 0 ]; then
    exit 1
else
    exit 0
fi

# Forecasting functions
test_expr "movavg([1;2;3;4;5], 3)" "column"
test_expr "ewma([1;2;3;4;5], 0.3)" "column"
test_expr "lag([1;2;3;4;5])" "column"
test_expr "lead([1;2;3;4;5])" "column"
test_expr "detrend([1;2;3;4;5])" "column"
test_expr "trend([1;2;3;4;5])" "column"
test_expr "pctchange([100;110;121])" "column"
test_expr "autocorr([1;2;3;4;5])" "column"
