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
    echo "$1
quit" | $SC 2>/dev/null | grep "^=" | head -1 | sed 's/^= //'
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
test_eq "isprime(17)" "1"
test_eq "isprime(18)" "0"
test_eq "isprime(2)" "1"

# ========== BITWISE ==========
echo ""
echo "=== Bitwise Operations ==="
test_eq "and(255,15)" "15"
test_eq "or(240,15)" "255"
test_eq "xor(255,15)" "240"
test_eq "not(0)" "-1"
test_eq "lsl(1,8)" "256"
test_eq "lsr(256,4)" "16"

# ========== HEX/BINARY ==========
echo ""
echo "=== Hexadecimal and Binary ==="
test_eq "0xFF" "255"
test_eq "0x100" "256"
test_eq "0xDEAD" "57005"
test_eq "0b1111" "15"
test_eq "0b10000000" "128"

# ========== ROMAN NUMERALS ==========
echo ""
echo "=== Roman Numerals ==="
test_eq "I" "1"
test_eq "IV" "4"
test_eq "IX" "9"
test_eq "XL" "40"
test_eq "MCMXCIX" "1999"
test_eq "MMXXV" "2025"

# ========== VARIABLES ==========
echo ""
echo "=== Variables ==="
test_contains "a=5
a*2" "= 10"
test_contains "x=3
y=4
sqrt(x^2+y^2)" "= 5"

# ========== USER FUNCTIONS ==========
echo ""
echo "=== User Functions ==="
test_contains "f(x)=x^2
f(5)" "= 25"
test_contains "g(x)=2*x+1
g(10)" "= 21"

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
test_contains "orb mu earth" "3.986004418"
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
