#!/bin/bash
# test_ieee754.sh - Comprehensive IEEE 754 Floating Point Compliance Tests
#
# Tests sc against IEEE 754-2008 standard requirements
# Usage: ./test_ieee754.sh [path-to-sc]
#
# Note: sc uses 128-bit extended precision, so some "failures" are actually
# features (larger exponent range, no overflow at IEEE double limits).

SC="${1:-../bin/sc}"
PASS=0
FAIL=0

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Test function: test_ieee "name" "expression" "expected" [exact]
test_ieee() {
    local name="$1"
    local expr="$2"
    local expected="$3"
    local exact="${4:-0}"
    
    local result=$(echo "$expr" | $SC 2>&1 | grep "= " | tail -1)
    local value="${result#*= }"
    
    local passed=0
    if [ "$exact" = "1" ]; then
        [ "$value" = "$expected" ] && passed=1
    else
        [[ "$result" == *"$expected"* ]] && passed=1
    fi
    
    if [ $passed -eq 1 ]; then
        echo -e "  ${GREEN}PASS${NC}: $name"
        ((PASS++))
    else
        echo -e "  ${RED}FAIL${NC}: $name"
        echo "        expr: $expr"
        echo "        expected: $expected"
        echo "        got: $value"
        ((FAIL++))
    fi
}

section() {
    echo ""
    echo -e "${YELLOW}=== $1 ===${NC}"
}

subsection() {
    echo -e "  ${CYAN}--- $1 ---${NC}"
}

echo "====================================================="
echo "  Comprehensive IEEE 754-2008 Compliance Tests"
echo "  Testing: $SC"
echo "====================================================="

# ============================================================
# IEEE 754 Section 5.4.1: Arithmetic Operations
# ============================================================
section "IEEE 754 Section 5: Basic Arithmetic"

subsection "Addition"
test_ieee "0 + 0 = 0" "0 + 0" "0" 1
test_ieee "1 + 0 = 1" "1 + 0" "1" 1
test_ieee "Inf + 1 = Inf" "(1/0) + 1" "Inf"
test_ieee "-Inf + 1 = -Inf" "(-1/0) + 1" "-Inf"
test_ieee "Inf + Inf = Inf" "(1/0) + (1/0)" "Inf"
test_ieee "-Inf + -Inf = -Inf" "(-1/0) + (-1/0)" "-Inf"
test_ieee "Inf + -Inf = NaN" "(1/0) + (-1/0)" "NaN"

subsection "Subtraction"
test_ieee "1 - 1 = 0" "1 - 1" "0" 1
test_ieee "Inf - 1 = Inf" "(1/0) - 1" "Inf"
test_ieee "Inf - Inf = NaN" "(1/0) - (1/0)" "NaN"
test_ieee "-Inf - (-Inf) = NaN" "(-1/0) - (-1/0)" "NaN"
test_ieee "Inf - (-Inf) = Inf" "(1/0) - (-1/0)" "Inf"

subsection "Multiplication"
test_ieee "2 * 3 = 6" "2 * 3" "6" 1
test_ieee "Inf * 2 = Inf" "(1/0) * 2" "Inf"
test_ieee "Inf * -2 = -Inf" "(1/0) * (-2)" "-Inf"
test_ieee "-Inf * -2 = Inf" "(-1/0) * (-2)" "Inf"
test_ieee "Inf * Inf = Inf" "(1/0) * (1/0)" "Inf"
test_ieee "-Inf * Inf = -Inf" "(-1/0) * (1/0)" "-Inf"
test_ieee "-Inf * -Inf = Inf" "(-1/0) * (-1/0)" "Inf"
test_ieee "Inf * 0 = NaN" "(1/0) * 0" "NaN"
test_ieee "0 * Inf = NaN" "0 * (1/0)" "NaN"
test_ieee "-Inf * 0 = NaN" "(-1/0) * 0" "NaN"

subsection "Division"
test_ieee "6 / 2 = 3" "6 / 2" "3" 1
test_ieee "1 / 0 = Inf" "1 / 0" "Inf"
test_ieee "-1 / 0 = -Inf" "-1 / 0" "-Inf"
test_ieee "0 / 0 = NaN" "0 / 0" "NaN"
test_ieee "Inf / 2 = Inf" "(1/0) / 2" "Inf"
test_ieee "Inf / -2 = -Inf" "(1/0) / (-2)" "-Inf"
test_ieee "1 / Inf = 0" "1 / (1/0)" "0"
test_ieee "-1 / Inf = 0" "-1 / (1/0)" "0"
test_ieee "1 / -Inf = 0" "1 / (-1/0)" "0"
test_ieee "Inf / Inf = NaN" "(1/0) / (1/0)" "NaN"
test_ieee "-Inf / Inf = NaN" "(-1/0) / (1/0)" "NaN"
test_ieee "Inf / -Inf = NaN" "(1/0) / (-1/0)" "NaN"

# ============================================================
# IEEE 754 Section 6.2: Operations with NaN
# ============================================================
section "IEEE 754 Section 6.2: NaN Propagation"

subsection "Arithmetic with NaN"
test_ieee "NaN + 1 = NaN" "(0/0) + 1" "NaN"
test_ieee "NaN - 1 = NaN" "(0/0) - 1" "NaN"
test_ieee "NaN * 2 = NaN" "(0/0) * 2" "NaN"
test_ieee "NaN / 2 = NaN" "(0/0) / 2" "NaN"
test_ieee "1 + NaN = NaN" "1 + (0/0)" "NaN"
test_ieee "1 - NaN = NaN" "1 - (0/0)" "NaN"
test_ieee "2 * NaN = NaN" "2 * (0/0)" "NaN"
test_ieee "2 / NaN = NaN" "2 / (0/0)" "NaN"
test_ieee "NaN + NaN = NaN" "(0/0) + (0/0)" "NaN"
test_ieee "NaN * NaN = NaN" "(0/0) * (0/0)" "NaN"
test_ieee "NaN / NaN = NaN" "(0/0) / (0/0)" "NaN"

subsection "Functions with NaN"
test_ieee "sqrt(NaN) = NaN" "sqrt(0/0)" "NaN"
test_ieee "sin(NaN) = NaN" "sin(0/0)" "NaN"
test_ieee "cos(NaN) = NaN" "cos(0/0)" "NaN"
test_ieee "tan(NaN) = NaN" "tan(0/0)" "NaN"
test_ieee "exp(NaN) = NaN" "exp(0/0)" "NaN"
test_ieee "ln(NaN) = NaN" "ln(0/0)" "NaN"
test_ieee "abs(NaN) = NaN" "abs(0/0)" "NaN"
test_ieee "atan(NaN) = NaN" "atan(0/0)" "NaN"

# ============================================================
# IEEE 754 Section 9.2.1: pow(x,y) Special Cases
# ============================================================
section "IEEE 754 Section 9.2.1: Power Function Special Cases"

subsection "x^0 = 1 for ANY x (including NaN)"
test_ieee "0^0 = 1" "0^0" "1"
test_ieee "1^0 = 1" "1^0" "1"
test_ieee "(-1)^0 = 1" "(-1)^0" "1"
test_ieee "2^0 = 1" "2^0" "1"
test_ieee "Inf^0 = 1" "(1/0)^0" "1"
test_ieee "-Inf^0 = 1" "(-1/0)^0" "1"
test_ieee "NaN^0 = 1" "(0/0)^0" "1"

subsection "1^y = 1 for ANY y (including NaN, Inf)"
test_ieee "1^1 = 1" "1^1" "1"
test_ieee "1^100 = 1" "1^100" "1"
test_ieee "1^(-100) = 1" "1^(-100)" "1"
test_ieee "1^Inf = 1" "1^(1/0)" "1"
test_ieee "1^(-Inf) = 1" "1^(-1/0)" "1"
test_ieee "1^NaN = 1" "1^(0/0)" "1"

subsection "(-1)^Inf = 1 (because |-1| = 1)"
test_ieee "(-1)^Inf = 1" "(-1)^(1/0)" "1"
test_ieee "(-1)^(-Inf) = 1" "(-1)^(-1/0)" "1"

subsection "0^y cases"
test_ieee "0^1 = 0" "0^1" "0"
test_ieee "0^2 = 0" "0^2" "0"
test_ieee "0^0.5 = 0" "0^0.5" "0"
test_ieee "0^Inf = 0" "0^(1/0)" "0"
test_ieee "0^(-1) = Inf" "0^(-1)" "Inf"
test_ieee "0^(-2) = Inf" "0^(-2)" "Inf"
test_ieee "0^(-0.5) = Inf" "0^(-0.5)" "Inf"

subsection "Inf^y cases"
test_ieee "Inf^1 = Inf" "(1/0)^1" "Inf"
test_ieee "Inf^2 = Inf" "(1/0)^2" "Inf"
test_ieee "Inf^0.5 = Inf" "(1/0)^0.5" "Inf"
test_ieee "Inf^(-1) = 0" "(1/0)^(-1)" "0"
test_ieee "Inf^(-2) = 0" "(1/0)^(-2)" "0"
test_ieee "Inf^(-0.5) = 0" "(1/0)^(-0.5)" "0"

subsection "(-Inf)^n cases (sign depends on odd/even)"
test_ieee "(-Inf)^1 = -Inf" "(-1/0)^1" "-Inf"
test_ieee "(-Inf)^2 = Inf" "(-1/0)^2" "Inf"
test_ieee "(-Inf)^3 = -Inf" "(-1/0)^3" "-Inf"
test_ieee "(-Inf)^4 = Inf" "(-1/0)^4" "Inf"
test_ieee "(-Inf)^(-1) = 0" "(-1/0)^(-1)" "0"
test_ieee "(-Inf)^(-2) = 0" "(-1/0)^(-2)" "0"

subsection "x^Inf cases (depends on |x| vs 1)"
test_ieee "2^Inf = Inf" "2^(1/0)" "Inf"
test_ieee "10^Inf = Inf" "10^(1/0)" "Inf"
test_ieee "0.5^Inf = 0" "0.5^(1/0)" "0"
test_ieee "0.1^Inf = 0" "0.1^(1/0)" "0"
test_ieee "2^(-Inf) = 0" "2^(-1/0)" "0"
test_ieee "0.5^(-Inf) = Inf" "0.5^(-1/0)" "Inf"

subsection "Negative base with non-integer exponent (complex)"
test_ieee "(-1)^0.5 = i" "(-1)^0.5" "i"
test_ieee "sqrt(-1) = i" "sqrt(-1)" "i"
test_ieee "sqrt(-4) = 2i" "sqrt(-4)" "2i"

subsection "NaN propagation in pow"
test_ieee "NaN^2 = NaN" "(0/0)^2" "NaN"
test_ieee "2^NaN = NaN" "2^(0/0)" "NaN"
test_ieee "NaN^NaN = NaN" "(0/0)^(0/0)" "NaN"

# ============================================================
# IEEE 754 Section 9.2.1: sqrt Special Cases
# ============================================================
section "IEEE 754 Section 9.2.1: Square Root"

test_ieee "sqrt(0) = 0" "sqrt(0)" "0" 1
test_ieee "sqrt(1) = 1" "sqrt(1)" "1" 1
test_ieee "sqrt(4) = 2" "sqrt(4)" "2" 1
test_ieee "sqrt(9) = 3" "sqrt(9)" "3" 1
test_ieee "sqrt(2) ≈ 1.414" "sqrt(2)" "1.41421356"
test_ieee "sqrt(Inf) = Inf" "sqrt(1/0)" "Inf"
test_ieee "sqrt(NaN) = NaN" "sqrt(0/0)" "NaN"
test_ieee "sqrt(-1) = i" "sqrt(-1)" "i"
test_ieee "sqrt(-4) = 2i" "sqrt(-4)" "2i"
test_ieee "sqrt(-9) = 3i" "sqrt(-9)" "3i"

# ============================================================
# IEEE 754 Section 9.2.1: exp/log Special Cases
# ============================================================
section "IEEE 754 Section 9.2.1: Exponential and Logarithm"

subsection "exp special cases"
test_ieee "exp(0) = 1" "exp(0)" "1" 1
test_ieee "exp(1) = e" "exp(1)" "2.71828"
test_ieee "exp(-1) = 1/e" "exp(-1)" "0.36787"
test_ieee "exp(Inf) = Inf" "exp(1/0)" "Inf"
test_ieee "exp(-Inf) = 0" "exp(-1/0)" "0"
test_ieee "exp(NaN) = NaN" "exp(0/0)" "NaN"

subsection "ln special cases"
test_ieee "ln(1) = 0" "ln(1)" "0" 1
test_ieee "ln(e) = 1" "ln(e)" "1"
test_ieee "ln(e^2) = 2" "ln(e^2)" "2"
test_ieee "ln(0) = -Inf" "ln(0)" "-Inf"
test_ieee "ln(Inf) = Inf" "ln(1/0)" "Inf"
test_ieee "ln(NaN) = NaN" "ln(0/0)" "NaN"
test_ieee "ln(-1) = πi (complex)" "ln(-1)" "3.14159"

subsection "exp/ln identities"
test_ieee "ln(exp(2)) = 2" "ln(exp(2))" "2"
test_ieee "exp(ln(2)) = 2" "exp(ln(2))" "2"
test_ieee "exp(ln(10)) = 10" "exp(ln(10))" "10"

# ============================================================
# IEEE 754 Section 9.2.1: Trigonometric Functions
# ============================================================
section "IEEE 754 Section 9.2.1: Trigonometry"

subsection "sin special cases"
test_ieee "sin(0) = 0" "sin(0)" "0" 1
test_ieee "sin(pi/6) = 0.5" "sin(pi/6)" "0.5"
test_ieee "sin(pi/2) = 1" "sin(pi/2)" "1"
test_ieee "sin(pi) ≈ 0" "sin(pi)" "0"
test_ieee "sin(Inf) = NaN" "sin(1/0)" "NaN"
test_ieee "sin(-Inf) = NaN" "sin(-1/0)" "NaN"
test_ieee "sin(NaN) = NaN" "sin(0/0)" "NaN"

subsection "cos special cases"
test_ieee "cos(0) = 1" "cos(0)" "1" 1
test_ieee "cos(pi/3) = 0.5" "cos(pi/3)" "0.5"
test_ieee "cos(pi) = -1" "cos(pi)" "-1"
test_ieee "cos(Inf) = NaN" "cos(1/0)" "NaN"
test_ieee "cos(-Inf) = NaN" "cos(-1/0)" "NaN"
test_ieee "cos(NaN) = NaN" "cos(0/0)" "NaN"

subsection "tan special cases"
test_ieee "tan(0) = 0" "tan(0)" "0" 1
test_ieee "tan(pi/4) = 1" "tan(pi/4)" "1"
test_ieee "tan(Inf) = NaN" "tan(1/0)" "NaN"
test_ieee "tan(NaN) = NaN" "tan(0/0)" "NaN"

subsection "atan special cases"
test_ieee "atan(0) = 0" "atan(0)" "0" 1
test_ieee "atan(1) = pi/4" "atan(1)" "0.78539816"
test_ieee "atan(Inf) = pi/2" "atan(1/0)" "1.5707963"
test_ieee "atan(-Inf) = -pi/2" "atan(-1/0)" "-1.5707963"
test_ieee "atan(NaN) = NaN" "atan(0/0)" "NaN"

subsection "Trigonometric identities"
test_ieee "sin²(1) + cos²(1) = 1" "sin(1)^2 + cos(1)^2" "1"
test_ieee "sin(pi/6)² = 0.25" "sin(pi/6)^2" "0.25"

# ============================================================
# IEEE 754 Section 9.2.1: Hyperbolic Functions
# ============================================================
section "IEEE 754 Section 9.2.1: Hyperbolic Functions"

test_ieee "sinh(0) = 0" "sinh(0)" "0" 1
test_ieee "cosh(0) = 1" "cosh(0)" "1" 1
test_ieee "tanh(0) = 0" "tanh(0)" "0" 1
test_ieee "sinh(1) ≈ 1.175" "sinh(1)" "1.175201"
test_ieee "cosh(1) ≈ 1.543" "cosh(1)" "1.543080"
test_ieee "sinh(Inf) = Inf" "sinh(1/0)" "Inf"
test_ieee "sinh(-Inf) = -Inf" "sinh(-1/0)" "-Inf"
test_ieee "cosh(Inf) = Inf" "cosh(1/0)" "Inf"
test_ieee "cosh(-Inf) = Inf" "cosh(-1/0)" "Inf"
test_ieee "tanh(Inf) = 1" "tanh(1/0)" "1"
test_ieee "tanh(-Inf) = -1" "tanh(-1/0)" "-1"
test_ieee "tanh(100) = 1" "tanh(100)" "1"
test_ieee "tanh(-100) = -1" "tanh(-100)" "-1"

subsection "Hyperbolic identity"
test_ieee "cosh²(1) - sinh²(1) = 1" "cosh(1)^2 - sinh(1)^2" "1"

# ============================================================
# IEEE 754 Section 9.2.1: Inverse Trig Domain
# ============================================================
section "IEEE 754: Inverse Trig Domain Errors"

test_ieee "asin(0) = 0" "asin(0)" "0" 1
test_ieee "asin(0.5) = pi/6" "asin(0.5)" "0.5235987"
test_ieee "asin(1) = pi/2" "asin(1)" "1.5707963"
test_ieee "asin(2) = NaN (domain)" "asin(2)" "NaN"
test_ieee "asin(-2) = NaN (domain)" "asin(-2)" "NaN"
test_ieee "acos(0) = pi/2" "acos(0)" "1.5707963"
test_ieee "acos(0.5) = pi/3" "acos(0.5)" "1.047197"
test_ieee "acos(1) = 0" "acos(1)" "0"
test_ieee "acos(2) = NaN (domain)" "acos(2)" "NaN"

# ============================================================
# Complex Number Operations
# ============================================================
section "Complex Number Extensions"

subsection "Basic complex arithmetic"
test_ieee "i * i = -1" "i * i" "-1"
test_ieee "i^2 = -1" "i^2" "-1"
test_ieee "i^3 = -i" "i^3" "-i"
test_ieee "i^4 = 1" "i^4" "1"
test_ieee "(1+i)*(1-i) = 2" "(1+i)*(1-i)" "2"
test_ieee "(2+3i)+(4+5i) = 6+8i" "(2+3i)+(4+5i)" "6 + 8i"

subsection "Complex magnitude and conjugate"
test_ieee "|0| = 0" "abs(0)" "0"
test_ieee "|3+4i| = 5" "abs(3+4i)" "5"
test_ieee "|i| = 1" "abs(i)" "1"
test_ieee "|1+i| = sqrt(2)" "abs(1+i)" "1.41421"
test_ieee "conj(3+4i) = 3-4i" "conj(3+4i)" "3 - 4i"
test_ieee "re(3+4i) = 3" "re(3+4i)" "3"
test_ieee "im(3+4i) = 4" "im(3+4i)" "4"

subsection "Euler's identity"
test_ieee "e^(i*pi) + 1 ≈ 0" "e^(i*pi) + 1" "0"

# ============================================================
# Precision and Accuracy
# ============================================================
section "Precision and Accuracy"

subsection "Repeating decimals"
test_ieee "1/3 precision" "1/3" "0.333333"
test_ieee "2/3 precision" "2/3" "0.666666"
test_ieee "1/7 precision" "1/7" "0.142857"
test_ieee "1/9 precision" "1/9" "0.111111"

subsection "Near-zero results"
test_ieee "0.1 + 0.2 ≈ 0.3" "0.1 + 0.2" "0.3"

subsection "Large numbers"
test_ieee "2^100 scientific" "2^100" "e+30"
test_ieee "2^1000 scientific" "2^1000" "e+301"
test_ieee "10^100 = 1e100" "10^100" "1e+100"

subsection "Small numbers"
test_ieee "1e-100" "1e-100" "1e-100"
test_ieee "1e-200 * 1e-200 underflow" "1e-200 * 1e-200" "0"

subsection "Integer precision"
test_ieee "20! exact" "20!" "2432902008176640000"
test_ieee "2^53 (IEEE double max int)" "2^53" "9007199254740992"

# ============================================================
# Inverse function identities
# ============================================================
section "Inverse Function Round-Trip"

test_ieee "asin(sin(0.5)) = 0.5" "asin(sin(0.5))" "0.5"
test_ieee "acos(cos(0.5)) = 0.5" "acos(cos(0.5))" "0.5"
test_ieee "atan(tan(0.5)) = 0.5" "atan(tan(0.5))" "0.5"

# Summary
echo ""
echo "====================================================="
TOTAL=$((PASS + FAIL))
PERCENT=$(echo "scale=1; $PASS * 100 / $TOTAL" | bc)
echo -e "  RESULTS: ${GREEN}$PASS passed${NC}, ${RED}$FAIL failed${NC}"
echo "  IEEE 754 Compliance: $PERCENT%"
echo "====================================================="

exit $FAIL
