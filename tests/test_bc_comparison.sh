#!/bin/bash
# test_bc_comparison.sh - Compare sc against bc
# Tests arithmetic, functions, precision, and edge cases

SC="${1:-../bin/sc}"
BC="bc -l"  # -l loads math library (sin, cos, etc.)

PASS=0
FAIL=0
SKIP=0
TOTAL=0

# Colors
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    CYAN='\033[0;36m'
    NC='\033[0m'
else
    GREEN=''
    RED=''
    YELLOW=''
    CYAN=''
    NC=''
fi

# Get result from sc (quiet pipe mode - no "= " prefix)
run_sc() {
    echo "$1" | $SC 2>/dev/null | tail -1
}

# Get result from bc
run_bc() {
    echo "scale=16; $1" | $BC 2>/dev/null | tr -d '\\\n' | head -1
}

# Compare two floating point numbers (within tolerance)
float_eq() {
    local a="$1"
    local b="$2"
    local tol="${3:-0.0000000001}"
    
    # Handle special cases
    [[ "$a" == "$b" ]] && return 0
    [[ -z "$a" || -z "$b" ]] && return 1
    
    # Handle infinity/nan
    [[ "$a" == "Inf" && "$b" == "Inf" ]] && return 0
    [[ "$a" == "-Inf" && "$b" == "-Inf" ]] && return 0
    [[ "$a" == "NaN" && "$b" == "NaN" ]] && return 0
    
    # Numeric comparison using awk
    awk -v a="$a" -v b="$b" -v tol="$tol" 'BEGIN {
        if (a == "" || b == "") exit 1
        diff = a - b
        if (diff < 0) diff = -diff
        mag = (a < 0 ? -a : a)
        if (mag < 1) mag = 1
        if (diff / mag < tol) exit 0
        exit 1
    }' 2>/dev/null
}

# Test and compare
test_compare() {
    local sc_expr="$1"
    local bc_expr="${2:-$1}"
    local desc="${3:-$sc_expr}"
    local tol="${4:-0.0000000001}"
    
    ((TOTAL++))
    
    local sc_result=$(run_sc "$sc_expr")
    local bc_result=$(run_bc "$bc_expr")
    
    # Clean up bc result (remove trailing zeros, leading zeros on decimals)
    bc_result=$(echo "$bc_result" | sed 's/^\./0./' | sed 's/^-\./-0./')
    
    if float_eq "$sc_result" "$bc_result" "$tol"; then
        echo -e "  ${GREEN}PASS${NC}: $desc"
        echo "        sc: $sc_result"
        echo "        bc:    $bc_result"
        ((PASS++))
    else
        echo -e "  ${RED}FAIL${NC}: $desc"
        echo "        sc: $sc_result"
        echo "        bc:    $bc_result"
        ((FAIL++))
    fi
}

# Test sc-only (bc doesn't support)
test_sc_only() {
    local expr="$1"
    local expected="$2"
    local desc="${3:-$expr}"
    
    ((TOTAL++))
    
    local result=$(run_sc "$expr")
    
    if float_eq "$result" "$expected" "0.0000001"; then
        echo -e "  ${GREEN}PASS${NC}: $desc (sc-only)"
        echo "        result: $result"
        ((PASS++))
    else
        echo -e "  ${RED}FAIL${NC}: $desc (sc-only)"
        echo "        got:      $result"
        echo "        expected: $expected"
        ((FAIL++))
    fi
}

# Skip test (feature not comparable)
test_skip() {
    local desc="$1"
    local reason="$2"
    ((TOTAL++))
    ((SKIP++))
    echo -e "  ${YELLOW}SKIP${NC}: $desc - $reason"
}

echo "=============================================="
echo "  sc vs bc Comparison Test Suite"
echo "=============================================="
echo ""

# ============================================
echo -e "${CYAN}=== Basic Arithmetic ===${NC}"
# ============================================

test_compare "2+3" "2+3" "addition: 2+3"
test_compare "100-37" "100-37" "subtraction: 100-37"
test_compare "6*7" "6*7" "multiplication: 6*7"
test_compare "355/113" "355/113" "division: 355/113"
test_compare "2^10" "2^10" "power: 2^10"
test_compare "2^0.5" "e(0.5*l(2))" "power: 2^0.5 (sqrt 2)"
test_compare "10^3" "10^3" "power: 10^3"

echo ""
echo -e "${CYAN}=== Negative Numbers ===${NC}"

test_compare "-5+3" "-5+3" "negative: -5+3"
test_compare "5*(-3)" "5*(-3)" "negative: 5*(-3)"
test_compare "-10/-2" "-10/-2" "negative: -10/-2"
test_compare "(-2)^3" "(-2)^3" "negative power: (-2)^3"
test_compare "(-2)^4" "(-2)^4" "negative power: (-2)^4"

echo ""
echo -e "${CYAN}=== Order of Operations ===${NC}"

test_compare "2+3*4" "2+3*4" "precedence: 2+3*4"
test_compare "(2+3)*4" "(2+3)*4" "parens: (2+3)*4"
test_compare "2^3^2" "2^(3^2)" "right-assoc power: 2^3^2"
test_compare "100/10/2" "100/10/2" "left-assoc division: 100/10/2"
test_compare "2+3*4^2" "2+3*4^2" "mixed: 2+3*4^2"

echo ""
echo -e "${CYAN}=== Fractions & Decimals ===${NC}"

test_compare "1/3" "1/3" "fraction: 1/3"
test_compare "1/7" "1/7" "fraction: 1/7"
test_compare "22/7" "22/7" "fraction: 22/7"
test_compare "0.1+0.2" "0.1+0.2" "decimal: 0.1+0.2"
test_compare "1.5*2.5" "1.5*2.5" "decimal: 1.5*2.5"
test_compare "3.14159*2" "3.14159*2" "decimal: 3.14159*2"

echo ""
echo -e "${CYAN}=== Large Numbers ===${NC}"

test_compare "2^20" "2^20" "large: 2^20"
test_compare "2^30" "2^30" "large: 2^30"
test_compare "10^10" "10^10" "large: 10^10"
test_compare "123456789*987654321" "123456789*987654321" "large multiply"
test_compare "2^50" "2^50" "large: 2^50"

echo ""
echo -e "${CYAN}=== Small Numbers ===${NC}"

test_compare "1/1000000" "1/1000000" "small: 1/1000000"
test_compare "0.000001*0.000001" "0.000001*0.000001" "small multiply"
test_compare "1e-10+1e-10" "0.0000000001+0.0000000001" "small add"

echo ""
echo -e "${CYAN}=== Square Root ===${NC}"

test_compare "sqrt(2)" "sqrt(2)" "sqrt(2)"
test_compare "sqrt(3)" "sqrt(3)" "sqrt(3)"
test_compare "sqrt(4)" "sqrt(4)" "sqrt(4)"
test_compare "sqrt(100)" "sqrt(100)" "sqrt(100)"
test_compare "sqrt(2)*sqrt(2)" "sqrt(2)*sqrt(2)" "sqrt(2)*sqrt(2)"
test_compare "sqrt(0.5)" "sqrt(0.5)" "sqrt(0.5)"

echo ""
echo -e "${CYAN}=== Exponential & Logarithm ===${NC}"

test_compare "exp(1)" "e(1)" "exp(1) = e"
test_compare "exp(0)" "e(0)" "exp(0) = 1"
test_compare "exp(2)" "e(2)" "exp(2)"
test_compare "ln(1)" "l(1)" "ln(1) = 0"
test_compare "ln(e)" "l(e(1))" "ln(e) = 1"
# Note: sc's log() is natural log (same as ln), not log10
test_compare "ln(10)" "l(10)" "ln(10)"
test_compare "ln(100)" "l(100)" "ln(100)"
test_compare "exp(ln(5))" "e(l(5))" "exp(ln(5)) = 5"

echo ""
echo -e "${CYAN}=== Trigonometry (radians) ===${NC}"

test_compare "sin(0)" "s(0)" "sin(0)"
test_compare "cos(0)" "c(0)" "cos(0)"
test_compare "sin(1)" "s(1)" "sin(1)"
test_compare "cos(1)" "c(1)" "cos(1)"
test_compare "sin(0.5)" "s(0.5)" "sin(0.5)"
test_compare "cos(0.5)" "c(0.5)" "cos(0.5)"
test_compare "atan(1)" "a(1)" "atan(1) = pi/4"
test_compare "sin(3.14159265/2)" "s(3.14159265/2)" "sin(pi/2) ≈ 1"
test_compare "cos(3.14159265)" "c(3.14159265)" "cos(pi) ≈ -1"

echo ""
echo -e "${CYAN}=== Trigonometric Identities ===${NC}"

test_compare "sin(0.5)^2+cos(0.5)^2" "s(0.5)^2+c(0.5)^2" "sin²+cos² = 1"
test_compare "sin(1)^2+cos(1)^2" "s(1)^2+c(1)^2" "sin²+cos² = 1 (x=1)"

echo ""
echo -e "${CYAN}=== Constants ===${NC}"

# bc doesn't have pi directly, use 4*atan(1)
test_compare "pi" "4*a(1)" "pi"
test_compare "e" "e(1)" "e"
test_compare "pi*2" "4*a(1)*2" "2*pi"
test_compare "e^2" "e(1)^2" "e^2"

echo ""
echo -e "${CYAN}=== Complex Expressions ===${NC}"

test_compare "sqrt(2)+sqrt(3)" "sqrt(2)+sqrt(3)" "sqrt(2)+sqrt(3)"
test_compare "(1+sqrt(5))/2" "(1+sqrt(5))/2" "golden ratio"
# bc's ^ only works with integers, use e(x*l(y)) for y^x
test_compare "exp(1)^pi" "e((4*a(1))*l(e(1)))" "e^pi" "0.000001"
test_compare "pi^exp(1)" "e(e(1)*l(4*a(1)))" "pi^e" "0.000001"
test_compare "ln(2)*ln(2)" "l(2)*l(2)" "ln(2)^2"

echo ""
echo -e "${CYAN}=== Precision Tests ===${NC}"

test_compare "1/3*3" "1/3*3" "1/3*3 (should be ~1)"
test_compare "sqrt(2)^2" "sqrt(2)^2" "sqrt(2)^2 (should be 2)"
test_compare "(1/49)*49" "(1/49)*49" "(1/49)*49"
test_compare "sin(pi)" "s(4*a(1))" "sin(pi) ≈ 0" "0.0000001"

echo ""
echo -e "${CYAN}=== sc-only Features ===${NC}"

# Complex numbers (bc doesn't support)
test_sc_only "sqrt(-1)" "i" "sqrt(-1) = i"
test_sc_only "i*i" "-1" "i*i = -1"
test_sc_only "i^2" "-1" "i^2 = -1"
test_sc_only "(1+i)*(1-i)" "2" "(1+i)(1-i) = 2"
test_sc_only "abs(3+4i)" "5" "|3+4i| = 5"
test_sc_only "exp(i*pi)+1" "0" "Euler's identity" 

# Factorial
test_sc_only "5!" "120" "5! = 120"
test_sc_only "10!" "3628800" "10!"
test_sc_only "0!" "1" "0! = 1"

# Hyperbolic functions
test_sc_only "sinh(0)" "0" "sinh(0)"
test_sc_only "cosh(0)" "1" "cosh(0)"
test_sc_only "tanh(0)" "0" "tanh(0)"

# Hex and binary input
test_sc_only "0xFF" "255" "hex: 0xFF"
test_sc_only "0b1010" "10" "binary: 0b1010"
test_sc_only "0x100+0b100" "260" "hex + binary"

# Roman numerals
test_sc_only "MMXXV" "2025" "Roman: MMXXV"
test_sc_only "MCMXCIX" "1999" "Roman: MCMXCIX"

# Previous answer
test_sc_only "5+5" "10" "setup ans"
# Can't easily test ans in this framework

echo ""
echo -e "${CYAN}=== Edge Cases ===${NC}"

test_compare "0+0" "0+0" "zero: 0+0"
test_compare "0*1000" "0*1000" "zero: 0*1000"
test_compare "1^1000" "1^1000" "one: 1^1000"
test_compare "0^0" "0^0" "zero: 0^0 (convention: 1)"
# bc gives 1 for 0^0, sc should too

test_sc_only "1/0" "Inf" "division by zero"
test_sc_only "-1/0" "-Inf" "negative div by zero"
test_sc_only "0/0" "NaN" "0/0 = NaN"
test_sc_only "Inf+1" "Inf" "Inf+1"
test_sc_only "Inf-Inf" "NaN" "Inf-Inf = NaN"
test_sc_only "Inf*0" "NaN" "Inf*0 = NaN"

echo ""
echo "=============================================="
echo -e "  RESULTS: ${GREEN}$PASS passed${NC}, ${RED}$FAIL failed${NC}, ${YELLOW}$SKIP skipped${NC}"
echo "  Total: $TOTAL tests"
echo "=============================================="

# Exit with failure if any tests failed
[ $FAIL -eq 0 ]
