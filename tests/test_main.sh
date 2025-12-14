#!/bin/bash
# test_main.sh - Main test runner for sc
# Runs comprehensive test suite

SC="${1:-../bin/sc}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=============================================="
echo "  sc - Scientific Calculator Test Suite"
echo "=============================================="
echo ""
echo "Calculator: $SC"
echo ""

# Colors
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    CYAN='\033[0;36m'
    NC='\033[0m'
else
    GREEN=''
    RED=''
    CYAN=''
    NC=''
fi

TOTAL_PASS=0
TOTAL_FAIL=0

# Function to run a test suite and capture results
run_suite() {
    local name="$1"
    local script="$2"
    local output
    local pass fail
    
    echo -e "${CYAN}Running: $name${NC}"
    output=$("$DIR/$script" "$SC" 2>&1)
    
    # Extract pass/fail counts - look for "N passed" pattern
    pass=$(echo "$output" | grep -oE '[0-9]+ passed' | tail -1 | grep -oE '[0-9]+')
    fail=$(echo "$output" | grep -oE '[0-9]+ failed' | tail -1 | grep -oE '[0-9]+')
    
    if [ -z "$pass" ]; then pass=0; fi
    if [ -z "$fail" ]; then fail=0; fi
    
    TOTAL_PASS=$((TOTAL_PASS + pass))
    TOTAL_FAIL=$((TOTAL_FAIL + fail))
    
    if [ "$fail" -eq 0 ]; then
        echo -e "  ${GREEN}✓ $pass tests passed${NC}"
    else
        echo -e "  ${RED}✗ $fail failures${NC} ($pass passed)"
        # Show failed tests
        echo "$output" | grep -E "^\s*FAIL" | head -5
    fi
    echo ""
}

# Run all test suites
run_suite "BEDMAS Order of Operations" "test_bedmas.sh"
run_suite "Core Functionality" "test_core.sh"
run_suite "IEEE 754 Compliance" "test_ieee754.sh"
run_suite "BC Comparison" "test_bc_comparison.sh"

# Summary
echo "=============================================="
if [ $TOTAL_FAIL -eq 0 ]; then
    echo -e "  ${GREEN}ALL TESTS PASSED${NC}"
else
    echo -e "  ${RED}SOME TESTS FAILED${NC}"
fi
echo "  Total: $TOTAL_PASS passed, $TOTAL_FAIL failed"
echo "=============================================="

# Exit with failure if any tests failed
[ $TOTAL_FAIL -eq 0 ]
