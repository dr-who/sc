#!/bin/bash
# test_cli.sh - Command-line interface tests
# Tests command-line usage of sc

SC="${SC:-./bin/sc}"
pass=0
fail=0

# Test command-line expression (with eval for proper quote handling)
test_cli() {
    local expected="$1"
    shift
    local args="$*"
    
    local result=$(eval "$SC $args" 2>&1)
    
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ sc $args"
        ((pass++))
    else
        echo "  ✗ sc $args"
        echo "    Expected: $expected"
        echo "    Got: $(echo "$result" | head -3)"
        ((fail++))
    fi
}

echo "Command-Line Interface Tests"
echo "============================="
echo ""

echo "Basic Expressions:"
test_cli "3" "1+2"
test_cli "8388608" "2^23"
test_cli "false" "1==2"
test_cli "true" "2==2"
test_cli "3.141592653589793" "pi"
test_cli "2.718281828459045" "e"
echo ""

echo "Math Functions (use double quotes for parentheses):"
test_cli "3628800" '"fact(10)"'
test_cli "1.414213562373095" '"sqrt(2)"'
test_cli "2" '"1+1"'
test_cli "120" '"5!"'
echo ""

echo "Demo Commands:"
test_cli "svd" "demo svd"
test_cli "Singular value" "demo svd"
test_cli "QR decomposition" "demo qr"
test_cli "eigenvalues" "demo eig"
echo ""

echo "Help Commands:"
test_cli "qr - Linear Algebra" "help qr"
test_cli "svd - Linear Algebra" "help svd"
test_cli "kmeans - Data Science" "help kmeans"
echo ""

echo "Info Commands:"
test_cli "Built-in Functions" "functions"
test_cli "Feature Flags" "features"
test_cli "Mathematical Constants" "constants"
echo ""

echo "Boolean Results:"
test_cli "true" '"5>3"'
test_cli "false" '"3>5"'
test_cli "true" '"abs(-5)==5"'
echo ""

echo "Matrix Operations:"
test_cli "3" '"rows([1,2,3;4,5,6;7,8,9])"'
test_cli "-2" '"det([1,2;3,4])"'
echo ""

echo "===================================="
echo "Results: $pass passed, $fail failed"

if [ $fail -gt 0 ]; then
    exit 1
fi
exit 0
