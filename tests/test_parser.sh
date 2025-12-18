#!/bin/bash
# test_parser.sh - Parser edge cases and syntax tests
# Tests the parser with various complex and edge case expressions

SC="${SC:-./bin/sc}"
pass=0
fail=0

test_expr() {
    local expected="$1"
    shift
    local expr="$*"
    
    local result=$(echo -e "$expr\nquit" | $SC 2>&1)
    
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ $expr"
        ((pass++))
    else
        echo "  ✗ $expr"
        echo "    Expected: $expected"
        echo "    Got: $(echo "$result" | tail -3)"
        ((fail++))
    fi
}

# Test that should produce an error
test_error() {
    local expected_err="$1"
    shift
    local expr="$*"
    
    local result=$(echo -e "$expr\nquit" | $SC 2>&1)
    
    if echo "$result" | grep -qi "error\|Error\|ERROR"; then
        echo "  ✓ $expr (correctly errored)"
        ((pass++))
    else
        echo "  ✗ $expr (should have errored)"
        echo "    Got: $(echo "$result" | tail -3)"
        ((fail++))
    fi
}

echo "Parser Edge Cases and Syntax Tests"
echo "==================================="
echo ""

echo "Basic Arithmetic:"
test_expr "5" "2+3"
test_expr "6" "2*3"
test_expr "8" "2^3"
test_expr "-1" "2-3"
test_expr "0.5" "1/2"
echo ""

echo "Unary Minus:"
test_expr "-5" "-5"
test_expr "-2" "-5+3"
test_expr "8" "3--5"
test_expr "-10" "-5*2"
test_expr "-4" "-2^2"
test_expr "4" "(-2)^2"
test_expr "5" "--5"
test_expr "-5" "---5"
echo ""

echo "Operator Precedence:"
test_expr "14" "2+3*4"
test_expr "20" "(2+3)*4"
test_expr "512" "2^3^2"
test_expr "64" "(2^3)^2"
test_expr "11" "2+3^2"
test_expr "25" "(2+3)^2"
echo ""

echo "Nested Parentheses:"
test_expr "5" "((2+3))"
test_expr "5" "(((2+3)))"
test_expr "45" "((2+3)*(4+5))"
test_expr "58" "2*((3+(4*5))+6)"
echo ""

echo "Function Calls:"
test_expr "2" "sqrt(4)"
test_expr "1" "sin(pi/2)"
test_expr "3628800" "fact(10)"
test_expr "120" "5!"
test_expr "2" "abs(-2)"
test_expr "3" "max([1,2,3])"
test_expr "1" "min([1,2,3])"
echo ""

echo "Nested Function Calls:"
test_expr "2" "sqrt(sqrt(16))"
test_expr "0" "abs(sin(pi))"
test_expr "5" "sqrt(abs(-25))"
test_expr "2" "log10(sqrt(10000))"
echo ""

echo "Constants:"
test_expr "3.141592653589793" "pi"
test_expr "2.718281828459045" "e"
test_expr "1.618033988749895" "golden"
test_expr "Inf" "Inf"
test_expr "NaN" "NaN"
echo ""

echo "Matrix Literals:"
test_expr "1  2" "[1,2]"
test_expr "1  2" "[1,2;3,4]"
test_expr "3  4" "[1,2;3,4]"
test_expr "2" "rows([1,2;3,4])"
test_expr "2" "cols([1,2;3,4])"
echo ""

echo "Matrix Indexing (via variable):"
test_expr "1" "A = [1,2;3,4]\nA(1,1)"
test_expr "4" "A = [1,2;3,4]\nA(2,2)"
test_expr "1  2" "A = [1,2;3,4]\nA(1,:)"
test_expr "1" "A = [1,2;3,4]\nA(:,1)"
echo ""

echo "Range Expressions:"
test_expr "1" "1:5"
test_expr "5" "1:5"
test_expr "1" "1:2:10"
test_expr "9" "1:2:10"
test_expr "10" "10:-1:1"
test_expr "1" "10:-1:1"
echo ""

echo "Variable Assignment:"
test_expr "5" "x = 5"
test_expr "10" "x = 5\nx * 2"
test_expr "25" "x = 5\ny = x^2\ny"
echo ""

echo "Multi-Statement (Semicolon):"
test_expr "10" "5;10"
test_expr "6" "x=2;y=3;x*y"
test_expr "1" "2+3;ans/5"
echo ""

echo "Comparison Operators:"
test_expr "true" "5>3"
test_expr "false" "3>5"
test_expr "true" "5>=5"
test_expr "true" "3<5"
test_expr "true" "5<=5"
test_expr "true" "5==5"
test_expr "false" "5==6"
test_expr "true" "5~=6"
test_expr "true" "5<>6"
echo ""

echo "Logical Operators:"
test_expr "true" "true && true"
test_expr "false" "true && false"
test_expr "true" "true || false"
test_expr "false" "false || false"
test_expr "false" "not(true)"
test_expr "true" "not(false)"
echo ""

echo "Element-wise Operations:"
test_expr "1  4" "[1,2] .* [1,2]"
test_expr "1  4" "[1,2] .^ 2"
test_expr "1  1" "[2,4] ./ [2,4]"
echo ""

echo "String Arguments (if supported):"
test_expr "150" "load fisheriris\nrows(meas)"
echo ""

echo "Complex Numbers:"
test_expr "1" "real(1+2i)"
test_expr "2" "imag(1+2i)"
test_expr "2.236067977" "abs(1+2i)"
test_expr "1 - 2i" "conj(1+2i)"
test_expr "5" "abs(3+4i)"
echo ""

echo "Special Values:"
test_expr "Inf" "1/0"
test_expr "-Inf" "-1/0"
test_expr "NaN" "0/0"
test_expr "true" "isnan(NaN)"
test_expr "true" "isinf(Inf)"
test_expr "false" "isfinite(Inf)"
echo ""

echo "Unicode Operators:"
test_expr "12" "3 × 4"
test_expr "3" "12 ÷ 4"
echo ""

echo "Comments:"
test_expr "8" "5 + 3 % this is a comment"
test_expr "10" "10 % everything after percent is ignored"
echo ""

echo "Edge Cases - Empty/Whitespace:"
test_expr "" ""
test_expr "" "   "
echo ""

echo "Edge Cases - Large Numbers:"
test_expr "1000000000000" "10^12"
test_expr "1e+100" "10^100"
echo ""

echo "Edge Cases - Small Numbers:"
test_expr "0.001" "10^-3"
test_expr "1e-10" "10^-10"
echo ""

echo "Multi-output Functions:"
test_expr "0.316227" "[Q,R] = qr([1,2;3,4])\nQ"
test_expr "5.464985" "[U,S,V] = svd([1,2;3,4])\nS"
test_expr "150" "load fisheriris\n[idx,C] = kmeans(meas, 3)\nrows(idx)"
echo ""

echo "Error Handling:"
test_error "expected '('" "undefined_variable"
test_error "error" "[1,2] + [1,2,3]"
echo ""

echo "===================================="
echo "Results: $pass passed, $fail failed"

if [ $fail -gt 0 ]; then
    exit 1
fi
exit 0
