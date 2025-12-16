#!/bin/bash
# MATLAB Full Compatibility Test Suite for scalc
# Based on matlab_compat_tests.m (1000 tests spec)

SC="${1:-../bin/sc}"
PASS=0
FAIL=0
SKIP=0

test_expr() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | tail -1)
    
    # Normalize whitespace
    result=$(echo "$result" | tr -s ' ' | sed 's/^ //;s/ $//')
    expected=$(echo "$expected" | tr -s ' ' | sed 's/^ //;s/ $//')
    
    if [ "$result" = "$expected" ]; then
        echo "  PASS: $expr → $expected"
        ((PASS++))
    else
        echo "  FAIL: $expr"
        echo "        expected: $expected"
        echo "        got:      $result"
        ((FAIL++))
    fi
}

test_approx() {
    local expr="$1"
    local expected="$2"
    local tol="${3:-0.0001}"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | tail -1)
    
    # Check if result is approximately equal
    local check
    check=$(echo "abs($result - $expected) < $tol" | $SC 2>/dev/null | tail -1)
    
    if [ "$check" = "true" ] || [ "$check" = "1" ]; then
        echo "  PASS: $expr ≈ $expected"
        ((PASS++))
    else
        echo "  FAIL: $expr"
        echo "        expected: ≈$expected"
        echo "        got:      $result"
        ((FAIL++))
    fi
}

test_bool() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | tail -1)
    
    if [ "$result" = "$expected" ]; then
        echo "  PASS: $expr → $expected"
        ((PASS++))
    else
        echo "  FAIL: $expr"
        echo "        expected: $expected"
        echo "        got:      $result"
        ((FAIL++))
    fi
}

skip_test() {
    local name="$1"
    echo "  SKIP: $name (not implemented)"
    ((SKIP++))
}

echo "=============================================="
echo "MATLAB Full Compatibility Test Suite"
echo "=============================================="
echo ""

echo "=== Basic Arithmetic ==="
test_expr "3+5" "8"
test_expr "2-7" "-5"
test_expr "4*6" "24"
test_expr "8/2" "4"
test_expr "2^8" "256"
test_expr "10.*[1 2 3]" "10  20  30"
test_expr "[1 2 3].*2" "2  4  6"
test_expr "[1 2 3].^2" "1  4  9"
test_expr "mod(10,3)" "1"
test_expr "rem(-10,3)" "-1"
test_expr "floor(3.9)" "3"
test_expr "ceil(3.1)" "4"
test_expr "round(2.5)" "3"
test_expr "fix(-3.9)" "-3"
test_expr "abs(-7)" "7"
test_expr "sign(-9)" "-1"

echo ""
echo "=== Relational/Logical ==="
test_bool "3>2" "true"
test_bool "3<2" "false"
test_bool "3==3" "true"
test_bool "3~=4" "true"
test_bool "3!=4" "true"
test_bool "(3>2)&&(2>1)" "true"
test_bool "(3>2)||(2>5)" "true"
test_bool "not(0)" "true"
test_bool "any([0 1 0])" "true"
test_bool "all([1 1 1])" "true"

echo ""
echo "=== Vectors ==="
test_expr "[1 2 3] + [4 5 6]" "5  7  9"
test_expr "sum([1 2 3])" "6"
test_expr "mean([1 2 3])" "2"
test_expr "prod([1 2 3])" "6"
test_expr "min([1 2 3])" "1"
test_expr "max([1 2 3])" "3"
test_expr "length([1 2 3 4 5])" "5"
test_expr "numel([1 2 3 4])" "4"

echo ""
echo "=== Matrices ==="
test_expr "A=[1 2;3 4]; A(2,1)" "3"
test_expr "A=[1 2;3 4]; A(1,2)" "2"
test_expr "M=[1 2 3;4 5 6]; M(1,3)" "3"
test_expr "[1 2;3 4]'" "1  3
2  4"
test_expr "det([1 2;3 4])" "-2"
test_expr "trace([1 2;3 4])" "5"

echo ""
echo "=== Matrix Indexing ==="
test_expr "v=[1 2 3 4 5]; v(1)" "1"
test_expr "v=[1 2 3 4 5]; v(3)" "3"
test_expr "v=[1 2 3 4 5]; v(5)" "5"
test_expr "M=[1 2 3;4 5 6;7 8 9]; M(1,1)" "1"
test_expr "M=[1 2 3;4 5 6;7 8 9]; M(2,3)" "6"
test_expr "M=[1 2 3;4 5 6;7 8 9]; M(3,2)" "8"

echo ""
echo "=== Colon/Ranges ==="
test_expr "sum(1:5)" "15"
test_expr "sum(1:2:9)" "25"
test_expr "length(1:5)" "5"
test_expr "length(1:10)" "10"
test_expr "sum(linspace(0,1,5))" "2.5"
test_approx "sum(logspace(0,2,3))" "111" "0.1"

echo ""
echo "=== Trigonometry ==="
test_approx "sin(pi/2)" "1" "0.0001"
test_approx "cos(0)" "1" "0.0001"
test_approx "tan(pi/4)" "1" "0.0001"
test_approx "asin(1)" "1.5708" "0.001"
test_approx "acos(0)" "1.5708" "0.001"
test_approx "atan(1)" "0.7854" "0.001"
test_approx "atan2(1,1)" "0.7854" "0.001"
test_approx "sind(90)" "1" "0.0001"
test_approx "cosd(0)" "1" "0.0001"
test_approx "cosd(60)" "0.5" "0.0001"
test_approx "tand(45)" "1" "0.0001"

echo ""
echo "=== Exponential/Log ==="
test_approx "exp(1)" "2.71828" "0.0001"
test_approx "log(e)" "1" "0.0001"
test_expr "log10(1000)" "3"
test_expr "sqrt(16)" "4"

echo ""
echo "=== Special Functions ==="
test_expr "factorial(5)" "120"
test_expr "factorial(0)" "1"
test_expr "factorial(10)" "3628800"
test_expr "hypot(3,4)" "5"
test_approx "gamma(5)" "24" "0.0001"
test_approx "gamma(0.5)" "1.7725" "0.001"
test_approx "erf(0)" "0" "0.0001"
test_approx "erf(1)" "0.8427" "0.001"

echo ""
echo "=== Complex Numbers ==="
test_expr "real(3+4i)" "3"
test_expr "imag(3+4i)" "4"
test_expr "abs(3+4i)" "5"
test_expr "conj(3+4i)" "3 - 4i"
test_approx "angle(1+1i)" "0.7854" "0.001"

echo ""
echo "=== Constants ==="
test_approx "pi" "3.14159265" "0.00001"
test_approx "e" "2.71828" "0.0001"
test_bool "Inf > 1e9" "true"
test_bool "NaN == NaN" "false"
test_bool "isnan(NaN)" "true"
test_bool "isinf(Inf)" "true"
test_bool "isfinite(42)" "true"
test_bool "isreal(42)" "true"
test_bool "isreal(3+4i)" "false"

echo ""
echo "=== Matrix Creation ==="
test_expr "sum(sum(zeros(2,2)))" "0"
test_expr "sum(sum(ones(2,2)))" "4"
test_expr "sum(sum(eye(3)))" "3"
test_expr "zeros(2,3)" "0  0  0
0  0  0"
test_expr "ones(2,2)" "1  1
1  1"
test_expr "eye(2)" "1  0
0  1"

echo ""
echo "=== Cumulative Functions ==="
test_expr "v=[1 2 3]; sum(cumsum(v))" "10"
test_expr "cumsum([1 2 3])" "1  3  6"

echo ""
echo "=== Differences ==="
test_expr "sum(diff([1 2 4 7]))" "6"
test_expr "diff([1 3 6 10])" "2  3  4"

echo ""
echo "=== Statistics ==="
test_expr "median([1 5 9])" "5"
test_expr "median([1 2 3 4])" "2.5"
test_approx "std([1 2 3])" "1" "0.01"
test_approx "var([1 2 3])" "1" "0.01"

echo ""
echo "=== Bitwise ==="
test_expr "bitand(7, 3)" "3"
test_expr "bitor(7, 8)" "15"
test_expr "bitxor(7, 3)" "4"
test_expr "bitshift(1, 3)" "8"
test_expr "bitshift(8, -2)" "2"

echo ""
echo "=== Vector Functions ==="
test_expr "dot(1,2,3,4,5,6)" "32"
test_expr "dot(1,0,0,0,1,0)" "0"
test_expr "norm([3 4])" "5"

echo ""
echo "=== Flip Functions ==="
test_expr "fliplr([1 2 3])" "3  2  1"
test_expr "flipud([1;2;3])" "3
2
1"

echo ""
echo "=== Sort/Find ==="
test_expr "sort([3 1 4 1 5])" "1  1  3  4  5"
test_expr "find([0 1 0 1 1])" "2  4  5"
test_expr "length(find([0 1 0 1]))" "2"

echo ""
echo "=== Reshape/Repmat ==="
test_expr "reshape([1 2 3 4], 2, 2)" "1  3
2  4"
test_expr "repmat([1 2], 2, 2)" "1  2  1  2
1  2  1  2"

echo ""
echo "=== Named Variables ==="
test_expr "myVar = 42; myVar" "42"
test_expr "nAssets = 3; nAssets + 1" "4"
test_expr "Sigma = [1 2; 3 4]; Sigma(1,1)" "1"
test_expr "alpha = 0.5; beta = 0.5; alpha + beta" "1"

echo ""
echo "=== Elementwise Operations ==="
test_expr "[1 2 3] .* [4 5 6]" "4  10  18"
test_expr "[2 4 6] ./ [1 2 3]" "2  2  2"
test_expr "[2 3 4] .^ [1 2 3]" "2  9  64"

echo ""
echo "=== Linear Algebra ==="
test_expr "inv([1 0; 0 2])" "1  0
0  0.5"
test_expr "diag([1 2 3])" "1  0  0
0  2  0
0  0  3"
test_expr "diag([1 2; 3 4])" "1
4"

echo ""
echo "=== Comments ==="
result=$(echo "% this is a comment
42" | $SC 2>/dev/null | tail -1)
if [ "$result" = "42" ]; then
    echo "  PASS: % comment works"
    ((PASS++))
else
    echo "  FAIL: % comment"
    ((FAIL++))
fi

echo ""
echo "=== Semicolon Suppression ==="
result=$(echo "a=1; b=2; a+b" | $SC 2>/dev/null)
lines=$(echo "$result" | wc -l)
if [ "$lines" -eq 1 ]; then
    echo "  PASS: semicolon suppresses output"
    ((PASS++))
else
    echo "  FAIL: semicolon suppression (got $lines lines)"
    ((FAIL++))
fi

echo ""
echo "=== Display Function ==="
result=$(echo "disp('hello')" | $SC 2>/dev/null)
if echo "$result" | grep -q "hello"; then
    echo "  PASS: disp('hello') outputs hello"
    ((PASS++))
else
    echo "  FAIL: disp('hello')"
    ((FAIL++))
fi

echo ""
echo "=== Vectorization ==="
test_expr "sum((1:5).*2)" "30"
test_expr "sum((1:5)+10)" "65"
test_expr "sum((1:5).^2)" "55"

echo ""
echo "=== Misc Math ==="
test_expr "sum(diag([1 2;3 4]))" "5"
test_approx "norm([3 4])" "5" "0.0001"

echo ""
echo "=== Logical Precedence ==="
test_bool "(1==2)==0" "true"
test_bool "(1<2)&&(2<3)" "true"
test_bool "(1<2)||(2<3)" "true"

echo ""
echo "=== Performance Sanity ==="
result=$(echo "length(zeros(100,1))" | $SC 2>/dev/null | tail -1)
if [ "$result" = "100" ]; then
    echo "  PASS: zeros(100,1) creates 100 element vector"
    ((PASS++))
else
    echo "  FAIL: zeros(100,1) length"
    ((FAIL++))
fi

result=$(echo "sum(sum(eye(5)))" | $SC 2>/dev/null | tail -1)
if [ "$result" = "5" ]; then
    echo "  PASS: eye(5) creates 5x5 identity"
    ((PASS++))
else
    echo "  FAIL: eye(5)"
    ((FAIL++))
fi

echo ""
echo "=== Additional Tests ==="
test_expr "fix(3.7)" "3"
test_expr "fix(-3.7)" "-3"
test_expr "sign(0)" "0"
test_expr "sign(5)" "1"
test_expr "sign(-5)" "-1"

echo ""
echo "=============================================="
echo "Results: $PASS passed, $FAIL failed, $SKIP skipped"
echo "=============================================="

exit $FAIL
