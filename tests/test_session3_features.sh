#!/bin/bash
# Test session 3 features

SC="${1:-../bin/sc}"

PASS=0
FAIL=0

test_expr() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>/dev/null | tr -d '\n' | sed 's/^ *//;s/ *$//')
    if [ "$result" = "$expected" ]; then
        echo "  PASS: $expr = $expected"
        ((PASS++))
    else
        echo "  FAIL: $expr = $result (expected $expected)"
        ((FAIL++))
    fi
}

echo "=== Testing Session 3 Features ==="
echo

echo "Testing signal functions..."
test_expr "sinc(0)" "1"
test_expr "sigmoid(0)" "0.5"
test_expr "step(1)" "1"
test_expr "step(-1)" "0"
test_expr "rect(0)" "1"
test_expr "tri(0)" "1"
test_expr "softplus(0)" "0.6931471805599453"

echo "Testing wrap functions..."
test_expr "wrap360(370)" "10"
test_expr "wrap180(270)" "-90"

echo "Testing log1p/expm1..."
test_expr "round(log1p(0.001)*1000000)" "1000"
test_expr "round(expm1(0.001)*1000000)" "1001"

echo "Testing advanced statistics..."
test_expr "round(geomean([1, 2, 4, 8])*100)" "283"
test_expr "round(harmmean([1, 2, 4])*100)" "171"
test_expr "skewness([1, 2, 3, 4, 5])" "0"
test_expr "mad([1, 2, 3, 4, 5])" "1.2"
test_expr "iqr([1, 2, 3, 4, 5, 6, 7, 8])" "3.5"

echo "Testing interp1..."
test_expr "interp1([0, 1, 2], [0, 1, 4], 1.5)" "2.5"
test_expr "interp1([0, 1], [0, 10], 0.5)" "5"

echo "Testing conv/deconv..."
test_expr "conv([1, 2], [1, 1])" "1  3  2"
test_expr "deconv([1, 3, 2], [1, 1])" "1  2"

echo
echo "========================================"
echo "Results: $PASS passed, $FAIL failed"
echo "========================================"
