#!/bin/bash
# Test script for Session 2 new features

SC="${1:-../bin/sc}"
PASS=0
FAIL=0

test_expr() {
    local expr="$1"
    local expected="$2"
    local result
    result=$(echo "$expr" | $SC 2>&1 | tail -1 | tr -d ' ')
    expected_clean=$(echo "$expected" | tr -d ' ')
    
    if [ "$result" = "$expected_clean" ]; then
        echo "  PASS: $expr = $expected"
        ((PASS++))
    else
        echo "  FAIL: $expr"
        echo "    Expected: $expected"
        echo "    Got:      $result"
        ((FAIL++))
    fi
}

echo "=== Testing Session 2 Features ==="
echo ""

echo "Testing cummax/cummin/rescale..."
test_expr "cummax([3 1 4])" "3  3  4"
test_expr "cummin([3 1 4])" "3  1  1"
test_expr "rescale([0 5 10])" "0  0.5  1"

echo "Testing eq/ne/lt/le/gt/ge..."
test_expr "eq([1 2], [1 3])" "1  0"
test_expr "ne([1 2], [1 3])" "0  1"
test_expr "lt([1 2], [2 2])" "1  0"
test_expr "le([1 2], [2 2])" "1  1"
test_expr "gt([1 2], [0 1])" "1  1"
test_expr "ge([1 2], [1 3])" "1  0"

echo "Testing times/rdivide/ldivide/plus/minus..."
test_expr "times([2 3], [4 5])" "8  15"
test_expr "rdivide([8 12], [2 3])" "4  4"
test_expr "plus([1 2], [3 4])" "4  6"
test_expr "minus([5 6], [1 2])" "4  4"

echo "Testing conv..."
test_expr "conv([1 2], [1 1])" "1  3  2"
test_expr "conv([1 2 3], [1 1])" "1  3  5  3"

echo "Testing corrcoef..."
test_expr "corrcoef([1 2 3], [2 4 6])" "1"
test_expr "corrcoef([1 2 3], [3 2 1])" "-1"

echo "Testing movmean/movsum..."
# movmean([1 2 3 4 5], 3) should give centered moving average
# movsum([1 2 3], 2) should give [1+2, 2+3] style sums

echo "Testing kron..."
test_expr "kron([1 2], [1 1])" "1  1  2  2"
test_expr "kron(2, 3)" "6"

echo "Testing trapz..."
test_expr "trapz([0 1])" "0.5"
test_expr "trapz([0 1 2])" "2"

echo "Testing gradient..."
test_expr "gradient([1 2 3])" "1  1  1"

echo "Testing vander..."
# vander([1 2], 2) = [1 1; 2 1]

echo "Testing toeplitz..."
test_expr "toeplitz([1])" "1"

echo "Testing hankel..."
test_expr "hankel([1])" "1"

echo ""
echo "========================================"
echo "Results: $PASS passed, $FAIL failed"
echo "========================================"

exit $FAIL

echo "Testing blkdiag..."
test_expr "sum(sum(blkdiag([1], [1])))" "2"

echo "Testing sub2ind/ind2sub..."
test_expr "sub2ind([3, 4], 2, 3)" "8"
test_expr "ind2sub([3, 4], 8)" "2  3"

echo "Testing clamp..."
test_expr "clamp(5, 0, 10)" "5"
test_expr "clamp(-5, 0, 10)" "0"
test_expr "clamp(15, 0, 10)" "10"

echo "Testing isapprox..."
test_expr "isapprox(1, 1)" "1"
test_expr "isapprox(1, 2)" "0"

echo "Testing colon..."
test_expr "colon(1, 3)" "1  2  3"
test_expr "colon(1, 2, 5)" "1  3  5"

echo "Testing cat..."
test_expr "cat(1, [1], [2])" "1"
test_expr "cat(2, [1], [2])" "1  2"

echo "Testing compan..."
# compan([1 -1]) should give a 1x1 matrix with value 1
test_expr "compan([1, 0])" "0"
