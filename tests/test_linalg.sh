#!/bin/bash
# test_linalg.sh - Linear Algebra tests
# Tests matrix operations, decompositions, and solvers

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

echo "Linear Algebra Tests"
echo "===================="
echo ""

echo "Matrix Construction:"
test_expr "2" "rows([1,2;3,4])"
test_expr "3" "cols([1,2,3;4,5,6])"
test_expr "3" "rows(eye(3))"
test_expr "3" "rows(zeros(3,3))"
test_expr "3" "rows(ones(3,3))"
test_expr "3" "rows(magic(3))"
echo ""

echo "Matrix Arithmetic:"
test_expr "5" "[1,2;3,4] + [1,1;1,1]"
test_expr "6" "[1,2;3,4] .* [2,2;2,2]"
test_expr "7" "[1,2;3,4] * [1;1]"
echo ""

echo "Determinant:"
test_expr "-2" "det([1,2;3,4])"
test_expr "-360" "det(magic(3))"
test_expr "1" "det(eye(5))"
test_expr "0" "det(zeros(3,3))"
echo ""

echo "Trace:"
test_expr "5" "trace([1,2;3,4])"
test_expr "15" "trace(magic(3))"
test_expr "5" "trace(eye(5))"
echo ""

echo "Matrix Inverse:"
test_expr "-2" "inv([1,2;3,4])"
test_expr "1" "inv(eye(3))"
test_expr "0.5" "inv([2,0;0,2])"
echo ""

echo "Transpose:"
test_expr "1  3" "transpose([1,2;3,4])"
test_expr "2  4" "transpose([1,2;3,4])"
echo ""

echo "Matrix Norms:"
test_expr "5.477225575" "norm([1,2;3,4])"
test_expr "1.732050807" "norm(eye(3))"
echo ""

echo "Eigenvalues (2x2):"
test_expr "5.372281" "eig([1,2;3,4])"
test_expr "-0.372281" "eig([1,2;3,4])"
echo ""

echo "QR Decomposition:"
test_expr "0.316227" "[Q,R] = qr([1,2;3,4])\nQ"
test_expr "1  0  0" "[Q,R] = qr(eye(3))\nQ"
echo ""

echo "SVD:"
test_expr "5.464985" "svd([1,2;3,4])"
test_expr "1" "svd(eye(3))"
echo ""

echo "LU Decomposition:"
test_expr "1  0" "[L,U] = lu([1,2;3,4])\nL"
echo ""

echo "Linear Solve:"
test_expr "1" "linsolve([2,0;0,2], [2;4])"
test_expr "2" "linsolve([2,0;0,2], [2;4])"
test_expr "0.05" "linsolve(magic(3), [1;2;3])"
echo ""

echo "Rank:"
test_expr "2" "rank([1,2;3,4])"
test_expr "3" "rank(magic(3))"
test_expr "1" "rank([1,1;2,2])"
echo ""

echo "Diagonal Operations:"
test_expr "1" "diag([1,2;3,4])"
test_expr "4" "diag([1,2;3,4])"
test_expr "1  0" "diag([1,2,3])"
echo ""

echo "Special Matrices:"
test_expr "45" "sum(magic(3))"
test_expr "1  0  0" "eye(3)"
echo ""

echo "===================================="
echo "Results: $pass passed, $fail failed"

if [ $fail -gt 0 ]; then
    exit 1
fi
exit 0
