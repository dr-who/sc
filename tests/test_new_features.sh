#!/bin/bash
# Test 40 new MATLAB-compatible features
# At least 5 tests per feature

SC="${1:-../bin/sc}"
PASS=0
FAIL=0

test_expr() {
    local expr="$1"
    local expected="$2"
    local result
    result=$("$SC" "$expr" 2>&1)
    # Use grep -F for fixed string matching, and -- to stop option processing
    if echo "$result" | grep -F -- "$expected" >/dev/null 2>&1; then
        ((PASS++))
    else
        echo "FAIL: $expr"
        echo "  Expected: $expected"
        echo "  Got: $result"
        ((FAIL++))
    fi
}

echo "=== Testing 40 New MATLAB-Compatible Features ==="
echo ""

# 1. isrow (5 tests)
echo "Testing isrow..."
test_expr "isrow([1 2 3])" "true"
test_expr "isrow([1; 2; 3])" "false"
test_expr "isrow([1 2; 3 4])" "false"
test_expr "isrow([5])" "true"
test_expr "isrow(5)" "false"

# 2. iscolumn (5 tests)
echo "Testing iscolumn..."
test_expr "iscolumn([1; 2; 3])" "true"
test_expr "iscolumn([1 2 3])" "false"
test_expr "iscolumn([1 2; 3 4])" "false"
test_expr "iscolumn([5])" "true"
test_expr "iscolumn(5)" "false"

# 3. issquare (5 tests)
echo "Testing issquare..."
test_expr "issquare([1 2; 3 4])" "true"
test_expr "issquare([1 2 3])" "false"
test_expr "issquare([1; 2; 3])" "false"
test_expr "issquare(eye(3))" "true"
test_expr "issquare(5)" "true"

# 4. issymmetric (5 tests)
echo "Testing issymmetric..."
test_expr "issymmetric([1 2; 2 1])" "true"
test_expr "issymmetric([1 2; 3 4])" "false"
test_expr "issymmetric(eye(3))" "true"
test_expr "issymmetric([1])" "true"
test_expr "issymmetric(5)" "true"

# 5. isdiag (5 tests)
echo "Testing isdiag..."
test_expr "isdiag(eye(3))" "true"
test_expr "isdiag([1 0; 0 2])" "true"
test_expr "isdiag([1 2; 3 4])" "false"
test_expr "isdiag([1])" "true"
test_expr "isdiag(5)" "true"

# 6. istriu (5 tests)
echo "Testing istriu..."
test_expr "istriu([1 2 3; 0 4 5; 0 0 6])" "true"
test_expr "istriu([1 2; 3 4])" "false"
test_expr "istriu(eye(3))" "true"
test_expr "istriu([1])" "true"
test_expr "istriu(5)" "true"

# 7. istril (5 tests)
echo "Testing istril..."
test_expr "istril([1 0 0; 2 3 0; 4 5 6])" "true"
test_expr "istril([1 2; 3 4])" "false"
test_expr "istril(eye(3))" "true"
test_expr "istril([1])" "true"
test_expr "istril(5)" "true"

# 8. issorted (5 tests)
echo "Testing issorted..."
test_expr "issorted([1 2 3 4 5])" "true"
test_expr "issorted([5 4 3 2 1])" "false"
test_expr "issorted([1 1 2 2 3])" "true"
test_expr "issorted([1])" "true"
test_expr "issorted(5)" "true"

# 9. ndims (5 tests)
echo "Testing ndims..."
test_expr "ndims([1 2 3])" "2"
test_expr "ndims([1 2; 3 4])" "2"
test_expr "ndims(eye(5))" "2"
test_expr "ndims([1])" "2"
test_expr "ndims(5)" "2"

# 10. triu (5 tests)
echo "Testing triu..."
test_expr "sum(sum(triu([1 2; 3 4])))" "7"
test_expr "triu([1])" "1"
test_expr "sum(sum(triu(eye(3))))" "3"
test_expr "triu([1 2 3 4])" "1  2  3  4"
test_expr "sum(sum(triu([1 2 3; 4 5 6; 7 8 9])))" "26"

# 11. tril (5 tests)
echo "Testing tril..."
test_expr "sum(sum(tril([1 2; 3 4])))" "8"
test_expr "tril([1])" "1"
test_expr "sum(sum(tril(eye(3))))" "3"
test_expr "sum(sum(tril([1 2 3; 4 5 6; 7 8 9])))" "34"
test_expr "sum(tril([1; 2; 3; 4]))" "10"

# 12. horzcat (5 tests)
echo "Testing horzcat..."
test_expr "horzcat([1 2], [3 4])" "1  2  3  4"
test_expr "sum(horzcat([1; 2], [3; 4]))" "10"
test_expr "horzcat([1], [2])" "1  2"
test_expr "horzcat(1, 2)" "1  2"
test_expr "numel(horzcat([1 2], [3 4]))" "4"

# 13. vertcat (5 tests)
echo "Testing vertcat..."
test_expr "sum(vertcat([1 2], [3 4]))" "10"
test_expr "rows(vertcat([1; 2], [3; 4]))" "4"
test_expr "sum(vertcat([1], [2]))" "3"
test_expr "sum(vertcat(1, 2))" "3"
test_expr "rows(vertcat([1 2 3], [4 5 6]))" "2"

# 14. squeeze (5 tests)
echo "Testing squeeze..."
test_expr "sum(squeeze([1 2 3]))" "6"
test_expr "sum(squeeze([1; 2; 3]))" "6"
test_expr "squeeze([1])" "1"
test_expr "squeeze(5)" "5"
test_expr "numel(squeeze([1 2; 3 4]))" "4"

# 15. zscore (5 tests)
echo "Testing zscore..."
test_expr "round(sum(zscore([1 2 3 4 5])))" "0"
test_expr "round(mean(zscore([10 20 30])))" "0"
test_expr "round(std(zscore([1 2 3 4 5])))" "1"
test_expr "numel(zscore([1 2 3]))" "3"
test_expr "zscore([5])" "0"

# 16. sortrows (5 tests)
echo "Testing sortrows..."
test_expr "sum(sum(sortrows([3 1; 1 2; 2 3])))" "12"
test_expr "sortrows([1])" "1"
test_expr "sum(sortrows([5]))" "5"
test_expr "rows(sortrows([3 1; 1 2; 2 3]))" "3"
test_expr "cols(sortrows([3 1; 1 2]))" "2"

# 17. union (5 tests)
echo "Testing union..."
test_expr "union([1 2 3], [2 3 4])" "1  2  3  4"
test_expr "union([1 2], [3 4])" "1  2  3  4"
test_expr "union([5], [5])" "5"
test_expr "union(1, 2)" "1  2"
test_expr "numel(union([1 1], [2 2]))" "2"

# 18. intersect (5 tests)
echo "Testing intersect..."
test_expr "intersect([1 2 3], [2 3 4])" "2  3"
test_expr "numel(intersect([1 2], [3 4]))" "0"
test_expr "intersect([5], [5])" "5"
test_expr "intersect(1, 1)" "1"
test_expr "intersect([1 1 2], [2 2 3])" "2"

# 19. setdiff (5 tests)
echo "Testing setdiff..."
test_expr "setdiff([1 2 3], [2])" "1  3"
test_expr "setdiff([1 2 3], [4 5])" "1  2  3"
test_expr "numel(setdiff([1 2 3], [1 2 3]))" "0"
test_expr "setdiff([5], [6])" "5"
test_expr "setdiff(1, 2)" "1"

# 20. setxor (5 tests)
echo "Testing setxor..."
test_expr "setxor([1 2 3], [2 3 4])" "1  4"
test_expr "setxor([1 2], [3 4])" "1  2  3  4"
test_expr "numel(setxor([1 2 3], [1 2 3]))" "0"
test_expr "setxor([5], [6])" "5  6"
test_expr "numel(setxor(1, 1))" "0"

# 21. maxk (5 tests)
echo "Testing maxk..."
test_expr "maxk([1 5 3 8 2], 3)" "8  5  3"
test_expr "maxk([10 20 30], 2)" "30  20"
test_expr "numel(maxk([1 2 3 4 5], 5))" "5"
test_expr "maxk([1 2 3], 1)" "3"
test_expr "maxk(5, 1)" "5"

# 22. mink (5 tests)
echo "Testing mink..."
test_expr "mink([1 5 3 8 2], 3)" "1  2  3"
test_expr "mink([10 20 30], 2)" "10  20"
test_expr "numel(mink([1 2 3 4 5], 5))" "5"
test_expr "mink([1 2 3], 1)" "1"
test_expr "mink(5, 1)" "5"

# 23. range (5 tests)
echo "Testing range..."
test_expr "range([1 5 3 8 2])" "7"
test_expr "range([10 20 30])" "20"
test_expr "range([5 5 5])" "0"
test_expr "range([1])" "0"
test_expr "range(5)" "0"

# 24. iqr (5 tests)
echo "Testing iqr..."
test_expr "iqr([1 2 3 4 5 6 7 8])" "3.5"
test_expr "iqr([1 2 3 4 5])" "2"
test_expr "iqr([10 20 30 40 50])" "20"
test_expr "iqr([1 1 1 1])" "0"
test_expr "iqr(5)" "0"

# 25. prctile (5 tests)
echo "Testing prctile..."
test_expr "prctile([1 2 3 4 5], 50)" "3"
test_expr "prctile([1 2 3 4 5], 0)" "1"
test_expr "prctile([1 2 3 4 5], 100)" "5"
test_expr "prctile([1 2 3 4 5], 25)" "2"
test_expr "prctile(5, 50)" "5"

# 26. cov (5 tests)
echo "Testing cov..."
test_expr "cov([1 2 3 4 5])" "2.5"
test_expr "cov([10 20 30])" "100"
test_expr "cov([5 5 5])" "0"
test_expr "cov([1])" "0"
test_expr "cov(5)" "0"

# 27. secd (5 tests)
echo "Testing secd..."
test_expr "secd(0)" "1"
test_expr "secd(60)" "2"
test_expr "round(secd(45)*1000)/1000" "1.414"
test_expr "round(secd(30)*1000)/1000" "1.155"
test_expr "secd(180)" "-1"

# 28. cscd (5 tests)
echo "Testing cscd..."
test_expr "cscd(90)" "1"
test_expr "cscd(30)" "2"
test_expr "round(cscd(45)*1000)/1000" "1.414"
test_expr "round(cscd(60)*1000)/1000" "1.155"
test_expr "cscd(270)" "-1"

# 29. cotd (5 tests)
echo "Testing cotd..."
test_expr "cotd(45)" "1"
test_expr "round(cotd(30)*1000)/1000" "1.732"
test_expr "round(cotd(60)*1000)/1000" "0.577"
test_expr "cotd(135)" "-1"
test_expr "cotd(225)" "1"

# 30. asecd (5 tests)
echo "Testing asecd..."
test_expr "asecd(1)" "0"
test_expr "asecd(2)" "60"
test_expr "asecd(-1)" "180"
test_expr "round(asecd(1.5))" "48"
test_expr "round(asecd(-2))" "120"

# 31. acscd (5 tests)
echo "Testing acscd..."
test_expr "acscd(1)" "90"
test_expr "acscd(2)" "30"
test_expr "acscd(-1)" "-90"
test_expr "round(acscd(1.5))" "42"
test_expr "round(acscd(-2))" "-30"

# 32. acotd (5 tests)
echo "Testing acotd..."
test_expr "acotd(0)" "90"
test_expr "acotd(1)" "45"
test_expr "acotd(-1)" "-45"
test_expr "round(acotd(0.5))" "63"
test_expr "round(acotd(2))" "27"

# 33. nthroot (5 tests)
echo "Testing nthroot..."
test_expr "nthroot(8, 3)" "2"
test_expr "nthroot(16, 4)" "2"
test_expr "nthroot(-27, 3)" "-3"
test_expr "nthroot(32, 5)" "2"
test_expr "nthroot(1, 10)" "1"

# 34. realsqrt (5 tests)
echo "Testing realsqrt..."
test_expr "realsqrt(4)" "2"
test_expr "realsqrt(9)" "3"
test_expr "round(realsqrt(2)*1000)/1000" "1.414"
test_expr "realsqrt(0)" "0"
test_expr "realsqrt(100)" "10"

# 35. reallog (5 tests)
echo "Testing reallog..."
test_expr "reallog(1)" "0"
test_expr "reallog(e)" "1"
test_expr "round(reallog(10)*1000)/1000" "2.303"
test_expr "round(reallog(100)*1000)/1000" "4.605"
test_expr "round(abs(reallog(0.5))*1000)/1000" "0.693"

# 36. realpow (5 tests)
echo "Testing realpow..."
test_expr "realpow(2, 3)" "8"
test_expr "realpow(3, 2)" "9"
test_expr "realpow(10, 0)" "1"
test_expr "realpow(-2, 3)" "-8"
test_expr "realpow(4, 0.5)" "2"

# ADDITIONAL TESTS: floor/ceil/trunc bug fix for large numbers
echo "Testing floor/ceil/trunc with large numbers..."
test_expr "floor(333333333333.3333)" "333333333333"
test_expr "ceil(333333333333.3333)" "333333333334"
test_expr "trunc(333333333333.3333)" "333333333333"
test_expr "floor(-333333333333.3333)" "-333333333334"
test_expr "trunc(-333333333333.3333)" "-333333333333"
test_expr "floor(1e15 + 0.5)" "1000000000000000"
test_expr "ceil(1e15 + 0.5)" "1000000000000001"

echo ""
echo "========================================"
echo "Results: $PASS passed, $FAIL failed"
echo "========================================"

exit $FAIL
