#!/bin/bash
# Test script for Fisher Iris functionality in scalc
# This checks specific values to ensure MATLAB compatibility

cd "$(dirname "$0")/.."
SC="./bin/sc"

pass=0
fail=0

test_contains() {
    local cmd="$1"
    local expected="$2"
    local result
    result=$(printf '%s\nquit\n' "$cmd" | $SC 2>&1)
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ $cmd"
        ((pass++))
    else
        echo "  ✗ $cmd (expected to contain: $expected)"
        echo "    got: $result"
        ((fail++))
    fi
}

test_exact() {
    local cmd="$1"
    local expected="$2"
    local result
    result=$(printf '%s\nquit\n' "$cmd" | $SC 2>&1 | tail -1)
    if [ "$result" = "$expected" ]; then
        echo "  ✓ $cmd"
        ((pass++))
    else
        echo "  ✗ $cmd"
        echo "    expected: $expected"
        echo "    got: $result"
        ((fail++))
    fi
}

# Helper for multi-line commands
test_multi() {
    local expected="$1"
    shift
    local result
    result=$(printf '%s\n' "$@" "quit" | $SC 2>&1 | tail -1)
    if [ "$result" = "$expected" ]; then
        echo "  ✓ $*"
        ((pass++))
    else
        echo "  ✗ $*"
        echo "    expected: $expected"
        echo "    got: $result"
        ((fail++))
    fi
}

test_multi_contains() {
    local expected="$1"
    shift
    local result
    result=$(printf '%s\n' "$@" "quit" | $SC 2>&1)
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ $*"
        ((pass++))
    else
        echo "  ✗ $* (expected to contain: $expected)"
        echo "    got: $result"
        ((fail++))
    fi
}

echo "Testing Fisher Iris Dataset Features"
echo "====================================="
echo ""

echo "Dataset Loading:"
test_contains "load fisheriris" "meas:    150x4"
test_contains "load fisheriris" "species: 150x1"

echo ""
echo "Matrix Size Functions:"
test_multi "150" "load fisheriris" "rows(meas)"
test_multi "4" "load fisheriris" "cols(meas)"
test_multi "600" "load fisheriris" "size(meas)"
test_multi "150" "load fisheriris" "size(meas, 1)"
test_multi "4" "load fisheriris" "size(meas, 2)"

echo ""
echo "Matrix Indexing:"
test_multi "5.1" "load fisheriris" "meas(1, 1)"
test_multi "1.5" "load fisheriris" "meas(10, 3)"
test_multi_contains "5.1  3.5" "load fisheriris" "meas(1:3, 1:2)"
test_multi_contains "4.9  3" "load fisheriris" "meas(1:3, 1:2)"

echo ""
echo "Named Variable Assignment:"
test_multi "50" "load fisheriris" "setosa = meas(1:50, :)" "rows(setosa)"
test_multi "50" "load fisheriris" "versicolor = meas(51:100, :)" "rows(versicolor)"
test_multi "50" "load fisheriris" "virginica = meas(101:150, :)" "rows(virginica)"

echo ""
echo "Statistical Functions:"
test_multi_contains "2.533" "load fisheriris" "setosa = meas(1:50, :)" "mean(setosa)"
test_multi_contains "3.573" "load fisheriris" "versicolor = meas(51:100, :)" "mean(versicolor)"
test_multi_contains "4.285" "load fisheriris" "virginica = meas(101:150, :)" "mean(virginica)"

echo ""
echo "Correlation Matrix:"
test_multi_contains "0.87175" "load fisheriris" "corr(meas)"
test_multi_contains "0.96275" "load fisheriris" "corr(meas)"

echo ""
echo "Covariance Matrix:"
test_multi_contains "0.6856" "load fisheriris" "cov(meas)"
test_multi_contains "3.113" "load fisheriris" "cov(meas)"

echo ""
echo "K-Means Clustering:"
test_multi "150" "load fisheriris" "[idx, C] = kmeans(meas, 3)" "rows(idx)"
test_multi_contains "50  0  0" "load fisheriris" "[idx, C] = kmeans(meas, 3)" "crosstab(idx, species)"

echo ""
echo "Linear Algebra:"
test_exact "det([1,2;3,4])" "-2"
test_multi_contains "-2  1" "inv([1,2;3,4])"
test_multi_contains "1.5  -0.5" "inv([1,2;3,4])"
test_contains "eig([1,2;3,4])" "5.3722"
test_exact "det(magic(3))" "-360"
test_contains "linsolve(magic(3), [1;2;3])" "0.05"

echo ""
echo "Machine Learning:"
test_multi_contains "0.96" "load fisheriris" "mdl = fitcknn(meas, species, 5)" "Ypred = predict(mdl, meas)" "accuracy(species, Ypred)"
test_multi_contains "50  0  0" "load fisheriris" "mdl = fitcknn(meas, species, 5)" "Ypred = predict(mdl, meas)" "confusionmat(species, Ypred)"
test_multi_contains "150" "load fisheriris" "mdl = fitcsvm(meas, species)" "Ypred = predict(mdl, meas)" "rows(Ypred)"
test_multi_contains "150" "load fisheriris" "mdl = fitctree(meas, species)" "Ypred = predict(mdl, meas)" "rows(Ypred)"
test_multi_contains "150" "load fisheriris" "mdl = fitcnb(meas, species)" "Ypred = predict(mdl, meas)" "rows(Ypred)"
test_multi_contains "0.96" "load fisheriris" "mdl = fitcnb(meas, species)" "Ypred = predict(mdl, meas)" "accuracy(species, Ypred)"
test_multi_contains "0" "load fisheriris" "Z = zscore(meas)" "mean(Z)"
test_multi_contains "105" "load fisheriris" "train_idx = cvpartition(species, 30)" "rows(train_idx)"

echo ""
echo "====================================="
echo "Results: $pass passed, $fail failed"

if [ $fail -gt 0 ]; then
    exit 1
fi
exit 0
