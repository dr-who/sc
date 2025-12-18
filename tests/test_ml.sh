#!/bin/bash
# test_ml.sh - Comprehensive Machine Learning tests
# Tests all ML functions: fitcknn, fitcsvm, fitctree, fitcnb, predict, accuracy, etc.

SC="${SC:-./bin/sc}"
pass=0
fail=0

# Test function: run expression and check if output contains expected string
test_expr() {
    local expected="$1"
    shift
    local expr="$*"
    
    # Replace newlines with actual newlines in expression
    local result=$(echo -e "$expr\nquit" | $SC 2>&1)
    
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ $expr"
        ((pass++))
    else
        echo "  ✗ $expr"
        echo "    Expected: $expected"
        echo "    Got: $(echo "$result" | tail -5)"
        ((fail++))
    fi
}

# Test function: run multi-line commands
test_multi() {
    local expected="$1"
    shift
    local commands=""
    for cmd in "$@"; do
        commands+="$cmd"$'\n'
    done
    commands+="quit"
    
    local result=$(echo "$commands" | $SC 2>&1)
    local desc=$(echo "$@" | tr '\n' ' ')
    
    if echo "$result" | grep -qF -- "$expected"; then
        echo "  ✓ $desc"
        ((pass++))
    else
        echo "  ✗ $desc"
        echo "    Expected: $expected"
        echo "    Got: $(echo "$result" | tail -5)"
        ((fail++))
    fi
}

echo "Comprehensive Machine Learning Tests"
echo "====================================="
echo ""

echo "KNN Classifier (fitcknn):"
test_multi "k-NN" "load fisheriris" "mdl = fitcknn(meas, species, 5)"
test_multi "0.96" "load fisheriris" "mdl = fitcknn(meas, species, 5)" "pred = predict(mdl, meas)" "accuracy(species, pred)"
test_multi "k=3" "load fisheriris" "mdl = fitcknn(meas, species, 3)"
test_multi "50  0  0" "load fisheriris" "mdl = fitcknn(meas, species, 5)" "pred = predict(mdl, meas)" "confusionmat(species, pred)"
echo ""

echo "Naive Bayes Classifier (fitcnb):"
test_multi "Naive Bayes" "load fisheriris" "mdl = fitcnb(meas, species)"
test_multi "0.96" "load fisheriris" "mdl = fitcnb(meas, species)" "pred = predict(mdl, meas)" "accuracy(species, pred)"
test_multi "150" "load fisheriris" "mdl = fitcnb(meas, species)" "pred = predict(mdl, meas)" "rows(pred)"
echo ""

echo "Decision Tree Classifier (fitctree):"
test_multi "Decision Tree" "load fisheriris" "mdl = fitctree(meas, species)"
test_multi "feature" "load fisheriris" "mdl = fitctree(meas, species)"
test_multi "150" "load fisheriris" "mdl = fitctree(meas, species)" "pred = predict(mdl, meas)" "rows(pred)"
echo ""

echo "SVM Classifier (fitcsvm):"
test_multi "linear SVM" "load fisheriris" "mdl = fitcsvm(meas, species)"
test_multi "150" "load fisheriris" "mdl = fitcsvm(meas, species)" "pred = predict(mdl, meas)" "rows(pred)"
echo ""

echo "Multiple Models:"
test_multi "0.96" "load fisheriris" "m1 = fitcknn(meas, species, 5)" "m2 = fitcnb(meas, species)" "p1 = predict(m1, meas)" "p2 = predict(m2, meas)" "accuracy(species, p2)"
echo ""

echo "Data Preprocessing (zscore):"
test_multi "0" "load fisheriris" "Z = zscore(meas)" "mean(Z)"
test_multi "1" "load fisheriris" "Z = zscore(meas)" "std(Z)"
echo ""

echo "Cross-Validation (cvpartition):"
test_multi "105" "load fisheriris" "train = cvpartition(species, 30)" "rows(train)"
test_multi "70" "load fisheriris" "train = cvpartition(species, 50)" "rows(train)"
echo ""

echo "Model Evaluation:"
test_multi "3" "load fisheriris" "mdl = fitcknn(meas, species, 5)" "pred = predict(mdl, meas)" "C = confusionmat(species, pred)" "rows(C)"
echo ""

echo "Edge Cases:"
# Test with simple 2-class data
test_multi "2 classes" "A = [1;2;3;4;5;6;7;8;9;10]" "B = [1;1;1;1;1;2;2;2;2;2]" "m = fitcknn(A, B, 3)"
echo ""

echo "===================================="
echo "Results: $pass passed, $fail failed"

if [ $fail -gt 0 ]; then
    exit 1
fi
exit 0
