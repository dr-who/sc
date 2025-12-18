#!/bin/bash
# run_all_tests.sh - Run all test suites

set -e  # Exit on first failure

SC="${SC:-./bin/sc}"
export SC

echo "============================================="
echo "SCALC COMPREHENSIVE TEST SUITE"
echo "============================================="
echo ""

echo ">>> Core Tests (make test)"
make test
echo ""

echo ">>> Fisher Iris Tests"
./tests/test_iris.sh
echo ""

echo ">>> Machine Learning Tests"
./tests/test_ml.sh
echo ""

echo ">>> Linear Algebra Tests"
./tests/test_linalg.sh
echo ""

echo ">>> Command-Line Interface Tests"
./tests/test_cli.sh
echo ""

echo ">>> Parser Edge Case Tests"
./tests/test_parser.sh
echo ""

echo "============================================="
echo "ALL TEST SUITES PASSED"
echo "============================================="
