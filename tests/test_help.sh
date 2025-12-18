#!/bin/bash
# Test that help works for all registered functions

SCALC=${1:-./bin/sc}
FAILED=0
PASSED=0

echo "Testing help system..."

# Extract function names from registry - match only the function name field
FUNCS=$(grep -E '^ *\{"[a-zA-Z_][a-zA-Z0-9_]*", *FH_' src/func_registry.c | sed 's/.*{"\([^"]*\)".*/\1/' | sort -u)

for func in $FUNCS; do
    result=$(echo "help $func
quit" | $SCALC 2>&1 | head -3)
    if echo "$result" | grep -q "Unknown function"; then
        echo "FAIL: help $func"
        FAILED=$((FAILED + 1))
    else
        PASSED=$((PASSED + 1))
    fi
done

echo ""
echo "Help test results: $PASSED passed, $FAILED failed"

if [ $FAILED -gt 0 ]; then
    exit 1
fi
exit 0
