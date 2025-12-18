#!/bin/bash
# Test that all functions have help entries

SCALC="${1:-bin/sc}"
FAILED=0
PASSED=0
FAILED_LIST=""

# Get all function names from registry (only valid identifiers followed by handler)
FUNCS=$(grep '{"[a-zA-Z_][a-zA-Z0-9_]*", FH_' src/func_registry.c | sed 's/.*{"\([^"]*\)".*/\1/' | sort -u)

for func in $FUNCS; do
    # Skip NULL and empty
    [ -z "$func" ] && continue
    
    # Test help
    result=$(echo "help $func" | $SCALC 2>&1 | grep -c "Unknown function")
    if [ "$result" -gt 0 ]; then
        FAILED_LIST="$FAILED_LIST $func"
        FAILED=$((FAILED + 1))
    else
        PASSED=$((PASSED + 1))
    fi
done

echo ""
echo "Help Coverage Test:"
echo "  Passed: $PASSED"
echo "  Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo "  Missing help:"
    for f in $FAILED_LIST; do
        echo "    $f"
    done
    exit 1
fi
exit 0
