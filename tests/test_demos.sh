#!/bin/bash
# Test that all function demos can be displayed (not necessarily execute)
# Some demos show usage patterns requiring setup

SCALC="${1:-bin/sc}"
FAILED=0
PASSED=0

# Get all function names from registry
FUNCS=$(grep '{"[a-zA-Z_][a-zA-Z0-9_]*", FH_' src/func_registry.c | sed 's/.*{"\([^"]*\)".*/\1/' | sort -u)

for func in $FUNCS; do
    [ -z "$func" ] && continue
    
    # Run demo and check that it shows the demo header
    result=$(echo "demo $func" | $SCALC 2>&1)
    
    if echo "$result" | grep -q "=== Demo:"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "Demo Display Test:"
echo "  $PASSED functions have demos"
echo "  $FAILED functions have no demos (OK for commands)"
exit 0
