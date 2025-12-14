#!/bin/bash
# download_ieee754_suite.sh - Download IBM FPgen IEEE 754 test suite
# From: https://github.com/sergev/ieee754-test-suite

set -e

SUITE_DIR="$(dirname "$0")/ieee754-test-suite"
BASE_URL="https://raw.githubusercontent.com/sergev/ieee754-test-suite/master"

# List of test files (33 fptest files)
FILES=(
    "Add-Cancellation-And-Subnorm-Result.fptest"
    "Add-Cancellation.fptest"
    "Add-Shift-And-Special-Significands.fptest"
    "Add-Shift.fptest"
    "Basic-Types-Inputs.fptest"
    "Basic-Types-Intermediate.fptest"
    "Compare-Different-Input-Field-Relations.fptest"
    "Corner-Rounding.fptest"
    "Decimal-Basic-Types-Inputs.fptest"
    "Decimal-Basic-Types-Intermediate.fptest"
    "Decimal-Clamping.fptest"
    "Decimal-Mul-Trailing-Zeros.fptest"
    "Decimal-Overflow.fptest"
    "Decimal-Rounding.fptest"
    "Decimal-Trailing-And-Leading-Zeros-Input.fptest"
    "Decimal-Trailing-And-Leading-Zeros-Result.fptest"
    "Decimal-Underflow.fptest"
    "Divide-Divide-By-Zero-Exception.fptest"
    "Divide-Trailing-Zeros.fptest"
    "Hamming-Distance.fptest"
    "Input-Special-Significand.fptest"
    "MultiplyAdd-Cancellation-And-Subnorm-Result.fptest"
    "MultiplyAdd-Cancellation.fptest"
    "MultiplyAdd-Shift-And-Special-Significands.fptest"
    "MultiplyAdd-Shift.fptest"
    "MultiplyAdd-Special-Events-Inexact.fptest"
    "MultiplyAdd-Special-Events-Overflow.fptest"
    "MultiplyAdd-Special-Events-Underflow.fptest"
    "Overflow.fptest"
    "Rounding.fptest"
    "Sticky-Bit-Calculation.fptest"
    "Underflow.fptest"
    "Vicinity-Of-Rounding-Boundaries.fptest"
    "syntax.txt"
    "README.md"
)

echo "Downloading IEEE 754 test suite from GitHub..."
echo "Source: https://github.com/sergev/ieee754-test-suite"
echo ""

mkdir -p "$SUITE_DIR"

# Try git clone first (faster if available)
if command -v git &> /dev/null; then
    echo "Attempting git clone..."
    if git clone --depth 1 https://github.com/sergev/ieee754-test-suite.git "$SUITE_DIR.tmp" 2>/dev/null; then
        rm -rf "$SUITE_DIR"
        mv "$SUITE_DIR.tmp" "$SUITE_DIR"
        echo "Successfully cloned repository"
        echo "Downloaded $(ls "$SUITE_DIR"/*.fptest 2>/dev/null | wc -l) test files"
        exit 0
    fi
    rm -rf "$SUITE_DIR.tmp"
    echo "Git clone failed, falling back to wget/curl..."
fi

# Fall back to downloading individual files
downloaded=0
failed=0

for file in "${FILES[@]}"; do
    echo -n "  Downloading $file... "
    if command -v curl &> /dev/null; then
        if curl -sLf "$BASE_URL/$file" -o "$SUITE_DIR/$file" 2>/dev/null; then
            echo "OK"
            ((downloaded++))
        else
            echo "FAILED"
            ((failed++))
        fi
    elif command -v wget &> /dev/null; then
        if wget -q "$BASE_URL/$file" -O "$SUITE_DIR/$file" 2>/dev/null; then
            echo "OK"
            ((downloaded++))
        else
            echo "FAILED"
            ((failed++))
        fi
    else
        echo "ERROR: Neither curl nor wget available"
        exit 1
    fi
done

echo ""
echo "Download complete: $downloaded files downloaded, $failed failed"

if [ $downloaded -gt 0 ]; then
    echo "Test files are in: $SUITE_DIR/"
    echo ""
    echo "To run tests:"
    echo "  python3 parse_fptest.py ../bin/sc ieee754-test-suite"
fi
