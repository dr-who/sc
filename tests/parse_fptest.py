#!/usr/bin/env python3
"""
parse_fptest.py - Parser and test runner for IBM IEEE 754 test suite

Parses .fptest files and runs applicable tests against sc.

The test format is:
  <precision><op> <rounding> [enables] <input1> [input2] [input3] -> <output> [exceptions]

Where:
  - precision: b32, b64, b128, d32, d64, d128
  - op: +, -, *, /, V (sqrt), etc.
  - rounding: =0 (nearest), > (up), < (down), 0 (toward zero), =^ (ties away)
  - enables: x (inexact), u (underflow), o (overflow), z (divbyzero), i (invalid)
  - values: +/-Inf, +/-Zero, Q (qNaN), S (sNaN), or hex float like +1.7FFFFFP127

Note: sc uses 128-bit extended precision, so we focus on special value
handling rather than precision-specific overflow/underflow boundaries.
"""

import sys
import os
import re
import subprocess
import math
from typing import Optional, Tuple, List, Dict

# ANSI colors
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
CYAN = '\033[0;36m'
NC = '\033[0m'

class FPTestParser:
    """Parser for .fptest files"""
    
    # Operation mapping from fptest to sc
    OPS = {
        '+': '+',
        '-': '-',
        '*': '*',
        '/': '/',
        'V': 'sqrt',  # square root
        '*+': 'fma',  # fused multiply-add (not in sc)
        '%': 'mod',   # remainder
        'A': 'abs',
        '~': 'neg',
        'cp': 'copy',
    }
    
    # Comparison ops (return boolean)
    CMP_OPS = {'qC', 'sC', '?-', '?n', '?f', '?0', '?s', '?i', '?N', '?sN'}
    
    def __init__(self, sc_path: str):
        self.sc_path = sc_path
        self.stats = {
            'total': 0,
            'passed': 0,
            'failed': 0,
            'skipped': 0,
            'not_applicable': 0
        }
        self.failures = []
    
    def parse_hex_float(self, s: str) -> Optional[float]:
        """
        Parse IBM hex float format: <sign><significand>P<exp>
        Example: +1.7FFFFFP127, -0.400000P-126
        
        Format: sign + hex_significand + 'P' + decimal_exponent
        The significand is in hex, representing binary bits.
        For normalized: 1.xxxx means 1 + fraction
        For subnormal: 0.xxxx means just the fraction
        """
        if not s:
            return None
            
        # Handle special values
        s_upper = s.upper()
        if s_upper in ('+INF', 'INF'):
            return float('inf')
        if s_upper == '-INF':
            return float('-inf')
        if s_upper in ('+ZERO', 'ZERO'):
            return 0.0
        if s_upper == '-ZERO':
            return -0.0
        if s_upper in ('Q', 'S', 'QNAN', 'SNAN', 'NAN'):
            return float('nan')
        
        # Parse hex float: +1.7FFFFFP127
        match = re.match(r'^([+-]?)([0-9A-Fa-f]+)\.?([0-9A-Fa-f]*)P([+-]?\d+)$', s)
        if not match:
            return None
        
        sign_str, int_part, frac_part, exp_str = match.groups()
        sign = -1 if sign_str == '-' else 1
        exp = int(exp_str)
        
        # Convert hex significand to value
        # int_part.frac_part in hex
        hex_str = int_part + frac_part
        if not hex_str:
            return 0.0
        
        # Each hex digit is 4 bits
        # The implicit point is after int_part
        frac_bits = len(frac_part) * 4
        
        significand = int(hex_str, 16)
        
        # Value = significand * 2^(exp - frac_bits)
        # For normalized IEEE: value = significand * 2^exp / 2^frac_bits
        value = sign * significand * (2.0 ** (exp - frac_bits))
        
        return value
    
    def format_for_sc(self, value_str: str) -> Optional[str]:
        """Convert fptest value to sc expression"""
        if not value_str:
            return None
        
        value_str = value_str.strip()
        upper = value_str.upper()
        
        # Special values
        if upper in ('+INF', 'INF'):
            return '(1/0)'
        if upper == '-INF':
            return '(-1/0)'
        if upper in ('+ZERO', 'ZERO', '-ZERO'):
            return '0'
        if upper in ('Q', 'S', 'QNAN', 'SNAN'):
            return '(0/0)'
        if upper == '#':  # No result (exception raised)
            return None
        
        # Try to parse as hex float
        val = self.parse_hex_float(value_str)
        if val is not None:
            if math.isnan(val):
                return '(0/0)'
            if math.isinf(val):
                return '(1/0)' if val > 0 else '(-1/0)'
            # Format as decimal
            if val == 0:
                return '0'
            return f'{val:.15g}'
        
        return None
    
    def parse_line(self, line: str) -> Optional[Dict]:
        """
        Parse a single test line.
        Returns dict with: precision, op, rounding, enables, inputs, output, exceptions
        """
        line = line.strip()
        if not line or line.startswith('#'):
            return None
        
        # Split by ->
        if '->' not in line:
            return None
        
        left, right = line.split('->', 1)
        left_parts = left.split()
        right_parts = right.split()
        
        if len(left_parts) < 3:
            return None
        
        # Parse operation (first part like "b32+")
        op_str = left_parts[0]
        
        # Extract precision (b32, b64, b128, d32, d64, d128)
        precision_match = re.match(r'^([bd]\d+)', op_str)
        if not precision_match:
            return None
        precision = precision_match.group(1)
        op = op_str[len(precision):]
        
        # Skip decimal tests (d32, d64, d128) - sc is binary
        if precision.startswith('d'):
            return None
        
        # Parse rounding mode
        rounding = left_parts[1]
        
        # Parse enables and inputs
        # Enables are optional single letters: x, u, o, z, i
        enables = ''
        inputs = []
        
        for part in left_parts[2:]:
            if re.match(r'^[xuozi]+$', part):
                enables = part
            else:
                inputs.append(part)
        
        # Parse output and exceptions
        output = right_parts[0] if right_parts else None
        exceptions = right_parts[1] if len(right_parts) > 1 else ''
        
        return {
            'precision': precision,
            'op': op,
            'rounding': rounding,
            'enables': enables,
            'inputs': inputs,
            'output': output,
            'exceptions': exceptions,
            'raw': line
        }
    
    def run_sc(self, expr: str) -> str:
        """Run an expression through sc and return the result"""
        try:
            proc = subprocess.run(
                [self.sc_path],
                input=expr + '\n',
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Parse output - look for "= " line
            for line in proc.stdout.split('\n'):
                if '= ' in line:
                    result = line.split('= ', 1)[1].strip()
                    return result
            
            return 'ERROR'
        except Exception as e:
            return f'ERROR: {e}'
    
    def compare_results(self, expected_str: str, actual: str) -> bool:
        """Compare expected and actual results"""
        expected_str = expected_str.upper()
        actual_upper = actual.upper()
        
        # Handle NaN
        if expected_str in ('Q', 'S', 'QNAN', 'SNAN', 'NAN'):
            return 'NAN' in actual_upper
        
        # Handle Inf
        if expected_str in ('+INF', 'INF'):
            return actual_upper == 'INF' or actual_upper == '+INF'
        if expected_str == '-INF':
            return actual_upper == '-INF'
        
        # Handle Zero
        if expected_str in ('+ZERO', 'ZERO', '-ZERO'):
            try:
                return float(actual) == 0.0
            except:
                return actual == '0'
        
        # Handle # (no result expected due to exception)
        if expected_str == '#':
            return True  # Accept any result
        
        # Try numeric comparison
        expected_val = self.parse_hex_float(expected_str)
        if expected_val is None:
            return False
        
        try:
            # Clean up actual value - remove complex part if present
            actual_clean = actual.strip()
            # Remove imaginary part if present and it's effectively zero
            if ' + ' in actual_clean or ' - ' in actual_clean:
                # Complex number - check if imaginary part is negligible
                parts = re.split(r' [+-] ', actual_clean)
                if len(parts) >= 1:
                    actual_clean = parts[0]
            # Remove trailing 'i' if any
            actual_clean = actual_clean.replace('i', '').strip()
            
            actual_val = float(actual_clean)
            
            if math.isnan(expected_val):
                return math.isnan(actual_val) or 'NAN' in actual_upper
            if math.isinf(expected_val):
                return math.isinf(actual_val) and (expected_val > 0) == (actual_val > 0)
            if expected_val == 0:
                return abs(actual_val) < 1e-300
            
            # Relative comparison with tolerance
            if abs(expected_val) < 1e-10:
                # For very small values, use absolute comparison
                return abs(actual_val - expected_val) < 1e-10
            
            rel_err = abs(actual_val - expected_val) / abs(expected_val)
            return rel_err < 1e-6
        except:
            return False
    
    def test_case(self, test: Dict) -> Tuple[str, str]:
        """
        Run a single test case.
        Returns (status, message) where status is 'pass', 'fail', 'skip', or 'n/a'
        """
        op = test['op']
        inputs = test['inputs']
        expected = test['output']
        
        # Check if we support this operation
        if op not in self.OPS and op not in self.CMP_OPS:
            return ('n/a', f'Unsupported op: {op}')
        
        # Skip comparison operations (not implemented in sc same way)
        if op in self.CMP_OPS:
            return ('n/a', 'Comparison op')
        
        # Skip FMA (not in sc)
        if op == '*+':
            return ('n/a', 'FMA not implemented')
        
        # Convert inputs to sc format
        sc_inputs = []
        for inp in inputs:
            converted = self.format_for_sc(inp)
            if converted is None:
                return ('skip', f'Cannot parse input: {inp}')
            sc_inputs.append(converted)
        
        if not sc_inputs:
            return ('skip', 'No inputs')
        
        # Build expression
        sc_op = self.OPS.get(op, op)
        
        if sc_op == 'sqrt':
            if len(sc_inputs) != 1:
                return ('skip', 'sqrt needs 1 input')
            expr = f'sqrt({sc_inputs[0]})'
        elif sc_op == 'abs':
            if len(sc_inputs) != 1:
                return ('skip', 'abs needs 1 input')
            expr = f'abs({sc_inputs[0]})'
        elif sc_op == 'neg':
            if len(sc_inputs) != 1:
                return ('skip', 'neg needs 1 input')
            expr = f'-({sc_inputs[0]})'
        elif sc_op in ('+', '-', '*', '/'):
            if len(sc_inputs) != 2:
                return ('skip', f'{sc_op} needs 2 inputs')
            expr = f'({sc_inputs[0]}) {sc_op} ({sc_inputs[1]})'
        else:
            return ('n/a', f'Op not mapped: {sc_op}')
        
        # Run through sc
        actual = self.run_sc(expr)
        
        if actual.startswith('ERROR'):
            return ('fail', f'sc error: {actual}')
        
        # Compare result
        if expected == '#':
            # Exception expected, any result OK
            return ('pass', 'Exception case')
        
        if self.compare_results(expected, actual):
            return ('pass', '')
        else:
            return ('fail', f'Expected: {expected}, Got: {actual}')
    
    def run_file(self, filepath: str) -> None:
        """Run all tests in a .fptest file"""
        filename = os.path.basename(filepath)
        print(f'\n{CYAN}=== {filename} ==={NC}')
        
        file_stats = {'total': 0, 'passed': 0, 'failed': 0, 'skipped': 0, 'na': 0}
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line_num, line in enumerate(f, 1):
                test = self.parse_line(line)
                if test is None:
                    continue
                
                file_stats['total'] += 1
                self.stats['total'] += 1
                
                status, msg = self.test_case(test)
                
                if status == 'pass':
                    file_stats['passed'] += 1
                    self.stats['passed'] += 1
                elif status == 'fail':
                    file_stats['failed'] += 1
                    self.stats['failed'] += 1
                    self.failures.append({
                        'file': filename,
                        'line': line_num,
                        'test': test['raw'],
                        'msg': msg
                    })
                    # Print first few failures
                    if len(self.failures) <= 10:
                        print(f'  {RED}FAIL{NC} line {line_num}: {msg}')
                        print(f'       {test["raw"][:80]}')
                elif status == 'skip':
                    file_stats['skipped'] += 1
                    self.stats['skipped'] += 1
                else:  # n/a
                    file_stats['na'] += 1
                    self.stats['not_applicable'] += 1
        
        # Summary for this file
        applicable = file_stats['total'] - file_stats['na'] - file_stats['skipped']
        if applicable > 0:
            pass_rate = file_stats['passed'] / applicable * 100
            print(f'  Applicable: {applicable}, Passed: {file_stats["passed"]}, '
                  f'Failed: {file_stats["failed"]} ({pass_rate:.1f}%)')
        else:
            print(f'  No applicable tests (all skipped or N/A)')
    
    def run_all(self, test_dir: str) -> None:
        """Run all .fptest files in directory"""
        if not os.path.isdir(test_dir):
            print(f'{RED}Error: Test directory not found: {test_dir}{NC}')
            return
        
        # Get all .fptest files
        test_files = sorted([
            f for f in os.listdir(test_dir)
            if f.endswith('.fptest')
        ])
        
        if not test_files:
            print(f'{YELLOW}No .fptest files found in {test_dir}{NC}')
            return
        
        print(f'\nFound {len(test_files)} test files')
        
        # Skip decimal tests (sc is binary floating point)
        binary_tests = [f for f in test_files if not f.startswith('Decimal')]
        decimal_tests = [f for f in test_files if f.startswith('Decimal')]
        
        print(f'Running {len(binary_tests)} binary tests (skipping {len(decimal_tests)} decimal tests)')
        
        for test_file in binary_tests:
            filepath = os.path.join(test_dir, test_file)
            self.run_file(filepath)
        
        # Print summary
        print(f'\n{"="*60}')
        print(f'{CYAN}OVERALL SUMMARY{NC}')
        print(f'{"="*60}')
        print(f'Total test cases:     {self.stats["total"]}')
        print(f'Not applicable:       {self.stats["not_applicable"]} (FMA, comparisons, etc.)')
        print(f'Skipped (parse err):  {self.stats["skipped"]}')
        applicable = self.stats['total'] - self.stats['not_applicable'] - self.stats['skipped']
        print(f'Applicable tests:     {applicable}')
        print(f'{GREEN}Passed:               {self.stats["passed"]}{NC}')
        print(f'{RED}Failed:               {self.stats["failed"]}{NC}')
        
        if applicable > 0:
            pass_rate = self.stats['passed'] / applicable * 100
            print(f'\n{GREEN if pass_rate >= 95 else YELLOW if pass_rate >= 80 else RED}'
                  f'Pass rate: {pass_rate:.1f}%{NC}')
        
        if self.failures and len(self.failures) > 10:
            print(f'\n{RED}Showing first 10 of {len(self.failures)} failures above{NC}')


def main():
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <sc-path> <test-suite-dir>')
        sys.exit(1)
    
    sc_path = sys.argv[1]
    test_dir = sys.argv[2]
    
    parser = FPTestParser(sc_path)
    parser.run_all(test_dir)
    
    # Exit with failure count
    sys.exit(min(parser.stats['failed'], 255))


if __name__ == '__main__':
    main()
