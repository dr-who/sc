# sc - Scientific Calculator

A high-precision scientific calculator with complex number support, matrix operations, and IEEE 754-2008 compliance.

Join us as we run 64/128/256/512/1024 and more floating point math on DOS, VIC-20s, MacOS, Linux, iPhones.

## Features

- **Compile time precision**: e.g. in 128-bit more it's ~34 significant decimal digits
- **Complex Numbers**: Full arithmetic with `i` as imaginary unit
- **Matrices**: Create, manipulate, and solve linear systems
- **Scientific Functions**: Trig, hyperbolic, exponential, logarithmic
- **Multiple Formats**: Decimal, hex (0x), binary (0b), octal (0o), Roman numerals
- **IEEE 754-2008**: Proper Inf, NaN, signed zero handling
- **Interactive REPL**: With command history and line editing
- **RPN Mode**: HP calculator-style reverse Polish notation
- **TVM**: Time Value of Money calculations
- **Statistical Functions**: Mean, stddev, regression

## Quick Start

```bash
# Build
make

# Run interactively
./bin/sc

# Command line
./bin/sc "sqrt(2)"
./bin/sc "exp(i*pi)+1"    # Euler's identity = 0

# Pipe mode
echo "2^100" | ./bin/sc
```

## Usage

```
>>> 2 + 3 * 4
= 14

>>> sqrt(-1)
= i

>>> (1+2i) * (3-4i)
= 11+2i

>>> exp(i*pi) + 1
= 0

>>> 1/3; ans*3
= 1

>>> A = [1,2;3,4]
>>> det(A)
= -2

>>> 0xFF + 0b1010
= 265
```

## Operators

| Operator | Description |
|----------|-------------|
| `+` `-` `*` `/` | Basic arithmetic |
| `^` | Exponentiation (right associative) |
| `%` | Modulo |
| `!` | Factorial (postfix) |
| `( )` | Grouping |

## Functions

### Trigonometric (radians)
`sin` `cos` `tan` `asin` `acos` `atan` `atan2`

### Hyperbolic
`sinh` `cosh` `tanh` `asinh` `acosh` `atanh`

### Exponential & Logarithmic
- `exp(x)` - e^x
- `log(x)` `ln(x)` - Natural logarithm
- `log10(x)` - Base-10 logarithm
- `log2(x)` - Base-2 logarithm

### Power & Root
`sqrt` `cbrt` `pow`

### Complex
`abs` `arg` `conj` `real` `imag`

### Other
`floor` `ceil` `round` `trunc` `frac` `min` `max` `gcd` `lcm`

## Constants

| Constant | Value |
|----------|-------|
| `pi` | 3.14159265358979323846... |
| `e` | 2.71828182845904523536... |
| `i` | √(-1) |
| `Inf` | Positive infinity |
| `NaN` | Not a number |
| `ans` | Previous result |

## Commands

| Command | Description |
|---------|-------------|
| `help` | Show help |
| `digits N` | Set display precision (default 16, max 37) |
| `round MODE` | Set rounding: nearest/ceil/floor/trunc/away |
| `vars` | Show variables |
| `clear` | Clear variables |
| `test` | Run test suite |
| `rpn` | Enter RPN mode |
| `tvm` | Time Value of Money mode |
| `stat` | Statistics mode |
| `quit` | Exit |

## Installation

### From Source
```bash
make
sudo make install
```

### RPM Package (Fedora/RHEL/CentOS)
```bash
# Build RPM
make rpm

# Install
sudo rpm -i ~/rpmbuild/RPMS/x86_64/sc-1.0.0-1.*.rpm

# Or use dnf/yum
sudo dnf install ~/rpmbuild/RPMS/x86_64/sc-1.0.0-1.*.rpm
```

## Testing

```bash
# Run all tests
make test

# Or run individual test suites:
./tests/test_bedmas.sh ./bin/sc      # 127 order-of-operations tests
./tests/test_main.sh ./bin/sc        # Comprehensive test runner

# IEEE 754 test suite (requires download)
cd tests
./download_ieee754_suite.sh          # Download IBM FPgen test suite
python3 parse_fptest.py ../bin/sc ieee754-test-suite
```

### Test Suites

| Suite | Tests | Description |
|-------|-------|-------------|
| BEDMAS | 127 | Order of operations, unary minus, brackets |
| Main | 117 | Core functionality |
| IEEE 754 | 201 | IEEE 754-2008 compliance |
| BC Compare | 97 | Comparison with GNU bc |
| IBM FPgen | 33 files | Full IEEE 754R test vectors |

## Directory Structure

```
sc/
├── src/           # Source code (C89)
├── bin/           # Built binary
├── obj/           # Object files
├── tests/         # Test suites
│   ├── test_bedmas.sh
│   ├── test_main.sh
│   ├── parse_fptest.py
│   └── ieee754-test-suite/
├── packaging/
│   ├── sc.spec    # RPM spec file
│   └── sc.1       # Man page
├── Makefile
├── README.md
└── LICENSE
```

## Notes

- `log()` and `ln()` both compute natural logarithm (use `log10()` for base-10)
- Exponentiation is right-associative: `2^3^2` = `2^(3^2)` = 512
- `ans` updates after each semicolon-separated expression

## Man Page

After installation, view the man page:
```bash
man sc
```

## License

MIT License

## Author

Stu
