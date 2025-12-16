# sc - Stu's Scientific Calculator

A high-precision scientific calculator designed to be the ultimate portable calculation tool. Runs on everything from 8-bit VIC-20 to modern Linux, with pure C89 code and no external dependencies.

## Highlights

- **128-bit Soft-Float Precision**: ~34 significant decimal digits, IEEE 754-2008 compliant
- **Pure C89**: No float/double in core library - runs on 8-bit through 64-bit platforms
- **Complex Numbers**: Full arithmetic with `i` as imaginary unit
- **Matrices**: Up to 8×8 with determinant, inverse, eigenvalues, QR decomposition
- **MATLAB Compatible**: 365+ functions - colon notation, indexing, set operations, statistics, predicates, signal processing, coordinate conversions
- **RPN Mode**: HP-style Reverse Polish Notation with 4/16 level stack and LASTX
- **Special Functions**: Gamma, Bessel, elliptic integrals, error function, Lambert W
- **Statistical Distributions**: Normal, t, chi-squared, F, binomial, Poisson
- **Multiple Platforms**: Linux, DOS (Watcom C), VIC-20 (cc65), Windows

## Quick Start

```bash
# Build for Linux (full features)
make

# Build for DOS 16-bit (medium features)
make CFLAGS="-DSCALC_MEDIUM"

# Build minimal for embedded/8-bit
make CFLAGS="-DSCALC_TINY"

# Run interactively
./bin/sc

# Command line evaluation
./bin/sc "sqrt(2)"
./bin/sc "exp(i*pi)+1"    # Euler's identity = 0
./bin/sc "gamma(5.5)"     # = 52.342777...

# Pipe mode
echo "2^100" | ./bin/sc
```

## Platform Builds

| Platform | Preset | Features | Precision |
|----------|--------|----------|-----------|
| Linux/macOS | (default) | All Phase 1-4 | 128-bit (8 limbs) |
| DOS 16-bit | `SCALC_MEDIUM` | Phase 1-2 | 128-bit (8 limbs) |
| DOS minimal | `SCALC_TINY` | Phase 1 core | 128-bit (8 limbs) |
| VIC-20 | `SCALC_VIC20` | Core + CORDIC trig | 64-bit (4 limbs) |
| Embedded | `SCALC_MINIMAL` | Basic + sqrt/exp | 64-bit (4 limbs) |

## Features by Phase

### Phase 1: Core Foundation
- Basic arithmetic: `+` `-` `*` `/` `^` `%`
- Square root, nth root: `sqrt(x)`, `cbrt(x)`, `x^(1/n)`
- Trigonometry: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`
- Hyperbolic: `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`
- Exponential/logarithm: `exp`, `log`, `ln`, `log10`, `log2`, `log(x,base)`
- Power: `pow(x,y)`, `x^y`
- RPN mode with LASTX register
- Display modes: FIX, SCI, ENG, ALL

### Phase 2: Scientific Essentials
- Gamma function: `gamma(x)`, `lgamma(x)` (log-gamma)
- Error function: `erf(x)`, `erfc(x)`
- Normal distribution: `normcdf(x)`, `norminv(p)`
- Bessel functions: `besselj0`, `besselj1`, `bessely0`, `bessely1`
- Modified Bessel: `besseli0`, `besseli1`, `besselk0`, `besselk1`
- Elliptic integrals: `ellipk(m)`, `ellipe(m)` (AGM method)
- Complex numbers: `abs`, `arg`, `conj`, `real`, `imag`, `sqrt`, `exp`, `log`, `pow`
- Matrices: `det`, `inv`, `transpose`, `trace`, `solve`, `eig`, `qr`

### Phase 3: Professional Features
- Student's t-distribution: `tcdf(t,df)`, `tinv(p,df)`
- Chi-squared: `chi2cdf(x,df)`, `chi2inv(p,df)`
- F-distribution: `fcdf(x,df1,df2)`
- Binomial: `binompdf(k,n,p)`, `binomcdf(k,n,p)`
- Poisson: `poissonpdf(k,lambda)`, `poissoncdf(k,lambda)`
- Numerical integration: `integrate(f,a,b)` (Gauss-Legendre 15-point)
- ODE solver: `ode(f,t0,y0,t1)` (RK4 method)
- Polynomial roots: `polyroots(coeffs)`
- TVM (Time Value of Money): PV, FV, PMT, NPER, RATE

### Phase 4: Advanced
- Lambert W function: `lambertw(x)` - solves W·e^W = x
- Eigenvalues: `eig(A)` - QR iteration method
- QR decomposition: `qr(A)` - Gram-Schmidt
- Polynomial operations: evaluate, derivative, roots

## Operators

| Operator | Description | Precedence |
|----------|-------------|------------|
| `+` `-` | Addition, subtraction | Low |
| `*` `/` `%` | Multiply, divide, modulo | Medium |
| `^` | Exponentiation (right-associative) | High |
| `!` | Factorial (postfix) | Highest |
| `( )` | Grouping | - |
| `=` `==` `!=` `<` `>` `<=` `>=` | Comparison | Lowest |
| `~=` | Approximately equal (within 5%) | Lowest |

## Unicode Operators

```
>>> 2 × 3 ÷ 2
= 3
>>> π × r²
= 3.14159... × r²
>>> √2
= 1.41421356...
```

| Unicode | ASCII | Meaning |
|---------|-------|---------|
| `×` | `*` | Multiply |
| `÷` | `/` | Divide |
| `√` | `sqrt()` | Square root |
| `∛` | `cbrt()` | Cube root |
| `π` | `pi` | Pi constant |

## Functions Reference

### Basic Math
```
sqrt(x)      Square root
cbrt(x)      Cube root  
abs(x)       Absolute value
floor(x)     Floor (round down)
ceil(x)      Ceiling (round up)
round(x)     Round to nearest
trunc(x)     Truncate toward zero
frac(x)      Fractional part
int(x)       Integer part
sign(x)      Sign (-1, 0, or 1)
min(a,b)     Minimum
max(a,b)     Maximum
gcd(a,b)     Greatest common divisor
lcm(a,b)     Least common multiple
```

### Trigonometric (radians by default)
```
sin(x)       Sine
cos(x)       Cosine
tan(x)       Tangent
asin(x)      Arc sine
acos(x)      Arc cosine
atan(x)      Arc tangent
atan2(y,x)   Two-argument arc tangent
sec(x)       Secant
csc(x)       Cosecant
cot(x)       Cotangent
```

### Hyperbolic
```
sinh(x)      Hyperbolic sine
cosh(x)      Hyperbolic cosine
tanh(x)      Hyperbolic tangent
asinh(x)     Inverse hyperbolic sine
acosh(x)     Inverse hyperbolic cosine
atanh(x)     Inverse hyperbolic tangent
```

### Exponential & Logarithmic
```
exp(x)       e^x
log(x)       Natural logarithm (ln)
ln(x)        Natural logarithm
log10(x)     Base-10 logarithm
log2(x)      Base-2 logarithm
log(x,b)     Logarithm base b
pow(x,y)     x raised to power y
```

### Complex Numbers
```
real(z)      Real part
imag(z)      Imaginary part
abs(z)       Magnitude |z|
arg(z)       Argument (phase angle)
conj(z)      Complex conjugate
```

### Special Functions
```
gamma(x)     Gamma function Γ(x)
lgamma(x)    Log-gamma ln|Γ(x)|
erf(x)       Error function
erfc(x)      Complementary error function (1 - erf)
normcdf(x)   Standard normal CDF Φ(x)
norminv(p)   Inverse normal CDF (probit)
besselj0(x)  Bessel function J₀
besselj1(x)  Bessel function J₁
bessely0(x)  Bessel function Y₀
bessely1(x)  Bessel function Y₁
besseli0(x)  Modified Bessel I₀
besseli1(x)  Modified Bessel I₁
besselk0(x)  Modified Bessel K₀
besselk1(x)  Modified Bessel K₁
ellipk(m)    Complete elliptic integral K(m)
ellipe(m)    Complete elliptic integral E(m)
lambertw(x)  Lambert W function W₀(x)
```

### Statistical Distributions
```
normcdf(x)         Standard normal CDF
norminv(p)         Normal inverse CDF
tcdf(t, df)        Student's t CDF
fcdf(x, df1, df2)  F-distribution CDF
chi2cdf(x, df)     Chi-squared CDF
binompdf(k,n,p)    Binomial PMF
binomcdf(k,n,p)    Binomial CDF
poissonpdf(k,λ)    Poisson PMF
poissoncdf(k,λ)    Poisson CDF
```

### Combinatorics
```
n!           Factorial (up to 100000!)
nPr(n,r)     Permutations
nCr(n,r)     Combinations (binomial coefficient)
```

## Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `pi`, `π` | 3.14159265358979323846... | Pi |
| `e` | 2.71828182845904523536... | Euler's number |
| `i` | √(-1) | Imaginary unit |
| `phi` | 1.61803398874989484820... | Golden ratio |
| `Inf` | ∞ | Positive infinity |
| `NaN` | - | Not a number |
| `ans` | - | Previous result |

## Commands

| Command | Description |
|---------|-------------|
| `help` | Show help |
| `digits N` | Set display precision (1-37) |
| `round MODE` | Set rounding: nearest/ceil/floor/trunc/away |
| `fix N` | Fixed decimal places display |
| `sci N` | Scientific notation display |
| `eng N` | Engineering notation display |
| `all` | Show all significant digits |
| `deg` | Angle mode: degrees |
| `rad` | Angle mode: radians (default) |
| `grad` | Angle mode: gradians |
| `vars` | Show variables |
| `funcs` | Show user functions |
| `functions` | Show all available functions |
| `clear` | Clear variables |
| `rpn` | Enter RPN mode |
| `alg` | Enter algebraic mode (default) |
| `tvm` | Time Value of Money mode |
| `stat` | Statistics mode |
| `test` | Run test suite |
| `quit`, `exit` | Exit calculator |

## RPN Mode

Enter RPN mode with the `rpn` command. Uses HP-style 4-level stack (X, Y, Z, T) with LASTX register.

```
rpn> 3 4 + 5 *
= 35

rpn> 2 3 ^
= 8

rpn> lastx
= 3        # Retrieved last X value

rpn> stack
T: 0
Z: 0  
Y: 8
X: 3
```

### RPN Operations
```
enter    Push X, duplicate X into Y
drop     Pop X
swap     Exchange X and Y
clear    Clear stack
dup      Duplicate X
r↓       Roll stack down
r↑       Roll stack up
lastx    Recall last X (before last operation)
stack    Display full stack
alg      Return to algebraic mode
```

## Matrix Operations

```
>>> A = [1,2;3,4]      # 2x2 matrix
>>> det(A)
= -2
>>> inv(A)
>>> A * inv(A)         # Identity matrix
>>> eig(A)             # Eigenvalues
>>> qr(A)              # QR decomposition
>>> solve(A, b)        # Solve Ax = b
```

## Variables and User Functions

```
>>> x = 42
>>> y = x^2 + 1
= 1765

>>> f(x) = x^2 + 2*x + 1
>>> f(3)
= 16

>>> g(x) = sin(x)/x
>>> g(0.1)
= 0.99833416646828...
```

## Plotting

```
>>> plot sin           # Plot sin(x) from -2π to 2π
>>> plot x^2 -5:5      # Plot x² from -5 to 5
>>> plot sin(x)*exp(-x/5) 0:20
```

## Statistics Mode

```
stat> clear            # Clear data
stat> 10 15 20 25 30   # Enter data points
stat> mean             # = 20
stat> sdev             # Standard deviation
stat> var              # Variance
stat> sum              # Sum of values
stat> n                # Count
```

## Installation

### From Source
```bash
make
sudo make install      # Installs to /usr/local/bin
```

### RPM Package (Fedora/RHEL/CentOS)
```bash
make rpm
sudo dnf install ~/rpmbuild/RPMS/x86_64/sc-*.rpm
```

### DOS Build (Watcom C)
```bash
# Cross-compile or use Watcom C on DOS
wmake -f Makefile.wat CFLAGS=-DSCALC_MEDIUM
```

## Testing

```bash
# Run all tests (624 tests)
make test

# Individual test suites
./tests/test_bedmas.sh ./bin/sc      # BEDMAS order of operations
./tests/test_core.sh ./bin/sc        # Core functionality
./tests/test_ieee754.sh ./bin/sc     # IEEE 754-2008 compliance

# IEEE 754 test vectors (requires download)
cd tests && ./download_ieee754_suite.sh
python3 parse_fptest.py ../bin/sc ieee754-test-suite
```

## Architecture

```
sc/
├── src/
│   ├── apf.c/h         # Arbitrary Precision Float (pure C89, no deps)
│   ├── apfc.c/h        # Complex number extension
│   ├── apfx.c/h        # Extended math functions (trig, exp, special)
│   ├── special_funcs.c/h  # Phase 2-4 special functions
│   ├── matrix*.c/h     # Matrix operations and linear algebra
│   ├── rpn.c/h         # RPN stack engine
│   ├── parser.c        # Expression parser
│   ├── lexer.c         # Token lexer
│   ├── stats.c/h       # Statistical functions
│   ├── tvm.c/h         # Time Value of Money
│   ├── config.h        # Build configuration and feature flags
│   └── ...
├── tests/              # Test suites
├── packaging/
│   ├── sc.1            # Man page
│   └── sc.spec         # RPM spec
└── Makefile
```

### Core Library (Zero Dependencies)

The APF (Arbitrary Precision Float) library is completely standalone:
- No stdio, stdlib, string.h, float, or double
- Pure integer arithmetic using 16-bit limbs
- Works on any platform with 16-bit and 32-bit integer support
- IEEE 754-2008 compliant rounding modes

## Technical Notes

### Precision
- 128-bit mantissa = ~38 decimal digits internal precision
- Display precision configurable 1-37 digits
- IEEE 754-2008 quadruple precision compatible

### Algorithms
- **Trig functions**: Taylor series with argument reduction (16-bit+), CORDIC (8-bit)
- **Logarithm**: ln((1+z)/(1-z)) series with range reduction
- **Square root**: Newton-Raphson (Heron's method)
- **Gamma**: Lanczos approximation (g=7, 9 coefficients)
- **Bessel**: Power series for small x, asymptotic for large x
- **Elliptic integrals**: Arithmetic-Geometric Mean (AGM)
- **Eigenvalues**: QR iteration with Wilkinson shift
- **Polynomial roots**: Durand-Kerner iteration

### Special Values
```
1/0 = Inf           # Positive infinity
-1/0 = -Inf         # Negative infinity
0/0 = NaN           # Not a number
Inf - Inf = NaN
0^0 = 1             # IEEE 754-2008
1^Inf = 1           # IEEE 754-2008
Inf^0 = 1           # IEEE 754-2008
```

## Comparison with Other Calculators

See [COMPARISON.md](COMPARISON.md) for detailed comparison with:
- HP-15C, HP-42S, HP-48/49/50
- TI-89, TI-Nspire
- Casio fx-991EX
- Free42, Qalculate!, SpeedCrunch

## License

MIT License - see [LICENSE](LICENSE)

## Author

Stu - https://github.com/dr-who/sc

## Acknowledgments

- HP calculator community for RPN inspiration
- Fabrice Bellard's LibBF for soft-float design patterns
- NIST Digital Library of Mathematical Functions for algorithm references
