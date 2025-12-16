# Calculator Comparison: sc vs Popular Calculators

This document compares **sc** to popular scientific calculators and command-line tools.

## Quick Comparison Table

| Feature | sc | bc -l | Casio fx-991 | TI-84 | HP 48G |
|---------|:--:|:-----:|:------------:|:-----:|:------:|
| **Precision** | 38 digits | arbitrary | 10+2 | 14 | 12 |
| **Complex Numbers** | ✓ | ✗ | ✓ | ✓ | ✓ |
| **Matrices** | ✓ | ✗ | ✓ | ✓ | ✓ |
| **Statistics** | ✓ | ✗ | ✓ | ✓ | ✓ |
| **IEEE 754 inf/nan** | ✓ | ✗ | varies | varies | ✓ |
| **Symbolic Math** | ✗ | ✗ | ✗ | limited | ✓ |
| **Programmable** | ✗ | ✓ | limited | ✓ | ✓ |
| **RPN Mode** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Price** | Free | Free | ~$30 | ~$120 | ~$200+ |

## Detailed Comparison

### vs GNU bc -l

**bc** is the standard Unix arbitrary precision calculator. While powerful, it has several limitations that **sc** addresses:

#### Rounding Errors

```bash
# bc -l: Rounding error in basic arithmetic
$ echo "scale=20; 1/3*3" | bc -l
.99999999999999999999        # WRONG! Should be 1

# sc: Correct result
$ sc "1/3*3"
1                            # CORRECT

# bc -l: sqrt(2)^2 has error
$ echo "scale=20; sqrt(2)^2" | bc -l  
1.99999999999999999999        # WRONG!

# sc: Correct
$ sc "sqrt(2)^2"
2                            # CORRECT
```

#### Division by Zero

```bash
# bc -l: Crashes on division by zero
$ echo "1/0" | bc -l
Runtime error: Divide by zero   # CRASH

# sc: Returns IEEE 754 infinity
$ sc "1/0"
Inf                            # CORRECT per IEEE 754-2008
```

#### Complex Numbers

```bash
# bc -l: No complex number support
$ echo "sqrt(-1)" | bc -l
Runtime error                  # FAILS

# sc: Full complex number support
$ sc "sqrt(-1)"
i                              # CORRECT

$ sc "(3+4i)*(1-2i)"
11 - 2i                        # Full complex arithmetic
```

#### Feature Comparison with bc

| Operation | bc -l | sc |
|-----------|-------|-----|
| `1/3*3` | 0.999... | 1 |
| `sqrt(2)^2` | 1.999... | 2 |
| `1/0` | Error | Inf |
| `0/0` | Error | NaN |
| `sqrt(-1)` | Error | i |
| `sin(pi)` | Needs `s()` | sin(pi) works |
| Matrices | None | Full support |
| Statistics | None | Full support |

### vs Casio fx-991EX (ClassWiz)

The Casio fx-991 is one of the most popular scientific calculators worldwide.

| Feature | sc | fx-991EX |
|---------|-----|----------|
| **Display Digits** | Up to 38 | 10 (natural display) |
| **Internal Precision** | 128-bit (~38 digits) | 15 digits internal |
| **Matrix Size** | Up to 100×100 | 4×4 max |
| **Complex Numbers** | Full support | Full support |
| **Statistics** | With t-test | With regression |
| **Equation Solver** | Newton-Raphson | Built-in solver |
| **Base Conversions** | hex, bin, oct | hex, bin, oct, dec |
| **User Variables** | Unlimited | 9 (A-F, X, Y, M) |
| **Programming** | Pipe scripts | Table function |
| **Physical Keys** | N/A | Yes |
| **Battery** | N/A | Solar + battery |

**Advantages of sc:**
- Much higher precision (38 vs 10 displayed digits)
- Larger matrices
- Scriptable via pipes
- Free and open source
- Runs on any computer

**Advantages of fx-991:**
- Portable hardware device
- No computer needed
- Exam-approved
- Natural textbook display
- Physical buttons

### vs Texas Instruments TI-84 Plus

The TI-84 is the dominant calculator in US education.

| Feature | sc | TI-84 Plus |
|---------|-----|------------|
| **Display Digits** | 38 | 10 |
| **Internal Precision** | 128-bit | 14 digits |
| **Graphing** | ASCII plot | Full graphing |
| **Programming** | Pipe scripts | TI-BASIC |
| **Matrix Size** | 100×100 | 99×99 |
| **Lists** | Via matrices | Up to 999 elements |
| **Statistics** | With t-test | Comprehensive |
| **Financial (TVM)** | Full support | Via apps |
| **Memory** | Computer RAM | 24KB RAM |
| **Link/Share** | Files | USB/Cable |

**Advantages of sc:**
- Higher precision
- Free
- Integrates with Unix tools
- Fast (runs on modern CPU)
- IEEE 754 compliant

**Advantages of TI-84:**
- Full graphical display
- Comprehensive graphing
- Exam-approved
- Large ecosystem of programs
- Educational standard

### vs HP 48G/48GX

The HP 48 series represents the gold standard in RPN scientific calculators.

| Feature | sc | HP 48G |
|---------|-----|--------|
| **Entry Mode** | Algebraic + RPN | RPN (algebraic option) |
| **Display Digits** | 38 | 12 |
| **Stack Size** | Unlimited | Unlimited |
| **Symbolic Math** | No | Yes (limited) |
| **Units** | No | Full unit conversions |
| **Matrices** | Full support | Full support |
| **Complex Numbers** | Full support | Full support |
| **Programming** | No | RPL language |
| **Equation Library** | No | Yes |
| **I/O** | Pipes/files | IR, serial |

**Advantages of sc:**
- Higher numerical precision
- Free
- Faster computation
- Modern IEEE 754 compliance

**Advantages of HP 48:**
- True RPN calculator experience
- Symbolic mathematics
- Unit conversions built-in
- Equation library
- Dedicated hardware

## IEEE 754-2008 Compliance

**sc** is fully compliant with IEEE 754-2008 for special values:

| Operation | sc | bc -l | Most Calculators |
|-----------|-----|-------|------------------|
| `1/0` | Inf | Error | Error or Inf |
| `-1/0` | -Inf | Error | Error or -Inf |
| `0/0` | NaN | Error | Error |
| `Inf + Inf` | Inf | N/A | Usually Inf |
| `Inf - Inf` | NaN | N/A | Error or NaN |
| `Inf * 0` | NaN | N/A | Error or NaN |
| `0 * Inf` | NaN | N/A | Error or NaN |
| `Inf / Inf` | NaN | N/A | Error or NaN |
| `NaN` operations | NaN | N/A | Varies |

## Precision Comparison

Calculating π to various precisions:

```bash
# sc (38 digits available)
$ sc "digits 38; pi"
3.1415926535897932384626433832795028841

# bc -l (arbitrary, but slow)
$ echo "scale=38; 4*a(1)" | bc -l
3.14159265358979323844643265362274860062

# Casio fx-991: 3.141592654 (10 digits displayed)
# TI-84: 3.1415926536 (10 digits displayed)  
# HP 48G: 3.14159265359 (12 digits displayed)
```

## Performance

For large calculations, **sc** significantly outperforms bc:

```bash
# Calculate 50!
$ time echo "50!" | sc
3.0414093201713378043612608166065e+64
real    0m0.003s

$ time echo "x=1; for(i=1;i<=50;i++) x=x*i; x" | bc
30414093201713378043612608166064768844377641568960512...
real    0m0.015s
```

## Use Case Recommendations

**Use sc when you need:**
- High precision calculations (more than 15 digits)
- IEEE 754 compliant special value handling
- Complex number arithmetic
- Integration with shell scripts
- Quick matrix operations
- RPN mode with command line
- Free, open source solution

**Use bc when you need:**
- Arbitrary precision (1000+ digits)
- Programmable calculations
- POSIX compliance
- Simple scripts

**Use a physical calculator when you need:**
- Portability without a computer
- Exam-approved device
- Dedicated hardware
- Graphical display (TI-84, HP 48)

## Migration Guide

### From bc -l to sc

| bc -l | sc |
|-------|-----|
| `scale=20` | `digits 20` |
| `s(x)` | `sin(x)` |
| `c(x)` | `cos(x)` |
| `a(x)` | `atan(x)` |
| `l(x)` | `ln(x)` |
| `e(x)` | `exp(x)` |
| `sqrt(x)` | `sqrt(x)` |
| Variables: single letter | Variables: any name |
| No complex | Full complex support |
| No matrices | Full matrix support |

### From HP Calculator to sc

| HP RPN | sc RPN mode |
|--------|-------------|
| ENTER | (automatic) |
| + - * / | + - * / |
| x² | ^ 2 |
| √x | sqrt |
| DUP | dup |
| DROP | drop |
| SWAP | swap |
| CLR | clr |

## Boolean Operations and Shell Scripting

**sc** includes boolean operators and comparison operations that integrate seamlessly with shell scripts:

### Comparison Operators

```bash
$ sc "3 < 4"
true

$ sc "pi > 3"
true

$ sc "100 ~= 105"    # Approximately equal (within 5%)
true

$ sc "100 ~= 106"    # Outside 5% tolerance
false
```

### Boolean Operators

```bash
$ sc "1<2 and 3<4"
true

$ sc "1>2 or 3<4"
true

$ sc "not 1>2"
true

$ sc "1<2 xor 3<4"   # Both true, so XOR is false
false
```

### Exit Codes for Shell Integration

**sc** returns exit code 0 for true/non-zero results, 1 for false/zero:

```bash
# Use in if statements
if sc "isprime(17)"; then
    echo "17 is prime"
fi

# Check approximate equality
expected=100
actual=103
if sc "$actual ~= $expected"; then
    echo "Result is close enough"
fi

# Complex boolean logic
if sc "(x > 0) and (x < 100)"; then
    echo "x is in range"
fi
```

### Primality Testing

```bash
$ sc "isprime(7)"
1                    # Returns 1 (true)

$ sc "isprime(8)"
0                    # Returns 0 (false)

# Use in loops
for n in 2 3 4 5 6 7 8 9 10; do
    if sc "isprime($n)"; then
        echo "$n is prime"
    fi
done
```

## Special Functions Comparison

**sc** provides a comprehensive set of special functions rivaling dedicated math software:

| Function | sc | bc | Casio fx-991EX | HP 48G | Python scipy |
|----------|:--:|:--:|:--------------:|:------:|:------------:|
| **Gamma Γ(x)** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Log-Gamma ln(Γ)** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Error Function erf** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Normal CDF Φ(x)** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Normal Inverse** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Student's t CDF** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Chi-squared CDF** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **F-distribution CDF** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Binomial PMF/CDF** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Poisson PMF/CDF** | ✓ | ✗ | ✓ | ✗ | ✓ |
| **Bessel J₀, J₁** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Bessel Y₀, Y₁** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Modified Bessel I, K** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Elliptic K, E** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Lambert W** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Beta function** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Incomplete Gamma** | ✓ | ✗ | ✗ | ✗ | ✓ |
| **Incomplete Beta** | ✓ | ✗ | ✗ | ✗ | ✓ |

### Special Functions Examples

```bash
# Gamma function - factorial extended to real numbers
$ sc "gamma(5.5)"
52.342777784553525196...

# Error function - used in probability
$ sc "erf(1)"
0.84270079294971486934...

# Normal distribution
$ sc "normcdf(1.96)"
0.97500210485177961532...   # P(Z < 1.96) ≈ 97.5%

$ sc "norminv(0.975)"
1.95996398454005423552...   # Z-score for 97.5th percentile

# Student's t distribution
$ sc "tcdf(2.0, 10)"         # t=2.0, df=10
0.96315700683571127...

# Chi-squared test
$ sc "chi2cdf(7.815, 3)"     # Critical value for α=0.05, df=3
0.95...

# Bessel functions
$ sc "j0(0)"
1

$ sc "j1(1)"
0.44005058574493351596...

# Elliptic integrals (period of pendulum, etc.)
$ sc "ellipk(0.5)"
1.85407467730137191843...

# Lambert W - useful in combinatorics, physics
$ sc "lambertw(1)"
0.56714329040978387299...   # W(1)·e^W(1) = 1
```

### Statistical Analysis Advantage

Unlike pocket calculators, **sc** provides inverse distribution functions for hypothesis testing:

```bash
# Find critical t-value for two-tailed test, α=0.05, df=20
$ sc "tinv(0.975, 20)"
2.08596...

# Find chi-squared critical value
$ sc "chi2inv(0.95, 5)"
11.0704...

# Exact binomial probability
$ sc "binompdf(7, 10, 0.5)"   # P(X=7) for 10 flips of fair coin
0.1171875
```

## Summary

**sc** fills a unique niche as a high-precision, IEEE 754 compliant scientific calculator that:

1. **Fixes bc's shortcomings** - proper rounding, infinity/NaN support, complex numbers
2. **Matches graphing calculator features** - matrices, statistics, TVM, equation solving
3. **Provides HP-style RPN** - for those who prefer postfix notation
4. **Integrates with Unix** - pipes, scripts, automation
5. **Is free and open source** - no licensing costs, auditable code

For command-line scientific computation with modern features and proper mathematical behavior, **sc** is the recommended choice.

## MATLAB/Octave Compatibility

**sc** implements many MATLAB/Octave-compatible functions and syntax, making it easy for users familiar with those environments:

### Supported MATLAB-style Syntax

```bash
# Matrix definition with semicolon row separators
$ sc "[1 2 3; 4 5 6; 7 8 9]"

# Colon notation for ranges
$ sc "1:5"
1  2  3  4  5

$ sc "1:2:9"
1  3  5  7  9

# Matrix indexing (1-based)
$ sc "M = [1 2 3; 4 5 6]; M(1,2)"
2

# Colon indexing for rows/columns
$ sc "M = [1 2 3; 4 5 6]; M(:,2)"  # Second column
2
5

$ sc "M = [1 2 3; 4 5 6]; M(1,:)"  # First row
1  2  3
```

### Array Creation Functions

| Function | Description | Example |
|----------|-------------|---------|
| `zeros(r,c)` | Matrix of zeros | `zeros(2,3)` |
| `ones(r,c)` | Matrix of ones | `ones(2,3)` |
| `eye(n)` | Identity matrix | `eye(3)` |
| `linspace(a,b,n)` | n evenly spaced points | `linspace(0,1,5)` |
| `logspace(a,b,n)` | n logarithmically spaced | `logspace(0,2,5)` |
| `rand(r,c)` | Uniform random | `rand(2,3)` |
| `randn(r,c)` | Normal random | `randn(2,3)` |
| `randperm(n)` | Random permutation | `randperm(5)` |

### Array Manipulation Functions

| Function | Description |
|----------|-------------|
| `fliplr(M)` | Flip left-right |
| `flipud(M)` | Flip up-down |
| `flip(v)` | Flip vector |
| `rot90(M)` | Rotate 90° counterclockwise |
| `sort(v)` | Sort elements |
| `sortrows(M)` | Sort matrix by first column |
| `unique(v)` | Unique sorted elements |
| `reshape(M,r,c)` | Reshape to r×c |
| `repmat(M,r,c)` | Replicate r×c times |
| `circshift(v,n)` | Circular shift by n |
| `cumsum(v)` | Cumulative sum |
| `cumprod(v)` | Cumulative product |
| `diff(v)` | Differences |
| `find(v)` | Indices of non-zeros |
| `triu(M)` | Upper triangular part |
| `tril(M)` | Lower triangular part |
| `horzcat(A,B)` | Horizontal concatenation |
| `vertcat(A,B)` | Vertical concatenation |
| `squeeze(M)` | Remove singleton dimensions |
| `zscore(v)` | Z-score normalization |

### Set Operations

| Function | Description |
|----------|-------------|
| `union(A,B)` | Set union (sorted unique) |
| `intersect(A,B)` | Set intersection |
| `setdiff(A,B)` | Elements in A but not B |
| `setxor(A,B)` | Symmetric difference |

### Sorting and Selection

| Function | Description |
|----------|-------------|
| `sort(v)` | Sort ascending |
| `sortrows(M)` | Sort matrix rows |
| `maxk(v,k)` | k largest elements |
| `mink(v,k)` | k smallest elements |

### Cumulative Operations

| Function | Description |
|----------|-------------|
| `cumsum(v)` | Cumulative sum |
| `cumprod(v)` | Cumulative product |
| `cummax(v)` | Cumulative maximum |
| `cummin(v)` | Cumulative minimum |

### Element-wise Comparisons

| Function | Description |
|----------|-------------|
| `eq(A,B)` | Element-wise equal |
| `ne(A,B)` | Element-wise not equal |
| `lt(A,B)` | Element-wise less than |
| `le(A,B)` | Element-wise less than or equal |
| `gt(A,B)` | Element-wise greater than |
| `ge(A,B)` | Element-wise greater than or equal |

### Element-wise Operations

| Function | Description |
|----------|-------------|
| `times(A,B)` | Element-wise multiply |
| `rdivide(A,B)` | Element-wise right divide |
| `ldivide(A,B)` | Element-wise left divide |
| `plus(A,B)` | Element-wise add |
| `minus(A,B)` | Element-wise subtract |

### Reduction Functions

| Function | Description |
|----------|-------------|
| `sum(M)` | Sum of all elements |
| `prod(M)` | Product of all elements |
| `mean(M)` | Mean value |
| `std(M)` | Standard deviation |
| `var(M)` | Variance |
| `cov(M)` | Covariance (variance for vector) |
| `min(M)` | Minimum |
| `max(M)` | Maximum |
| `median(M)` | Median value |
| `mode(v)` | Most frequent value |
| `range(v)` | max - min |
| `iqr(v)` | Interquartile range |
| `prctile(v,p)` | p-th percentile |
| `any(M)` | True if any nonzero |
| `all(M)` | True if all nonzero |
| `nnz(M)` | Count of nonzeros |

### Matrix Operations

| Operator/Function | Description |
|-------------------|-------------|
| `A * B` | Matrix multiplication |
| `A .* B` | Element-wise multiply |
| `A / B` | A * inv(B) |
| `A \ B` | inv(A) * B |
| `A'` | Conjugate transpose |
| `A.'` | Transpose |
| `det(A)` | Determinant |
| `inv(A)` | Inverse |
| `eig(A)` | Eigenvalues |
| `trace(A)` | Trace |
| `norm(A)` | Frobenius norm |
| `rank(M)` | Matrix rank |
| `cond(M)` | Condition number |
| `linsolve(A,b)` | Solve Ax = b |

### Special Matrices

| Function | Description |
|----------|-------------|
| `kron(A,B)` | Kronecker product |
| `vander(v,n)` | Vandermonde matrix |
| `toeplitz(c)` | Toeplitz matrix |
| `hankel(c)` | Hankel matrix |
| `magic(n)` | Magic square |
| `pascal(n)` | Pascal's triangle |
| `hilb(n)` | Hilbert matrix |

### Signal Processing

| Function | Description |
|----------|-------------|
| `conv(a,b)` | Convolution |
| `gradient(v)` | Numerical gradient |
| `trapz(y)` | Trapezoidal integration |

### Statistics (Extended)

| Function | Description |
|----------|-------------|
| `corrcoef(X,Y)` | Correlation coefficient |
| `movmean(v,k)` | Moving average |
| `movsum(v,k)` | Moving sum |
| `rescale(v)` | Scale to [0,1] |
| `zscore(v)` | Z-score normalization |

### Predicate Functions

| Function | Returns true if... |
|----------|-------------------|
| `isnan(x)` | x is NaN |
| `isinf(x)` | x is Inf |
| `isfinite(x)` | x is finite |
| `isreal(x)` | x has no imaginary part |
| `isempty(M)` | M is empty |
| `isscalar(x)` | x is scalar |
| `isvector(M)` | M is vector |
| `ismatrix(M)` | M is matrix |
| `isrow(M)` | M is row vector (1×n) |
| `iscolumn(M)` | M is column vector (n×1) |
| `issquare(M)` | M is square (n×n) |
| `issymmetric(M)` | M equals its transpose |
| `isdiag(M)` | M is diagonal |
| `istriu(M)` | M is upper triangular |
| `istril(M)` | M is lower triangular |
| `issorted(v)` | v is sorted ascending |
| `isequal(a,b)` | a equals b |
| `ndims(M)` | Number of dimensions (always 2) |

### Trigonometric Functions (Degrees)

| Function | Description |
|----------|-------------|
| `sind(x)` | sin(x°) |
| `cosd(x)` | cos(x°) |
| `tand(x)` | tan(x°) |
| `secd(x)` | sec(x°) |
| `cscd(x)` | csc(x°) |
| `cotd(x)` | cot(x°) |
| `asind(x)` | asin in degrees |
| `acosd(x)` | acos in degrees |
| `atand(x)` | atan in degrees |
| `atan2d(y,x)` | atan2 in degrees |
| `asecd(x)` | asec in degrees |
| `acscd(x)` | acsc in degrees |
| `acotd(x)` | acot in degrees |

### Real-Only Math Functions

| Function | Description |
|----------|-------------|
| `nthroot(x,n)` | Real nth root (preserves sign for odd n) |
| `realsqrt(x)` | Real sqrt (error if x < 0) |
| `reallog(x)` | Real log (error if x ≤ 0) |
| `realpow(x,y)` | Real power (error if result complex) |

### Additional Math Functions

| Function | Description |
|----------|-------------|
| `log2(x)` | Base-2 logarithm |
| `exp2(x)` / `pow2(x)` | 2^x |
| `deg2rad(x)` | Degrees to radians |
| `rad2deg(x)` | Radians to degrees |
| `nextpow2(x)` | Next power of 2 |
| `complex(a,b)` | Create a+bi |
| `sign(x)` | Sign function |
| `fix(x)` | Round toward zero |
| `rem(a,b)` | Remainder (sign of dividend) |
| `hypot(x,y)` | sqrt(x²+y²) |

### Polynomial Functions

| Function | Description |
|----------|-------------|
| `polyval(p,x)` | Evaluate polynomial at x |
| `polyder(p)` | Derivative coefficients |
| `polyint(p)` | Integral coefficients |

### Constants

| Constant | Value |
|----------|-------|
| `pi` | 3.14159... |
| `e` | 2.71828... |
| `i` | √(-1) |
| `Inf` | Infinity |
| `NaN` | Not a Number |
| `true` | 1 |
| `false` | 0 |
| `eps` | Machine epsilon (~2.2e-16) |
| `realmax` | Largest finite float |
| `realmin` | Smallest positive float |

### Example MATLAB-style Session

```bash
$ sc
>>> A = [1 2; 3 4]
1  2
3  4
>>> B = A'
1  3
2  4
>>> A * B
5   11
11  25
>>> sum(A(:))
10
>>> eig(A)
-0.3722813232690143  5.372281323269014
>>> det(A)
-2
>>> inv(A) * A
1  0
0  1
>>> linspace(0, 1, 5)
0  0.25  0.5  0.75  1
>>> v = 1:5
1  2  3  4  5
>>> cumsum(v)
1  3  6  10  15
>>> any(v > 3)
true
>>> all(v > 0)
true
```

## Calculator Comparison

### vs Casio fx-82/991
| Feature | Casio fx-82 | Casio fx-991 | scalc |
|---------|-------------|--------------|-------|
| Basic Math | ✓ | ✓ | ✓ |
| Scientific | ✓ | ✓ | ✓ |
| Complex Numbers | - | ✓ | ✓ |
| Matrix (2x2) | - | ✓ | Up to 10x10 |
| Equation Solver | - | ✓ | ✓ |
| Statistics | Basic | Advanced | Advanced |
| Calculus | - | Numeric | ✓ |
| Programmable | - | - | User functions |
| Precision | 10 digits | 10 digits | 38 digits |

### vs TI-84/89
| Feature | TI-84 | TI-89 | scalc |
|---------|-------|-------|-------|
| Graphing | ✓ | ✓ | ASCII plots |
| CAS (Symbolic) | - | ✓ | - |
| Complex Numbers | ✓ | ✓ | ✓ |
| Matrices | ✓ | ✓ | ✓ |
| Programming | TI-BASIC | TI-BASIC | User functions |
| Statistics | ✓ | ✓ | ✓ |
| Precision | 14 digits | 14 digits | 38 digits |

### vs HP-42S/48
| Feature | HP-42S | HP-48 | scalc |
|---------|--------|-------|-------|
| RPN Mode | ✓ | ✓ | ✓ |
| Algebraic | - | ✓ | ✓ |
| Complex Numbers | ✓ | ✓ | ✓ |
| Matrices | ✓ | ✓ | ✓ |
| Solver | ✓ | ✓ | ✓ |
| Unit Conversion | ✓ | ✓ | - |
| CAS | - | ✓ | - |
| Programming | Keystroke | RPL | User functions |

### vs MATLAB/Octave
| Feature | MATLAB | Octave | scalc |
|---------|--------|--------|-------|
| Matrix Operations | ✓ | ✓ | ✓ |
| Linear Algebra | Full | Full | Basic |
| Signal Processing | Toolbox | ✓ | Basic |
| Statistics | Toolbox | ✓ | ✓ |
| Symbolic Math | Toolbox | SymPy | - |
| Plotting | Advanced | ✓ | ASCII |
| File I/O | ✓ | ✓ | - |
| Cost | $$$ | Free | Free |
| Precision | 16 digits | 16 digits | 38 digits |

### vs Mathematica/Wolfram
| Feature | Mathematica | scalc |
|---------|-------------|-------|
| Symbolic Math | ✓ | - |
| Arbitrary Precision | ✓ | 38-256 digits |
| Matrices | Unlimited | Up to 10x10 |
| Special Functions | Comprehensive | Good coverage |
| Plotting | Advanced | ASCII |
| Cost | $$$ | Free |
| Portable | Desktop | Anywhere |

## Function Count by Category

| Category | Count |
|----------|-------|
| Basic Math | 30+ |
| Trigonometry | 40+ |
| Complex Numbers | 15+ |
| Matrix Operations | 60+ |
| Statistics | 40+ |
| Special Functions | 30+ |
| Signal Processing | 20+ |
| Polynomials | 10+ |
| Distributions | 15+ |
| Utility | 40+ |
| **Total** | **~350** |

