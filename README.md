# scalc - Scientific Calculator

A MATLAB-compatible scientific calculator with 128-bit precision, 640+ built-in functions, matrix operations, complex numbers, ML classifiers, and time series forecasting. Runs on modern systems and DOS 2.11+.

## Features

- **640+ built-in functions** across 25+ categories
- **128-bit precision** (37 significant digits)
- **Complex number support** with `i` notation
- **Matrix operations** with MATLAB-style syntax
- **Machine learning** classifiers (k-NN, Naive Bayes, SVM, Decision Tree)
- **Time series** forecasting (movavg, ewma, trend, forecast, autocorr)
- **Tables and timetables** with groupsummary, head, tail
- **RPN mode** for HP calculator enthusiasts
- **DOS compatible** - runs on 8088 with 256KB RAM

## Quick Start

```bash
make && ./bin/sc
```

```matlab
>>> 2+2
= 4
>>> sin(pi/2)
= 1
>>> A = [1,2;3,4]; inv(A)
>>> exp(i*pi) + 1      % Euler's identity = 0
>>> help sin           % Function help
>>> demo matrix        % Interactive demo
```

## Syntax Notes

- `%` starts a comment (MATLAB-style). Use `mod(a,b)` for modulo.
- Matrices: `[1,2;3,4]` creates 2x2 matrix
- Complex: `3+4i` or `complex(3,4)`

## Function Reference (640+ functions)

### Constants (12)

| Function | Description | Example |
|----------|-------------|---------|
| `pi` | π = 3.14159... | `pi` → 3.14159... |
| `e` | Euler's number | `e` → 2.71828... |
| `i` | Imaginary unit | `i^2` → -1 |
| `Inf` | Infinity | `1/0` → Inf |
| `NaN` | Not a Number | `0/0` → NaN |
| `true` | Boolean true (1) | `true` → 1 |
| `false` | Boolean false (0) | `false` → 0 |
| `eps` | Machine epsilon | `eps` → ~1e-37 |
| `realmax` | Largest number | `realmax` |
| `realmin` | Smallest positive | `realmin` |
| `phi` / `golden` | Golden ratio | `phi` → 1.618... |
| `tau` | 2π | `tau` → 6.283... |

### Trigonometry (54)

| Function | Description | Example |
|----------|-------------|---------|
| `sin(x)` | Sine (radians) | `sin(pi/2)` → 1 |
| `cos(x)` | Cosine | `cos(0)` → 1 |
| `tan(x)` | Tangent | `tan(pi/4)` → 1 |
| `asin(x)` | Arcsine | `asin(1)` → π/2 |
| `acos(x)` | Arccosine | `acos(0)` → π/2 |
| `atan(x)` | Arctangent | `atan(1)` → π/4 |
| `atan2(y,x)` | Four-quadrant arctan | `atan2(1,1)` → π/4 |
| `sinh(x)` | Hyperbolic sine | `sinh(0)` → 0 |
| `cosh(x)` | Hyperbolic cosine | `cosh(0)` → 1 |
| `tanh(x)` | Hyperbolic tangent | `tanh(0)` → 0 |
| `asinh(x)` | Inverse hyperbolic sine | `asinh(0)` → 0 |
| `acosh(x)` | Inverse hyperbolic cosine | `acosh(1)` → 0 |
| `atanh(x)` | Inverse hyperbolic tangent | `atanh(0)` → 0 |
| `sec(x)` | Secant | `sec(0)` → 1 |
| `csc(x)` | Cosecant | `csc(pi/2)` → 1 |
| `cot(x)` | Cotangent | `cot(pi/4)` → 1 |
| `asec(x)` | Arcsecant | `asec(1)` → 0 |
| `acsc(x)` | Arccosecant | `acsc(1)` → π/2 |
| `acot(x)` | Arccotangent | `acot(1)` → π/4 |
| `sech(x)` | Hyperbolic secant | `sech(0)` → 1 |
| `csch(x)` | Hyperbolic cosecant | `csch(1)` |
| `coth(x)` | Hyperbolic cotangent | `coth(1)` |
| `asech(x)` | Inverse hyperbolic secant | `asech(1)` → 0 |
| `acsch(x)` | Inverse hyperbolic cosecant | `acsch(1)` |
| `acoth(x)` | Inverse hyperbolic cotangent | `acoth(2)` |
| `sind(x)` | Sine (degrees) | `sind(90)` → 1 |
| `cosd(x)` | Cosine (degrees) | `cosd(0)` → 1 |
| `tand(x)` | Tangent (degrees) | `tand(45)` → 1 |
| `asind(x)` | Arcsine (degrees) | `asind(1)` → 90 |
| `acosd(x)` | Arccosine (degrees) | `acosd(0)` → 90 |
| `atand(x)` | Arctangent (degrees) | `atand(1)` → 45 |
| `atan2d(y,x)` | Four-quadrant arctan (deg) | `atan2d(1,1)` → 45 |
| `sinc(x)` | Sinc function sin(πx)/(πx) | `sinc(0)` → 1 |
| `sinpi(x)` | sin(πx) | `sinpi(0.5)` → 1 |
| `cospi(x)` | cos(πx) | `cospi(0.5)` → 0 |
| Aliases: `arcsin`, `arccos`, `arctan`, `arcsinh`, `arccosh`, `arctanh`, `arcsec`, `arccsc`, `arccot`, `cosec`, `cotan`, `cotd`, `cscd`, `secd`, `acotd`, `acscd`, `asecd` |

### Exponential (10)

| Function | Description | Example |
|----------|-------------|---------|
| `exp(x)` | e^x | `exp(1)` → e |
| `exp2(x)` | 2^x | `exp2(3)` → 8 |
| `expm1(x)` | e^x - 1 (accurate for small x) | `expm1(0)` → 0 |
| `pow2(x)` | 2^x | `pow2(10)` → 1024 |
| `sqrt(x)` | Square root | `sqrt(4)` → 2 |
| `cbrt(x)` | Cube root | `cbrt(27)` → 3 |
| `nthroot(x,n)` | nth root | `nthroot(16,4)` → 2 |
| `hypot(a,b)` | √(a²+b²) | `hypot(3,4)` → 5 |
| `realsqrt(x)` | Real square root | `realsqrt(4)` → 2 |
| `realpow(x,y)` | Real power | `realpow(2,3)` → 8 |

### Logarithmic (7)

| Function | Description | Example |
|----------|-------------|---------|
| `log(x)` / `ln(x)` | Natural logarithm | `ln(e)` → 1 |
| `log10(x)` | Base-10 logarithm | `log10(100)` → 2 |
| `log2(x)` | Base-2 logarithm | `log2(8)` → 3 |
| `log1p(x)` | ln(1+x) (accurate for small x) | `log1p(0)` → 0 |
| `logb(x,b)` | Base-b logarithm | `logb(8,2)` → 3 |
| `reallog(x)` | Real logarithm | `reallog(e)` → 1 |

### Complex Numbers (12)

| Function | Description | Example |
|----------|-------------|---------|
| `real(z)` / `re(z)` | Real part | `real(3+4i)` → 3 |
| `imag(z)` / `im(z)` | Imaginary part | `imag(3+4i)` → 4 |
| `conj(z)` | Complex conjugate | `conj(3+4i)` → 3-4i |
| `abs(z)` | Magnitude | `abs(3+4i)` → 5 |
| `arg(z)` / `angle(z)` / `phase(z)` | Argument (angle) | `arg(1+i)` → π/4 |
| `complex(a,b)` | Create complex number | `complex(3,4)` → 3+4i |
| `cart2pol(x,y)` | Cartesian to polar | `cart2pol(1,1)` |
| `pol2cart(r,θ)` | Polar to Cartesian | `pol2cart(1,pi/4)` |

### Rounding (13)

| Function | Description | Example |
|----------|-------------|---------|
| `floor(x)` | Round toward -∞ | `floor(3.7)` → 3 |
| `ceil(x)` | Round toward +∞ | `ceil(3.2)` → 4 |
| `round(x)` | Round to nearest | `round(3.5)` → 4 |
| `trunc(x)` / `fix(x)` | Round toward zero | `trunc(-3.7)` → -3 |
| `frac(x)` | Fractional part | `frac(3.7)` → 0.7 |
| `sign(x)` / `signum(x)` | Sign (-1, 0, 1) | `sign(-5)` → -1 |
| `mod(a,b)` | Modulo (sign of b) | `mod(17,5)` → 2 |
| `rem(a,b)` | Remainder (sign of a) | `rem(-17,5)` → -2 |
| `clamp(x,lo,hi)` / `clip` | Clamp to range | `clamp(5,0,3)` → 3 |
| `int(x)` | Integer part | `int(3.7)` → 3 |

### Number Theory (52)

| Function | Description | Example |
|----------|-------------|---------|
| `gcd(a,b)` | Greatest common divisor | `gcd(12,18)` → 6 |
| `lcm(a,b)` | Least common multiple | `lcm(4,6)` → 12 |
| `isprime(n)` | Primality test | `isprime(17)` → 1 |
| `nextprime(n)` | Next prime after n | `nextprime(10)` → 11 |
| `prevprime(n)` | Previous prime | `prevprime(10)` → 7 |
| `fact(n)` / `factorial(n)` | n! | `fact(5)` → 120 |
| `ncr(n,r)` / `comb` / `choose` | Combinations C(n,r) | `ncr(5,2)` → 10 |
| `npr(n,r)` / `perm` | Permutations P(n,r) | `npr(5,2)` → 20 |
| `fib(n)` / `fibonacci(n)` | Fibonacci number | `fib(10)` → 55 |
| `lucas(n)` | Lucas number | `lucas(10)` → 123 |
| `totient(n)` / `eulerphi(n)` | Euler's totient | `totient(12)` → 4 |
| `even(n)` | Is n even? | `even(4)` → 1 |
| `odd(n)` | Is n odd? | `odd(3)` → 1 |
| `catalan(n)` | Catalan number | `catalan(5)` → 42 |
| `bell(n)` | Bell number | `bell(5)` → 52 |
| `stirling2(n,k)` | Stirling number 2nd kind | `stirling2(5,3)` |
| `derangements(n)` / `subfactorial(n)` | Derangements | `derangements(4)` → 9 |
| `divisorsum(n)` / `sigma(n)` | Sum of divisors | `sigma(12)` → 28 |
| `divisors(n)` | List of divisors | `divisors(12)` |
| `factor(n)` | Prime factorization | `factor(12)` |
| `mobius(n)` / `mu(n)` | Möbius function | `mobius(6)` → 1 |
| `omega(n)` | Distinct prime factors | `omega(12)` → 2 |
| `bigomega(n)` / `Omega(n)` | Prime factors with mult | `bigomega(12)` → 3 |
| `ndigits(n)` / `numdigits(n)` | Number of digits | `ndigits(1234)` → 4 |
| `digitsum(n)` / `digsum(n)` | Sum of digits | `digitsum(1234)` → 10 |
| `digitroot(n)` / `digroot(n)` | Digital root | `digitroot(1234)` → 1 |
| `isperfect(n)` | Is perfect number? | `isperfect(28)` → 1 |
| `isabundant(n)` | Is abundant? | `isabundant(12)` → 1 |
| `isdeficient(n)` | Is deficient? | `isdeficient(8)` → 1 |
| `issquarefree(n)` | Is squarefree? | `issquarefree(15)` → 1 |
| `fact2(n)` / `factorial2(n)` | Double factorial n!! | `fact2(5)` → 15 |
| `pochhammer(x,n)` / `rising(x,n)` | Rising factorial | `pochhammer(3,4)` |
| `falling(x,n)` | Falling factorial | `falling(5,3)` |
| `primepi(n)` | Count primes ≤ n | `primepi(100)` → 25 |
| `nthprime(n)` | nth prime number | `nthprime(10)` → 29 |
| `partition(n)` / `partitions(n)` | Partition number | `partition(10)` → 42 |

### Special Functions (27)

| Function | Description | Example |
|----------|-------------|---------|
| `gamma(x)` / `tgamma(x)` | Gamma function | `gamma(5)` → 24 |
| `lgamma(x)` | Log-gamma function | `lgamma(5)` |
| `beta(a,b)` | Beta function | `beta(2,3)` → 0.0833 |
| `erf(x)` | Error function | `erf(1)` → 0.8427 |
| `erfc(x)` | Complementary error func | `erfc(0)` → 1 |
| `erfinv(x)` | Inverse error function | `erfinv(0.5)` |
| `besselj(n,x)` / `besselJ` | Bessel J | `besselj(0,1)` |
| `bessely(n,x)` / `besselY` | Bessel Y | `bessely(0,1)` |
| `zeta(s)` / `riemann_zeta` | Riemann zeta | `zeta(2)` → π²/6 |
| `sigmoid(x)` / `logistic(x)` | 1/(1+e^-x) | `sigmoid(0)` → 0.5 |
| `heaviside(x)` / `step(x)` | Step function | `heaviside(1)` → 1 |
| `digamma(x)` / `psi(x)` | Digamma function | `digamma(1)` |
| `harmonic(n)` | Harmonic number H_n | `harmonic(10)` |
| `genharmonic(n,m)` / `harmonic2` | Generalized harmonic | `genharmonic(10,2)` |
| `softplus(x)` | ln(1+e^x) | `softplus(0)` |
| `rect(x)` / `rectangle(x)` | Rectangle function | `rect(0)` → 1 |
| `tri(x)` / `triangle(x)` | Triangle function | `tri(0)` → 1 |

### Statistics (25)

| Function | Description | Example |
|----------|-------------|---------|
| `sum(v)` | Sum of elements | `sum([1,2,3])` → 6 |
| `prod(v)` | Product of elements | `prod([1,2,3])` → 6 |
| `mean(v)` | Arithmetic mean | `mean([1,2,3])` → 2 |
| `median(v)` | Median | `median([1,2,3])` → 2 |
| `mode(v)` | Mode | `mode([1,1,2])` → 1 |
| `std(v)` / `sd(v)` | Standard deviation | `std([1,2,3])` |
| `var(v)` | Variance | `var([1,2,3])` |
| `min(v)` | Minimum | `min([3,1,2])` → 1 |
| `max(v)` | Maximum | `max([1,3,2])` → 3 |
| `range(v)` | Range (max-min) | `range([1,5])` → 4 |
| `geomean(v)` | Geometric mean | `geomean([1,2,4])` → 2 |
| `harmmean(v)` | Harmonic mean | `harmmean([1,2,4])` |
| `skewness(v)` | Skewness | `skewness(v)` |
| `kurtosis(v)` | Kurtosis | `kurtosis(v)` |
| `mad(v)` | Median absolute deviation | `mad(v)` |
| `iqr(v)` | Interquartile range | `iqr(v)` |
| `rms(v)` | Root mean square | `rms([1,2,3])` |
| `sumsq(v)` | Sum of squares | `sumsq([1,2,3])` → 14 |
| `meansq(v)` | Mean of squares | `meansq([1,2,3])` |
| `prctile(v,p)` / `pct` | Percentile | `prctile(v,50)` |
| `zscore(v)` | Z-score normalization | `zscore(v)` |
| `cov(x,y)` | Covariance | `cov(x,y)` |
| `corrcoef(x,y)` | Correlation coefficient | `corrcoef(x,y)` |

### Forecasting / Time Series (12)

| Function | Description | Example |
|----------|-------------|---------|
| `movavg(x,k)` | Simple moving average | `movavg([1;2;3;4;5], 3)` |
| `ewma(x,alpha)` | Exponentially weighted MA | `ewma([1;2;3;4;5], 0.3)` |
| `trend(x)` | Extract linear trend | `trend([1;3;5;7;9])` |
| `detrend(x)` | Remove linear trend | `detrend([1;3;5;7;9])` |
| `lag(x,k)` | Lag by k periods | `lag([1;2;3;4;5])` |
| `lead(x,k)` | Lead by k periods | `lead([1;2;3;4;5])` |
| `pctchange(x)` | Percentage change | `pctchange([100;110;121])` |
| `autocorr(x)` | Autocorrelation function | `autocorr([1;2;3;4;5])` |
| `forecast(x,n)` | Linear forecast | `forecast([1;2;3;4;5], 3)` |
| `expsmooth(x,a)` | Exponential smoothing | `expsmooth([1;2;3], 0.3)` |
| `movmean(x,k)` | Centered moving mean | `movmean([1;2;3;4;5], 3)` |
| `movsum(x,k)` | Centered moving sum | `movsum([1;2;3;4;5], 3)` |

### Probability (15)

| Function | Description | Example |
|----------|-------------|---------|
| `normpdf(x)` / `normalpdf` | Normal PDF | `normpdf(0)` → 0.399 |
| `normcdf(x)` / `normalcdf` | Normal CDF | `normcdf(0)` → 0.5 |
| `norminv(p)` / `normalinv` / `probit` | Inverse normal CDF | `norminv(0.5)` → 0 |
| `rand()` | Uniform random [0,1) | `rand()` |
| `randn()` | Standard normal random | `randn()` |
| `randi(n)` | Random integer 1 to n | `randi(6)` |
| `chi2cdf(x,k)` / `chicdf` | Chi-squared CDF | `chi2cdf(3.84,1)` |
| `tcdf(x,df)` / `student_cdf` | Student's t CDF | `tcdf(2,10)` |
| `randperm(n)` | Random permutation | `randperm(5)` |

### Linear Algebra (24)

| Function | Description | Example |
|----------|-------------|---------|
| `det(A)` | Determinant | `det([1,2;3,4])` → -2 |
| `inv(A)` / `inverse(A)` | Inverse | `inv([1,2;3,4])` |
| `trace(A)` / `tr(A)` | Trace (diagonal sum) | `trace(eye(3))` → 3 |
| `transpose(A)` / `trans(A)` | Transpose | `transpose([1,2;3,4])` |
| `norm(v)` | Vector/matrix norm | `norm([3,4])` → 5 |
| `rank(A)` | Matrix rank | `rank([1,2;2,4])` → 1 |
| `eig(A)` / `eigenvalues(A)` | Eigenvalues | `eig([1,2;3,4])` |
| `svd(A)` | Singular values | `svd([1,2;3,4])` |
| `qr(A)` | QR decomposition | `qr([1,2;3,4])` |
| `lu(A)` | LU decomposition | `lu([1,2;3,4])` |
| `chol(A)` / `cholesky(A)` | Cholesky decomposition | `chol(A)` |
| `pinv(A)` | Pseudoinverse | `pinv(A)` |
| `linsolve(A,b)` | Solve Ax=b | `linsolve(A,b)` |
| `cond(A)` | Condition number | `cond(A)` |
| `null(A)` / `nullspace(A)` | Null space | `null(A)` |
| `schur(A)` | Schur decomposition | `schur(A)` |
| `mldivide(A,b)` / `ldivide` | Matrix left divide A\\b | `mldivide(A,b)` |

### Matrix Creation (35)

| Function | Description | Example |
|----------|-------------|---------|
| `zeros(m,n)` | Zero matrix | `zeros(2,3)` |
| `ones(m,n)` | Matrix of ones | `ones(2,3)` |
| `eye(n)` / `identity(n)` | Identity matrix | `eye(3)` |
| `diag(v)` | Diagonal matrix | `diag([1,2,3])` |
| `linspace(a,b,n)` | Linear spacing | `linspace(0,1,5)` |
| `logspace(a,b,n)` | Logarithmic spacing | `logspace(0,2,3)` |
| `reshape(A,m,n)` | Reshape matrix | `reshape(1:6,2,3)` |
| `sort(v)` | Sort elements | `sort([3,1,2])` |
| `unique(v)` | Unique elements | `unique([1,1,2])` |
| `flip(v)` | Flip vector | `flip([1,2,3])` |
| `fliplr(A)` | Flip left-right | `fliplr(A)` |
| `flipud(A)` | Flip up-down | `flipud(A)` |
| `rot90(A)` | Rotate 90° | `rot90(A)` |
| `dot(a,b)` | Dot product | `dot([1,2],[3,4])` → 11 |
| `cross(a,b)` | Cross product | `cross([1,0,0],[0,1,0])` |
| `kron(A,B)` | Kronecker product | `kron(A,B)` |
| `magic(n)` | Magic square | `magic(3)` |
| `pascal(n)` | Pascal matrix | `pascal(4)` |
| `hilb(n)` | Hilbert matrix | `hilb(3)` |
| `invhilb(n)` | Inverse Hilbert | `invhilb(3)` |
| `vander(v)` | Vandermonde matrix | `vander([1,2,3])` |
| `toeplitz(c,r)` | Toeplitz matrix | `toeplitz([1,2,3])` |
| `hankel(c)` | Hankel matrix | `hankel([1,2,3])` |
| `compan(p)` | Companion matrix | `compan([1,0,-1])` |
| `blkdiag(A,B,...)` | Block diagonal | `blkdiag(A,B)` |
| `cat(dim,A,B)` | Concatenate | `cat(1,A,B)` |
| `horzcat(A,B)` | Horizontal concat | `horzcat(A,B)` |
| `vertcat(A,B)` | Vertical concat | `vertcat(A,B)` |
| `circshift(A,k)` | Circular shift | `circshift([1,2,3],1)` |
| `repmat(A,m,n)` | Replicate matrix | `repmat([1,2],2,3)` |
| `triu(A)` | Upper triangular | `triu(A)` |
| `tril(A)` | Lower triangular | `tril(A)` |
| `sortrows(A)` | Sort rows | `sortrows(A)` |
| `find(v)` | Find nonzero indices | `find([0,1,2])` |

### Matrix Query (16)

| Function | Description | Example |
|----------|-------------|---------|
| `size(A)` | Matrix dimensions | `size([1,2;3,4])` → [2,2] |
| `length(v)` | Length of vector | `length([1,2,3])` → 3 |
| `rows(A)` | Number of rows | `rows([1,2;3,4])` → 2 |
| `cols(A)` | Number of columns | `cols([1,2;3,4])` → 2 |
| `numel(A)` | Number of elements | `numel([1,2;3,4])` → 4 |
| `isempty(A)` | Is empty? | `isempty([])` → 1 |
| `isscalar(x)` | Is scalar? | `isscalar(5)` → 1 |
| `isvector(v)` | Is vector? | `isvector([1,2])` → 1 |
| `issquare(A)` | Is square? | `issquare(eye(3))` → 1 |
| `nnz(A)` | Count nonzeros | `nnz([0,1,2])` → 2 |
| `ndims(A)` | Number of dimensions | `ndims(A)` |
| `isrow(v)` | Is row vector? | `isrow([1,2])` → 1 |
| `iscolumn(v)` | Is column vector? | `iscolumn([1;2])` → 1 |
| `ismatrix(A)` | Is matrix? | `ismatrix([1,2;3,4])` → 1 |
| `issorted(v)` | Is sorted? | `issorted([1,2,3])` → 1 |
| `issymmetric(A)` | Is symmetric? | `issymmetric(eye(3))` → 1 |

### Signal Processing (12)

| Function | Description | Example |
|----------|-------------|---------|
| `diff(v)` | Differences | `diff([1,3,6])` → [2,3] |
| `cumsum(v)` | Cumulative sum | `cumsum([1,2,3])` → [1,3,6] |
| `cumprod(v)` | Cumulative product | `cumprod([1,2,3])` → [1,2,6] |
| `conv(a,b)` | Convolution | `conv([1,2],[1,1])` |
| `deconv(a,b)` | Deconvolution | `deconv([1,3,2],[1,1])` |
| `trapz(v)` | Trapezoidal integration | `trapz([0,1,4])` → 2.5 |
| `gradient(v)` | Numerical gradient | `gradient([1,4,9])` |
| `cummin(v)` | Cumulative minimum | `cummin([3,1,2])` |
| `cummax(v)` | Cumulative maximum | `cummax([1,3,2])` |
| `movmean(v,k)` | Moving average | `movmean([1,2,3,4],2)` |
| `movsum(v,k)` | Moving sum | `movsum([1,2,3,4],2)` |
| `interp1(x,y,xi)` | 1D interpolation | `interp1([0,1],[0,1],0.5)` |

### Angle Conversion (8)

| Function | Description | Example |
|----------|-------------|---------|
| `deg2rad(x)` | Degrees to radians | `deg2rad(180)` → π |
| `rad2deg(x)` | Radians to degrees | `rad2deg(pi)` → 180 |
| `wrapToPi(x)` / `wrap(x)` | Wrap to [-π, π] | `wrapToPi(4)` |
| `wrapTo180(x)` / `wrap180(x)` | Wrap to [-180, 180] | `wrapTo180(270)` → -90 |
| `wrapTo360(x)` / `wrap360(x)` | Wrap to [0, 360) | `wrapTo360(370)` → 10 |

### Logical (28)

| Function | Description | Example |
|----------|-------------|---------|
| `any(v)` | Any nonzero? | `any([0,1,0])` → 1 |
| `all(v)` | All nonzero? | `all([1,1,0])` → 0 |
| `isnan(x)` | Is NaN? | `isnan(0/0)` → 1 |
| `isinf(x)` | Is infinite? | `isinf(1/0)` → 1 |
| `isfinite(x)` | Is finite? | `isfinite(1)` → 1 |
| `isreal(x)` | Is real? | `isreal(3+4i)` → 0 |
| `and(a,b)` / `band` | Logical AND | `and(1,0)` → 0 |
| `or(a,b)` / `bor` | Logical OR | `or(1,0)` → 1 |
| `not(x)` / `bnot` | Logical NOT | `not(0)` → 1 |
| `xor(a,b)` / `bxor` | Logical XOR | `xor(1,0)` → 1 |
| `nand(a,b)` | Logical NAND | `nand(1,1)` → 0 |
| `nor(a,b)` | Logical NOR | `nor(0,0)` → 1 |
| `implies(a,b)` / `imply` | Logical implication | `implies(1,0)` → 0 |
| `eq(a,b)` | Equal | `eq(1,1)` → 1 |
| `ne(a,b)` | Not equal | `ne(1,2)` → 1 |
| `lt(a,b)` | Less than | `lt(1,2)` → 1 |
| `le(a,b)` | Less or equal | `le(1,1)` → 1 |
| `gt(a,b)` | Greater than | `gt(2,1)` → 1 |
| `ge(a,b)` | Greater or equal | `ge(1,1)` → 1 |
| `approxeq(a,b)` / `isapprox` | Approximately equal | `approxeq(1,1.0001)` |
| `isequal(a,b)` | Exactly equal | `isequal(1,1)` → 1 |
| `isinteger(x)` | Is integer? | `isinteger(3)` → 1 |
| `isnumeric(x)` | Is numeric? | `isnumeric(3)` → 1 |
| `islogical(x)` | Is logical? | `islogical(true)` → 1 |
| `isnegative(x)` | Is negative? | `isnegative(-1)` → 1 |
| `ispositive(x)` | Is positive? | `ispositive(1)` → 1 |

### Bitwise (6)

| Function | Description | Example |
|----------|-------------|---------|
| `bitand(a,b)` | Bitwise AND | `bitand(12,10)` → 8 |
| `bitor(a,b)` | Bitwise OR | `bitor(12,10)` → 14 |
| `bitxor(a,b)` | Bitwise XOR | `bitxor(12,10)` → 6 |
| `bitshift(a,n)` | Bit shift | `bitshift(4,2)` → 16 |
| `shl(a,n)` | Shift left | `shl(1,3)` → 8 |
| `shr(a,n)` | Shift right | `shr(8,2)` → 2 |

### Polynomials (6)

| Function | Description | Example |
|----------|-------------|---------|
| `polyval(p,x)` | Evaluate polynomial | `polyval([1,0,-1],2)` → 3 |
| `polyder(p)` | Derivative of polynomial | `polyder([1,0,-1])` |
| `polyint(p)` | Integral of polynomial | `polyint([2,1])` |
| `roots2(a,b,c)` / `quadroots` / `quadratic` | Quadratic roots | `roots2(1,0,-1)` → [1,-1] |

### Set Operations (4)

| Function | Description | Example |
|----------|-------------|---------|
| `union(a,b)` | Set union | `union([1,2],[2,3])` |
| `intersect(a,b)` | Set intersection | `intersect([1,2],[2,3])` |
| `setdiff(a,b)` | Set difference | `setdiff([1,2,3],[2])` |
| `setxor(a,b)` | Symmetric difference | `setxor([1,2],[2,3])` |

### Coordinates (4)

| Function | Description | Example |
|----------|-------------|---------|
| `cart2pol(x,y)` | Cartesian to polar | `cart2pol(1,1)` |
| `pol2cart(r,θ)` | Polar to Cartesian | `pol2cart(1,pi/4)` |
| `cart2sph(x,y,z)` | Cartesian to spherical | `cart2sph(1,1,1)` |
| `sph2cart(r,θ,φ)` | Spherical to Cartesian | `sph2cart(1,0,0)` |

### Utility (5)

| Function | Description | Example |
|----------|-------------|---------|
| `nextpow2(n)` | Next power of 2 | `nextpow2(5)` → 3 |
| `normalize(v)` | Normalize vector | `normalize([3,4])` |
| `rescale(v)` | Rescale to [0,1] | `rescale([1,2,3])` |
| `disp(x)` | Display value | `disp(pi)` |
| `lerp(a,b,t)` | Linear interpolation | `lerp(0,10,0.5)` → 5 |

### Data Science (8)

| Function | Description | Example |
|----------|-------------|---------|
| `pca(X)` | Principal components | `pca(X)` |
| `pcareduce(X,k)` | PCA dimensionality reduction | `pcareduce(X,2)` |
| `kmeans(X,k)` | K-means clustering | `kmeans(X,3)` |
| `silhouette(X,idx)` | Silhouette scores | `silhouette(X,idx)` |
| `pdist(X)` / `dist` / `distance` | Pairwise distances | `pdist(X)` |
| `manhattan(a,b)` | Manhattan distance | `manhattan([0,0],[3,4])` → 7 |

### Machine Learning Classifiers

| Function | Description | Example |
|----------|-------------|---------|
| `fitcknn(X,Y,k)` | k-NN classifier | `mdl = fitcknn(X,Y,5)` |
| `fitcnb(X,Y)` | Naive Bayes | `mdl = fitcnb(X,Y)` |
| `fitctree(X,Y)` | Decision tree | `mdl = fitctree(X,Y)` |
| `fitcsvm(X,Y)` | SVM classifier | `mdl = fitcsvm(X,Y)` |
| `predict(mdl,X)` | Predict with model | `Ypred = predict(mdl,Xtest)` |
| `accuracy(Y,Ypred)` | Classification accuracy | `accuracy(Ytrue,Ypred)` |
| `confusionmat(Y,Ypred)` | Confusion matrix | `C = confusionmat(Y,Ypred)` |
| `cvpartition(Y,pct)` | Train/test split | `idx = cvpartition(Y,30)` |
| `crosstab(A,B)` | Cross-tabulation | `crosstab(pred,actual)` |

### Tables and Timetables

Create and manipulate tabular data with time-series support:

```matlab
% Create timetable with timestamp, customer ID, and usage data
times = [1704067200; 1704153600; 1704240000; 1704326400]   % Unix timestamps
custid = [1; 1; 2; 2]
usage = [100; 150; 200; 180]
timetable TT times custid usage

% View table contents
TT                        % Print full table
head TT 3                 % First 3 rows
tail TT 2                 % Last 2 rows
height(TT)                % Number of rows

% Group and aggregate
groupsummary MonthlySum TT times month sum usage       % Sum by month
groupsummary CustMean TT custid mean usage             % Mean by customer
groupcounts CustCounts TT custid                        % Count by customer

% Sort and manipulate
sortrows MonthlySum       % Sort by first column
```

| Command | Description |
|---------|-------------|
| `timetable Name cols...` | Create timetable from column vectors |
| `table Name cols...` | Create regular table |
| `head TT [n]` | Show first n rows (default 5) |
| `tail TT [n]` | Show last n rows (default 5) |
| `height(TT)` | Number of rows |
| `groupsummary Res T col method data` | Group and aggregate |
| `groupsummary Res T col bin method data` | Group with time binning |
| `groupcounts Res T col` | Count by group |
| `sortrows TT` | Sort by first column |

Time binning options for groupsummary: `year`, `month`, `day`, `hour`, `minute`, `week`, `quarter`

Aggregation methods: `sum`, `mean`, `min`, `max`, `numel`

---

## Commands

| Command | Description |
|---------|-------------|
| `help` | Show help |
| `help <func>` | Help for specific function |
| `demo <func>` | Interactive demo |
| `functions` | List all 640+ functions |
| `constants` | Show mathematical constants |
| `features` | Show enabled features |
| `test` | Run internal test suite |
| `bench` | Run benchmarks |
| `digits N` | Set display precision (1-37) |
| `format short` / `format long` | MATLAB-style format (4 or 15 digits) |
| `fix N` / `sci N` / `eng N` | Display format |
| `deg` / `rad` / `grad` | Angle mode |
| `rpn` / `alg` | RPN or algebraic mode |
| `vars` | Show variables |
| `clear` | Clear variables |
| `quit` | Exit |

## Matrix Operations

```matlab
>>> A = [1,2;3,4]        % Create 2x2 matrix
>>> B = eye(2)           % Identity matrix
>>> A * B                % Matrix multiply
>>> A .* B               % Element-wise multiply
>>> A(1,:)               % First row
>>> A(:,2)               % Second column
>>> A(1:2,1:2)           % Submatrix
>>> det(A)               % Determinant = -2
>>> inv(A)               % Inverse
>>> eig(A)               % Eigenvalues
>>> [Q,R] = qr(A)        % QR decomposition
>>> svd(A)               % Singular values
```

## Machine Learning Example

```matlab
>>> load fisheriris                    % Load iris dataset
>>> mdl = fitcknn(meas, species, 5)    % Train k-NN (k=5)
>>> Ypred = predict(mdl, meas)         % Predict
>>> accuracy(species, Ypred)           % ~96% accuracy
>>> confusionmat(species, Ypred)       % Confusion matrix
```

## Built-in Datasets

| Dataset | Command | Variables | Description |
|---------|---------|-----------|-------------|
| Fisher's Iris | `load fisheriris` | `meas` (150×4), `species` (150×1) | Flower classification |
| Hald Cement | `load hald` | `X` (13×4), `y` (13×1) | Regression benchmark |
| Imports-85 | `load imports85` | `auto` (98×15), `price` (98×1) | Automobile prices |

```matlab
>>> load hald           % Load cement dataset
>>> mean(y)             % Average heat: 95.42
>>> corrcoef(X(:,1), y) % Correlation with calcium aluminate

>>> load imports85      % Load automobile dataset  
>>> mean(price)         % Average price: $12,677
>>> max(price)          % Most expensive: $36,880
```

## Building

```bash
make                # Full build (Linux/macOS)
make test           # Run 627 tests
./bin/sc            # Run calculator

# Cross-compile for DOS
make CFLAGS="-DPLATFORM_DOS"
```

## License

MIT License - see [LICENSE](LICENSE)

## New in This Release

### Forecasting & Time Series (50+ functions)
- **Moving averages**: `movavg(x, k)`, `movmean(x, k)`, `movsum(x, k)`
- **Exponential smoothing**: `ewma(x, alpha)`, `expsmooth(x, alpha)`
- **Trend analysis**: `trend(x)`, `detrend(x)`
- **Time series operations**: `lag(x, k)`, `lead(x, k)`, `pctchange(x)`
- **Correlation**: `autocorr(x, maxlag)`, `crosscorr(x, y, maxlag)`

### Date/Time Functions
- **Date creation**: `datenum(y,m,d,h,mi,s)`, `datetime({...})`
- **Date parts**: `year(t)`, `month(t)`, `day(t)`, `hour(t)`, `minute(t)`, `second(t)`
- **Duration conversion**: `days(n)`, `hours(n)`, `minutes(n)`, `seconds(n)`
- **Period boundaries**: `startofmonth(t)`, `startofyear(t)`, `startofday(t)` (aliases: `som`, `soy`, `sod`)
- **Time binning**: `dateshift(t, 'start', unit)` - shift to period start

### Table Operations
- **Creation**: `timetable name col1 col2 ...`, `table name col1 col2 ...`
- **Aggregation**: `groupsummary Result T groupcol method datacol` with time binning support
- **Counting**: `groupcounts Result T groupcol`
- **Inspection**: `head T n`, `tail T n`, `height(T)`, `width(T)`
- **Sorting**: `sortrows tablename`

### Time Binning in groupsummary
Group data by time periods in a single step:
```matlab
% Instead of:
months = dateshift(TT.times, 'start', 'month')
timetable MM months custid usage
groupsummary MonthlyUsage MM months sum usage

% Now do this:
groupsummary MonthlyUsage TT times month sum usage
```
Supported bins: `year`, `month`, `week`, `day`, `hour`, `minute`, `dayname`

### Number Theory Additions
- **Figurate numbers**: `triangular(n)`, `pentagonal(n)`, `hexagonal(n)`
- **Tests**: `istriangular(n)`, `isperfectsquare(n)`, `ispower(n)`
- **Digit operations**: `reversedigits(n)`, `ispalindrome(n)`
- **Combinatorics**: `bell(n)`, `stirling2(n,k)`, `partition(n)`, `derangements(n)`

### Bug Fixes
- **Critical**: Fixed exponential notation parsing hang for large exponents like `1e1000000`
- Fixed case-insensitive help lookup (e.g., `help SIN` now works)
- Added precision-specific constants: `epsp` (native epsilon), `precision` (bit width)

### Help System
- 640+ functions with comprehensive help documentation
- All functions have working `help function` and `demo function`
- Run `help categories` to see function categories

### Data Cleaning Functions
- **Missing data**: `fillmissing(x, method)` - linear interpolation, mean, or previous value
- **Remove NaN**: `rmmissing(x)` - remove all NaN values from vector
- **Smoothing**: `smoothdata(x, method, window)` - moving average or Gaussian smoothing

### Interactive Visualization Commands

#### `iplot` - Interactive Function Plotter
Full-screen function plotting with pan and zoom.

```matlab
>>> iplot sin(x)           % Plot with default range
>>> iplot x^2 -5 5         % Plot with custom x range
>>> iplot sin(x)*cos(2*x)  % Complex expressions
```

**Controls:**
- **Arrow keys**: Pan left/right (X) and up/down (Y)
- **+/-**: Zoom in/out
- **a**: Toggle auto-Y scaling
- **r**: Reset to original view
- **q/ESC**: Exit

#### `mandelbrot` - Mandelbrot Set Viewer
Interactive fractal explorer with 256-color rendering.

```matlab
>>> mandelbrot    % Launch viewer
>>> mandel        % Alias
```

**Controls:**
- **Arrow keys**: Pan around
- **+/-**: Zoom in/out (auto-adjusts iterations)
- **i/I**: Decrease/increase max iterations
- **c**: Toggle color/ASCII mode
- **w/e/s**: Jump to cardioid/elephant valley/seahorse valley
- **r**: Reset view
- **q/ESC**: Exit

#### `tscatter` / `planets` / `solar` - Solar System Viewer
Interactive time-based solar system simulation.

```matlab
>>> tscatter           % All planets (Mercury through Saturn)
>>> tscatter inner     % Inner planets (Mercury, Venus, Earth, Mars)
>>> tscatter outer     % Outer planets (Jupiter, Saturn, Uranus, Neptune)
>>> tscatter moon      % Earth-Moon system (follows Earth+Moon barycenter)
>>> planets all        % Alias for tscatter
>>> solar              % Alias for tscatter
```

**Controls:**
- **←/→**: Step time backward/forward (×10 of dt)
- **↑/↓**: Increase/decrease time step
- **Space**: Play/pause animation
- **+/-**: Zoom in/out
- **f/F**: Cycle follow target (sun → earth → earth+moon → each planet)
- **1**: Inner view (2 AU, follow sun)
- **2**: Outer view (35 AU, follow sun)
- **3**: Moon view (0.01 AU, follow earth+moon barycenter)
- **4**: Follow Earth
- **r**: Reset to t=0
- **q/ESC**: Exit

### Astronomy Functions

#### `planet("name", t_days)` - Planet Positions
Returns planet position as complex number (x + iy) in AU (Astronomical Units).

```matlab
>>> planet("earth", 0)           % Earth at t=0
= 0.98 + 0.199i

>>> planet("mars", 365)          % Mars after 1 year
= -1.49 - 0.32i

>>> abs(planet("jupiter", 0))    % Jupiter's distance from sun
= 5.2

>>> real(planet("venus", 100))   % Venus x-coordinate
= 0.203
```

**Available bodies:** sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, moon

#### `mandelbrot_iter(x, y, [max_iter])` - Mandelbrot Iteration Count
Returns the number of iterations before escape (or max_iter if in set).

```matlab
>>> mandelbrot_iter(0, 0, 100)      % Origin (in set) = 100
>>> mandelbrot_iter(-2.5, 0, 100)   % Outside = 1
>>> mandelbrot_iter(0.5, 0.5, 100)  % Escapes after 5 iterations
```
