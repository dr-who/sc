# New MATLAB/Octave/Mathematica Compatible Features (Session 2)

## Matrix Analysis (4 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `rank(M)` | Matrix rank | `rank([1 2; 2 4])` → 1 |
| `cond(M)` | Condition number | `cond([1 2; 3 4])` → 14.9 |
| `linsolve(A,b)` | Solve Ax = b | `linsolve([1 2;3 4],[5;11])` |
| `gradient(v)` | Numerical gradient | `gradient([1 4 9 16 25])` |

## Element-wise Comparisons (6 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `eq(A,B)` | Element-wise equal | `eq([1 2],[1 3])` → [1 0] |
| `ne(A,B)` | Element-wise not equal | `ne([1 2],[1 3])` → [0 1] |
| `lt(A,B)` | Element-wise less than | `lt([1 2],[2 2])` → [1 0] |
| `le(A,B)` | Element-wise ≤ | `le([1 2],[2 2])` → [1 1] |
| `gt(A,B)` | Element-wise greater | `gt([1 2],[0 1])` → [1 1] |
| `ge(A,B)` | Element-wise ≥ | `ge([1 2],[1 3])` → [1 0] |

## Element-wise Operations (5 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `times(A,B)` | Element-wise multiply | `times([1 2],[3 4])` → [3 8] |
| `rdivide(A,B)` | Element-wise right divide | `rdivide([4 6],[2 3])` → [2 2] |
| `ldivide(A,B)` | Element-wise left divide | `ldivide([2 3],[4 6])` → [2 2] |
| `plus(A,B)` | Element-wise add | `plus([1 2],[3 4])` → [4 6] |
| `minus(A,B)` | Element-wise subtract | `minus([3 4],[1 2])` → [2 2] |

## Cumulative Operations (3 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `cummax(v)` | Cumulative maximum | `cummax([3 1 4])` → [3 3 4] |
| `cummin(v)` | Cumulative minimum | `cummin([3 1 4])` → [3 1 1] |
| `rescale(v)` | Scale to [0,1] | `rescale([1 2 3])` → [0 0.5 1] |

## Statistics (3 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `corrcoef(X,Y)` | Correlation coefficient | `corrcoef([1 2 3],[2 4 6])` → 1 |
| `movmean(v,k)` | Moving average | `movmean([1:7], 3)` |
| `movsum(v,k)` | Moving sum | `movsum([1:7], 3)` |

## Signal Processing (2 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `conv(a,b)` | Convolution | `conv([1 2],[1 1])` → [1 3 2] |
| `trapz(y)` | Trapezoidal integration | `trapz([1 2 3])` → 4 |

## Special Matrices (4 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `kron(A,B)` | Kronecker product | `kron([1 2],[1;2])` |
| `vander(v,n)` | Vandermonde matrix | `vander([1 2 3])` |
| `toeplitz(c)` | Toeplitz matrix | `toeplitz([1 2 3])` |
| `hankel(c)` | Hankel matrix | `hankel([1 2 3])` |

## Total: 27 New Functions

Combined with previous session:
- Session 1: 40 functions
- Session 2: 27 functions
- **Total new: 67 MATLAB-compatible functions**

## Examples

```matlab
% Element-wise comparisons
eq([1 2 3], [1 2 4])       % [1 1 0]
lt([1 2 3], [2 2 2])       % [1 0 0]

% Cumulative operations
cummax([3 1 4 1 5])        % [3 3 4 4 5]
rescale([0 50 100])        % [0 0.5 1]

% Statistics
corrcoef([1 2 3], [3 2 1]) % -1 (negative correlation)
movmean([1:10], 3)         % 3-element moving average

% Signal processing
conv([1 2 3], [1 1])       % [1 3 5 3]
trapz([0 1 4 9])           % 10.5 (trapezoidal integral)

% Special matrices
kron([1 2], [1; 2])        % [1 2; 2 4]
vander([1 2 3], 3)         % Vandermonde matrix
toeplitz([1 2 3])          % Symmetric Toeplitz
hankel([1 2 3 4])          % Hankel matrix

% Linear algebra
linsolve([1 2; 3 4], [5; 11])  % [1; 2]
```

## Additional Functions (Session 2 Part 2)

### Index Conversion (2 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `sub2ind([m,n],r,c)` | Subscripts to linear index | `sub2ind([3,4],2,3)` → 8 |
| `ind2sub([m,n],idx)` | Linear index to subscripts | `ind2sub([3,4],8)` → [2 3] |

### Matrix Construction (3 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `blkdiag(A,B)` | Block diagonal | `blkdiag([1 2],[3 4])` |
| `cat(dim,A,B)` | Concatenate along dimension | `cat(1,[1],[2])` |
| `compan(p)` | Companion matrix | `compan([1 -3 2])` |

### Utility Functions (4 functions)

| Function | Description | Example |
|----------|-------------|---------|
| `clamp(x,lo,hi)` | Limit to range | `clamp(15, 0, 10)` → 10 |
| `isapprox(a,b)` | Approximate equality | `isapprox(0.1+0.2, 0.3)` → 1 |
| `colon(a,b)` | Explicit colon | `colon(1,5)` → [1 2 3 4 5] |
| `linspace2(a,b,n)` | Explicit linspace | `linspace2(0,1,5)` |

## Summary

**Session 2 Total: 37 new functions**

Categories:
- Matrix Analysis: rank, cond, linsolve, gradient (4)
- Element-wise Comparisons: eq, ne, lt, le, gt, ge (6)
- Element-wise Operations: times, rdivide, ldivide, plus, minus (5)
- Cumulative: cummax, cummin, rescale (3)
- Statistics: corrcoef, movmean, movsum (3)
- Signal Processing: conv, trapz (2)
- Special Matrices: kron, vander, toeplitz, hankel (4)
- Index Conversion: sub2ind, ind2sub (2)
- Matrix Construction: blkdiag, cat, compan (3)
- Utility: clamp, isapprox, colon, linspace2 (4)

**Combined Total: 284 functions**
