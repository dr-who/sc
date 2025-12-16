/*
 * func_misc.c - Signal processing, polynomials, angles, sets, logic, bitwise,
 *               coordinates, utility, constants, data science, text
 * C89 compliant
 */
#include "func_defs.h"

/* Signal processing functions */
const FuncDoc func_signal_docs[] = {
    {"conv", CAT_SIGNAL, "conv(a, b)",
     "Convolution of vectors.",
     {"conv([1,2,3], [1,1])", "conv([1,0,1], [1,2,1])", NULL, NULL, NULL, NULL},
     "deconv"},
    
    {"deconv", CAT_SIGNAL, "deconv(a, b)",
     "Deconvolution. Polynomial division.",
     {"deconv([1,3,3,1], [1,1])", "deconv([1,4,4,1], [1,2])", "deconv([1,2,1], [1,1])", NULL, NULL, NULL},
     "conv"},
    
    {"diff", CAT_SIGNAL, "diff(v)",
     "Differences between adjacent elements. Length reduces by 1.",
     {"diff([1,4,9,16]) = [3,5,7]", "diff([1,2,4,7,11])", "diff([1;2;3])", NULL, NULL, NULL},
     "cumsum, gradient"},
    
    {"gradient", CAT_SIGNAL, "gradient(v)",
     "Numerical gradient (central differences).",
     {"gradient([1,2,4,7,11])", "gradient([0,1,4,9,16])", NULL, NULL, NULL, NULL},
     "diff"},
    
    {"cumsum", CAT_SIGNAL, "cumsum(v)",
     "Cumulative sum.",
     {"cumsum([1,2,3,4]) = [1,3,6,10]", "cumsum(1:5)", "cumsum([1;2;3;4])", NULL, NULL, NULL},
     "cumprod, sum, diff"},
    
    {"cumprod", CAT_SIGNAL, "cumprod(v)",
     "Cumulative product.",
     {"cumprod([1,2,3,4]) = [1,2,6,24]", "cumprod(1:5)", NULL, NULL, NULL, NULL},
     "cumsum, prod"},
    
    {"cummin", CAT_SIGNAL, "cummin(v)",
     "Cumulative minimum.",
     {"cummin([3,1,4,1,5]) = [3,1,1,1,1]", "cummin([5,4,3,2,1])", NULL, NULL, NULL, NULL},
     "cummax, min"},
    
    {"cummax", CAT_SIGNAL, "cummax(v)",
     "Cumulative maximum.",
     {"cummax([3,1,4,1,5]) = [3,3,4,4,5]", "cummax([1,2,3,4,5])", NULL, NULL, NULL, NULL},
     "cummin, max"},
    
    {"movmean", CAT_SIGNAL, "movmean(v, k)",
     "Moving average with window size k.",
     {"movmean([1,2,3,4,5], 3)", "movmean(1:10, 5)", NULL, NULL, NULL, NULL},
     "movsum, mean"},
    
    {"movsum", CAT_SIGNAL, "movsum(v, k)",
     "Moving sum with window size k.",
     {"movsum([1,2,3,4,5], 3)", "movsum(1:10, 3)", NULL, NULL, NULL, NULL},
     "movmean, sum"},
    
    {"interp1", CAT_SIGNAL, "interp1(x, y, xi)",
     "1D linear interpolation.",
     {"interp1([0,1,2], [0,1,4], 0.5)", "interp1([0,1,2], [0,1,4], 1.5)", NULL, NULL, NULL, NULL},
     "linspace"},
    
    {"trapz", CAT_SIGNAL, "trapz(y) or trapz(x, y)",
     "Trapezoidal numerical integration.",
     {"trapz([1,2,3,4]) = 7.5", "trapz([0,1,2,3], [0,1,4,9])", NULL, NULL, NULL, NULL},
     "cumsum"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Polynomial functions */
const FuncDoc func_poly_docs[] = {
    {"polyval", CAT_POLY, "polyval(p, x)",
     "Evaluate polynomial p at x. p = [a_n, ..., a_1, a_0].",
     {"polyval([1,0,-2], 2) = 2", "polyval([1,-6,11,-6], 1) = 0", "polyval([1,2,1], -1) = 0", NULL, NULL, NULL},
     "roots2, polyder"},
    
    {"polyder", CAT_POLY, "polyder(p)",
     "Polynomial derivative.",
     {"polyder([1,0,-2]) = [2,0]", "polyder([1,2,3,4])", "polyder([3,2,1])", NULL, NULL, NULL},
     "polyint, polyval"},
    
    {"polyint", CAT_POLY, "polyint(p) or polyint(p, c)",
     "Polynomial integral. c is integration constant (default 0).",
     {"polyint([2,0]) = [1,0,0]", "polyint([3,2,1])", "polyint([1,1], 5)", NULL, NULL, NULL},
     "polyder"},
    
    {"roots2", CAT_POLY, "roots2(a, b, c)",
     "Roots of quadratic ax^2 + bx + c = 0. Returns [x1, x2].",
     {"roots2(1, 0, -4) = [2,-2]", "roots2(1, -5, 6) = [3,2]", "roots2(1, 0, 1) = [i,-i]", "roots2(1, -2, 1) = [1,1]", NULL, NULL},
     "polyval"},
    
    {"quadroots", CAT_POLY, "quadroots(a, b, c)",
     "Alias for roots2.",
     {"quadroots(1, -3, 2)", "quadroots(1, 0, -9)", NULL, NULL, NULL, NULL},
     "roots2"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Angle conversion functions */
const FuncDoc func_angle_docs[] = {
    {"deg2rad", CAT_ANGLE, "deg2rad(x)",
     "Convert degrees to radians. x * pi / 180.",
     {"deg2rad(180) = 3.1416", "deg2rad(90) = 1.5708", "deg2rad(45) = 0.7854", "deg2rad(360) = 6.2832", NULL, NULL},
     "rad2deg"},
    
    {"rad2deg", CAT_ANGLE, "rad2deg(x)",
     "Convert radians to degrees. x * 180 / pi.",
     {"rad2deg(pi) = 180", "rad2deg(pi/2) = 90", "rad2deg(1) = 57.2958", "rad2deg(2*pi) = 360", NULL, NULL},
     "deg2rad"},
    
    {"wrapToPi", CAT_ANGLE, "wrapToPi(x)",
     "Wrap angle to [-pi, pi].",
     {"wrapToPi(4) = 0.8584", "wrapToPi(-4) = -0.8584", "wrapToPi(2*pi) = 0", NULL, NULL, NULL},
     "wrap180, wrapTo180"},
    
    {"wrap180", CAT_ANGLE, "wrap180(x)",
     "Wrap angle to [-180, 180] degrees.",
     {"wrap180(270) = -90", "wrap180(-270) = 90", "wrap180(360) = 0", "wrap180(450) = 90", NULL, NULL},
     "wrap360, wrapToPi"},
    
    {"wrap360", CAT_ANGLE, "wrap360(x)",
     "Wrap angle to [0, 360) degrees.",
     {"wrap360(-90) = 270", "wrap360(370) = 10", "wrap360(720) = 0", "wrap360(-10) = 350", NULL, NULL},
     "wrap180"},
    
    {"wrapTo180", CAT_ANGLE, "wrapTo180(x)",
     "Alias for wrap180.",
     {"wrapTo180(270) = -90", "wrapTo180(-270) = 90", NULL, NULL, NULL, NULL},
     "wrap180"},
    
    {"wrapTo360", CAT_ANGLE, "wrapTo360(x)",
     "Alias for wrap360.",
     {"wrapTo360(-90) = 270", "wrapTo360(370) = 10", NULL, NULL, NULL, NULL},
     "wrap360"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Set operations */
const FuncDoc func_set_docs[] = {
    {"union", CAT_SET, "union(a, b)",
     "Set union. Unique sorted elements in either a or b.",
     {"union([1,2,3], [2,3,4]) = [1,2,3,4]", "union([1,1,2], [2,3])", NULL, NULL, NULL, NULL},
     "intersect, setdiff"},
    
    {"intersect", CAT_SET, "intersect(a, b)",
     "Set intersection. Unique sorted elements in both a and b.",
     {"intersect([1,2,3], [2,3,4]) = [2,3]", "intersect([1,2], [3,4]) = []", NULL, NULL, NULL, NULL},
     "union, setdiff"},
    
    {"setdiff", CAT_SET, "setdiff(a, b)",
     "Set difference. Elements in a but not in b.",
     {"setdiff([1,2,3,4], [2,4]) = [1,3]", "setdiff([1,2,3], [1,2,3]) = []", NULL, NULL, NULL, NULL},
     "union, intersect, setxor"},
    
    {"setxor", CAT_SET, "setxor(a, b)",
     "Set exclusive or. Elements in a or b but not both.",
     {"setxor([1,2,3], [2,3,4]) = [1,4]", "setxor([1,2], [2,3])", NULL, NULL, NULL, NULL},
     "union, intersect, setdiff"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Logical functions */
const FuncDoc func_logic_docs[] = {
    {"any", CAT_LOGIC, "any(v)",
     "True if any element is nonzero.",
     {"any([0,0,1,0]) = 1", "any([0,0,0]) = 0", "any([1,2,3]) = 1", NULL, NULL, NULL},
     "all"},
    
    {"all", CAT_LOGIC, "all(v)",
     "True if all elements are nonzero.",
     {"all([1,2,3,4]) = 1", "all([1,0,1]) = 0", "all([]) = 1", NULL, NULL, NULL},
     "any"},
    
    {"eq", CAT_LOGIC, "eq(a, b) or a == b",
     "Element-wise equality.",
     {"eq(3, 3) = 1", "eq([1,2,3], [1,2,4])", "3 == 3 = true", NULL, NULL, NULL},
     "ne"},
    
    {"ne", CAT_LOGIC, "ne(a, b) or a ~= b",
     "Element-wise inequality.",
     {"ne(3, 4) = 1", "ne([1,2,3], [1,2,3])", "3 ~= 4 = true", NULL, NULL, NULL},
     "eq"},
    
    {"lt", CAT_LOGIC, "lt(a, b) or a < b",
     "Element-wise less than.",
     {"lt(3, 4) = 1", "3 < 4 = true", "lt([1,2,3], 2)", NULL, NULL, NULL},
     "le, gt, ge"},
    
    {"le", CAT_LOGIC, "le(a, b) or a <= b",
     "Element-wise less than or equal.",
     {"le(3, 3) = 1", "3 <= 4 = true", "le([1,2,3], 2)", NULL, NULL, NULL},
     "lt, gt, ge"},
    
    {"gt", CAT_LOGIC, "gt(a, b) or a > b",
     "Element-wise greater than.",
     {"gt(4, 3) = 1", "4 > 3 = true", "gt([1,2,3], 2)", NULL, NULL, NULL},
     "ge, lt, le"},
    
    {"ge", CAT_LOGIC, "ge(a, b) or a >= b",
     "Element-wise greater than or equal.",
     {"ge(3, 3) = 1", "4 >= 3 = true", "ge([1,2,3], 2)", NULL, NULL, NULL},
     "gt, lt, le"},
    
    {"isnan", CAT_LOGIC, "isnan(x)",
     "Test for NaN (Not a Number).",
     {"isnan(0/0) = 1", "isnan(1/0) = 0", "isnan(5) = 0", "isnan(NaN) = 1", NULL, NULL},
     "isinf, isfinite"},
    
    {"isinf", CAT_LOGIC, "isinf(x)",
     "Test for infinity.",
     {"isinf(1/0) = 1", "isinf(-1/0) = 1", "isinf(5) = 0", "isinf(Inf) = 1", NULL, NULL},
     "isnan, isfinite"},
    
    {"isfinite", CAT_LOGIC, "isfinite(x)",
     "Test for finite value (not NaN or Inf).",
     {"isfinite(5) = 1", "isfinite(Inf) = 0", "isfinite(NaN) = 0", "isfinite(0) = 1", NULL, NULL},
     "isnan, isinf"},
    
    {"isreal", CAT_LOGIC, "isreal(x)",
     "Test if imaginary part is zero.",
     {"isreal(5) = 1", "isreal(3+4i) = 0", "isreal(2i) = 0", "isreal(0) = 1", NULL, NULL},
     "im"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Bitwise functions */
const FuncDoc func_bitwise_docs[] = {
    {"bitand", CAT_BITWISE, "bitand(a, b)",
     "Bitwise AND.",
     {"bitand(12, 10) = 8", "bitand(255, 15) = 15", "bitand(7, 3) = 3", NULL, NULL, NULL},
     "bitor, bitxor"},
    
    {"bitor", CAT_BITWISE, "bitor(a, b)",
     "Bitwise OR.",
     {"bitor(12, 10) = 14", "bitor(8, 4) = 12", "bitor(5, 3) = 7", NULL, NULL, NULL},
     "bitand, bitxor"},
    
    {"bitxor", CAT_BITWISE, "bitxor(a, b)",
     "Bitwise XOR.",
     {"bitxor(12, 10) = 6", "bitxor(255, 15) = 240", "bitxor(5, 3) = 6", NULL, NULL, NULL},
     "bitand, bitor"},
    
    {"bnot", CAT_BITWISE, "bnot(a)",
     "Bitwise NOT (one's complement).",
     {"bnot(0) = -1", "bnot(255) = -256", "bnot(-1) = 0", NULL, NULL, NULL},
     "bitand, bitor"},
    
    {"bitshift", CAT_BITWISE, "bitshift(a, k)",
     "Bit shift. Positive k shifts left, negative shifts right.",
     {"bitshift(1, 4) = 16", "bitshift(16, -2) = 4", "bitshift(7, 1) = 14", NULL, NULL, NULL},
     "shl, shr"},
    
    {"shl", CAT_BITWISE, "shl(a, k)",
     "Shift left k bits. Same as bitshift(a, k).",
     {"shl(1, 8) = 256", "shl(5, 2) = 20", "shl(1, 10) = 1024", NULL, NULL, NULL},
     "shr, bitshift"},
    
    {"shr", CAT_BITWISE, "shr(a, k)",
     "Shift right k bits. Same as bitshift(a, -k).",
     {"shr(256, 8) = 1", "shr(20, 2) = 5", "shr(1024, 5) = 32", NULL, NULL, NULL},
     "shl, bitshift"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Coordinate transforms */
const FuncDoc func_coord_docs[] = {
    {"cart2sph", CAT_COORD, "cart2sph(x, y, z)",
     "Cartesian to spherical coordinates. Returns [azimuth, elevation, r].",
     {"cart2sph(1, 0, 0)", "cart2sph(0, 0, 1)", "cart2sph(1, 1, 1)", NULL, NULL, NULL},
     "sph2cart"},
    
    {"sph2cart", CAT_COORD, "sph2cart(az, el, r)",
     "Spherical to Cartesian coordinates. Returns [x, y, z].",
     {"sph2cart(0, 0, 1)", "sph2cart(pi/2, 0, 1)", "sph2cart(0, pi/2, 1)", NULL, NULL, NULL},
     "cart2sph"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Utility functions */
const FuncDoc func_util_docs[] = {
    {"eps", CAT_UTIL, "eps",
     "Machine epsilon. Smallest x where 1 + x > 1.",
     {"eps", "1 + eps > 1 = true", "1 + eps/2 > 1 = false", NULL, NULL, NULL},
     "realmax, realmin"},
    
    {"realmax", CAT_UTIL, "realmax",
     "Largest finite floating-point number.",
     {"realmax", "realmax > 1e300 = true", NULL, NULL, NULL, NULL},
     "realmin, eps"},
    
    {"realmin", CAT_UTIL, "realmin",
     "Smallest positive normalized floating-point number.",
     {"realmin", "realmin > 0 = true", "realmin < 1e-300 = true", NULL, NULL, NULL},
     "realmax, eps"},
    
    {"nextpow2", CAT_UTIL, "nextpow2(n)",
     "Exponent of next power of 2. Smallest k where 2^k >= |n|.",
     {"nextpow2(5) = 3", "nextpow2(9) = 4", "nextpow2(1) = 0", "nextpow2(1000) = 10", NULL, NULL},
     "log2"},
    
    {"normalize", CAT_UTIL, "normalize(v)",
     "Normalize vector to unit length.",
     {"normalize([3,4]) = [0.6,0.8]", "normalize([1,0,0])", "norm(normalize([3,4])) = 1", NULL, NULL, NULL},
     "norm"},
    
    {"rescale", CAT_UTIL, "rescale(v)",
     "Rescale to [0, 1] range.",
     {"rescale([1,2,3,4,5])", "rescale([0,10,100])", NULL, NULL, NULL, NULL},
     "normalize"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Constants */
const FuncDoc func_const_docs[] = {
    {"pi", CAT_CONST, "pi",
     "Circle constant. Ratio of circumference to diameter. 3.14159265...",
     {"pi = 3.1416", "2*pi = 6.2832", "sin(pi) = 0", "cos(pi) = -1", NULL, NULL},
     "e, tau"},
    
    {"e", CAT_CONST, "e",
     "Euler's number. Base of natural logarithm. 2.71828182...",
     {"e = 2.7183", "ln(e) = 1", "exp(1) = e", "e^2 = 7.3891", NULL, NULL},
     "pi, exp"},
    
    {"i", CAT_CONST, "i",
     "Imaginary unit. sqrt(-1).",
     {"i^2 = -1", "i*i = -1", "sqrt(-1) = i", "e^(i*pi) = -1", NULL, NULL},
     "complex, re, im"},
    
    {"phi", CAT_CONST, "phi",
     "Golden ratio. (1 + sqrt(5)) / 2 = 1.61803398...",
     {"phi = 1.618", "phi^2 - phi = 1", "1/phi = phi - 1", NULL, NULL, NULL},
     "golden"},
    
    {"golden", CAT_CONST, "golden",
     "Golden ratio. Alias for phi.",
     {"golden = 1.618", "golden = phi", NULL, NULL, NULL, NULL},
     "phi"},
    
    {"tau", CAT_CONST, "tau",
     "Tau = 2*pi. Full circle in radians.",
     {"tau = 6.2832", "tau/2 = pi", "sin(tau) = 0", NULL, NULL, NULL},
     "pi"},
    
    {"Inf", CAT_CONST, "Inf",
     "Positive infinity.",
     {"1/0 = Inf", "Inf + 1 = Inf", "Inf > 1e308 = true", "-Inf < 0 = true", NULL, NULL},
     "NaN, isinf"},
    
    {"NaN", CAT_CONST, "NaN",
     "Not a Number. Result of undefined operations.",
     {"0/0 = NaN", "Inf - Inf = NaN", "NaN == NaN = false", "isnan(NaN) = 1", NULL, NULL},
     "Inf, isnan"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Data science functions */
const FuncDoc func_datasci_docs[] = {
    {"pca", CAT_DATASCI, "[coeff,score,latent] = pca(X)",
     "Principal component analysis. coeff=loadings, score=transformed data, latent=eigenvalues.",
     {"pca([1,2;3,4;5,6;7,8])", "pca([1,0;0,1;1,1])", NULL, NULL, NULL, NULL},
     "svd, cov"},
    
    {"pcareduce", CAT_DATASCI, "pcareduce(X, k)",
     "Reduce dimensionality to k principal components.",
     {"pcareduce([1,2,3;4,5,6;7,8,9], 2)", "pcareduce([1,2;3,4;5,6], 1)", NULL, NULL, NULL, NULL},
     "pca"},
    
    {"kmeans", CAT_DATASCI, "[idx,C] = kmeans(X, k)",
     "K-means clustering. idx=cluster assignments, C=centroids.",
     {"kmeans([1,1;1,2;8,8;8,9], 2)", "kmeans([0,0;1,0;0,1;1,1], 2)", NULL, NULL, NULL, NULL},
     "silhouette, pdist"},
    
    {"silhouette", CAT_DATASCI, "silhouette(X, idx)",
     "Silhouette values for cluster quality. Values in [-1, 1], higher is better.",
     {"silhouette([1,1;2,2;8,8;9,9], [1;1;2;2])", NULL, NULL, NULL, NULL, NULL},
     "kmeans"},
    
    {"pdist", CAT_DATASCI, "pdist(X)",
     "Pairwise distances between rows of X.",
     {"pdist([0,0;3,4;0,5])", "pdist([1,0;0,1;1,1])", NULL, NULL, NULL, NULL},
     "kmeans"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Text analysis functions */
const FuncDoc func_text_docs[] = {
    {"bow", CAT_TEXT, "bow(doc1, doc2, ...)",
     "Bag of words. Creates term-document matrix from text documents.",
     {"bow(\"hello world\", \"world peace\")", NULL, NULL, NULL, NULL, NULL},
     "tfidf"},
    
    {"tfidf", CAT_TEXT, "tfidf(doc1, doc2, ...)",
     "TF-IDF weighted term-document matrix.",
     {"tfidf(\"hello world\", \"world peace\")", NULL, NULL, NULL, NULL, NULL},
     "bow"},
    
    {"wordvec", CAT_TEXT, "wordvec(word)",
     "Simple word vector representation.",
     {"wordvec(\"hello\")", NULL, NULL, NULL, NULL, NULL},
     "bow, tfidf"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
