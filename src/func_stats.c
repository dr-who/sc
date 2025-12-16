/*
 * func_stats.c - Statistics function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_stats_docs[] = {
    {"sum", CAT_STATS, "sum(v)",
     "Sum of elements. For matrix, returns column sums.",
     {"sum([1,2,3,4,5]) = 15", "sum([1;2;3]) = 6", "sum(1:100) = 5050", "sum([]) = 0", NULL, NULL},
     "prod, mean, cumsum"},
    
    {"prod", CAT_STATS, "prod(v)",
     "Product of elements. For matrix, returns column products.",
     {"prod([1,2,3,4,5]) = 120", "prod([2;2;2]) = 8", "prod(1:5) = 120", NULL, NULL, NULL},
     "sum, cumprod"},
    
    {"mean", CAT_STATS, "mean(v)",
     "Arithmetic mean. mean(v) = sum(v) / length(v).",
     {"mean([1,2,3,4,5]) = 3", "mean([10,20,30]) = 20", "mean(1:100) = 50.5", "mean([2,4,6,8]) = 5", NULL, NULL},
     "median, geomean, harmmean"},
    
    {"median", CAT_STATS, "median(v)",
     "Median value. Middle value when sorted (average of two middle values if even count).",
     {"median([1,2,3,4,5]) = 3", "median([1,2,3,4]) = 2.5", "median([5,1,3]) = 3", "median([1]) = 1", NULL, NULL},
     "mean, mode, prctile"},
    
    {"mode", CAT_STATS, "mode(v)",
     "Most frequent value. Returns smallest if tie.",
     {"mode([1,2,2,3]) = 2", "mode([1,1,2,2,3]) = 1", "mode([1,2,3]) = 1", NULL, NULL, NULL},
     "mean, median"},
    
    {"std", CAT_STATS, "std(v)",
     "Standard deviation (sample, N-1 denominator).",
     {"std([1,2,3,4,5]) = 1.5811", "std([2,4,4,4,5,5,7,9]) = 2.1381", "std([1,1,1]) = 0", NULL, NULL, NULL},
     "var, mean"},
    
    {"var", CAT_STATS, "var(v)",
     "Variance (sample, N-1 denominator). var = std^2.",
     {"var([1,2,3,4,5]) = 2.5", "var([2,4,4,4,5,5,7,9]) = 4.5714", "var([1,1,1]) = 0", NULL, NULL, NULL},
     "std, mean"},
    
    {"min", CAT_STATS, "min(v) or min(a, b)",
     "Minimum value. For two args, returns smaller.",
     {"min([3,1,4,1,5]) = 1", "min(5, 3) = 3", "min([10,20,5,15]) = 5", "min(-5, 5) = -5", NULL, NULL},
     "max, range"},
    
    {"max", CAT_STATS, "max(v) or max(a, b)",
     "Maximum value. For two args, returns larger.",
     {"max([3,1,4,1,5]) = 5", "max(5, 3) = 5", "max([10,20,5,15]) = 20", "max(-5, 5) = 5", NULL, NULL},
     "min, range"},
    
    {"range", CAT_STATS, "range(v)",
     "Range: max - min.",
     {"range([1,2,3,4,5]) = 4", "range([10,20,30]) = 20", "range([5,5,5]) = 0", NULL, NULL, NULL},
     "min, max"},
    
    {"geomean", CAT_STATS, "geomean(v)",
     "Geometric mean. (product of elements)^(1/n). For positive values only.",
     {"geomean([1,2,4,8]) = 2.8284", "geomean([2,8]) = 4", "geomean([1,10,100]) = 10", NULL, NULL, NULL},
     "mean, harmmean"},
    
    {"harmmean", CAT_STATS, "harmmean(v)",
     "Harmonic mean. n / sum(1/x_i). For positive values only.",
     {"harmmean([1,2,4]) = 1.7143", "harmmean([2,3,6]) = 3", "harmmean([1,1,1]) = 1", NULL, NULL, NULL},
     "mean, geomean"},
    
    {"skewness", CAT_STATS, "skewness(v)",
     "Skewness. Measure of asymmetry. 0 for symmetric distributions.",
     {"skewness([1,2,3,4,5]) = 0", "skewness([1,2,2,3,10]) = 1.329", NULL, NULL, NULL, NULL},
     "kurtosis, std"},
    
    {"kurtosis", CAT_STATS, "kurtosis(v)",
     "Kurtosis. Measure of tail heaviness. 3 for normal distribution.",
     {"kurtosis([1,2,3,4,5]) = 1.7", "kurtosis([1,1,5,5]) = 1", NULL, NULL, NULL, NULL},
     "skewness, std"},
    
    {"mad", CAT_STATS, "mad(v)",
     "Mean absolute deviation from median.",
     {"mad([1,2,3,4,5]) = 1.2", "mad([1,1,1,1,10]) = 1.44", NULL, NULL, NULL, NULL},
     "std, iqr"},
    
    {"iqr", CAT_STATS, "iqr(v)",
     "Interquartile range. 75th percentile - 25th percentile.",
     {"iqr([1,2,3,4,5,6,7]) = 4", "iqr(1:100) = 50", NULL, NULL, NULL, NULL},
     "prctile, mad"},
    
    {"rms", CAT_STATS, "rms(v)",
     "Root mean square. sqrt(mean(v.^2)).",
     {"rms([1,2,3]) = 2.1602", "rms([3,4]) = 3.5355", "rms([1,1,1]) = 1", NULL, NULL, NULL},
     "mean, meansq"},
    
    {"sumsq", CAT_STATS, "sumsq(v)",
     "Sum of squares. sum(v.^2).",
     {"sumsq([1,2,3]) = 14", "sumsq([3,4]) = 25", "sumsq(1:5) = 55", NULL, NULL, NULL},
     "meansq, rms"},
    
    {"meansq", CAT_STATS, "meansq(v)",
     "Mean of squares. mean(v.^2).",
     {"meansq([1,2,3]) = 4.6667", "meansq([3,4]) = 12.5", "meansq([2,2,2]) = 4", NULL, NULL, NULL},
     "sumsq, rms"},
    
    {"prctile", CAT_STATS, "prctile(v, p)",
     "p-th percentile of v. p in range [0, 100].",
     {"prctile([1:10], 50) = 5.5", "prctile([1:100], 25) = 25.75", "prctile([1,2,3,4], 75) = 3.25", NULL, NULL, NULL},
     "median, iqr"},
    
    {"zscore", CAT_STATS, "zscore(v)",
     "Standardize to zero mean and unit variance. (v - mean) / std.",
     {"zscore([1,2,3,4,5])", "zscore([10,20,30])", NULL, NULL, NULL, NULL},
     "std, mean, normalize"},
    
    {"cov", CAT_STATS, "cov(X)",
     "Covariance matrix. For vector, returns variance.",
     {"cov([1,2,3,4,5]) = 2.5", "cov([1,2;3,4;5,6])", NULL, NULL, NULL, NULL},
     "corrcoef, var"},
    
    {"corrcoef", CAT_STATS, "corrcoef(X)",
     "Correlation coefficient matrix. Normalized covariance.",
     {"corrcoef([1,2;3,4;5,6])", "corrcoef([1,2,3],[1,2,3])", NULL, NULL, NULL, NULL},
     "cov"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Probability distributions */
const FuncDoc func_prob_docs[] = {
    {"normpdf", CAT_PROB, "normpdf(x) or normpdf(x, mu, sigma)",
     "Normal (Gaussian) probability density function.",
     {"normpdf(0) = 0.3989", "normpdf(0, 0, 1) = 0.3989", "normpdf(1, 0, 1) = 0.2420", "normpdf(0, 0, 2) = 0.1995", NULL, NULL},
     "normcdf, norminv"},
    
    {"normcdf", CAT_PROB, "normcdf(x) or normcdf(x, mu, sigma)",
     "Normal cumulative distribution function.",
     {"normcdf(0) = 0.5", "normcdf(1) = 0.8413", "normcdf(1.96) = 0.975", "normcdf(-1) = 0.1587", "normcdf(2) = 0.9772", NULL},
     "normpdf, norminv"},
    
    {"norminv", CAT_PROB, "norminv(p) or norminv(p, mu, sigma)",
     "Inverse normal CDF (quantile function).",
     {"norminv(0.5) = 0", "norminv(0.975) = 1.96", "norminv(0.8413) = 1", "norminv(0.025) = -1.96", NULL, NULL},
     "normcdf, normpdf"},
    
    {"tcdf", CAT_PROB, "tcdf(t, df)",
     "Student's t cumulative distribution function.",
     {"tcdf(0, 10) = 0.5", "tcdf(2.228, 10) = 0.975", "tcdf(-2, 5) = 0.0509", NULL, NULL, NULL},
     "tinv"},
    
    {"chi2cdf", CAT_PROB, "chi2cdf(x, df)",
     "Chi-squared cumulative distribution function.",
     {"chi2cdf(0, 5) = 0", "chi2cdf(5, 5) = 0.5841", "chi2cdf(11.07, 5) = 0.95", NULL, NULL, NULL},
     "tcdf"},
    
    {"rand", CAT_PROB, "rand or rand(n) or rand(m, n)",
     "Uniform random numbers in [0, 1). rand(n) returns n values.",
     {"rand", "rand(5)", "rand(3, 3)", "rand(2, 4)", NULL, NULL},
     "randn, randi"},
    
    {"randn", CAT_PROB, "randn or randn(n) or randn(m, n)",
     "Standard normal random numbers (mean 0, std 1).",
     {"randn", "randn(5)", "randn(3, 3)", "mean(randn(1000))", NULL, NULL},
     "rand, randi"},
    
    {"randi", CAT_PROB, "randi(imax) or randi(imax, n)",
     "Random integers from 1 to imax.",
     {"randi(6)", "randi(10, 5)", "randi(100, 3, 3)", NULL, NULL, NULL},
     "rand, randn"},
    
    {"randperm", CAT_PROB, "randperm(n)",
     "Random permutation of 1:n.",
     {"randperm(5)", "randperm(10)", "sort(randperm(6))", NULL, NULL, NULL},
     "rand"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
