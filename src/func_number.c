/*
 * func_number.c - Number theory function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_number_docs[] = {
    {"gcd", CAT_NUMBER, "gcd(a, b)",
     "Greatest common divisor. Largest positive integer dividing both a and b.",
     {"gcd(12, 8) = 4", "gcd(100, 35) = 5", "gcd(17, 13) = 1", "gcd(48, 18) = 6", "gcd(0, 5) = 5", NULL},
     "lcm"},
    
    {"lcm", CAT_NUMBER, "lcm(a, b)",
     "Least common multiple. Smallest positive integer divisible by both a and b.",
     {"lcm(4, 6) = 12", "lcm(3, 5) = 15", "lcm(12, 8) = 24", "lcm(7, 11) = 77", NULL, NULL},
     "gcd"},
    
    {"isprime", CAT_NUMBER, "isprime(n)",
     "Test if n is prime. Returns 1 if prime, 0 otherwise.",
     {"isprime(2) = 1", "isprime(17) = 1", "isprime(4) = 0", "isprime(97) = 1", "isprime(1) = 0", NULL},
     "nextprime, prevprime, primepi"},
    
    {"nextprime", CAT_NUMBER, "nextprime(n)",
     "Smallest prime greater than n.",
     {"nextprime(10) = 11", "nextprime(100) = 101", "nextprime(2) = 3", "nextprime(20) = 23", NULL, NULL},
     "prevprime, isprime"},
    
    {"prevprime", CAT_NUMBER, "prevprime(n)",
     "Largest prime less than n. Returns 0 if n <= 2.",
     {"prevprime(10) = 7", "prevprime(100) = 97", "prevprime(20) = 19", "prevprime(3) = 2", NULL, NULL},
     "nextprime, isprime"},
    
    {"primepi", CAT_NUMBER, "primepi(n)",
     "Prime counting function. Number of primes <= n.",
     {"primepi(10) = 4", "primepi(100) = 25", "primepi(2) = 1", "primepi(1000) = 168", NULL, NULL},
     "nthprime, isprime"},
    
    {"nthprime", CAT_NUMBER, "nthprime(n)",
     "The n-th prime number. nthprime(1) = 2, nthprime(2) = 3, etc.",
     {"nthprime(1) = 2", "nthprime(10) = 29", "nthprime(100) = 541", "nthprime(5) = 11", NULL, NULL},
     "primepi, isprime"},
    
    {"fact", CAT_NUMBER, "fact(n)",
     "Factorial. n! = 1 * 2 * 3 * ... * n. fact(0) = 1.",
     {"fact(0) = 1", "fact(5) = 120", "fact(10) = 3628800", "fact(3) = 6", "fact(7) = 5040", NULL},
     "factorial, ncr, npr"},
    
    {"factorial", CAT_NUMBER, "factorial(n)",
     "Factorial. Alias for fact.",
     {"factorial(5) = 120", "factorial(0) = 1", "factorial(6) = 720", NULL, NULL, NULL},
     "fact"},
    
    {"factorial2", CAT_NUMBER, "factorial2(n)",
     "Double factorial. n!! = n * (n-2) * (n-4) * ... down to 1 or 2.",
     {"factorial2(5) = 15", "factorial2(6) = 48", "factorial2(7) = 105", "factorial2(0) = 1", NULL, NULL},
     "fact"},
    
    {"ncr", CAT_NUMBER, "ncr(n, r)",
     "Binomial coefficient. n choose r = n! / (r! * (n-r)!).",
     {"ncr(5, 2) = 10", "ncr(10, 3) = 120", "ncr(6, 0) = 1", "ncr(6, 6) = 1", "ncr(10, 5) = 252", NULL},
     "npr, fact, comb"},
    
    {"npr", CAT_NUMBER, "npr(n, r)",
     "Permutations. n permute r = n! / (n-r)!.",
     {"npr(5, 2) = 20", "npr(10, 3) = 720", "npr(4, 4) = 24", "npr(6, 2) = 30", NULL, NULL},
     "ncr, fact, perm"},
    
    {"comb", CAT_NUMBER, "comb(n, r)",
     "Combinations. Alias for ncr.",
     {"comb(5, 2) = 10", "comb(10, 5) = 252", "comb(8, 3) = 56", NULL, NULL, NULL},
     "ncr, perm"},
    
    {"perm", CAT_NUMBER, "perm(n, r)",
     "Permutations. Alias for npr.",
     {"perm(5, 2) = 20", "perm(4, 4) = 24", "perm(7, 3) = 210", NULL, NULL, NULL},
     "npr, comb"},
    
    {"fibonacci", CAT_NUMBER, "fibonacci(n)",
     "Fibonacci number. F(0)=0, F(1)=1, F(n)=F(n-1)+F(n-2).",
     {"fibonacci(0) = 0", "fibonacci(1) = 1", "fibonacci(10) = 55", "fibonacci(20) = 6765", "fibonacci(7) = 13", NULL},
     "lucas"},
    
    {"lucas", CAT_NUMBER, "lucas(n)",
     "Lucas number. L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2).",
     {"lucas(0) = 2", "lucas(1) = 1", "lucas(10) = 123", "lucas(5) = 11", NULL, NULL},
     "fibonacci"},
    
    {"catalan", CAT_NUMBER, "catalan(n)",
     "Catalan number. C(n) = (2n)! / ((n+1)! * n!).",
     {"catalan(0) = 1", "catalan(3) = 5", "catalan(5) = 42", "catalan(10) = 16796", NULL, NULL},
     "ncr"},
    
    {"divisors", CAT_NUMBER, "divisors(n)",
     "Number of positive divisors of n.",
     {"divisors(12) = 6", "divisors(1) = 1", "divisors(100) = 9", "divisors(17) = 2", "divisors(36) = 9", NULL},
     "divisorsum, totient"},
    
    {"divisorsum", CAT_NUMBER, "divisorsum(n)",
     "Sum of positive divisors of n.",
     {"divisorsum(12) = 28", "divisorsum(1) = 1", "divisorsum(6) = 12", "divisorsum(28) = 56", NULL, NULL},
     "divisors, totient"},
    
    {"totient", CAT_NUMBER, "totient(n)",
     "Euler's totient. Count of integers 1 to n coprime to n.",
     {"totient(1) = 1", "totient(10) = 4", "totient(12) = 4", "totient(7) = 6", "totient(100) = 40", NULL},
     "eulerphi, gcd"},
    
    {"eulerphi", CAT_NUMBER, "eulerphi(n)",
     "Euler's phi function. Alias for totient.",
     {"eulerphi(10) = 4", "eulerphi(12) = 4", "eulerphi(17) = 16", NULL, NULL, NULL},
     "totient"},
    
    {"mobius", CAT_NUMBER, "mobius(n)",
     "Mobius function. Returns 0 if n has squared prime factor, (-1)^k otherwise.",
     {"mobius(1) = 1", "mobius(2) = -1", "mobius(6) = 1", "mobius(4) = 0", "mobius(30) = -1", NULL},
     "totient"},
    
    {"even", CAT_NUMBER, "even(n)",
     "Test if n is even. Returns 1 if even, 0 if odd.",
     {"even(4) = 1", "even(7) = 0", "even(0) = 1", "even(-2) = 1", NULL, NULL},
     "odd"},
    
    {"odd", CAT_NUMBER, "odd(n)",
     "Test if n is odd. Returns 1 if odd, 0 if even.",
     {"odd(7) = 1", "odd(4) = 0", "odd(1) = 1", "odd(-3) = 1", NULL, NULL},
     "even"},
    
    {"ndigits", CAT_NUMBER, "ndigits(n)",
     "Number of decimal digits in |n|. ndigits(0) = 1.",
     {"ndigits(123) = 3", "ndigits(1000) = 4", "ndigits(0) = 1", "ndigits(-42) = 2", NULL, NULL},
     "digitsum"},
    
    {"digitsum", CAT_NUMBER, "digitsum(n)",
     "Sum of decimal digits of |n|.",
     {"digitsum(123) = 6", "digitsum(999) = 27", "digitsum(1000) = 1", "digitsum(-42) = 6", NULL, NULL},
     "digitroot, ndigits"},
    
    {"digitroot", CAT_NUMBER, "digitroot(n)",
     "Repeated digit sum until single digit. digitroot(n) = 1 + (n-1) mod 9.",
     {"digitroot(123) = 6", "digitroot(999) = 9", "digitroot(12345) = 6", "digitroot(0) = 0", NULL, NULL},
     "digitsum"},
    
    {"isperfect", CAT_NUMBER, "isperfect(n)",
     "Test if n is a perfect number (equals sum of proper divisors).",
     {"isperfect(6) = 1", "isperfect(28) = 1", "isperfect(12) = 0", "isperfect(496) = 1", NULL, NULL},
     "isabundant, isdeficient"},
    
    {"isabundant", CAT_NUMBER, "isabundant(n)",
     "Test if n is abundant (proper divisor sum > n).",
     {"isabundant(12) = 1", "isabundant(6) = 0", "isabundant(18) = 1", "isabundant(20) = 1", NULL, NULL},
     "isperfect, isdeficient"},
    
    {"isdeficient", CAT_NUMBER, "isdeficient(n)",
     "Test if n is deficient (proper divisor sum < n).",
     {"isdeficient(8) = 1", "isdeficient(6) = 0", "isdeficient(10) = 1", "isdeficient(21) = 1", NULL, NULL},
     "isperfect, isabundant"},
    
    {"issquarefree", CAT_NUMBER, "issquarefree(n)",
     "Test if n has no squared prime factors.",
     {"issquarefree(6) = 1", "issquarefree(4) = 0", "issquarefree(30) = 1", "issquarefree(18) = 0", NULL, NULL},
     "mobius, isprime"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
