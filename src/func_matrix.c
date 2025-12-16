/*
 * func_matrix.c - Matrix creation and manipulation function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_matrix_docs[] = {
    {"zeros", CAT_MATRIX, "zeros(n) or zeros(m, n)",
     "Create matrix of zeros.",
     {"zeros(3)", "zeros(2, 4)", "zeros(1, 5)", "zeros(3, 1)", NULL, NULL},
     "ones, eye"},
    
    {"ones", CAT_MATRIX, "ones(n) or ones(m, n)",
     "Create matrix of ones.",
     {"ones(3)", "ones(2, 4)", "sum(ones(5,5)) = 25", NULL, NULL, NULL},
     "zeros, eye"},
    
    {"eye", CAT_MATRIX, "eye(n) or eye(m, n)",
     "Identity matrix. Ones on diagonal, zeros elsewhere.",
     {"eye(3)", "eye(2, 4)", "A*eye(n) = A", "det(eye(5)) = 1", NULL, NULL},
     "zeros, ones, diag"},
    
    {"diag", CAT_MATRIX, "diag(v) or diag(A) or diag(A, k)",
     "Create diagonal matrix from vector, or extract diagonal from matrix.",
     {"diag([1,2,3])", "diag([1,2;3,4]) = [1;4]", "diag(eye(3)) = [1;1;1]", "diag([1,2], 1)", NULL, NULL},
     "eye, blkdiag"},
    
    {"linspace", CAT_MATRIX, "linspace(a, b, n)",
     "Create n linearly spaced points from a to b (inclusive).",
     {"linspace(0, 1, 5)", "linspace(0, 10, 11)", "linspace(-1, 1, 3)", "linspace(0, 2*pi, 100)", NULL, NULL},
     "logspace"},
    
    {"logspace", CAT_MATRIX, "logspace(a, b, n)",
     "Create n logarithmically spaced points from 10^a to 10^b.",
     {"logspace(0, 2, 3)", "logspace(-1, 1, 3)", "logspace(0, 3, 4)", NULL, NULL, NULL},
     "linspace"},
    
    {"magic", CAT_MATRIX, "magic(n)",
     "Magic square. All rows, columns, and diagonals sum to same value.",
     {"magic(3)", "sum(magic(4))", "magic(5)", NULL, NULL, NULL},
     "pascal, hilb"},
    
    {"pascal", CAT_MATRIX, "pascal(n)",
     "Pascal matrix. Contains binomial coefficients.",
     {"pascal(4)", "pascal(5)", "det(pascal(5)) = 1", NULL, NULL, NULL},
     "magic, hilb"},
    
    {"hilb", CAT_MATRIX, "hilb(n)",
     "Hilbert matrix. H(i,j) = 1/(i+j-1). Notoriously ill-conditioned.",
     {"hilb(3)", "hilb(4)", "cond(hilb(5))", "inv(hilb(3))", NULL, NULL},
     "invhilb, pascal"},
    
    {"invhilb", CAT_MATRIX, "invhilb(n)",
     "Inverse Hilbert matrix. Exact integer entries.",
     {"invhilb(3)", "hilb(3)*invhilb(3)", NULL, NULL, NULL, NULL},
     "hilb"},
    
    {"vander", CAT_MATRIX, "vander(v)",
     "Vandermonde matrix. V(i,j) = v(i)^(n-j).",
     {"vander([1,2,3])", "vander([1,2,3,4])", NULL, NULL, NULL, NULL},
     "compan"},
    
    {"toeplitz", CAT_MATRIX, "toeplitz(c) or toeplitz(c, r)",
     "Toeplitz matrix. Constant along diagonals.",
     {"toeplitz([1,2,3])", "toeplitz([1,2,3], [1,4,5])", NULL, NULL, NULL, NULL},
     "hankel"},
    
    {"hankel", CAT_MATRIX, "hankel(c) or hankel(c, r)",
     "Hankel matrix. Constant along anti-diagonals.",
     {"hankel([1,2,3])", "hankel([1,2,3], [3,4,5])", NULL, NULL, NULL, NULL},
     "toeplitz"},
    
    {"compan", CAT_MATRIX, "compan(p)",
     "Companion matrix for polynomial coefficients.",
     {"compan([1,0,-2])", "eig(compan([1,-6,11,-6]))", NULL, NULL, NULL, NULL},
     "vander"},
    
    {"reshape", CAT_MATRIX, "reshape(A, m, n)",
     "Reshape matrix to m-by-n. Total elements must match.",
     {"reshape([1,2,3,4,5,6], 2, 3)", "reshape([1;2;3;4], 2, 2)", "reshape(1:12, 3, 4)", NULL, NULL, NULL},
     "size"},
    
    {"repmat", CAT_MATRIX, "repmat(A, m, n)",
     "Repeat matrix m times vertically and n times horizontally.",
     {"repmat([1,2], 2, 3)", "repmat(eye(2), 2, 2)", "repmat([1;2], 3, 1)", NULL, NULL, NULL},
     "kron"},
    
    {"rot90", CAT_MATRIX, "rot90(A) or rot90(A, k)",
     "Rotate matrix 90 degrees counterclockwise k times.",
     {"rot90([1,2;3,4])", "rot90([1,2,3], 2)", "rot90(eye(3))", NULL, NULL, NULL},
     "flip, fliplr, flipud"},
    
    {"flip", CAT_MATRIX, "flip(A) or flip(A, dim)",
     "Flip matrix along dimension. Default flips rows.",
     {"flip([1,2,3])", "flip([1;2;3])", "flip([1,2;3,4])", "flip([1,2;3,4], 2)", NULL, NULL},
     "fliplr, flipud"},
    
    {"fliplr", CAT_MATRIX, "fliplr(A)",
     "Flip matrix left-right.",
     {"fliplr([1,2,3])", "fliplr([1,2;3,4])", "fliplr(eye(3))", NULL, NULL, NULL},
     "flipud, rot90"},
    
    {"flipud", CAT_MATRIX, "flipud(A)",
     "Flip matrix up-down.",
     {"flipud([1;2;3])", "flipud([1,2;3,4])", "flipud(eye(3))", NULL, NULL, NULL},
     "fliplr, rot90"},
    
    {"triu", CAT_MATRIX, "triu(A) or triu(A, k)",
     "Upper triangular part. k=0 main diagonal, k>0 above, k<0 below.",
     {"triu([1,2,3;4,5,6;7,8,9])", "triu(ones(3,3), 1)", "triu(magic(4))", NULL, NULL, NULL},
     "tril, diag"},
    
    {"tril", CAT_MATRIX, "tril(A) or tril(A, k)",
     "Lower triangular part. k=0 main diagonal, k<0 below, k>0 above.",
     {"tril([1,2,3;4,5,6;7,8,9])", "tril(ones(3,3), -1)", "tril(magic(4))", NULL, NULL, NULL},
     "triu, diag"},
    
    {"blkdiag", CAT_MATRIX, "blkdiag(A, B, ...)",
     "Block diagonal matrix from input matrices.",
     {"blkdiag([1,2;3,4], [5])", "blkdiag(eye(2), eye(3))", NULL, NULL, NULL, NULL},
     "diag"},
    
    {"cat", CAT_MATRIX, "cat(dim, A, B, ...)",
     "Concatenate arrays along dimension dim.",
     {"cat(1, [1,2], [3,4])", "cat(2, [1;2], [3;4])", NULL, NULL, NULL, NULL},
     "horzcat, vertcat"},
    
    {"horzcat", CAT_MATRIX, "[A, B]",
     "Horizontal concatenation. Same as cat(2, A, B).",
     {"[1,2,3]", "[[1;2], [3;4]]", "[eye(2), ones(2,3)]", NULL, NULL, NULL},
     "vertcat, cat"},
    
    {"vertcat", CAT_MATRIX, "[A; B]",
     "Vertical concatenation. Same as cat(1, A, B).",
     {"[1;2;3]", "[[1,2]; [3,4]]", "[eye(2); ones(3,2)]", NULL, NULL, NULL},
     "horzcat, cat"},
    
    {"circshift", CAT_MATRIX, "circshift(A, k)",
     "Circular shift of elements. Positive k shifts right/down.",
     {"circshift([1,2,3,4], 1)", "circshift([1;2;3], -1)", "circshift([1,2,3,4], 2)", NULL, NULL, NULL},
     "rot90"},
    
    {"sort", CAT_MATRIX, "sort(v) or sort(A, dim)",
     "Sort elements in ascending order.",
     {"sort([3,1,4,1,5]) = [1,1,3,4,5]", "sort([3;1;2])", "sort([3,1;4,2])", NULL, NULL, NULL},
     "sortrows, unique"},
    
    {"sortrows", CAT_MATRIX, "sortrows(A) or sortrows(A, col)",
     "Sort rows by specified column(s).",
     {"sortrows([3,1;1,2;2,3])", "sortrows([3,1;1,2;2,3], 2)", NULL, NULL, NULL, NULL},
     "sort"},
    
    {"unique", CAT_MATRIX, "unique(v)",
     "Unique elements in sorted order.",
     {"unique([1,2,1,3,2]) = [1,2,3]", "unique([3,1,2,1])", NULL, NULL, NULL, NULL},
     "sort"},
    
    {"kron", CAT_MATRIX, "kron(A, B)",
     "Kronecker tensor product.",
     {"kron([1,2], [3,4])", "kron(eye(2), [1,2;3,4])", NULL, NULL, NULL, NULL},
     "repmat"},
    
    {"dot", CAT_MATRIX, "dot(a, b)",
     "Dot product of vectors.",
     {"dot([1,2,3], [4,5,6]) = 32", "dot([1,0,0], [0,1,0]) = 0", NULL, NULL, NULL, NULL},
     "cross, norm"},
    
    {"cross", CAT_MATRIX, "cross(a, b)",
     "Cross product of 3D vectors.",
     {"cross([1,0,0], [0,1,0]) = [0,0,1]", "cross([1,2,3], [4,5,6])", NULL, NULL, NULL, NULL},
     "dot"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Matrix query functions */
const FuncDoc func_matquery_docs[] = {
    {"size", CAT_MATQUERY, "size(A) or size(A, dim)",
     "Matrix dimensions. Returns [rows, cols] or specific dimension.",
     {"size([1,2,3]) = [1,3]", "size([1;2;3]) = [3,1]", "size(eye(4)) = [4,4]", "size([1,2;3,4], 1) = 2", NULL, NULL},
     "length, numel, rows, cols"},
    
    {"length", CAT_MATQUERY, "length(A)",
     "Length of largest dimension. max(size(A)).",
     {"length([1,2,3]) = 3", "length([1;2;3;4;5]) = 5", "length([1,2;3,4]) = 2", NULL, NULL, NULL},
     "size, numel"},
    
    {"rows", CAT_MATQUERY, "rows(A)",
     "Number of rows. Same as size(A, 1).",
     {"rows([1,2,3]) = 1", "rows([1;2;3]) = 3", "rows(eye(5)) = 5", NULL, NULL, NULL},
     "cols, size"},
    
    {"cols", CAT_MATQUERY, "cols(A)",
     "Number of columns. Same as size(A, 2).",
     {"cols([1,2,3]) = 3", "cols([1;2;3]) = 1", "cols(eye(5)) = 5", NULL, NULL, NULL},
     "rows, size"},
    
    {"numel", CAT_MATQUERY, "numel(A)",
     "Number of elements. rows * cols.",
     {"numel([1,2,3]) = 3", "numel([1,2;3,4]) = 4", "numel(eye(5)) = 25", NULL, NULL, NULL},
     "size, length"},
    
    {"ndims", CAT_MATQUERY, "ndims(A)",
     "Number of dimensions. Always 2 for matrices.",
     {"ndims([1,2,3]) = 2", "ndims(5) = 2", NULL, NULL, NULL, NULL},
     "size"},
    
    {"isempty", CAT_MATQUERY, "isempty(A)",
     "Test if matrix is empty (any dimension is 0).",
     {"isempty([]) = 1", "isempty([1]) = 0", "isempty(zeros(0,5)) = 1", NULL, NULL, NULL},
     "isscalar, isvector"},
    
    {"isscalar", CAT_MATQUERY, "isscalar(A)",
     "Test if A is 1x1.",
     {"isscalar(5) = 1", "isscalar([1,2]) = 0", "isscalar([5]) = 1", NULL, NULL, NULL},
     "isvector, ismatrix"},
    
    {"isvector", CAT_MATQUERY, "isvector(A)",
     "Test if A is row or column vector.",
     {"isvector([1,2,3]) = 1", "isvector([1;2;3]) = 1", "isvector([1,2;3,4]) = 0", NULL, NULL, NULL},
     "isrow, iscolumn, isscalar"},
    
    {"isrow", CAT_MATQUERY, "isrow(A)",
     "Test if A is row vector (1xn).",
     {"isrow([1,2,3]) = 1", "isrow([1;2;3]) = 0", "isrow(5) = 1", NULL, NULL, NULL},
     "iscolumn, isvector"},
    
    {"iscolumn", CAT_MATQUERY, "iscolumn(A)",
     "Test if A is column vector (mx1).",
     {"iscolumn([1;2;3]) = 1", "iscolumn([1,2,3]) = 0", "iscolumn(5) = 1", NULL, NULL, NULL},
     "isrow, isvector"},
    
    {"ismatrix", CAT_MATQUERY, "ismatrix(A)",
     "Test if A is 2D matrix. Always true.",
     {"ismatrix([1,2;3,4]) = 1", "ismatrix(5) = 1", "ismatrix([]) = 1", NULL, NULL, NULL},
     "isvector, issquare"},
    
    {"issquare", CAT_MATQUERY, "issquare(A)",
     "Test if A is square (rows == cols).",
     {"issquare(eye(3)) = 1", "issquare([1,2;3,4]) = 1", "issquare([1,2,3]) = 0", NULL, NULL, NULL},
     "ismatrix"},
    
    {"issorted", CAT_MATQUERY, "issorted(v)",
     "Test if vector is sorted ascending.",
     {"issorted([1,2,3,4]) = 1", "issorted([1,3,2]) = 0", "issorted([5,5,5]) = 1", NULL, NULL, NULL},
     "sort"},
    
    {"issymmetric", CAT_MATQUERY, "issymmetric(A)",
     "Test if A equals its transpose.",
     {"issymmetric([1,2;2,3]) = 1", "issymmetric([1,2;3,4]) = 0", "issymmetric(eye(3)) = 1", NULL, NULL, NULL},
     "trans"},
    
    {"nnz", CAT_MATQUERY, "nnz(A)",
     "Number of nonzero elements.",
     {"nnz([1,0,2,0,3]) = 3", "nnz(eye(5)) = 5", "nnz(zeros(3,3)) = 0", NULL, NULL, NULL},
     "find"},
    
    {"find", CAT_MATQUERY, "find(A)",
     "Indices of nonzero elements.",
     {"find([0,1,0,2,0]) = [2,4]", "find([1,2;0,3])", "find(eye(3))", NULL, NULL, NULL},
     "nnz"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
