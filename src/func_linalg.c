/*
 * func_linalg.c - Linear algebra function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_linalg_docs[] = {
    {"det", CAT_LINALG, "det(A)",
     "Matrix determinant. A must be square.",
     {"det([1,2;3,4]) = -2", "det(eye(3)) = 1", "det([1,0;0,1]) = 1", "det([2,0;0,3]) = 6", NULL, NULL},
     "inv, trace, rank"},
    
    {"inv", CAT_LINALG, "inv(A)",
     "Matrix inverse. A must be square and non-singular.",
     {"inv([1,2;3,4])", "inv(eye(3))", "inv([2,0;0,3])", NULL, NULL, NULL},
     "det, mldivide, linsolve"},
    
    {"trace", CAT_LINALG, "trace(A)",
     "Sum of diagonal elements.",
     {"trace(eye(5)) = 5", "trace([1,2;3,4]) = 5", "trace([1,0,0;0,2,0;0,0,3]) = 6", NULL, NULL, NULL},
     "det, diag"},
    
    {"trans", CAT_LINALG, "trans(A)",
     "Matrix transpose. Rows become columns.",
     {"trans([1,2,3])", "trans([1,2;3,4])", NULL, NULL, NULL, NULL},
     "transpose"},
    
    {"transpose", CAT_LINALG, "transpose(A)",
     "Matrix transpose. Alias for trans.",
     {"transpose([1,2,3])", "transpose([1,2;3,4])", NULL, NULL, NULL, NULL},
     "trans"},
    
    {"rank", CAT_LINALG, "rank(A)",
     "Matrix rank. Number of linearly independent rows (or columns).",
     {"rank([1,0;0,1]) = 2", "rank([1,2;2,4]) = 1", "rank(eye(5)) = 5", "rank(zeros(3,3)) = 0", NULL, NULL},
     "det, null"},
    
    {"cond", CAT_LINALG, "cond(A)",
     "Condition number. Ratio of largest to smallest singular value.",
     {"cond(eye(3))", "cond([1,2;3,4])", "cond([1,0;0,2])", NULL, NULL, NULL},
     "svd, norm"},
    
    {"norm", CAT_LINALG, "norm(A) or norm(A, p)",
     "Matrix or vector norm. Default is 2-norm (Frobenius for matrix).",
     {"norm([3,4]) = 5", "norm([1,2,2]) = 3", "norm([1,0;0,1])", NULL, NULL, NULL},
     "cond"},
    
    {"linsolve", CAT_LINALG, "linsolve(A, b)",
     "Solve linear system Ax = b.",
     {"linsolve([1,0;0,1], [3;4])", "linsolve(eye(3), [1;2;3])", NULL, NULL, NULL, NULL},
     "mldivide, inv"},
    
    {"mldivide", CAT_LINALG, "A \\ b",
     "Matrix left division. Solves Ax = b. Same as linsolve(A, b).",
     {"[1,0;0,1] \\ [3;4]", "eye(2) \\ [5;6]", NULL, NULL, NULL, NULL},
     "linsolve, inv"},
    
    {"chol", CAT_LINALG, "chol(A)",
     "Cholesky decomposition. Returns upper triangular R where A = R'*R. A must be positive definite.",
     {"chol([4,2;2,2])", "chol(eye(3))", "chol([9,0;0,4])", NULL, NULL, NULL},
     "lu, qr, svd"},
    
    {"qr", CAT_LINALG, "[Q,R] = qr(A)",
     "QR decomposition. A = Q*R where Q is orthogonal and R is upper triangular.",
     {"qr([1,2;3,4])", "qr(eye(3))", "qr([1,0;0,1])", NULL, NULL, NULL},
     "lu, chol, svd"},
    
    {"svd", CAT_LINALG, "[U,S,V] = svd(A)",
     "Singular value decomposition. A = U*S*V'. S is diagonal with singular values.",
     {"svd([1,2;3,4])", "svd(eye(3))", "svd([1,0;0,2])", NULL, NULL, NULL},
     "eig, qr"},
    
    {"lu", CAT_LINALG, "[L,U] = lu(A)",
     "LU decomposition. A = L*U where L is lower triangular and U is upper triangular.",
     {"lu([1,2;3,4])", "lu(eye(3))", "lu([2,1;4,3])", NULL, NULL, NULL},
     "qr, chol"},
    
    {"eig", CAT_LINALG, "eig(A) or [V,D] = eig(A)",
     "Eigenvalues and eigenvectors. D is diagonal matrix of eigenvalues, V has eigenvectors.",
     {"eig([1,2;2,1])", "eig(eye(3))", "eig([2,0;0,3])", NULL, NULL, NULL},
     "svd"},
    
    {"null", CAT_LINALG, "null(A)",
     "Null space basis. Columns span the null space of A.",
     {"null([1,2,3;4,5,6])", "null([1,0;0,0])", NULL, NULL, NULL, NULL},
     "rank"},
    
    {"schur", CAT_LINALG, "[U,T] = schur(A)",
     "Schur decomposition. A = U*T*U' where U is unitary and T is upper triangular.",
     {"schur([1,2;0,3])", "schur(eye(3))", "schur([2,1;1,2])", NULL, NULL, NULL},
     "eig, qr"},
    
    {"pinv", CAT_LINALG, "pinv(A)",
     "Moore-Penrose pseudoinverse. For non-square or singular matrices.",
     {"pinv([1,2;3,4])", "pinv([1;2;3])", "pinv(eye(3))", NULL, NULL, NULL},
     "inv, svd"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
