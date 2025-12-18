#!/bin/sh
# Matrix stress tests for scalc
# Run with: ./tests/matrix_stress.sh [path_to_sc]

SC="${1:-./bin/sc}"

echo "Matrix stress tests for: $SC"
echo "=============================="
echo ""

# Test 1: SVD on identity matrices
echo "Test 1: SVD on identity matrices..."
result=$($SC << 'EOF'
svd(eye(3))
svd(eye(5))
svd(eye(10))
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: SVD error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 2: SVD on random matrices
echo "Test 2: SVD on general matrices..."
result=$($SC << 'EOF'
svd([1,2,3;4,5,6;7,8,9])
svd([1,0;0,1;1,1])
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: SVD error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 3: rand operations
echo "Test 3: Random matrix operations..."
result=$($SC << 'EOF'
A = rand(10,10)*10
B = A - 5
mean(B(:))
quit
EOF
)
# Check that result isn't all -5 (the bug)
if echo "$result" | grep -q "^-5$"; then
    echo "FAIL: rand subtraction bug"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 4: Matrix multiplication chains
echo "Test 4: Matrix multiplication chains..."
result=$($SC << 'EOF'
A = [1,2;3,4]
B = A * A
C = B * A
D = C * A
trace(D)
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: multiplication chain error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 5: Large matrix operations
echo "Test 5: Large matrix operations..."
result=$($SC << 'EOF'
A = rand(50,50)
B = A' * A
trace(B)
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: large matrix error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 6: PCA
echo "Test 6: PCA operations..."
result=$($SC << 'EOF'
load iris
[coeff, score, latent] = pca(data)
size(coeff)
quit
EOF
)
if echo "$result" | grep -q "Error.*PCA"; then
    echo "FAIL: PCA error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 7: QR decomposition
echo "Test 7: QR decomposition..."
result=$($SC << 'EOF'
qr(eye(3))
qr(eye(5))
qr([1,2,3;4,5,6;7,8,9])
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: QR error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

# Test 8: Repeated operations (memory leak test)
echo "Test 8: Repeated operations..."
result=$($SC << 'EOF'
A = eye(10); svd(A)
A = eye(10); svd(A)
A = eye(10); svd(A)
A = eye(10); svd(A)
A = eye(10); svd(A)
quit
EOF
)
if echo "$result" | grep -q "Error"; then
    echo "FAIL: repeated operations error"
    echo "$result"
    exit 1
else
    echo "PASS"
fi

echo ""
echo "=============================="
echo "All stress tests passed!"
