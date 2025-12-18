#!/bin/bash
# PCA and K-Means Clustering Demo

SC=${1:-bin/sc}

echo "=== PCA + K-Means Clustering Demo ==="
echo ""
echo "Step 1: Create 8 points in 2 natural clusters"
$SC << 'EOF'
X = [1.2, 0.8; 0.9, 1.1; 1.0, 1.0; 1.1, 0.9; 5.1, 4.9; 4.8, 5.2; 5.0, 5.0; 4.9, 5.1]
quit
EOF

echo ""
echo "Step 2: Run k-means with k=2"
$SC << 'EOF'
X = [1.2, 0.8; 0.9, 1.1; 1.0, 1.0; 1.1, 0.9; 5.1, 4.9; 4.8, 5.2; 5.0, 5.0; 4.9, 5.1]
[idx, centroids] = kmeans(X, 2)
idx
quit
EOF

echo ""
echo "Step 3: Show centroids"
$SC << 'EOF'
X = [1.2, 0.8; 0.9, 1.1; 1.0, 1.0; 1.1, 0.9; 5.1, 4.9; 4.8, 5.2; 5.0, 5.0; 4.9, 5.1]
[idx, centroids] = kmeans(X, 2)
centroids
quit
EOF

echo ""
echo "Step 4: Clustering quality (silhouette score)"
$SC << 'EOF'
X = [1.2, 0.8; 0.9, 1.1; 1.0, 1.0; 1.1, 0.9; 5.1, 4.9; 4.8, 5.2; 5.0, 5.0; 4.9, 5.1]
[idx, centroids] = kmeans(X, 2)
silhouette(X, idx)
quit
EOF

echo ""
echo "Step 5: PCA on 5D data reduced to 2D"
$SC << 'EOF'
Y = [1,2,3,0,0; 2,1,3,0,0; 1,1,2,0,0; 0,0,0,3,2; 0,0,0,2,3; 0,0,0,3,3]
[idx, c] = kmeans(Y, 2)
idx
pcareduce(Y, 2)
quit
EOF

echo ""
echo "=== Complete! ==="
