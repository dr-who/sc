#!/bin/bash
# Document Clustering Demo using sc
# Demonstrates: k-means clustering + PCA dimensionality reduction

SC=${SC:-./bin/sc}

echo "=== Document Clustering with PCA Demo ==="
echo ""
echo "8 documents with 4 features (simulating bag-of-words)"
echo "Docs 1-4: animal/pet topics"
echo "Docs 5-8: weather/outdoor topics"
echo ""

$SC << 'EOF'
X = [2,1,0,0; 1,2,0,0; 2,2,0,1; 1,1,0,0; 0,0,2,1; 0,0,1,2; 0,1,2,2; 0,0,1,1]
rows(X)
cols(X)
[I,C] = kmeans(X, 2)
I
C
P = pcareduce(X, 2)
P
silhouette(X, I)
quit
EOF

echo ""
echo "I = cluster assignments, C = centroids, P = 2D PCA projection"
