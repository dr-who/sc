% MATLAB Test Script: Comprehensive Matrix Operations
% This script demonstrates various matrix manipulation techniques in MATLAB.
%% 1. Initialization and Random Matrix Creation
% Clear workspace and command window
clear;
clc;
% Set random number generator seed for reproducibility
rng('default');
disp('--- 1. Matrix Creation ---');
% Create a 4x5 matrix of uniform random numbers in the range [0, 1]
A = rand(4, 5);
disp('Matrix A (4x5, random [0, 1]):');
disp(A);
% Create a 3x3 matrix of random integers in the range [1, 10]
B = randi(10, 3, 3);
disp('Matrix B (3x3, random integers [1, 10]):');
disp(B);
%% 2. Matrix Slicing and Indexing
disp('--- 2. Slicing and Indexing ---');
% Extract the first two rows and all columns of A
slice_A_rows = A(1:2, :);
disp('First two rows of A:');
disp(slice_A_rows);
% Extract the elements in the 2nd to 4th rows and 3rd to 5th columns of A
slice_A_sub = A(2:4, 3:5);
disp('Sub-matrix of A (rows 2-4, cols 3-5):');
disp(slice_A_sub);
% Use linear indexing to get specific elements (e.g., the 4th element)
element_4 = A(4);
disp(['The 4th element of A using linear indexing: ', num2str(element_4)]);
% Use logical indexing to find elements in A greater than 0.5
logical_A = A > 0.5;
disp('Logical matrix for elements > 0.5 in A:');
disp(logical_A);
%% 3. Updating Ranges and Elements
disp('--- 3. Updating Ranges and Elements ---');
% Update a specific element: set A(1, 1) to 999
A(1, 1) = 999;
disp('A after setting A(1, 1) to 999:');
disp(A(1,1));
% Update a range: set elements in the last column of A to 0
A(:, end) = 0;
disp('A after setting the last column to 0:');
disp(A);
% Update elements based on logical indexing: set all elements > 0.5 (excluding the 999) to 1
A(A > 0.5) = 1;
disp('A after setting all elements > 0.5 (except A(1,1)) to 1:');
disp(A);
%% 4. Joining (Concatenation)
disp('--- 4. Joining (Concatenation) ---');
% Horizontal concatenation: join two matrices (must have same number of rows)
C = rand(4, 2);
H_joined = [A, C]; % Using comma or space for horizontal concatenation
disp('Horizontally concatenated matrix [A, C] (4x7):');
disp(H_joined);
% Vertical concatenation: join two matrices (must have same number of columns)
D = ones(1, 7);
V_joined = [H_joined; D]; % Using semicolon for vertical concatenation
disp('Vertically concatenated matrix [H_joined; D] (5x7):');
disp(V_joined);
%% 5. Complicated Matrix Operations
disp('--- 5. Complicated Operations ---');
% Matrix multiplication (A*C is not valid because inner dimensions don't match (5 vs 4))
% We can multiply A (4x5) by a 5x2 matrix (like C in H_joined, but we use a new one for clarity)
E = rand(5, 2);
Matrix_Mult_Result = A * E;
disp('Matrix multiplication result (A * E, 4x2):');
disp(Matrix_Mult_Result);
% Element-wise multiplication (.*) (requires matrices of the same size)
% Example using two 4x5 matrices (A and A+1)
Element_Wise_Mult = A .* (A + 1);
disp('Element-wise multiplication (A .* (A+1), 4x5):');
disp(Element_Wise_Mult);
% Solving a linear system (using backslash operator)
% For an equation Ax = b, we can solve for x using x = A\b
% Note: A must be a square matrix for standard linear solving.
% Let's create a square matrix for this example.
F = magic(3); % A 3x3 magic square matrix
b = [1; 3; 5]; % A 3x1 column vector
x = F \ b;
disp('Solution to Fx = b using backslash operator (x = F\b):');
disp(x);
% Get the transpose of B
B_transpose = B';
disp('Transpose of B:');
disp(B_transpose);
% Find eigenvalues of F
eigenvalues_F = eig(F);
disp('Eigenvalues of F:');
disp(eigenvalues_F);
