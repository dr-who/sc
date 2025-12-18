% Simple MATLAB Matrix Test Script using the Iris Dataset
% 1. Load the Iris dataset
% The 'load fisheriris' command loads the data into two variables:
% 'meas' (a matrix of measurements) and 'species' (a cell array of species names).
load fisheriris;
disp('--- Iris Dataset Loaded ---');
% Display the size of the measurement matrix
fprintf('Size of measurement matrix (rows x columns): %d x %d\n', size(meas, 1), size(meas, 2));
disp('The matrix "meas" contains 4 features (columns) for 150 specimens (rows).');
% 2. Matrix Indexing Examples
% Access a specific element (row 10, column 3 - petal length of the 10th flower)
element = meas(10, 3);
fprintf('\n--- Matrix Indexing Examples ---\n');
fprintf('Petal length of the 10th specimen: %.2f cm\n', element);
% Access an entire column (all rows, column 1 - sepal length)
sepal_lengths = meas(:, 1);
fprintf('Number of sepal length measurements (size of column vector): %d\n', length(sepal_lengths));
fprintf('First 5 sepal lengths: %.1f, %.1f, %.1f, %.1f, %.1f\n', sepal_lengths(1:5));
% Access an entire row (row 50, all columns - all measurements of the 50th flower)
measurements_50 = meas(50, :);
disp('Measurements for the 50th specimen (sepal length, sepal width, petal length, petal width):');
disp(measurements_50);
% Access a sub-matrix (first 5 rows and first 2 columns)
sub_matrix = meas(1:5, 1:2);
disp('Sub-matrix of the first 5 rows and first 2 columns:');
disp(sub_matrix);
% Logical indexing to find all sepal lengths greater than 7.0 cm
long_sepal_lengths = meas(meas(:, 1) > 7.0, 1);
fprintf('Number of specimens with sepal length > 7.0 cm: %d\n', length(long_sepal_lengths));
% 3. Simple Matrix Test/Operation
% Calculate the average petal width for all specimens
average_petal_width = mean(meas(:, 4));
fprintf('\n--- Simple Matrix Operation ---\n');
fprintf('Average petal width for all specimens: %.2f cm\n', average_petal_width);
% Separate data for 'setosa' species using logical indexing
% The 'species' variable is a cell array, so we compare strings
is_setosa = strcmp(species, 'setosa');
setosa_measurements = meas(is_setosa, :);
fprintf('Number of setosa specimens found: %d\n', size(setosa_measurements, 1));
