% MATLAB script for Fisher Iris data classification with cross-validation
% 1. Load the Fisher Iris dataset
load fisheriris; % Loads 'meas' (measurements) and 'species' (labels)
X = meas;
Y = species;
% Set random seed for reproducibility of cross-validation partitions
rng(1);
% 2. Complex Array Actions: Data Preprocessing and Transformation
% Standardize features (mean 0, std 1) using the z-score function
X_scaled = zscore(X);
% Convert species names to numeric labels (1, 2, 3) for easier matrix operations
% [~, ~, labels] = unique(Y); % Alternative method
% Or use grp2idx for more control and consistent mapping
[labels, ~, ~] = grp2idx(Y);
num_classes = max(labels); % Get the number of unique classes
% Perform one-hot encoding for the target variable Y
% Create a matrix of size [num_samples, num_classes] with zeros and ones
num_samples = size(X_scaled, 1);
Y_onehot = zeros(num_samples, num_classes);
for i = 1:num_classes
    Y_onehot(:, i) = (labels == i);
end
% 3. Building a Model and Cross Validation
% Define the number of folds for k-fold cross-validation
K = 10;
% Create a stratified cross-validation partition object to ensure
% each fold has roughly the same proportion of each species
cv_partition = cvpartition(Y, 'KFold', K);
% Initialize an array to store classification errors for each fold
fold_errors = zeros(K, 1);
fprintf('Starting %d-fold cross-validation...\n', K);
% Loop through each fold
for i = 1:K
    % Get training and testing indices for the current fold
    train_idx = training(cv_partition, i);
    test_idx = test(cv_partition, i);
    % Extract training and testing data using logical indexing (array action)
    X_train = X_scaled(train_idx, :);
    Y_train = Y(train_idx);
    X_test = X_scaled(test_idx, :);
    Y_test = Y(test_idx);
    % Train a Support Vector Machine (SVM) model
    % The 'fitcsvm' function handles multi-class classification internally
    % with appropriate options (e.g., 'ClassNames' is handled by default).
    model = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', false);
    % Predict the species for the test set
    Y_pred = predict(model, X_test);
    % Calculate the error rate for this fold
    % Error is the proportion of misclassified observations
    error_rate = sum(~strcmp(Y_pred, Y_test)) / numel(Y_test);
    fold_errors(i) = error_rate;
    fprintf('  Fold %d error rate: %.4f\n', i, error_rate);
end
% 4. Evaluate Performance
% Calculate the average error rate across all folds
mean_error = mean(fold_errors);
mean_accuracy = 1 - mean_error;
fprintf('\nAverage cross-validation error rate: %.4f\n', mean_error);
fprintf('Average cross-validation accuracy: %.2f%%\n', mean_accuracy * 100);
% Note: Built-in functions like 'crossval(Mdl)' provide a simpler, automated
% way to get the generalized error for many standard models.
