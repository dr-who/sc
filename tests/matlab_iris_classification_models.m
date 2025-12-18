% Load the Fisher iris data
load fisheriris;
X = meas; % Features: sepal length, sepal width, petal length, petal width
Y = species; % Labels: 'setosa', 'versicolor', 'virginica'
% 1. Data Partitioning (stratified holdout for balanced classes)
rng(1); % For reproducibility
cvp = cvpartition(Y, 'Holdout', 0.3); % 70% training, 30% testing
XTrain = X(cvp.training, :);
YTrain = Y(cvp.training, :);
XTest = X(cvp.test, :);
YTest = Y(cvp.test, :);
% 2. Train various classification models
disp('Training models...');
% k-Nearest Neighbor (k=5)
Mdl_kNN = fitcknn(XTrain, YTrain, 'NumNeighbors', 5);
% Support Vector Machine (using a linear kernel as a starting point)
% For multi-class, fitcsvm uses a one-versus-all approach by default
Mdl_SVM = fitcsvm(XTrain, YTrain, 'KernelFunction', 'linear', 'Standardize', true);
% Decision Tree
Mdl_Tree = fitctree(XTrain, YTrain);
% Naive Bayes
Mdl_NB = fitcnb(XTrain, YTrain);
% 3. Evaluate models on the test set
disp('Evaluating models...');
% k-NN
YPred_kNN = predict(Mdl_kNN, XTest);
Acc_kNN = sum(strcmp(YPred_kNN, YTest)) / numel(YTest) * 100;
% SVM
YPred_SVM = predict(Mdl_SVM, XTest);
Acc_SVM = sum(strcmp(YPred_SVM, YTest)) / numel(YTest) * 100;
% Decision Tree
YPred_Tree = predict(Mdl_Tree, XTest);
Acc_Tree = sum(strcmp(YPred_Tree, YTest)) / numel(YTest) * 100;
% Naive Bayes
YPred_NB = predict(Mdl_NB, XTest);
Acc_NB = sum(strcmp(YPred_NB, YTest)) / numel(YTest) * 100;
% 4. Summarize results in a table for easy comparison
ModelNames = {'k-NN'; 'SVM'; 'Decision Tree'; 'Naive Bayes'};
Accuracy = [Acc_kNN; Acc_SVM; Acc_Tree; Acc_NB];
ComparisonTable = table(ModelNames, Accuracy);
% Sort table by accuracy in descending order
ComparisonTable = sortrows(ComparisonTable, 'Accuracy', 'descend');
disp('--- Model Comparison Results (Test Accuracy %) ---');
disp(ComparisonTable);
% 5. (Optional) Visualize the decision boundary for two features (Petal Length/Width are best)
figure;
gscatter(X(:,3), X(:,4), Y, 'rgb', 'osd');
xlabel('Petal Length (cm)');
ylabel('Petal Width (cm)');
title('Fisher Iris Data Scatter Plot by Species');
legend('Setosa', 'Versicolor', 'Virginica', 'Location', 'best');
