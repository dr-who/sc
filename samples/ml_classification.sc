load fisheriris
mdl_knn = fitcknn(meas, species, 5)
pred_knn = predict(mdl_knn, meas)
accuracy(species, pred_knn)
confusionmat(species, pred_knn)
mdl_nb = fitcnb(meas, species)
pred_nb = predict(mdl_nb, meas)
accuracy(species, pred_nb)
confusionmat(species, pred_nb)
mdl_tree = fitctree(meas, species)
pred_tree = predict(mdl_tree, meas)
accuracy(species, pred_tree)
quit
