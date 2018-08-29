1. Run "featurization.R". It featurizes the smoted dataset given in file "Balanced_equal_set_sorted.csv". This data is saved in featuredata_smote.csv file. (Already provided  for direct use)
  

2. Run "testingExperimentsModified.R" to make predictions for leave one gene out method. In this file
  line 48: change threshold value to .1295
  line 114: final predictions are saved in "final_prediction.csv" file


3. Run "calculateclassificationscore.R" to get the values of area under ROC and PR-curves. 
   Choose file "final_predictions.csv"
   The accuracy, area under ROC and PR-curves are saved on features variable. Print it to see the results. 


[ The feature contaninning file is already given named "featuredata_smote.csv".
  Balanced_equal_set.csv contains the gene sequences.
  Balanced_equal_set_sorted.csv contains the genes sorted chronologically. 
  best predictions are given in "pred_smote_1393.csv file".
 ]