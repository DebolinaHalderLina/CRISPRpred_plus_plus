library(e1071)
library(ROCR)
library(pracma)
#setwd('/media/trishna/Study/L4T1/Thesis/CRISPRpred-master/byUs')
#setwd('..')
#dir = getwd()
#print(dir)
#source(paste0(dir,"/R_src/datapartition/datapartition.R"))
#source(paste0(dir,"/R_src/featureformula/featureformula.R"))
#setwd(paste0(dir,"/data-raw/"))

datapartition = function(data, targetcolumn, leaveonegene = FALSE, genename, nfold = 10) {
  if (leaveonegene == TRUE) {
    if (!genename %in% data[targetcolumn][,1]) {
      cat("Please provide a valid Gene name\n")
    }else{
      training = data[!data[targetcolumn][,1] %in% genename,]
      testing = data[data[targetcolumn][,1] %in% genename,]
    }
  }else{
    len = length(data[,1]) / nfold
    training = 1:(len * (nfold - 1))
    testing = (length(training) + 1):length(data[,1])
    training = data[training,]
    testing = data[testing,]
  }
  return(list(training,testing))
}

featureformula = function(featurelist) {
  predictionfeature = featurelist[length(featurelist)]
  righthand =
    paste(featurelist[-length(featurelist)], collapse = "+")
  paste0(predictionfeature,"~", righthand)
}




#fetaures after running randForest on smoted dataset
randomforestfeatures = read.csv('featuresSMOTE_copy.csv')

#What is fcres? list of genes available in dataSet
fcres = read.csv('Ballanced_equal_set.csv')

features = c()
for(i in 1:length(randomforestfeatures[,5])){
  if(randomforestfeatures[i,5] >1.5){
    features = c(features,toString(randomforestfeatures[i,1]))
  }
}

features[length(features) + 1] = "geneThreshold"
cat('Reading features complete: total',  length(features), ' features ...\n')


#smoted data here
#ontargetdata = dataW01; #khalead bhai er data
ontargetdata = read.csv('Smote_data_1withGeneThresh_Mod.csv')
columns = names(ontargetdata)
for(i in 1:length(columns)){
  if(!columns[i] %in% features){
    ontargetdata[columns[i]] = NULL
  }
}
gene = fcres$Target.gene
ontargetdata = data.frame(ontargetdata, gene)
cat('Reading featured data complete ...\n')

genes = c(
  'CCDC101', 'CD13','CD15', 'CD28' ,'CD33' ,'CD43' ,'CD45' ,'CD5' ,'CUL3', 'H2-K', 'HPRT1' ,'MED12', 'NF1' ,'NF2', 'TADA1', 'TADA2B', 'THY1'
)
data = ontargetdata
predValue = c()
for (i in 1:length(genes)) {
  dataset = datapartition(
    data, targetcolumn = "gene", leaveonegene = TRUE, genename = genes[i]
  )
 
  training = dataset[[1]]
  testing = dataset[[2]]
  bformula = featureformula(features)
  svmmodel = svm(
    as.formula(bformula), training, cost = 1
  )

  svmpred = predict(svmmodel, testing)
  svmpred = as.vector(svmpred)
  predValue = c(predValue, svmpred)
  cat("Running SVM iteration no. = ",i,"  ...\n")
}

#here goes sorted smoted data file
sortedresults = read.csv('Ballanced_equal_set_sorted.csv')
svmpred = predValue
trueV= sortedresults$score_drug_gene_threshold
svmprediction = prediction(svmpred, trueV)
svmROC = performance(svmprediction,"tpr","fpr")
svmACC = performance(svmprediction, "acc")
svmAUC = performance(svmprediction, "auc")
svmPR = performance(svmprediction,"prec","rec")
perf  <- performance(svmprediction, "prec", "rec")
xy    <-
  data.frame(recall = perf@x.values[[1]], precision = perf@y.values[[1]])
xy <- subset(xy,!is.nan(xy$precision))
xy <- rbind(c(0, 0), xy)
svmAUPR  <- trapz(xy$recall, xy$precision)
par(mar = c(4, 5, 1.98, 1))
plot(
  svmROC, main = "ROC-Curve", cex.main = 1.7, cex.lab = 1.7, box.lty = 7, box.lwd = 4,xaxis.cex.axis =
    1.3,yaxis.cex.axis = 1.3, lwd = 4, yaxis.lwd = 4, xaxis.lwd = 4, yaxis.las = 1
)


write.csv(svmpred,"pred_smote_1.50_womfeheat.csv")