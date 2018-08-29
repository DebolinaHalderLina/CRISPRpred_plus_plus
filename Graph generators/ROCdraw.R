library(ROCR)
#data(ROCR.simple)
#newob1= read.csv(file.choose())
#newob2= read.csv(file.choose())
#newob3= read.csv(file.choose())
#newob4= read.csv(file.choose())
#newob5= read.csv(file.choose())

#observed = read.csv(file.choose())

preds <- cbind(p1 = newob1$x, p1 = newob2$x, p1 = newob3$x,p1 = newob4$x,p1 = newob5$x)
n <- 5 # you have n models
colors <- c('red', 'blue','green','orange','black') # 2 colors
for (i in 1:n) {
  plot(performance(prediction(preds[,i],observed$result),"mat",), 
       add=(i!=1),col=colors[i],lwd=2)
  
}

#for (i in 1:n) {
  #plot(performance(prediction(preds[,i],observed$score_drug_gene_threshold),"mat"), 
   #    add=(i!=1),col=colors[3],lwd=2)
  #plot(performance(prediction(preds[,i],observed$result),"sens","spec"), 
       #add=(i!=1),col=colors[i],lwd=2)
#}

#for (i in 1:n) {
  #plot(performance(prediction(preds[,i],observed$score_drug_gene_threshold),"mat"), 
       #add=(i!=1),col=colors[2],lwd=2)
  #plot(performance(prediction(preds[,i],observed$result),"prec","rec"), 
       #add=(i!=1),col=colors[i],lwd=2)
#}