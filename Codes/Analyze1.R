# load data
library(data.table)
library(e1071)
library(pROC)
library(ranger)
setwd("C:/Users/Subho/Dropbox/Blood-brain/")

data1 = data.table(read.csv("Data/All descriptors  - BBB.csv"))
data2 = fread("Data/BBB set 2.csv")
data3 = fread("Data/BBB-415.csv")

names1 = data1[,Source.File]
for(rep in c(".mol2"," ","-","_")){
  names1 = gsub(rep,"",names1)
}
data1 = data1[,Source.File := toupper(names1)]
data11 = data.frame(cbind(data1,data2))

# remove NA cols
NA.counts = as.numeric(apply(data11, 2, function(x) sum(is.na(x))))
data11 = data11[,which(NA.counts==0)]

## add to basak data
names2 = data3[,Name]
for(rep in c(".mol2"," ","-","_")){
  names2 = gsub(rep,"",names2)
}
data3 = data3[,Name := toupper(names2)]

# ******************************** #
# Model Diudea data
# ******************************** #
# map diudea names to basak ones
data4 = data.table(merge(data11, data3, by.x="Source.File",by.y="Name"))
data4[, Class01 := 1]
data4[Class=="n", Class01 := 0]

# external validation
set.seed(11212018)
test = sample(1:nrow(data4),50,replace=F)
rfmod = ranger(as.factor(Class01)~.-Class, data=data4[-test,-c(1,256,512,513)],
               importance="impurity",probability = TRUE, seed=083120181)
preds = predict(rfmod, data4[test,], type="response")$predictions[,2]
roc(data4$Class01[test],preds)

# check variable importances
imps = sort(rfmod$variable.importance, decreasing=TRUE)
round(imps[1:10],2)

# 10-fold cv
nfolds = 10
set.seed(11212018)
n = nrow(data4)
folds = caret::createFolds(1:n,nfolds)
pred.vec1 = rep(0, nrow(data4))
for(k in 1:nfolds){
  test.k = folds[[k]]
  rfmod = ranger(as.factor(Class01)~.-Class, data=data4[-test.k,-c(1,256,512,513)],
                 importance="impurity",probability = TRUE, seed=083120181)
  pred.vec1[test.k] = predict(rfmod, data4[test.k,], type="response")$predictions[,2]
}
roc(data4$Class01,pred.vec1)

# ******************************** #
# Diudea data: Repeated External CV
# ******************************** #
nrep = 1e2
set.seed(11212018)
metrics.mat = matrix(0, nrow=nrep, ncol=6)
model.list = vector("list",nrep)
test.list = model.list
pb = txtProgressBar(0,nrep)
for(i in 1:nrep){
  test.i = sample(1:n, ceiling(.25*n), replace=F)
  ntest = length(test.i)
  imod = ranger(as.factor(Class01)~.-Class, data=data4[-test.i,-c(1,256,512,513)],
                 importance="impurity", probability = TRUE, seed=083120181)
  preds = predict(imod, data4[test.i], type="response")$predictions[,2]
  
  # get AUC
  test.y = data4[test.i, Class01]
  metrics.mat[i,1] = roc(test.y, preds)$auc
  
  # get 20% lift
  pred.df = data.table(cbind(test.y,preds))
  pred.df = pred.df[order(preds, decreasing=T)]
  metrics.mat[i,2] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
  model.list[[i]] = imod
  test.list[[i]] = test.i
  setTxtProgressBar(pb,i)
  
  # get binary classification metrics
  imod = ranger(as.factor(Class01)~.-Class, data=data4[-test.i,-c(1,256,512,513)], seed=083120181)
  pred.y = as.numeric(paste(predict(imod, data4[test.i], type="response")$predictions))
  TP = sum(test.y*pred.y)
  FN = sum(test.y) - TP
  TN = sum((1-test.y)*(1-pred.y))
  FP = sum(1-test.y) - TN
  SE = TP/(TP+FN)
  SP = TN/(TN+FP)
  Q = (TP+TN)/ntest
  C = (TP*TN - FN*FP)/sqrt(TP+FN)/sqrt(TP+FP)/sqrt(TN+FN)/sqrt(TN+FP)
  metrics.mat[i,3:6] = c(SE,SP,Q,C)
}
close(pb)
Diudea.metrics = data.frame(cbind(1:nrep, metrics.mat))
names(Diudea.metrics) = c("Test.set","AUC","Lift","SE","SP","Q","C")

# ******************************** #
# Diudea data: top theta% preds
# ******************************** #
set.seed(11212018)
auc.vec1 = rep(0,20)
theta.vec = (0:19)/20
out.mat1 = matrix(0,nrow=20, ncol=2)

for(i in 1:20){
  
  cat("\nDoing theta =",theta.vec[i],"\n")
  # Initialize quantities
  auc.vec = rep(0, ntest)
  lift.vec = auc.vec
  
  pb = txtProgressBar(0,ntest)
  for(j in 1:ntest){
    test.j = test.list[[j]]
    imps = sort(model.list[[j]]$variable.importance, decreasing=TRUE)
    formula.x = paste(names(imps)[1:ceiling(theta.vec[i]*length(imps))],collapse = "+")
    my.formula = paste0("as.factor(Class01)~",formula.x)
    jmod = ranger(as.formula(my.formula), data=data4[-test.j],
                  importance="impurity",probability = TRUE, seed=083120181)
    preds = predict(jmod, data4[test.j], type="response")$predictions[,2]
    
    # get AUC
    test.y = data4[test.j, Class01]
    auc.vec[j] = roc(test.y, preds)$auc
    
    # get 20% lift
    pred.df = data.table(cbind(test.y,preds))
    pred.df = pred.df[order(preds, decreasing=T)]
    lift.vec[j] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
    setTxtProgressBar(pb,j)
  }
  close(pb)
  
  # save mean outputs
  out.mat1[i,] = c(mean(auc.vec),mean(lift.vec))
}
close(pb)
out.mat1

# ******************************** #
# Model Basak data
# ******************************** #
basak2 = fread("Data/BBB-413-Polly.csv")
names(basak2)[1] = "No."
basak3 = fread("Data/BBB-413-Triplet.csv")
names(basak3)[1] = "No."
basak = Reduce(merge, list(data3, basak2, basak3))
# remove NA
NA.rows = apply(basak,1,function(x) -99 %in% x)
basak = basak[NA.rows==0]
basak = basak[Name %in% data4$Source.File]
basak[, Class01 := 1]
basak[Class=="n", Class01 := 0]

# external validation
set.seed(11212018)
test = sample(1:nrow(basak),50,replace=F)
rfmod1 = ranger(as.factor(Class01)~., data=basak[-test,-(1:5)],
               importance="impurity",probability = TRUE, seed=083120181)
preds = predict(rfmod1, basak[test,], type="response")$predictions[,2]
roc(basak$Class01[test],preds)

# check variable importances
imps = sort(rfmod1$variable.importance, decreasing=TRUE)
round(imps[1:10],2)

# 10-fold cv
nfolds = 10
set.seed(11212018)
n = nrow(basak)
folds = caret::createFolds(1:n,nfolds)
pred.vec2 = rep(0, n)
for(k in 1:nfolds){
  test.k = folds[[k]]
  rfmod1 = ranger(as.factor(Class01)~., data=basak[-test.k,-(1:5)],
                  importance="impurity",probability = TRUE, seed=083120181)
  pred.vec2[test.k] = predict(rfmod1, basak[test.k,], type="response")$predictions[,2]
}
roc(basak$Class01,pred.vec2)

# ******************************** #
# Basak data: Repeated External CV
# ******************************** #
ntest = 1e2
set.seed(11212018)
auc.vec = rep(0 ,ntest)
lift.vec = auc.vec
model.list = vector("list",ntest)
test.list = model.list
pb = txtProgressBar(0,ntest)
for(i in 1:ntest){
  test.i = sample(1:n, ceiling(.25*n), replace=F)
  imod = ranger(as.factor(Class01)~.-Class, data=basak[-test.i,-(1:5)],
                importance="impurity",probability = TRUE, seed=083120181)
  preds = predict(imod, basak[test.i], type="response")$predictions[,2]
  
  # get AUC
  test.y = basak[test.i, Class01]
  auc.vec[i] = roc(test.y, preds)$auc
  
  # get 20% lift
  pred.df = data.table(cbind(test.y,preds))
  pred.df = pred.df[order(preds, decreasing=T)]
  lift.vec[i] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
  model.list[[i]] = imod
  test.list[[i]] = test.i
  setTxtProgressBar(pb,i)
}
close(pb)
Basak.metrics = data.frame(Test.set=1:ntest, AUC=auc.vec, Lift=lift.vec)

# ******************************** #
# Basak data: top theta% preds
# ******************************** #
set.seed(11212018)
out.mat2 = matrix(0,nrow=20, ncol=2)

for(i in 1:20){
  
  cat("\nDoing theta =",theta.vec[i],"\n")
  # Initialize quantities
  auc.vec = rep(0, ntest)
  lift.vec = auc.vec
  
  pb = txtProgressBar(0,ntest)
  for(j in 1:ntest){
    test.j = test.list[[j]]
    imps = sort(model.list[[j]]$variable.importance, decreasing=TRUE)
    formula.x = paste(names(imps)[1:ceiling(theta.vec[i]*length(imps))],collapse = "+")
    my.formula = paste0("as.factor(Class01)~",formula.x)
    jmod = ranger(as.formula(my.formula), data=basak[-test.j],
                  importance="impurity",probability = TRUE, seed=083120181)
    preds = predict(jmod, basak[test.j], type="response")$predictions[,2]
    
    # get AUC
    test.y = basak[test.j, Class01]
    auc.vec[j] = roc(test.y, preds)$auc
    
    # get 20% lift
    pred.df = data.table(cbind(test.y,preds))
    pred.df = pred.df[order(preds, decreasing=T)]
    lift.vec[j] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
    setTxtProgressBar(pb,j)
  }
  close(pb)
  
  # save mean outputs
  out.mat2[i,] = c(mean(auc.vec), mean(lift.vec))
}
close(pb)
out.mat2

# ******************************** #
# Combined model
# ******************************** #
combined = merge(data4[,-c(256,512,513)],basak[,-c(1,3:5,204)],by.x="Source.File",by.y="Name")

# external validation
set.seed(11212018)
test = sample(1:nrow(combined),50,replace=F)
rfmod1 = ranger(as.factor(Class01)~.-Class, data=combined[-test],
                importance="impurity",probability = TRUE, seed=083120181)
preds = predict(rfmod1, combined[test], type="response")$predictions[,2]
roc(combined$Class01[test],preds)

# check variable importances
imps = sort(rfmod1$variable.importance, decreasing=TRUE)
round(imps[1:10],2)

# 10-fold cv
set.seed(11212018)
folds = caret::createFolds(1:nrow(combined),nfolds)
pred.vec3 = rep(0, nrow(combined))
for(k in 1:nfolds){
  test.k = folds[[k]]
  rfmod1 = ranger(as.factor(Class01)~.-Class, data=combined[-test.k],
                  importance="impurity", probability = TRUE, seed=083120181)
  pred.vec3[test.k] = predict(rfmod1, combined[test.k], type="response")$predictions[,2]
}
roc(combined$Class01,pred.vec3)

# ******************************** #
# Combined model: Repeated External CV
# ******************************** #
ntest = 1e2
n = nrow(combined)
set.seed(11212018)
auc.vec = rep(0 ,ntest)
lift.vec = auc.vec
model.list = vector("list",ntest)
test.list = model.list
pb = txtProgressBar(0,ntest)
for(i in 1:ntest){
  test.i = sample(1:n, ceiling(.25*n), replace=F)
  imod = ranger(as.factor(Class01)~.-Class, data=combined[-test.i],
                importance="impurity",probability = TRUE, seed=083120181)
  preds = predict(imod, combined[test.i], type="response")$predictions[,2]
  
  # get AUC
  test.y = combined[test.i, Class01]
  auc.vec[i] = roc(test.y, preds)$auc
  
  # get 20% lift
  pred.df = data.table(cbind(test.y,preds))
  pred.df = pred.df[order(preds, decreasing=T)]
  lift.vec[i] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
  model.list[[i]] = imod
  test.list[[i]] = test.i
  setTxtProgressBar(pb,i)
}
close(pb)
Combined.metrics = data.frame(Test.set=1:ntest, AUC=auc.vec, Lift=lift.vec)

# ******************************** #
# Combined model: top theta% preds
# ******************************** #
set.seed(11212018)
out.mat3 = matrix(0,nrow=20, ncol=2)

for(i in 1:20){
  
  cat("\nDoing theta =",theta.vec[i],"\n")
  # Initialize quantities
  auc.vec = rep(0, ntest)
  lift.vec = auc.vec
  
  pb = txtProgressBar(0,ntest)
  for(j in 1:ntest){
    test.j = test.list[[j]]
    imps = sort(model.list[[j]]$variable.importance, decreasing=TRUE)
    formula.x = paste(names(imps)[1:ceiling(theta.vec[i]*length(imps))],collapse = "+")
    my.formula = paste0("as.factor(Class01)~",formula.x)
    jmod = ranger(as.formula(my.formula), data=combined[-test.j],
                  importance="impurity",probability = TRUE, seed=083120181)
    preds = predict(jmod, combined[test.j], type="response")$predictions[,2]
    
    # get AUC
    test.y = combined[test.j, Class01]
    auc.vec[j] = roc(test.y, preds)$auc
    
    # get 20% lift
    pred.df = data.table(cbind(test.y,preds))
    pred.df = pred.df[order(preds, decreasing=T)]
    lift.vec[j] = mean(pred.df[(1:ceiling(.2*nrow(pred.df))),test.y])/.2
    setTxtProgressBar(pb,j)
  }
  close(pb)
  
  # save mean outputs
  out.mat3[i,] = c(mean(auc.vec), mean(lift.vec))
}
close(pb)
out.mat3

# ******************************** #
# Plot final results
# ******************************** #

pdf("plot_AUC.pdf",height=5,width=5)
plot(theta.vec*100,out.mat1[,1], type='l',lwd=2,
     xlab="Percent of top predictors selected", ylab="AUC",
     ylim=c(.6,.9))
lines(theta.vec*100,out.mat2[,1], lwd=2, col="red")
lines(theta.vec*100,out.mat3[,1], lwd=2, col="blue")
legend("bottomright",c("Basak Lab","Diudea Lab","Combined"),
       lwd=2, col=c("red","black","blue"),
       title="Descriptor set")
dev.off()

pdf("plot_Lift.pdf",height=5,width=5)
plot(theta.vec*100,out.mat1[,2], type='l',lwd=2,
     xlab="Percent of top predictors selected", ylab="20% Lift",
     ylim=c(4,5))
lines(theta.vec*100,out.mat2[,2], lwd=2, col="red")
lines(theta.vec*100,out.mat3[,2], lwd=2, col="blue")
legend("bottomright",c("Basak Lab","Diudea Lab","Combined"),
       lwd=2, col=c("red","black","blue"),
       title="Descriptor set")
dev.off()

# ******************************** #
# Save outputs
# ******************************** #
# full model
Full.model.metrics = list(Basak.metrics, Diudea.metrics, Combined.metrics)
save(Full.model.metrics, file="Full_model_metrics.Rda")

# validation metrics
Validation.metrics = list(out.mat1, out.mat2, out.mat3)
save(Validation.metrics, file="Validation_metrics.Rda")
