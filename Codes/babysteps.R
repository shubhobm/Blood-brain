# load data 1
library(data.table)
library(e1071)
library(pROC)

data1 = data.table(read.csv("C:/Users/Subho/Dropbox/Blood-brain/All descriptors  - BBB.csv"))
data2 = data.table(read.csv("C:/Users/Subho/Dropbox/Blood-brain/BBB set 2.csv"))
data3 = data.table(read.csv("C:/Users/Subho/Dropbox/Blood-brain/BBB415.csv"))

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

View(sort(names1))
View(sort(names2))
sum(names2 %in% names1)

# find closest match
require(e1071)
closest_match = names1
distance = rep(0, length(names1))
for(i in 1:length(names1)){
  dist.i = as.numeric(lapply(names2, function(x)
    hamming.distance(strsplit(names1[i],"")[[1]], strsplit(x,"")[[1]])))
  closest_match[i] = names2[which.min(dist.i)]
  distance[i] = dist.i[which.min(dist.i)]
}

df = data.table::data.table(cbind(names1, closest_match, distance))
df$distance = as.integer(df$distance)
df = df[order(names1)]

View(df[distance>2])

data4 = data.table(merge(data11, data3, by.x="Source.File",by.y="Name")[,-c(1,256,512,513)])
data4[, Class01 := 1]
data4[Class=="n", Class01 := 0]

library(ranger)
set.seed(11212018)
test = sample(1:nrow(data4),50,replace=F)
rfmod = ranger(as.factor(Class01)~.-Class, data=data4[-test,],
               importance="impurity",probability = TRUE, seed=083120181)
preds = predict(rfmod, data4[test,], type="response")$predictions[,2]
roc(data4$Class01[test],preds)

# check variable importances
imps = sort(rfmod$variable.importance, decreasing=TRUE)
imps[1:50]


