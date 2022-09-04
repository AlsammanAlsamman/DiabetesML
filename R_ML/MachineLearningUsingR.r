#!/usr/bin/env Rscript
# studying gene expression in diabetes to find gene markers using machine learning methods

# Clear the environment
rm(list=ls())

# load libraries
library("caret")
library("e1071")
library("randomForest")
library("rpart")
library("rpart.plot")
library(dplyr)
library(preprocessCore)
library(rattle)
library(limma)
library(edgeR)

# install packages
# install.packages("caret")
# install.packages("e1071")
# install.packages("randomForest")
# install.packages("rpart")
# install.packages("rpart.plot")
# install.packages("ROCR")
# install.packages("ROSE")
# install.packages("dplyr")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#BiocManager::install("preprocessCore")
# install.packages("rattle")

path="/home/samman/Documents/publishing/Diabetes_Amira/GEO_DATA"
# load data gene expression data
# load data GSE15932_GPL570.tsv
GSE15932_GPL570 <- read.table(file.path(path,"GSE15932_GPL570.tsv"), header=TRUE, sep="\t", row.names=1)
# load data GSE30208_GPL6102.tsv
GSE30208_GPL6102 <- read.table(file.path(path,"GSE30208_GPL6102.tsv"), header=TRUE, sep="\t", row.names=1)
# load data GSE55098_GPL570.tsv
GSE55098_GPL570 <- read.table(file.path(path,"GSE55098_GPL570.tsv"), header=TRUE, sep="\t", row.names=1)   

# Tagerted gene markers load
TargetMarkers <- read.table(file.path(path,"Target_Markers.tsv"), header=TRUE, sep="\t", row.names=1)

# select only the genes that are in the TargetMarkers
GSE15932_GPL570 <- GSE15932_GPL570[rownames(GSE15932_GPL570) %in% rownames(TargetMarkers),]
GSE30208_GPL6102 <- GSE30208_GPL6102[rownames(GSE30208_GPL6102) %in% rownames(TargetMarkers),]
GSE55098_GPL570 <- GSE55098_GPL570[rownames(GSE55098_GPL570) %in% rownames(TargetMarkers),]
# check number of rows
nrow(GSE15932_GPL570)
nrow(GSE30208_GPL6102)
nrow(GSE55098_GPL570)

# read phenotype data
phenotype <- read.table(file.path(path,"PhenoTypes.txt"), header=TRUE, sep="\t", row.names=1)
# select only the samples that are in the phenotype data
GSE15932_GPL570<-GSE15932_GPL570[,colnames(GSE15932_GPL570) %in% rownames(phenotype)]
GSE30208_GPL6102<-GSE30208_GPL6102[,colnames(GSE30208_GPL6102) %in% rownames(phenotype)]
GSE55098_GPL570<-GSE55098_GPL570[,colnames(GSE55098_GPL570) %in% rownames(phenotype)]

# save to SELECTED folder
write.table(GSE15932_GPL570, file.path(path,"SELECTED/GSE15932_GPL570.tsv"), sep="\t", quote=FALSE)
write.table(GSE30208_GPL6102, file.path(path,"SELECTED/GSE30208_GPL6102.tsv"), sep="\t", quote=FALSE)
write.table(GSE55098_GPL570, file.path(path,"SELECTED/GSE55098_GPL570.tsv"), sep="\t", quote=FALSE)

# Normalization using limma and edgeR
#GSE55098_GPL570
dgList <- DGEList(counts=as.matrix(GSE55098_GPL570))
d0 <- calcNormFactors(dgList, method = "TMM")
GSE55098_GPL570 <- cpm(d0)
#GSE30208_GPL6102
dgList <- DGEList(counts=as.matrix(GSE30208_GPL6102))
d0 <- calcNormFactors(dgList, method = "TMM")
GSE30208_GPL6102 <- cpm(d0)
#GSE15932_GPL570
dgList <- DGEList(counts=as.matrix(GSE15932_GPL570))
d0 <- calcNormFactors(dgList, method = "TMM")
GSE15932_GPL570 <- cpm(d0)

# save to NORMALIZED folder
write.table(GSE15932_GPL570, file.path(path,"NORMALIZED/GSE15932_GPL570.tsv"), sep="\t", quote=FALSE)
write.table(GSE30208_GPL6102, file.path(path,"NORMALIZED/GSE30208_GPL6102.tsv"), sep="\t", quote=FALSE)
write.table(GSE55098_GPL570, file.path(path,"NORMALIZED/GSE55098_GPL570.tsv"), sep="\t", quote=FALSE)


# Merge GSE15932_GPL570 and GSE55098_GPL570 data by rownames
GSE15932_GPL570_GSE55098_GPL570 <- merge(GSE15932_GPL570, GSE55098_GPL570, by=0)
# replace rownames by first column
rownames(GSE15932_GPL570_GSE55098_GPL570) <- GSE15932_GPL570_GSE55098_GPL570[,1]
# remove first column
GSE15932_GPL570_GSE55098_GPL570 <- GSE15932_GPL570_GSE55098_GPL570[,-1]


# Transpose data
GSE15932_GPL570_GSE55098_GPL570 <- t(GSE15932_GPL570_GSE55098_GPL570)
GSE30208_GPL6102 <- t(GSE30208_GPL6102)
# change to data frame
GSE15932_GPL570_GSE55098_GPL570 <- as.data.frame(GSE15932_GPL570_GSE55098_GPL570)
GSE30208_GPL6102 <- as.data.frame(GSE30208_GPL6102)
# save column names
colnamesGSE15932_GPL570_GSE55098_GPL570 <- colnames(GSE15932_GPL570_GSE55098_GPL570)
colnamesGSE30208_GPL6102 <- colnames(GSE30208_GPL6102)
# change column names to M and number
colnames(GSE15932_GPL570_GSE55098_GPL570) <- paste("M",1:ncol(GSE15932_GPL570_GSE55098_GPL570),sep="")
colnames(GSE30208_GPL6102) <- paste("M",1:ncol(GSE30208_GPL6102),sep="")
# Add phenotype column to the data
GSE15932_GPL570_GSE55098_GPL570$PhenoType <- phenotype[rownames(GSE15932_GPL570_GSE55098_GPL570),]
GSE30208_GPL6102$PhenoType <- phenotype[rownames(GSE30208_GPL6102),]
# convert PhenoType to factor
GSE15932_GPL570_GSE55098_GPL570$PhenoType <- as.factor(GSE15932_GPL570_GSE55098_GPL570$PhenoType)
GSE30208_GPL6102$PhenoType <- as.factor(GSE30208_GPL6102$PhenoType)

# Save to MERGED folder
write.table(GSE15932_GPL570_GSE55098_GPL570, file.path(path,"MERGED/GSE15932_GPL570_GSE55098_GPL570.tsv"), sep="\t", quote=FALSE)
write.table(GSE30208_GPL6102, file.path(path,"MERGED/GSE30208_GPL6102.tsv"), sep="\t", quote=FALSE)
# save TargetMarkers to MERGED folder
write.table(TargetMarkers, file.path(path,"MERGED/TargetMarkers.tsv"), sep="\t", quote=FALSE)
#save column names to MERGED folder
write.table(colnamesGSE15932_GPL570_GSE55098_GPL570, file.path(path,"MERGED/colnamesGSE15932_GPL570_GSE55098_GPL570.tsv"), sep="\t", quote=FALSE)
write.table(colnamesGSE30208_GPL6102, file.path(path,"MERGED/colnamesGSE30208_GPL6102.tsv"), sep="\t", quote=FALSE)

# ML using Decision Tree
# GSE15932_GPL570_GSE55098_GPL570
# split into training and testing
set.seed(123)
inTrain <- createDataPartition(y=GSE15932_GPL570_GSE55098_GPL570$PhenoType, p=0.7, list=FALSE)
training <- GSE15932_GPL570_GSE55098_GPL570[inTrain,]
testing <- GSE15932_GPL570_GSE55098_GPL570[-inTrain,]

# Machine learning method Decision tree
set.seed(123)
dt_fit <- train(PhenoType ~ ., 
                data = training, 
                method = "rpart")
# predict on testing data
dt_pred <- predict(dt_fit, testing)

# confusion matrix
# create folder for GSE15932_GPL570_GSE55098_GPL570 using Decision Tree
dir.create(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/DecisionTree"), recursive = TRUE)

# sink to file
sink(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/DecisionTree/ConfusionMatrix.txt"))
confusionMatrix(dt_pred, testing$PhenoType)
# accuracy
mean(dt_pred == testing$PhenoType)
# sensitivity
mean(dt_pred[testing$PhenoType == "Case"] == "Case")
# specificity
mean(dt_pred[testing$PhenoType == "Control"] == "Control")
# plot the model
# sink off
sink()




# save plot as pdf and save to GSE15932_GPL570_GSE55098_GPL570 using Decision Tree folder
pdf(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/DecisionTree/DecisionTree.pdf"))
rpart.plot(dt_fit$finalModel, # middle graph
type=4,
extra=101,
box.palette="GnBu",
branch.lty=3,
shadow.col="gray",
nn=TRUE
)
dev.off()

# important genes
ImportantMarkers <- dt_fit$finalModel$variable.importance
# as data frame
ImportantMarkers <- as.data.frame(ImportantMarkers)
ImportantMarkers
# remove M from the gene names
ImportantMarkers$MarkerN<-gsub("M","",rownames(ImportantMarkers))
# convert to numeric
ImportantMarkers$MarkerName<-as.numeric(ImportantMarkers$MarkerN)
ImportantMarkers$MarkerName<-colnamesGSE15932_GPL570_GSE55098_GPL570[ImportantMarkers$MarkerName]
# get gene names from TargetMarkers
ImportantMarkers$GeneName<-TargetMarkers[ImportantMarkers$MarkerName,]
# save ImportantMarkers to GSE15932_GPL570_GSE55098_GPL570 using Decision Tree folder
write.table(ImportantMarkers, file.path(path,"GSE15932_GPL570_GSE55098_GPL570/DecisionTree/ImportantMarkers.tsv"), sep="\t", quote=FALSE)


# use random forest and repeat the analysis
# create folder for GSE15932_GPL570_GSE55098_GPL570 using Random Forest
dir.create(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/RandomForest"), recursive = TRUE)
set.seed(123)
rf_fit <- train(PhenoType ~ ., 
                data = training, 
                method = "rf")
# predict on testing data
rf_pred <- predict(rf_fit, testing)

# sink to GSE15932_GPL570_GSE55098_GPL570 using Random Forest folder
sink(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/RandomForest/ConfusionMatrix.txt"))
# confusion matrix
confusionMatrix(rf_pred, testing$PhenoType)
# accuracy
mean(rf_pred == testing$PhenoType)
# sensitivity
mean(rf_pred[testing$PhenoType == "Case"] == "Case")
# specificity
mean(rf_pred[testing$PhenoType == "Control"] == "Control")
# sink off
sink()

# plot the model and save to GSE15932_GPL570_GSE55098_GPL570 using Random Forest folder
pdf(file.path(path,"GSE15932_GPL570_GSE55098_GPL570/RandomForest/RandomForest.pdf"))
plot(rf_fit$finalModel)
dev.off()
# plot marker importance
MI<-rf_fit$finalModel$importance
# as data frame
MI<-as.data.frame(MI)
MI$MarkerNumber<-rownames(MI)
# sort by importance
MI<-MI[order(MI$MeanDecreaseGini,decreasing=TRUE),]
MI<-head(MI,10)
# remove M from the gene names
MI$MarkerN<-gsub("M","",MI$MarkerNumber)
# convert to numeric
MI$MarkerName<-as.numeric(MI$MarkerN)
MI$MarkerName<-colnamesGSE15932_GPL570_GSE55098_GPL570[MI$MarkerName]
# get gene names from TargetMarkers
MI$GeneName<-TargetMarkers[MI$MarkerName,]
# save MI to GSE15932_GPL570_GSE55098_GPL570 using Random Forest folder
write.table(MI, file.path(path,"GSE15932_GPL570_GSE55098_GPL570/RandomForest/MarkerImportance.tsv"), sep="\t", quote=FALSE)


# using GSE30208_GPL6102
# using decision tree
# create folder for GSE30208_GPL6102 using Decision Tree
dir.create(file.path(path,"GSE30208_GPL6102/DecisionTree"), recursive = TRUE)
# split into training and testing
set.seed(123)
inTrain <- createDataPartition(y=GSE30208_GPL6102$PhenoType, p=0.7, list=FALSE)
training <- GSE30208_GPL6102[inTrain,]
testing <- GSE30208_GPL6102[-inTrain,]
# Machine learning method Decision tree
set.seed(123)
dt_fit <- train(PhenoType ~ ., 
                data = training, 
                method = "rpart")
# predict on testing data
dt_pred <- predict(dt_fit, testing)
# sink to GSE30208_GPL6102 using Decision Tree folder
sink(file.path(path,"GSE30208_GPL6102/DecisionTree/ConfusionMatrix.txt"))
# confusion matrix
confusionMatrix(dt_pred, testing$PhenoType)
# sink off
sink()

# plot the model to GSE30208_GPL6102 using Decision Tree folder
pdf(file.path(path,"GSE30208_GPL6102/DecisionTree/DecisionTree.pdf"))
rpart.plot(dt_fit$finalModel, # middle graph
type=4,
extra=101,
box.palette="GnBu",
branch.lty=3,
shadow.col="gray",
nn=TRUE
)
dev.off()
# gene importance
ImportantMarkers <- dt_fit$finalModel$variable.importance
# as data frame
ImportantMarkers <- as.data.frame(ImportantMarkers)
ImportantMarkers
# remove M from the gene names
ImportantMarkers$MarkerN<-gsub("M","",rownames(ImportantMarkers))
# convert to numeric
ImportantMarkers$MarkerName<-as.numeric(ImportantMarkers$MarkerN)
ImportantMarkers$MarkerName<-colnamesGSE30208_GPL6102[ImportantMarkers$MarkerName]
# get gene names from TargetMarkers
ImportantMarkers$GeneName<-TargetMarkers[ImportantMarkers$MarkerName,]
# save ImportantMarkers to GSE30208_GPL6102 using Decision Tree folder
write.table(ImportantMarkers, file.path(path,"GSE30208_GPL6102/DecisionTree/ImportantMarkers.tsv"), sep="\t", quote=FALSE)


# using random forest
# create folder for GSE30208_GPL6102 using Random Forest
dir.create(file.path(path,"GSE30208_GPL6102/RandomForest"), recursive = TRUE)
set.seed(123)
rf_fit <- train(PhenoType ~ ., 
                data = training, 
                method = "rf")
# predict on testing data
rf_pred <- predict(rf_fit, testing)
# sink to GSE30208_GPL6102 using Random Forest folder
sink(file.path(path,"GSE30208_GPL6102/RandomForest/ConfusionMatrix.txt"))
# confusion matrix
confusionMatrix(rf_pred, testing$PhenoType)
# sink off
sink()

# important genes
MI<-rf_fit$finalModel$importance
# as data frame
MI<-as.data.frame(MI)
MI$MarkerNumber<-rownames(MI)
# sort by importance
MI<-MI[order(MI$MeanDecreaseGini,decreasing=TRUE),]
MI<-head(MI,10)
# remove M from the gene names
MI$MarkerN<-gsub("M","",MI$MarkerNumber)
# convert to numeric
MI$MarkerName<-as.numeric(MI$MarkerN)
MI$MarkerName<-colnamesGSE30208_GPL6102[MI$MarkerName]
# get gene names from TargetMarkers
MI$GeneName<-TargetMarkers[MI$MarkerName,]

# save MI to GSE30208_GPL6102 using Random Forest folder
write.table(MI, file.path(path,"GSE30208_GPL6102/RandomForest/MarkerImportance.tsv"), sep="\t", quote=FALSE)

# plot model to GSE30208_GPL6102 using Random Forest folder
pdf(file.path(path,"GSE30208_GPL6102/RandomForest/RandomForest.pdf"))
plot(rf_fit$finalModel)
dev.off()






