#install.packages("dbscan") #provide "dbscan()" and "hdbscan()" function
library(dbscan) 
#install.packages("diptest") #Compute Hartigan's dip test statistic for unimodality / multimodality and provide a test with simulation based p-values
library(diptest)

setwd("/home/ztt/Research/DBMAC-master/code/")
source("func.R")

###Generate example data (with single-density clusters)
dataFrame <- generateSingleDensityData(0.7)
#write.csv(dataFrame, file = "/home/ztt/Research/DBMAC-master/code/example/DBMAC/data/exampleData.csv",row.names = FALSE)

dataSet <- as.matrix(dataFrame[,1:(ncol(dataFrame)-1)]) 
dataSet <- scale(dataSet) 
trueLabels <- as.numeric(dataFrame[,ncol(dataFrame)])
plot(dataSet,pch=16,cex=0.5, xlab="", ylab="", xaxt="n", yaxt="n")

#**********************************DBMAC*****************************************#
data <- dataSet
num <- nrow(data)
d <- ncol(data)

epscut <- 0.04*(1.1)^(d-2)
delta <- 0.02*(1.1)^(d-2)
significanceLevel <- 0.05

###get noise and clusters partitioning features from the original data
partitionFeatures <- getPartitionFeatures(data, epscut, delta, significanceLevel)
num_of_neighbors_features <- partitionFeatures[1][[1]]

###labels the noisy objects and clustered objects by using k-means clustering
labels <- replicate(nrow(num_of_neighbors_features),0)
kmeansData <- as.matrix(scale(num_of_neighbors_features))

set.seed(1)
kmeansResult <- kmeans(kmeansData,2)
clusterLabel <- which((kmeansResult$center[,1])==max(kmeansResult$center[,1]),arr.ind=TRUE)
labels[kmeansResult$cluster== clusterLabel[1]] <- 1
set.seed(NULL)

clusterData <- data[labels!= 0,]
noiseData <- data[labels== 0,]
plot(clusterData,cex=0.4,pch=16,xlim=c(min(data[,1]),max(data[,1])),ylim=c(min(data[,2]),max(data[,2])), xlab="", ylab="", xaxt="n", yaxt="n")

###partition performance evaluation --- G_mean (Table 1 in paper)
#P_t = nrow(clusterData[trueLabels[labels != 0] != 0,])
#N_f = nrow(clusterData[trueLabels[labels != 0] == 0,])
#N_t = nrow(noiseData[trueLabels[labels==0] ==0,])
#P_f = nrow(noiseData[trueLabels[labels==0]!=0,])
#G_mean = sqrt(P_t/(P_t + P_f) * N_t/(N_t + N_f))
#print(P_t)
#print(N_f)
#print(N_t)
#print(P_f)
#print(G_mean)

###clustering for clustered objects with DBSCAN(ExampleData:Eps=0.030,minPts=3)
###determine parameters according to k-dist graph
#kNNdistplot(clusterData,k=d+1)
dbscanResult <- dbscan(clusterData,0.030,d+1)
predictedLabels <- replicate(nrow(data),0)
predictedLabels[labels != 0] <- dbscanResult$cluster
plot(clusterData,col=dbscanResult$cluster, cex=0.5,pch=16, xlim=c(min(data[,1]),max(data[,1])),ylim=c(min(data[,2]),max(data[,2])), xlab="", ylab="", xaxt="n", yaxt="n")

#finalResult <- cbind(dataFrame, matrix(predictedLabels,nrow(dataFrame),1))
#write.csv(finalResult, file = "/home/ztt/Research/DBMAC-master/code/example/DBMAC/results/DBMAC(0.030-0.904-0.980).csv",row.names = FALSE)

###Calculates AMI with regard to all objects (clustered objects + noisy objects) in data 
trueLabelsList <- list()
predictedLabelsList <- list()
predictedLabelsList <- c(predictedLabelsList, list(predictedLabels))
trueLabelsList <- c(trueLabelsList, list(trueLabels))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics1)

###alculates AMI with regard to clustered objects in data
trueLabelsList <- list()
predictedLabelsList <- list()
predictedLabelsList <- c(predictedLabelsList, list(predictedLabels[trueLabels != 0]))
trueLabelsList <- c(trueLabelsList, list(trueLabels[trueLabels != 0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

predictedLabels[predictedLabels==0]<-8
plot(dataSet,cex=0.5,pch=16, col= predictedLabels, xlab="", ylab="", xaxt="n", yaxt="n")
