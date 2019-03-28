#install.packages("dbscan") #provide "dbscan()" and "hdbscan()" function
library(dbscan) 
#install.packages("diptest") #Compute Hartigan's dip test statistic for unimodality / multimodality and provide a test with simulation based p-values
library(diptest)
#install.packages("devtools") #to help with installing custom R package from github
library(devtools)
#devtools::install_github('gui11aume/unf')  #provide "pfaff()" function to test uniformity for multi-dimensional data
library(unf)      

setwd("/home/ztt/Research/DBMAC-master/code/")
source("func.R")

###Generate example data D1 and D2(with varying-densities clusters)
dataFrame <- generateVaryingDensityData1(0.7) #D1:generateVaryingDensityData1(0.7); D2:generateVaryingDensityData2(0.7)
#write.csv(dataFrame, file = "/home/ztt/Research/DBMAC-master/code/example/DBMAC-II/data/D1.csv",row.names = FALSE)

###generate 2D synthetic data DS1,DS2,DS3,DS4 in experiments
#dataFrame <- generate2DsyntheticDataSet1(0.7)
#dataFrame <- generate2DsyntheticDataSet2(0.6)
#dataFrame <- generate2DsyntheticDataSet3(0.7)
#dataFrame <- generate2DsyntheticDataSet4(0.7)

#generate multi-dimensional synthetic data
#clusterSizeList <- c(600,600,600,600,500,500,500,500,500,500,500)
#dataFrame <- generateMultiDimensionalDataSet(11,3,0.8,clusterSizeList)

dataSet <- as.matrix(dataFrame[,1:(ncol(dataFrame)-1)])
dataSet <- scale(dataSet)
trueLabels <- as.numeric(dataFrame[,ncol(dataFrame)])
plot(dataSet,pch=16,cex=0.5, xlab="", ylab="", xaxt="n", yaxt="n")

#**********************************DBMAC-II*****************************************#
data <- dataSet
num <- nrow(data)
d <- ncol(data)

epscut <- 0.04*(1.1)^(d-2)
delta <- 0.02*(1.1)^(d-2)
significanceLevel <- 0.05

###the matrix for storaging clustering result
clusteringResult <- matrix(replicate((ncol(data)+2),0),nrow=1,ncol=(ncol(data)+2)) 

###determine whether the dataset follows the uniform distribution
### in this case, we employ Friedman-Rafsky's minimal spanning tree (MST) based test for uniformity testing of multi-dimensional data
set.seed(1)
p_value <- pfaff(t(data))
set.seed(NULL)

while(p_value < significanceLevel){
  
  ###get noise and clusters partitioning features from the original data
  partitionFeatures <- getPartitionFeatures(data, epscut, delta, significanceLevel)
  num_of_neighbors_features <- partitionFeatures[1][[1]]
  ModesInformationMatrix <- partitionFeatures[2][[1]]

  ###labels the noisy objects and clustered objects by using k-means clustering
  labels <- replicate(nrow(num_of_neighbors_features),0)
  kmeansData <- as.matrix(scale(num_of_neighbors_features))

  set.seed(1)
  kmeansResult <- kmeans(kmeansData,ModesInformationMatrix[1,1])
  clusterLabel <- which((kmeansResult$center[,1])==max(kmeansResult$center[,1]),arr.ind=TRUE)
  labels[kmeansResult$cluster== clusterLabel[1]] <- 1
  set.seed(NULL)

  clusterData <- data[labels!= 0,]
  noiseData <- data[labels== 0,]
  #plot(clusterData,cex=0.5,pch=16,xlim=c(min(data[,1]),max(data[,1])),ylim=c(min(data[,2]),max(data[,2])), xlab="", ylab="", xaxt="n", yaxt="n")
                                                                              
  ###clustering for clustered objects with DBSCAN(minPts=3)
  ###D1:Eps1=0.035(iteration 1),Eps2=0.060(iteration 2)
  ###D2:Eps1=0.035(iteration 1),Eps2=0.060(iteration 2)
  ###determine parameters according to k-dist graph
  #kNNdistplot(clusterData,k=d+1)
  dbscanResult <- dbscan(clusterData,0.035,d+1)
  predictedLabels <- replicate(nrow(data),0)
  predictedLabels[labels != 0] <- dbscanResult$cluster
  #plot(clusterData,col=dbscanResult$cluster, cex=0.5,pch=16, xlim=c(min(data[,1]),max(data[,1])),ylim=c(min(data[,2]),max(data[,2])), xlab="", ylab="", xaxt="n", yaxt="n")
  
  partitionResult <- cbind(data, matrix(trueLabels,nrow(data),1), matrix(predictedLabels, nrow(data),1))
  #plot(partitionResult[,1:2], col = partitionResult[,4], pch=16, cex=0.5)
  
  ###extract the clusters identified by DBSACN as the part of final resut
  ###the remaining "noise" data will be regarded as the raw data in the next iteration 
  noiseData <- partitionResult[1:num,][partitionResult[1:num,ncol(partitionResult)]==0,]
  dbscanData <- partitionResult[1:num,][partitionResult[1:num,ncol(partitionResult)]!=0,]
  dbscanData[,ncol(dbscanData)] <- dbscanData[,ncol(dbscanData)] + max(clusteringResult[,ncol(clusteringResult)])
  
  clusteringResult <- rbind(clusteringResult,dbscanData)
  num <- nrow(noiseData)
  
  ###uniformity test for remaining "noise" data 
  set.seed(1)
  data <- partitionResult[partitionResult[,ncol(partitionResult)]==0,1:(ncol(partitionResult)-1)]
  gapFillingData <- dbscanData[sample(c(1:nrow(dbscanData)),mean(ModesInformationMatrix[,3])*nrow(dbscanData)),1:(ncol(partitionResult)-1)]
  data <- rbind(data,gapFillingData)
  plot(data[,1:2], pch=16, cex=0.4, xlab="", ylab="", xaxt="n", yaxt="n")
  p_value <- pfaff(t(data[,1:(ncol(data)-1)]))
  set.seed(NULL)
  
  trueLabels <- data[,ncol(data)]
  data <- data[,1:(ncol(data)-1)]

}

clusteringResult <- clusteringResult[-1,]
finalResult <- rbind(clusteringResult, noiseData)

###Calculates AMI with regard to all objects (clustered objects + noisy objects) in data 
trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(finalResult[,ncol(finalResult)-1]))
predictedLabelsList <- c(predictedLabelsList, list(finalResult[,ncol(finalResult)]))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics1)

###alculates AMI with regard to clustered objects in data
trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(finalResult[,ncol(finalResult)-1][finalResult[,ncol(finalResult)-1]!=0]))
predictedLabelsList <- c(predictedLabelsList, list(finalResult[,ncol(finalResult)][finalResult[,ncol(finalResult)-1]!=0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

#write.csv(finalResult, file = "/home/ztt/Research/DBMAC-master/code/example/DBMAC-II/results/D1-DBMAC II(0.035-0.060-0.862-0.922).csv.csv",row.names = FALSE)

finalResult[,ncol(finalResult)][finalResult[,ncol(finalResult)] == 0] <- 8
plot(finalResult[,1:2],cex=0.5,pch=16, col= finalResult[,ncol(finalResult)], xlab="", ylab="", xaxt="n", yaxt="n")


# *********************************SKINNYDIP******************************************#
source("skinnydip.R")
trueLabels <- as.numeric(dataFrame[,ncol(dataFrame)])
predictedLabels <- skinnyDipClusteringFullSpace(dataSet)

SkinnyDip_Result <- cbind(dataFrame, matrix(predictedLabels,nrow(dataFrame),1))

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(SkinnyDip_Result[,ncol(SkinnyDip_Result)-1]))
predictedLabelsList <- c(predictedLabelsList, list(SkinnyDip_Result[,ncol(SkinnyDip_Result)]))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList) 
print(metrics1)

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(SkinnyDip_Result[,ncol(SkinnyDip_Result)-1][SkinnyDip_Result[,ncol(SkinnyDip_Result)-1]!=0]))
predictedLabelsList <- c(predictedLabelsList, list(SkinnyDip_Result[,ncol(SkinnyDip_Result)][SkinnyDip_Result[,ncol(SkinnyDip_Result)-1]!=0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

predictedLabels[predictedLabels==0]<-8
plot(dataSet,cex=0.4, pch=16, col=predictedLabels, xaxt="n", yaxt="n", xlab="", ylab="")


#*********************************DBSCAN***********************************************#
###determine parameters according to k-dist graph
dbscanData <- dataSet; 
#kNNdistplot(dbscanData,k=ncol(dataSet)+1)

dbscanResult <- dbscan(dbscanData,0.035,ncol(dataSet)+1)
predictedLabels <- dbscanResult$cluster

Dbscan_Result <- cbind(dataFrame, matrix(predictedLabels,nrow(dataFrame),1))

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(Dbscan_Result[,ncol(Dbscan_Result)-1]))
predictedLabelsList <- c(predictedLabelsList, list(Dbscan_Result[,ncol(Dbscan_Result)]))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics1)

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(Dbscan_Result[,ncol(Dbscan_Result)-1][Dbscan_Result[,ncol(Dbscan_Result)-1]!=0]))
predictedLabelsList <- c(predictedLabelsList, list(Dbscan_Result[,ncol(Dbscan_Result)][Dbscan_Result[,ncol(Dbscan_Result)-1]!=0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

predictedLabels[predictedLabels==0]<-8
plot(dataSet,cex=0.5,pch=16, col= predictedLabels, xlab="", ylab="", xaxt="n", yaxt="n")

###systematically varied the parameters of Eps and presented the best AMI results
trueLabels <- as.numeric(dataFrame[,ncol(dataFrame)])
d <- ncol(dataSet)

eps_cuts <- seq(0.005,0.5,0.005)
trueLabelsList <- list()
predictedLabelsList <- list()

for(epscut in eps_cuts){
  dbscanResult <- dbscan(dataSet, epscut, d+1)
  predictedLabelsList <- c(predictedLabelsList, list(dbscanResult$cluster))
  trueLabelsList <- c(trueLabelsList, list(trueLabels))
}

metrics <- calculateMetrics(trueLabelsList, predictedLabelsList)
dataSetResults <- cbind(matrix(eps_cuts, ncol=1), metrics)

print(dataSetResults)


# *******************************OPTICS******************************************#
res <- optics(dataSet, eps=1, minPts=15) # we used the same parameters (Eps, minPts) as DBSCAN in experiments

opticsResult <- extractXi(res, xi=0.030)
plot(opticsResult)
predictedLabels <- opticsResult$cluster

optics_Result <- cbind(dataFrame, matrix(predictedLabels,nrow(dataFrame),1))

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(optics_Result[,ncol(optics_Result)-1]))
predictedLabelsList <- c(predictedLabelsList, list(optics_Result[,ncol(optics_Result)]))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics1)

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(optics_Result[,ncol(optics_Result)-1][optics_Result[,ncol(optics_Result)-1]!=0]))
predictedLabelsList <- c(predictedLabelsList, list(optics_Result[,ncol(optics_Result)][optics_Result[,ncol(optics_Result)-1]!=0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

predictedLabels[predictedLabels==1] <- 8
plot(dataSet,col=predictedLabels,pch=16,cex=0.5,xlab="", ylab="", xaxt="n", yaxt="n")


###systematically varied the parameters of steepness threshold and presented the best AMI results
res <- optics(dataSet, eps=1, minPts=15) # we used the same parameters (Eps, minPts) as DBSCAN in experiments
trueLabels <- as.numeric(dataFrame[,ncol(dataFrame)])
d <- ncol(dataSet)

eps_cuts <- seq(0.001,0.1,0.001)
trueLabelsList <- list()
predictedLabelsList <- list()

for(epscut in eps_cuts){
  opticsResult <- extractXi(res, xi=epscut)
  predictedLabelsList <- c(predictedLabelsList, list(opticsResult$cluster))
  trueLabelsList <- c(trueLabelsList, list(trueLabels))
}

metrics <- calculateMetrics(trueLabelsList, predictedLabelsList)
dataSetResults <- cbind(matrix(eps_cuts, ncol=1), metrics)

print(dataSetResults)


# *********************************************HDBSCAN***********************************************************#
hdbscanResult <- hdbscan(dataSet, minPts=22)
predictedLabels <- hdbscanResult$cluster
plot(dataSet, col=predictedLabels, pch=16,cex=0.5)

hdbscan_Result <- cbind(dataFrame, matrix(predictedLabels,nrow(dataFrame),1))

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(hdbscan_Result[,ncol(hdbscan_Result)-1]))
predictedLabelsList <- c(predictedLabelsList, list(hdbscan_Result[,ncol(hdbscan_Result)]))
metrics1 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics1)

trueLabelsList <- list()
predictedLabelsList <- list()
trueLabelsList <- c(trueLabelsList, list(hdbscan_Result[,ncol(hdbscan_Result)-1][hdbscan_Result[,ncol(hdbscan_Result)-1]!=0]))
predictedLabelsList <- c(predictedLabelsList, list(hdbscan_Result[,ncol(hdbscan_Result)][hdbscan_Result[,ncol(hdbscan_Result)-1]!=0]))
metrics2 <- calculateMetrics(trueLabelsList, predictedLabelsList)
print(metrics2)

predictedLabels[predictedLabels==0]<-8
plot(dataSet,cex=0.5,pch=16, col= predictedLabels, xlab="", ylab="", xaxt="n", yaxt="n")


###systematically varied the parameters of minPts and presented the best AMI results
eps_cuts <- seq(10,40,1)

trueLabelsList <- list()
predictedLabelsList <- list()

for(epscut in eps_cuts){
  cl <- hdbscan(dataSet, minPts=epscut)
  predictedLabelsList <- c(predictedLabelsList, list(cl$cluster))
  trueLabelsList <- c(trueLabelsList, list(trueLabels))
}

metrics <- calculateMetrics(trueLabelsList, predictedLabelsList)
dataSetResults <- cbind(matrix(eps_cuts, ncol=1), metrics)

print(dataSetResults)

