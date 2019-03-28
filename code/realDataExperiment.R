#install.packages("dbscan") #provide "dbscan()" and "hdbscan()" function
library(dbscan) 
#install.packages("diptest") #Compute Hartigan's dip test statistic for unimodality / multimodality and provide a test with simulation based p-values
library(diptest)
#install.packages("png") #provide "readPNG()" function to generate a matrix for a raster image from a PNG image
library(png)

## real data
exData <- read.csv("/home/ztt/Research/DBMAC-master/experiments/real-world data/spatialnetwork.data", header=FALSE) 

set.seed(63)
subsetIndices <- sample(1:nrow(exData),12000) 
exDataSubset <- exData[subsetIndices,] 
testDataMatrix <- as.matrix(exDataSubset[,c(2:3)]) 
dataFrame <- testDataMatrix

dataSet <- as.matrix(dataFrame[,1:ncol(dataFrame)])
dataSet <- scale(dataSet)
plot(dataSet,pch=16,cex=0.5, xlab="", ylab="", xaxt="n", yaxt="n")

#**********************************DBMAC-II*****************************************#
data <- dataSet
num <- nrow(data)
d <- ncol(data)

epscut <- 0.04*(1.1)^(d-2)
delta <- 0.02*(1.1)^(d-2)
significanceLevel <- 0.05

###the matrix for storaging clustering result
clusteringResult <- matrix(replicate((ncol(data)+1),0),nrow=1,ncol=(ncol(data)+1)) 

###identify the clusters with the largest and the second largest densities
iteration_max = 2
i = 1
while( i <= iteration_max){
  
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
  ###Eps1=0.040(iteration 1),Eps2=0.080(iteration 2)
  ###determine parameters according to k-dist graph
  #kNNdistplot(clusterData,k=d+1)
  dbscanResult <- dbscan(clusterData,0.040,d+1)
  predictedLabels <- replicate(nrow(data),0)
  predictedLabels[labels != 0] <- dbscanResult$cluster
  plot(clusterData,col=dbscanResult$cluster, cex=0.5,pch=16, xlim=c(min(data[,1]),max(data[,1])),ylim=c(min(data[,2]),max(data[,2])), xlab="", ylab="", xaxt="n", yaxt="n")
  
  partitionResult <- cbind(data, matrix(predictedLabels, nrow(data),1))
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
  #p_value <- pfaff(t(data[,1:(ncol(data)-1)]))
  set.seed(NULL)
  
  data <- data[,1:(ncol(data))]

  i = i + 1

}

clusteringResult <- clusteringResult[-1,]
finalResult <- rbind(clusteringResult, noiseData)

lim <- par() 
ima <- readPNG("/home/ztt/Research/DBMAC-master/experiments/real-world data/jutland.png") 
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
finalResult[finalResult[,3]==0,3]<-8
points(finalResult[,1],finalResult[,2],col=finalResult[,3], xaxt="n",yaxt="n",xlab="",ylab="",pch=20,cex=0.8,cex.main=1.9,main="Skinny-dip clustering") 

text(0.9,-0.3,"Aalborg",cex=1.5)
text(1.9,1.2,"Frederikshavn",cex=1.2) 
text(0.5,1.05,"Hjørring",cex=1.2) 
text(-0.55,0.7,"Blokhus",cex=1.2)

