#********************************Algorithm Function************************************#

###multiscale analysis : get noise and clusters partitioning features from the original data  
getPartitionFeatures <- function(data, epscut, delta, significanceLevel){
  
  i <- 1
  num_of_neighbors_features <- matrix(nrow=0,ncol=0)
  
  set.seed(1)   
  while(epscut <= 0.4*(1.1^(ncol(data)-2))) {
    fr <- frNN(data, epscut)
    num_of_neighbors_0 <- sapply(fr$id,length)
    num_of_neighbors <- dataContinue(num_of_neighbors_0) ## make discrete num of neighbors data continue
    
    num_of_neighbors <- num_of_neighbors[order(num_of_neighbors)]
    p_value <- dip.test(num_of_neighbors)$p.value
    
    if(p_value < significanceLevel){

      ## extract the number of intervals contained in the univariate data and the starting coordinate of the last interval
      ModesInformation <- extractNumofModalIntervals(num_of_neighbors, significanceLevel, debug=FALSE)
      
      ## extract the features of continuous radius(epscut) with the same number of intervals as the first non-unimodal density curve
      if(ModesInformation[2] > 0){
        if(ncol(num_of_neighbors_features) < 1){

          ModesInformationMatrix <- matrix(ModesInformation,1,3) 
          num_of_neighbors_features <- matrix(num_of_neighbors_0,nrow=nrow(data),1)
        } else{
          
          if(ModesInformationMatrix[i-1,1] != ModesInformation[1]){
            if(ncol(num_of_neighbors_features) == 1){
              num_of_neighbors_features <- matrix(nrow=0,ncol=0)
              i <- 0
            } else{
              break
            }
          } else{
            ModesInformationMatrix <- rbind(ModesInformationMatrix, matrix(ModesInformation,1,3))
            num_of_neighbors_features <- cbind(num_of_neighbors_features, matrix(num_of_neighbors_0,nrow=nrow(data),1))
          }
        }
        i <- i + 1
      }     
    } else{
      
      if(ncol(num_of_neighbors_features) > 1){
        break  
      } else if(ncol(num_of_neighbors_features) == 1){
        num_of_neighbors_features <- matrix(nrow=0,ncol=0)
        i <- 1
      }
    }
    
    epscut <- epscut + delta  
  }

  set.seed(NULL) 
  return(list(num_of_neighbors_features, ModesInformationMatrix, epscut)) 
}

###make discrete num of neighbors data continue
dataContinue <- function(data){
  data_raw <- as.data.frame(table(data))
  data_raw[,1] <- as.character(data_raw[,1])
  data_raw[,1] <- as.double(data_raw[,1])
  
  data_contious <- runif(data_raw[1,2], (data_raw[1,1]-(data_raw[2,1]-data_raw[1,1])), data_raw[1,1])
  
  for(j in 2:nrow(data_raw)){
    data_contious <- c(data_contious,runif(data_raw[j,2],data_raw[j-1,1],data_raw[j,1])) 
  }
  
  data <- as.numeric(data_contious)
  return(data)
}

### extract the number of intervals contained in the univariate data and the starting coordinate of the last interval
extractNumofModalIntervals <- function(data,significanceLevel,debug){
    if(is.unsorted(data)){
        stop("data must be sorted")
    }

    ## Find the raw clusters using recursive approach in UNIDIP
    ## Note: here we're saying that we're not testing a modal interval. This means that, if the full distribution is not multimodal, it will only return its estimate of the mode (not the full distribution)
    clustersRaw <- getClustersInInterval(data,1,length(data),"----",FALSE,debug,significanceLevel) 

    ## Consolidation
    clusters <- mergeIntervals(data,clustersRaw, debug)
    
    clusterStarts <- clusters[seq(1,length(clusters),2)]
    clusterEnds <- clusters[seq(2,length(clusters),2)]
    
    numModes <- length(clusterStarts)
    
    clusterIntervalStart <- (data[clusterEnds[numModes-1]] + data[clusterStarts[numModes]])/2
    
    if(data[clusterStarts[1]] <= 1){    
      data[clusterStarts[1]] <- 5
    }
    
    density_ratio = (2*data[clusterStarts[1]])/(data[clusterStarts[numModes]] + data[clusterEnds[numModes]])
    
    return(c(numModes, clusterIntervalStart, density_ratio))
} 

### Find the raw clusters using recursive approach in UNIDIP
### Note that here the indexes are always passed/returned in a global sense
getClustersInInterval <- function(data,indexStart,indexEnd, prefix,testingModalIntervalOnly, debug, significanceLevel){
    
    if(debug){
        print(sprintf("%sChecking interval [%d,%d]",prefix,indexStart,indexEnd))
    }
    ## Subset the data...that is, we want to recursively look at only the data between indexStart and indexEnd and search for modes in that distribution
    dataSubset <- data[indexStart:indexEnd]

    ## Run the dip test on the data subset. Todo: here we run both the dip.test and the dip, which does the procedure twice in total, but a simpler solution doesn't seem easily achievable
    dipTestResult <- dip.test(dataSubset)
    dipResult <- dip(dataSubset,full.result = TRUE)
    modalIntervalLeft <- indexStart + dipResult$lo.hi[1] - 1
    modalIntervalRight <- indexStart + dipResult$lo.hi[2] - 1

    ## Check for non-significance using our significance threshold. If the result is non-significant, we'll assume we only have one cluster (unimodal)
    if(dipTestResult$p.value > significanceLevel){
        if(testingModalIntervalOnly){
            if(debug){
                print(sprintf("%sModal Interval [%d,%d] is unimodally distributed (p-value %f)...returning it as a cluster...",prefix, indexStart,indexEnd, dipTestResult$p.value))
            }            
            return(c(indexStart,indexEnd))
        } else{
            ## Here we know we're finding the "last" cluster. For the unimodal case where the mode is indeed just a point of a small interval, the dip test has the tendency to
            ## return a very small modal interval (makes sense). This is bad in our case because it means that our core cluster is typically going to be very small in relation to
            ## the others that are found. For this reason we need a mechanism for finding out what the "full" core cluster is
            ## Our mechanism: mirror the data and run the dip on it. We're sure that it's then multimodal, and the dip should find "fully" one of the modes as a larger interval
            if(debug){
                print(sprintf("%sInterval [%d,%d] is unimodally distributed (p-value %f)...the modal interval found in was [%d,%d]. Proceeding to mirror data to extract a fatter cluster here...",prefix,indexStart,indexEnd,dipTestResult$p.value, modalIntervalLeft,modalIntervalRight))
            }

            ## Get left and right-mirrored results
            leftMirroredDataSet <- mirrorDataSet(dataSubset, TRUE)
            leftMirroredDataSetResult <- dip(leftMirroredDataSet, full.result=TRUE)
            rightMirroredDataSet <- mirrorDataSet(dataSubset, FALSE)
            rightMirroredDataSetResult <- dip(rightMirroredDataSet, full.result=TRUE)
            if(leftMirroredDataSetResult$dip > rightMirroredDataSetResult$dip){
                clusterRange <- mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData(leftMirroredDataSetResult$lo.hi[1],leftMirroredDataSetResult$lo.hi[2], modalIntervalLeft, modalIntervalRight, length(dataSubset), TRUE, indexStart)
                if(debug){
                    print(sprintf("%sModal interval on the left-mirrored data was [%d,%d]...which corresponds to a cluser (which we'll return now) in the original data of [%d,%d].",prefix,leftMirroredDataSetResult$lo.hi[1],leftMirroredDataSetResult$lo.hi[2],clusterRange[1],clusterRange[2]))
                }
                return(clusterRange)
            } else{
                clusterRange <- mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData(rightMirroredDataSetResult$lo.hi[1],rightMirroredDataSetResult$lo.hi[2], modalIntervalLeft, modalIntervalRight, length(dataSubset), FALSE, indexStart)
                if(debug){
                    print(sprintf("%sModal interval on the right-mirrored data was [%d,%d]...which corresponds to a cluser (which we'll return now) in the original data of [%d,%d].",prefix,rightMirroredDataSetResult$lo.hi[1],rightMirroredDataSetResult$lo.hi[2],clusterRange[1],clusterRange[2]))
                }
                return(clusterRange)
            }            
        }            
    }

    if(debug){
        print(sprintf("%sModal interval [%d,%d], p=%g",prefix,modalIntervalLeft,modalIntervalRight, dipTestResult$p.value))
    }
    
    ## Otherwise, expand the modal interval to see if it has more than one cluster
    modalIntervalClusters <- Recall(data, modalIntervalLeft, modalIntervalRight,paste(prefix,"----",sep=""), TRUE, debug, significanceLevel)

    ## Now we need to look at the various cases.
    ## If we only have a left interval, we just need to proceed in it...there must be at least one cluster there
    ## If we only have a right interval, we just need to proceed in it...there must be at least one cluster there
    ## If we have both, we need to consider both...there COULD be one or more on either side
    leftIntervalExists <- (indexStart < modalIntervalLeft)
    rightIntervalExists <- (indexEnd > modalIntervalRight)
    if(!leftIntervalExists & !rightIntervalExists){
        stop("We found a statistical multimodality, but the modal interval is the full interval! This should never happen!")
    }
    if(!leftIntervalExists & rightIntervalExists){
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] leaves no other points to the left, so we can continue to the right...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        rightClusters <- Recall(data,modalIntervalRight+1,indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,rightClusters))
    }
    if(leftIntervalExists & !rightIntervalExists){
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] leaves no other points to the right, so we can continue to the left...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        leftClusters <- Recall(data,indexStart,modalIntervalLeft-1,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters))
    }

    ## Otherwise, we have the general case of both intervals (left and right)
    ## Here we need to check both, including the closest cluster from the modal interval    
    if(length(modalIntervalClusters)>2){        
        ## More than one cluster in the modal interval, so include the closest for the test of each extreme        
        leftClusters <- Recall(data,indexStart,modalIntervalClusters[2],paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        rightClusters <- Recall(data,modalIntervalClusters[length(modalIntervalClusters)-1],indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters,rightClusters)) 
    } else{
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] is unimodal with intervals left and right of it. Checking these neighbouring intervals with the modal interval...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        ## Single cluster in modal interval. We hence know that there exists cluster(s) outside the modal interval. Find (them) by just focusing on the extreme intervals
        leftClusters <- Recall(data, indexStart, modalIntervalRight,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)    
        rightClusters <- Recall(data,modalIntervalLeft,indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters,rightClusters)) 
    }
}

### Shifts to zero start, then mirrors, then returns the mirrored data
## E.g. input c(2,3,4), output c(-2,-1,0,1,2)
## Assumes an ordered input
mirrorDataSet <- function(data, mirrorLeft){
    mirroredDataSet <- data
    if(mirrorLeft){
        minValue <- min(data)
        dataShifted <- data-minValue
        dataShiftedGtZero <- dataShifted[dataShifted > 0]
        dataShiftedGtZeroMirrored <- -dataShiftedGtZero
        mirroredDataSet <- c(dataShiftedGtZeroMirrored, 0, dataShiftedGtZero)
    } else{
        maxValue <- max(data)
        dataShifted <- data-maxValue
        dataShiftedLtZero <- dataShifted[dataShifted < 0]
        dataShiftedLtZeroMirrored <- -dataShiftedLtZero
        mirroredDataSet <- c(dataShiftedLtZero, 0, dataShiftedLtZeroMirrored)
    }
    return(mirroredDataSet[order(mirroredDataSet)])
}

#mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData
mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData <- function(lowerIndexToMap, upperIndexToMap, lowerFallbackIndex, upperFallbackIndex, lengthOfOriginalData, mirrorLeft, offsetIndex){
    ## Let's say our original data had a length of 2
    ## Then the mirrored data will have a length of 3
    ## The mirrored data will always have an odd length
    ## The zero point will be at lengthOfOriginalData

    if((lowerIndexToMap<lengthOfOriginalData) & (upperIndexToMap>lengthOfOriginalData)){
        return(c(lowerFallbackIndex, upperFallbackIndex))
    }    

    lowerIndexMapped <- mapIndexToOrderedMirroredDataIndexInOriginalOrderedData(lowerIndexToMap, lengthOfOriginalData, mirrorLeft)
    upperIndexMapped <- mapIndexToOrderedMirroredDataIndexInOriginalOrderedData(upperIndexToMap, lengthOfOriginalData, mirrorLeft)
    if(lowerIndexMapped > upperIndexMapped){
        return(c(upperIndexMapped, lowerIndexMapped)+offsetIndex-1)
    }
    return(c(lowerIndexMapped, upperIndexMapped)+offsetIndex-1)    
}
mapIndexToOrderedMirroredDataIndexInOriginalOrderedData <- function(indexToMap, lengthOfOriginalData, mirrorLeft){
    indexMirrored <- indexToMap-lengthOfOriginalData
    if(indexMirrored<0){
        indexMirrored <- -indexMirrored
    }
    if(mirrorLeft){
        if(indexToMap>=lengthOfOriginalData){
            return(indexToMap-lengthOfOriginalData+1)
        }else{
            return(lengthOfOriginalData-indexToMap+1)
        }
        return(indexMirrored+1)
    } else{
        if(indexToMap>lengthOfOriginalData){
            return((2*lengthOfOriginalData) - indexToMap)
        } else{
            return(indexToMap)
        }
    }
}

# merge any overlapped intervals and intervals that aren't statistically-significant enough to cause multimodality
mergeIntervals <- function(orderedData,intervals,debug){
    ## We first need to find the merged clusters (any overlaps are merged), such that we only have mutually-exclusive clusters
    clusterStartingIndices <- intervals[seq(1,length(intervals),2)]
    clusterEndingIndices <- intervals[seq(2,length(intervals),2)]
    clusterStartingIndicesOrderMappings <- order(clusterStartingIndices)
    clusterStartingIndicesOrdered <- clusterStartingIndices[clusterStartingIndicesOrderMappings]
    clusterEndingIndicesOrdered <- clusterEndingIndices[clusterStartingIndicesOrderMappings]
    
    clusters <- c(clusterStartingIndicesOrdered[1],clusterEndingIndicesOrdered[1])
    for(i in tail(seq_along(clusterStartingIndicesOrdered), length(clusterStartingIndicesOrdered)-1)){
        ## Get the current interval
        startInQuestion <- clusters[length(clusters)-1]
        endInQuestion <- clusters[length(clusters)]

        if(endInQuestion < clusterStartingIndicesOrdered[i]){
            clusters <- c(clusters,clusterStartingIndicesOrdered[i],clusterEndingIndicesOrdered[i])
        } else if(endInQuestion < clusterEndingIndicesOrdered[i]){
            clusters[length(clusters)] <- clusterEndingIndicesOrdered[i]
        }
    }

    ## Now we do our "consolidation" step
    ## The idea is that we merge in any "tails", "fringes" or "artefacts" that aren't statistically-significant enough to cause multimodality.
    ## How? We know that our clusters are ordered.
    ## We iterate though our clusters and perform the dip test on the range defined by successfive pairs
    ## If a pair has a non-significant multimodality, we call the entire range defined by that successive pair a single cluster
    consolidatedClusters <- consolidateClusters(orderedData,clusters,1,debug)
    return(consolidatedClusters)
}

##consolidateClusters
consolidateClusters <- function(orderedData, clusters,index,debug){
    ## If index > length-1 done
    ## do dip
    ## If significant
    ## increment index and recurse with index++
    ## else
    ## merge and recurse with index

    numClustersLeft <- length(clusters)/2
    if(index > (numClustersLeft-1)){
        return(clusters)
    }
    startingIndex <- (index*2)-1
    endingIndex <- (index+1)*2
    startingPointIndex <- clusters[startingIndex]
    endingPointIndex <- clusters[endingIndex]
    dipTestResult <- dip.test(orderedData[startingPointIndex:endingPointIndex])
    if(dipTestResult$p.value < 0.05){
        ## significant multimodality...continue with the next index
        if(debug){
            print(sprintf("Range %d to %d is significant...we're happy with that cluster!",startingPointIndex,endingPointIndex))
        }
        return(Recall(orderedData,clusters,index+1, debug))
    } else{
        if(debug){
            print(sprintf("Range %d to %d is not significant: merging",startingPointIndex,endingPointIndex))
        }
        ## not significant...merge and repeat with the same index        
        clusters <- c(clusters[1:startingIndex],clusters[endingIndex:length(clusters)])
        return(Recall(orderedData,clusters,index, debug))        
    }      
}


#********************************Evaluation Index Calculation************************************#

## Calculates AMI for each pair of labels.
calculateMetrics <- function(truelabelslist,predictedlabelslist){
    numCases <- length(truelabelslist)
    outputfilename <- "/home/ztt/Research/DBMAC-master/code/tmp/metrics.csv"
    file.remove(c("/home/ztt/Research/DBMAC-master/code/tmp/truelabels.csv","/home/ztt/Research/DBMAC-master/code/tmp/predictedlabels.csv",outputfilename))
    for(i in 1:numCases){
        write.table(matrix(truelabelslist[[i]],1),"/home/ztt/Research/DBMAC-master/code/tmp/truelabels.csv",append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
        write.table(matrix(predictedlabelslist[[i]],1),"/home/ztt/Research/DBMAC-master/code/tmp/predictedlabels.csv",append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
    }
    
    system(sprintf("/usr/local/MATLAB/R2016b/bin/matlab -nodisplay -r \"cd '/home/ztt/Research/DBMAC-master/code/scripts/'; calcami('/home/ztt/Research/DBMAC-master/code/tmp/truelabels.csv','/home/ztt/Research/DBMAC-master/code/tmp/predictedlabels.csv','%s')\"",outputfilename),ignore.stdout= TRUE)

    res <- readLines(outputfilename)
    numResults <- length(res)
    resultsMatrix <- matrix(replicate(numResults*3,-1),numResults,3)
    for(i in 1:numResults){
        resultsMatrix[i,] <- as.numeric(unlist(strsplit(res[i],split=",")))
    }

    return(resultsMatrix)
}


#********************************Data Generation Function************************************#

###---example data

#generate single-density data (example data for DBMAC)
generateSingleDensityData <- function(noiseFraction0to1=0.8){
  
  if(noiseFraction0to1==0.8){
    set.seed(100);
  } else{
    set.seed(120);
  }
  
  clusterSize <- 400;
  xVals <- c();
  yVals <- c();
  
  c1x <- runif(clusterSize*0.8, 0.84, 0.94); c1y <- runif(clusterSize*0.8, 0.15, 0.23); 
  c2x <- runif(clusterSize*0.8, 0.7, 0.8); c2y <- runif(clusterSize*0.8, 0.15, 0.23);
  c3x <- runif(clusterSize*3.2, 0.1, 0.14); c3y <- runif(clusterSize*3.2, 0.1, 0.9); 
  c4x <- runif(clusterSize*0.5, 0.7, 0.8); c4y <- runif(clusterSize*0.5, 0.85, 0.9);
  
  c5x <- runif(clusterSize*2.25, 0.4, 0.55); c5y <- runif(clusterSize*2.25, 0.4, 0.55); 
  dataMatrix_circle <- matrix(c(c5x, c5y), length(c5x), 2); 
  dataMatrix_circle <- as.matrix(dataMatrix_circle[sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c5x), c(0.475,0.475)), length(c5x), 2, byrow=TRUE))^2))<=0.075,]); 
  c5x <- as.numeric(dataMatrix_circle[,1]); c5y <- as.numeric(dataMatrix_circle[,2]); 
  
  c6x <- runif(clusterSize*0.8, 0.7, 0.8); c6y <- runif(clusterSize*0.8, 0.27, 0.35); 
  c7x <- runif(clusterSize*0.8, 0.14, 0.34); c7y <- runif(clusterSize*0.8, 0.1, 0.14); 
  c8x <- runif(clusterSize*1, 0.8, 0.85); c8y <- runif(clusterSize*1, 0.7, 0.9);
  
  xVals <- c(xVals, c1x,c2x,c3x,c4x, c5x, c6x, c7x, c8x);
  yVals <- c(yVals, c1y,c2y,c3y,c4y, c5y, c6y, c7y, c8y);
  numNonNoisePoints <- length(xVals);
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1));
  xVals <- c(xVals, runif(numNoisePoints,0,1));
  yVals <- c(yVals, runif(numNoisePoints,0,1));
  dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2);
  
  groundTruthLabels <- replicate(length(xVals),0)
  
  groundTruthLabels[xVals>0.84 & xVals<0.94 & yVals>0.15 & yVals<0.23] <- 1
  groundTruthLabels[xVals>0.7 & xVals<0.8 & yVals>0.15 & yVals<0.23] <- 2
  groundTruthLabels[(xVals>0.1 & xVals<0.14 & yVals>0.1 & yVals<0.9) | (xVals>0.14 & xVals<0.34 & yVals>0.1 & yVals<0.14)] <- 3  
  groundTruthLabels[(xVals>0.7 & xVals<0.8 & yVals>0.85 & yVals<0.9) | (xVals>0.8 & xVals<0.85 & yVals>0.7 & yVals<0.9)] <- 4  
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.475,0.475)), length(xVals), 2, byrow=TRUE))^2))<0.075] <- 5
  groundTruthLabels[xVals>0.7 & xVals<0.8 & yVals>0.27 & yVals<0.35] <- 6
  
  dataMatrix <- cbind(dataMatrix, matrix(groundTruthLabels, length(xVals), 1))
  
  ## random permute
  dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
  set.seed(NULL);
  return(dataMatrix);    
}

#generate varying-densities data (D1 and D2 for DBMAC-II)
##D1
generateVaryingDensityData1 <- function(noiseFraction0to1=0.8){
  
  if(noiseFraction0to1==0.8){
    set.seed(190);
  } else{
    set.seed(34);
  }
  clusterSize <- 200;
  xVals <- c();
  yVals <- c();
  
  c1x <- rnorm(clusterSize, 0.9, 0.02); c1y <- rnorm(clusterSize, 0.2, 0.02);
  c2x <- rnorm(clusterSize, 0.8, 0.02); c2y <- rnorm(clusterSize, 0.2, 0.02);
  c3x <- runif(clusterSize*1.9, 0.14, 0.18); c3y <- runif(clusterSize*1.9, 0.1, 0.9); 
  c4x <- runif(clusterSize, 0.75, 0.85); c4y <- runif(clusterSize, 0.86, 0.9); 
  
  c5x <- runif(clusterSize*0.6, 0.45, 0.55); c5y <- runif(clusterSize*0.6, 0.45, 0.55); 
  dataMatrix_circle <- matrix(c(c5x, c5y), length(c5x), 2); 
  dataMatrix_circle <- as.matrix(dataMatrix_circle[sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c5x), c(0.5,0.5)), length(c5x), 2, byrow=TRUE))^2))<=0.05,]); 
  c5x <- as.numeric(dataMatrix_circle[,1]); c5y <- as.numeric(dataMatrix_circle[,2]); 
  
  c6x <- rnorm(clusterSize, 0.8, 0.02); c6y <- rnorm(clusterSize, 0.3, 0.02);
  c7x <- runif(clusterSize*1.4, 0.85, 0.89); c7y <- runif(clusterSize*1.4, 0.75, 0.9);
  
  xVals <- c(xVals, c1x,c2x,c3x,c4x, c5x, c6x, c7x);
  yVals <- c(yVals, c1y,c2y,c3y,c4y, c5y, c6y, c7y);
  numNonNoisePoints <- length(xVals);
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1));
  #xVals <- c(xVals, runif(round(numNoisePoints*0.4),0,0.6), runif(numNoisePoints-round(numNoisePoints*0.4),0.6,1));
  xVals <- c(xVals, runif(numNoisePoints,0,1));
  yVals <- c(yVals, runif(numNoisePoints,0,1));
  dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2);
  
  groundTruthLabels <- replicate(length(xVals),0)
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.9,0.2)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 1
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.8,0.2)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 2
  groundTruthLabels[xVals>0.14 & xVals<0.18 & yVals>0.1 & yVals<0.9] <- 3
  groundTruthLabels[(xVals>0.75 & xVals<0.85 & yVals>0.86 & yVals<0.9) | (xVals>0.85 & xVals<0.89 & yVals>0.75 & yVals<0.9)] <- 4
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.5,0.5)), length(xVals), 2, byrow=TRUE))^2))<0.05] <- 5
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.8,0.3)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 6
  
  dataMatrix <- cbind(dataMatrix, matrix(groundTruthLabels, length(xVals), 1))
  
  ## random permute
  dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
  set.seed(NULL);
  return(dataMatrix);    
}
##D2
generateVaryingDensityData2 <- function(noiseFraction0to1=0.8){

  if(noiseFraction0to1==0.8){
    set.seed(192)
  } else{
    set.seed(20)
  }
  clusterSize <- 200
  xVals <- c()
  yVals <- c()

  c1x <- runif(clusterSize*5, 0.05, 0.35); c1y <- runif(clusterSize*5, 0.7, 0.85);
  dataMatrix_circle <- matrix(c(c1x, c1y), length(c1x), 2);
  dataMatrix_circle <- as.matrix(dataMatrix_circle[sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c1x), c(0.2,0.7)), length(c1x), 2, byrow=TRUE))^2))<0.15 & sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c1x), c(0.2,0.7)), length(c1x), 2, byrow=TRUE))^2))>0.1,]);
  c1x <- as.numeric(dataMatrix_circle[,1]); c1y <- as.numeric(dataMatrix_circle[,2]);
  
  c2x <- runif(clusterSize*2, 0.15, 0.45); c2y <- runif(clusterSize*2, 0.2, 0.4);
  c3x <- runif(clusterSize*2, 0.6, 0.9); c3y <- runif(clusterSize*2, 0.65, 0.85);
  c4x <- runif(clusterSize*1.3, 0.7, 0.9); c4y <- runif(clusterSize*1.3, 0.21, 0.26);
  c5x <- runif(clusterSize*1.3, 0.7, 0.9); c5y <- runif(clusterSize*1.3, 0.29, 0.34);

  xVals <- c(xVals, c1x,c2x,c3x,c4x,c5x)
  yVals <- c(yVals, c1y,c2y,c3y,c4y,c5y)
  numNonNoisePoints <- length(xVals)
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
  xVals <- c(xVals, runif(numNoisePoints,0,1))
  yVals <- c(yVals, runif(numNoisePoints,0,1))
  dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2)

  groundTruthLabels <- replicate(length(xVals),0)
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.2,0.7)), length(xVals), 2, byrow=TRUE))^2))<0.15 & sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.2,0.7)), length(xVals), 2, byrow=TRUE))^2))>0.1 & dataMatrix[,2] >= 0.7] <- 1
  groundTruthLabels[xVals>0.15 & xVals<0.45 & yVals>0.2 & yVals<0.4] <- 2
  groundTruthLabels[xVals>0.6 & xVals<0.9 & yVals>0.65 & yVals<0.85] <- 3
  groundTruthLabels[xVals>0.7 & xVals<0.9 & yVals>0.21 & yVals<0.26] <- 4
  groundTruthLabels[xVals>0.7 & xVals<0.9 & yVals>0.29 & yVals<0.34] <- 5

  dataMatrix <- cbind(dataMatrix, matrix(groundTruthLabels, length(xVals), 1))

  ## random permute
  dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
  set.seed(NULL);
  return(dataMatrix);
}

###---2D synthetic data in experiments 

##DS1
generate2DsyntheticDataSet1 <- function(noiseFraction0to1){
  
  ##cluster centers
  clusterCenters1 <- matrix(c(0.35,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.65,0.5,0.4,0.6,0.35,0.65,0.5,0.4,0.6,0.5),9,2)
  clusterCenters2 <- matrix(c(0.1,0.2,0.4,0.6,0.8,0.9,0.6,0.2,0.9,0.1,0.8,0.4),6,2) 
  clusterCenters <- rbind(clusterCenters1,clusterCenters2)
  
  ##generate data points
  set.seed(1)
  clusterSizeList <- c(replicate(nrow(clusterCenters1),300),replicate(nrow(clusterCenters2),300))
  numNonNoisePoints <- sum(clusterSizeList)
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
  numTotalPoints <- numNonNoisePoints+numNoisePoints
  
  numClusters <- nrow(clusterCenters1) + nrow(clusterCenters2)
  
  ## Begin with a data set of just the noise points
  ## The first "clusteringDimensionality" columns are the features, the last column is the label (in this case, all labeled zero to represent noise)
  data <- matrix(c(runif(numNoisePoints*2),replicate(numNoisePoints,0)), numNoisePoints, 2+1) #noise
  
  ## Standard deviation for our Gaussians
  gaussianSdList <- c(replicate(nrow(clusterCenters1),0.01),replicate(nrow(clusterCenters2),0.03))
  
  for(i in 1:numClusters){
    
    clusterCenter <- as.numeric(clusterCenters[i,])
    
    clusterPoints <- matrix(c(replicate(clusterSizeList[i]*2,0),replicate(clusterSizeList[i],i)),clusterSizeList[i],2+1)
    noisePointsWithinBounds <- replicate(nrow(data),FALSE)
    noisePointsWithinBounds[1:numNoisePoints] <- TRUE
    
    upperBounds <- clusterCenter+2*gaussianSdList[i] ## Assign noise points based on 1sd bound
    lowerBounds <- clusterCenter-2*gaussianSdList[i]
    for(j in 1:2){
      clusterPoints[,j] <- rnorm(clusterSizeList[i], clusterCenter[j], gaussianSdList[i])
    }
    data[sqrt(rowSums((data[,1:2]-matrix(replicate(nrow(data), c(clusterCenters[i,])), nrow(data), 2, byrow=TRUE))^2)) < 2*gaussianSdList[i],2+1] <- i
    
    data <- rbind(data, clusterPoints)
  }
  
  ## random permute
  data <- data[sample(1:nrow(data),nrow(data)),] 

  set.seed(NULL)
  return(data)
}

##DS2
generate2DsyntheticDataSet2 <- function(noiseFraction0to1){
  
  set.seed(100)
  
  clusterSize <- 200
  xVals <- c()
  yVals <- c()
  
  c1x <- runif(clusterSize*500, 0.05, 0.95); c1y <- runif(clusterSize*500, 0.05, 0.95)
  dataMatrix_circle <- matrix(c(c1x, c1y), length(c1x), 2)
  dataMatrix_circle <- as.matrix(dataMatrix_circle[sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c1x), c(0.5,0.5)), length(c1x), 2, byrow=TRUE))^2))<0.45 & sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c1x), c(0.5,0.5)), length(c1x), 2, byrow=TRUE))^2))>0.435,])
  c1x <- as.numeric(dataMatrix_circle[,1]); c1y <- as.numeric(dataMatrix_circle[,2])
  
  c2x <- runif(clusterSize*1, 0.25, 0.35); c2y <- runif(clusterSize*1, 0.65, 0.75)
  c3x <- runif(clusterSize*1, 0.65, 0.75); c3y <- runif(clusterSize*1, 0.65, 0.75)
  
  c4x <- runif(clusterSize*18, 0.2, 0.8); c4y <- runif(clusterSize*18, 0.2, 0.5)
  dataMatrix_circle <- matrix(c(c4x, c4y), length(c4x), 2)
  dataMatrix_circle <- as.matrix(dataMatrix_circle[sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c4x), c(0.5,0.5)), length(c4x), 2, byrow=TRUE))^2))<0.3 & sqrt(rowSums((dataMatrix_circle-matrix(replicate(length(c4x), c(0.5,0.5)), length(c4x), 2, byrow=TRUE))^2))>0.25,])
  c4x <- as.numeric(dataMatrix_circle[,1]); c4y <- as.numeric(dataMatrix_circle[,2])
  
  xVals <- c(xVals, c1x,c2x,c3x,c4x)
  yVals <- c(yVals, c1y,c2y,c3y,c4y)
  numNonNoisePoints <- length(xVals)
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
  xVals <- c(xVals, runif(numNoisePoints,0,1))
  yVals <- c(yVals, runif(numNoisePoints,0,1))
  dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2)
  
  groundTruthLabels <- replicate(length(xVals),0)
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.5,0.5)), length(xVals), 2, byrow=TRUE))^2))<0.45 & sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.5,0.5)), length(xVals), 2, byrow=TRUE))^2))>0.435] <- 1
  groundTruthLabels[xVals>0.25 & xVals<0.35 & yVals>0.65 & yVals<0.75] <- 2
  groundTruthLabels[xVals>0.65 & xVals<0.75 & yVals>0.65 & yVals<0.75] <- 3
  groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.5,0.5)), length(xVals), 2, byrow=TRUE))^2))<0.3 & sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.5,0.5)), length(xVals), 2, byrow=TRUE))^2))>0.25 & dataMatrix[,2]<=0.5] <- 4
  
  dataMatrix <- cbind(dataMatrix, matrix(groundTruthLabels, length(xVals), 1))
  
  ## random permute
  dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
  set.seed(NULL);
  return(dataMatrix);
}

#DS3
generate2DsyntheticDataSet3 <- function(noiseFraction0to1){
  dataFrame <- read.csv("/home/ztt/Research/DBMAC-master/experiments/cure-t2.4k(4.76%).csv", header=FALSE)
  dataFrame[,ncol(dataFrame)] <- as.numeric(dataFrame[,ncol(dataFrame)])
    
  set.seed(10)
  data1 <- dataFrame[dataFrame[,ncol(dataFrame)]==3 | dataFrame[,ncol(dataFrame)]==4,]    
  data1 <- data1[sample(1:nrow(data1),900),]
  data2 <- dataFrame[dataFrame[,ncol(dataFrame)]==1 | dataFrame[,ncol(dataFrame)]==2,]
  data2 <- data2[sample(1:nrow(data2),600),]
  data3 <- dataFrame[dataFrame[,ncol(dataFrame)]==5,]
  data3 <- data3[sample(1:nrow(data3),180),]

  dataFrame <- dataFrame[dataFrame[,ncol(dataFrame)]==0,]
  dataFrame[dataFrame[,ncol(dataFrame)]==0,ncol(dataFrame)] <- 6

  dataFrame <- rbind(dataFrame,data1,data2,data3)
  set.seed(NULL)

  dataFrame <- AddNoisetoClusterData(dataFrame,noiseFraction0to1, -1.2, 1.2, -1.2, 1.2)

  dataFrame[sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(0.92,-0.04)), nrow(dataFrame), 2, byrow=TRUE))^2)) < 0.08,2+1] <- 1
  dataFrame[sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(0.92,-0.56)), nrow(dataFrame), 2, byrow=TRUE))^2)) < 0.08,2+1] <- 2
  dataFrame[sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(-0.3,-0.3)), nrow(dataFrame), 2, byrow=TRUE))^2)) < 0.7,2+1] <- 6
  dataFrame[sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(-0.9,0.8)), nrow(dataFrame), 2, byrow=TRUE))^2)) + sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(-0.2,0.8)), nrow(dataFrame), 2, byrow=TRUE))^2)) < 0.8,2+1] <- 3
  dataFrame[sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(0.9,0.8)), nrow(dataFrame), 2, byrow=TRUE))^2)) + sqrt(rowSums((dataFrame[,1:2]-matrix(replicate(nrow(dataFrame), c(0.2,0.8)), nrow(dataFrame), 2, byrow=TRUE))^2)) < 0.8,2+1] <- 4
  dataFrame[dataFrame[,1] > -0.1 & dataFrame[,1] < 0.1 & dataFrame[,2] > 0.79 & dataFrame[,2] < 0.81,2+1] <- 5

  ## random permute
  dataMatrix <- dataFrame[sample(1:nrow(dataFrame), nrow(dataFrame)),]
  set.seed(NULL);
  return(dataMatrix);

}

AddNoisetoClusterData <- function(dataFrame,noiseFraction0to1, x1, x2, y1, y2){
  if(noiseFraction0to1==0.8){
    set.seed(122);
  } else{
    set.seed(100);  
  }
  
  numNonNoisePoints <- nrow(dataFrame);
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1));
  
  xVals <- c();
  yVals <- c();
  xVals <- c(xVals, runif(numNoisePoints,x1,x2));
  yVals <- c(yVals, runif(numNoisePoints,y1,y2));
  NoisePointsMatrix <- matrix(c(xVals, yVals), length(xVals), 2);
  
  groundTruthLabels <- replicate(numNoisePoints,0);
  NoisePointsMatrix <- cbind(NoisePointsMatrix,matrix(groundTruthLabels, numNoisePoints, 1));
  
  dataMatrix <- rbind(dataFrame, NoisePointsMatrix);
  set.seed(NULL);
  return(dataMatrix);    
}

#DS4
generate2DsyntheticDataSet4 <- function(noiseFraction0to1){
  
  set.seed(100)
  
  clusterSize <- 300
  xVals <- c()
  yVals <- c()
  
  c1x <- runif(clusterSize*3, 0.03, 0.08); c1y <- runif(clusterSize*3, 0.2, 0.8)
  c2x <- runif(clusterSize*0.5, 0.08, 0.18); c2y <- runif(clusterSize*0.5, 0.2, 0.25)
  c3x <- runif(clusterSize*0.5, 0.08, 0.18); c3y <- runif(clusterSize*0.5, 0.75, 0.8)
  
  c4x <- runif(clusterSize*3, 0.23, 0.28); c4y <- runif(clusterSize*3, 0.2, 0.8)
  c5x <- runif(clusterSize*0.4, 0.28, 0.36); c5y <- runif(clusterSize*0.4, 0.2, 0.25)
  
  c6x <- runif(clusterSize*7.8, 0.41, 0.46); c6y <- runif(clusterSize*7.8, 0.2, 0.8)
  c7x <- runif(clusterSize*1.17, 0.46, 0.55); c7y <- runif(clusterSize*1.17, 0.2, 0.25)
  c8x <- runif(clusterSize*7.8, 0.55, 0.6); c8y <- runif(clusterSize*7.8, 0.2, 0.8)
  
  c9x <- runif(clusterSize*0.7, 0.65, 0.8); c9y <- runif(clusterSize*0.7, 0.75, 0.8)
  c10x <- runif(clusterSize*1.375, 0.65, 0.7); c10y <- runif(clusterSize*1.375, 0.475, 0.75)
  c11x <- runif(clusterSize*0.5, 0.7, 0.8); c11y <- runif(clusterSize*0.5, 0.475, 0.525)
  c12x <- runif(clusterSize*1.375, 0.75, 0.8); c12y <- runif(clusterSize*1.375, 0.2, 0.475)
  c13x <- runif(clusterSize*0.5, 0.65, 0.75); c13y <- runif(clusterSize*0.5, 0.2, 0.25)
  
  c14x <- runif(clusterSize*0.85, 0.83, 1); c14y <- runif(clusterSize*0.85, 0.75, 0.8)
  c15x <- runif(clusterSize*2.75, 0.885, 0.935); c15y <- runif(clusterSize*2.75, 0.2, 0.75)
  
  xVals <- c(xVals, c1x,c2x,c3x,c4x,c5x,c6x,c7x,c8x,c9x,c10x,c11x,c12x,c13x,c14x,c15x)
  yVals <- c(yVals, c1y,c2y,c3y,c4y,c5y,c6y,c7y,c8y,c9y,c10y,c11y,c12y,c13y,c14y,c15y)
  numNonNoisePoints <- length(xVals)
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
  xVals <- c(xVals, runif(numNoisePoints,-0.05,1.05))
  yVals <- c(yVals, runif(numNoisePoints,-0.05,1.05))
  dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2)
  
  groundTruthLabels <- replicate(length(xVals),0)
  groundTruthLabels[xVals>0.03 & xVals<0.08 & yVals>0.2 & yVals<0.8] <- 1
  groundTruthLabels[xVals>0.08 & xVals<0.18 & yVals>0.2 & yVals<0.25] <- 1
  groundTruthLabels[xVals>0.08 & xVals<0.18 & yVals>0.75 & yVals<0.8] <- 1
  
  groundTruthLabels[xVals>0.23 & xVals<0.28 & yVals>0.2 & yVals<0.8] <- 2
  groundTruthLabels[xVals>0.28 & xVals<0.36 & yVals>0.2 & yVals<0.25] <- 2
  
  groundTruthLabels[xVals>0.41 & xVals<0.46 & yVals>0.2 & yVals<0.8] <- 3
  groundTruthLabels[xVals>0.46 & xVals<0.55 & yVals>0.2 & yVals<0.25] <- 3
  groundTruthLabels[xVals>0.55 & xVals<0.6 & yVals>0.2 & yVals<0.8] <- 3
  
  groundTruthLabels[xVals>0.65 & xVals<0.8 & yVals>0.75 & yVals<0.8] <- 4
  groundTruthLabels[xVals>0.65 & xVals<0.7 & yVals>0.475 & yVals<0.75] <- 4
  groundTruthLabels[xVals>0.7 & xVals<0.8 & yVals>0.475 & yVals<0.525] <- 4
  groundTruthLabels[xVals>0.75 & xVals<0.8 & yVals>0.2 & yVals<0.475] <- 4
  groundTruthLabels[xVals>0.65 & xVals<0.75 & yVals>0.2 & yVals<0.25] <- 4
  
  groundTruthLabels[xVals>0.83 & xVals<1 & yVals>0.75 & yVals<0.8] <- 5
  groundTruthLabels[xVals>0.885 & xVals<0.935 & yVals>0.2 & yVals<0.75] <- 5
  
  dataMatrix <- cbind(dataMatrix, matrix(groundTruthLabels, length(xVals), 1))
  
  ## random permute
  dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
  set.seed(NULL);
  return(dataMatrix);
}


###---multi-dimensional synthetic data in experiments 
generateMultiDimensionalDataSet <- function(numClusters, clusteringDimensionality, noiseFraction0to1, clusterSizeList){
  
  ##CLUSTER CENTERS
  set.seed(10)
  clusterCenter1 <- matrix(c(0.45,0.5,runif(clusteringDimensionality-2, 0.1, 0.9)),1,clusteringDimensionality)
  set.seed(NULL)
  clusterCenter2 <- clusterCenter1
  clusterCenter2[1,1] <- 0.45+0.08*1.1^0  ##D=2(0),3(0),4(1.7),5(1.9),6(3.1),7(4),8(5)
  clusterCenter3 <- clusterCenter1
  clusterCenter3[1,1] <- 0.45+0.04*1.1^0
  clusterCenter3[1,2] <- 0.5+0.07*1.1^0
  clusterCenter4 <- clusterCenter3
  clusterCenter4[1,1] <- 0.45+0.04*1.1^0 + 0.08*1.1^0
  
  clusterCenters2D <- matrix(c(0.15,0.15,0.3,0.3,0.7,0.7,0.85,0.85,0.7,0.3,0.85,0.15,0.85,0.15,0.7,0.3),8,2)
  
  set.seed(100)
  clusterCenters <- clusterCenters2D[sample(1:nrow(clusterCenters2D),numClusters-4),]
  clusterCenters <- cbind(clusterCenters, matrix(runif((numClusters-4)*(clusteringDimensionality-2), 0.1, 0.9), numClusters-4, clusteringDimensionality-2))
  set.seed(NULL)
  
  clusterCenters <- rbind(clusterCenter1,clusterCenter2,clusterCenter3,clusterCenter4, clusterCenters)
  
  ###GENERATE DATA
  set.seed(1)
  numNonNoisePoints <- sum(clusterSizeList)
  numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
  numTotalPoints <- numNonNoisePoints+numNoisePoints
  
  ## Begin with a data set of just the noise points
  ## The first "clusteringDimensionality" columns are the features, the last column is the label (in this case, all labeled zero to represent noise)
  data <- matrix(c(runif(numNoisePoints*clusteringDimensionality),replicate(numNoisePoints,0)), numNoisePoints, clusteringDimensionality+1) #noise
  
  ## Standard deviation for our Gaussians
  gaussianSdList <- c(0.015,0.015,0.015,0.015,replicate(numClusters-4,0.03))
  
  for(i in 1:numClusters){
    
    clusterCenter <- as.numeric(clusterCenters[i,])
    
    clusterPoints <- matrix(c(replicate(clusterSizeList[i]*clusteringDimensionality,0),replicate(clusterSizeList[i],i)),clusterSizeList[i],clusteringDimensionality+1)
    noisePointsWithinBounds <- replicate(nrow(data),FALSE)
    noisePointsWithinBounds[1:numNoisePoints] <- TRUE
    
    upperBounds <- clusterCenter+2*gaussianSdList[i] ## Assign noise points based on 1sd bound
    lowerBounds <- clusterCenter-2*gaussianSdList[i]
    for(j in 1:clusteringDimensionality){
      clusterPoints[,j] <- rnorm(clusterSizeList[i], clusterCenter[j], gaussianSdList[i])
      noisePointsWithinBounds <- noisePointsWithinBounds & (data[,j]>lowerBounds[j]) & (data[,j]<upperBounds[j])
    }
    
    data[noisePointsWithinBounds,clusteringDimensionality+1] <- i ## Assign noise points falling in cluster range to that cluster
    data <- rbind(data, clusterPoints)
  }
  
  data <- data[sample(1:nrow(data),nrow(data)),] ## random permute
  set.seed(NULL)
  return(data)
}

