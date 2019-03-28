skinnyDipClusteringFullSpace <- function(data, significanceLevel=0.05,debug=FALSE){
    ## Hypercube-detection
    hypercubes <- findClusterHypercubes(diag(ncol(data)),matrix(nrow=0,ncol=2),data,significanceLevel,debug)
    ## Object-assignment
    labels <- replicate(nrow(data),0)
    for(i in 1:length(hypercubes)){
        hypercube <- hypercubes[[i]]
        ## Get the objects that fall within this hypercube
        objectsInHypercube <- replicate(nrow(data),TRUE)
        for(j in 1:nrow(hypercube)){
            objectsInHypercube <- (objectsInHypercube & data[,j]>=hypercube[j,1] & data[,j]<=hypercube[j,2])
        }
        labels[objectsInHypercube] <- i
    }
    print(sprintf("Found %d clusters (plus noise)", max(labels)))
    return(labels) 
}

## Recursive method for finding the bounding hypercubes for clusters in the provided data, which is assumed to be a matrix where rows are observations and columns are attributes
## Subspace is an orthonormal matrix. The columns represent the directions for our data space. It is in these directions which we take univariate projections.
## ExistingHypercube is an (m x 2) matrix (m rows, representing the dimensions of the subspace. The two entries on each row are the bounds for the edge of the hypercube in the corresponding dimension.
findClusterHypercubes <- function(subspace, existingHypercube, filteredData, significanceLevel, debug){
    subspaceDimension <- ncol(subspace)
    if(nrow(existingHypercube)>=subspaceDimension){
        ## Our hypercube is complete...return it
        #print(existingHypercube)
        return(list(existingHypercube))
    }
    if(nrow(filteredData)==0){
        ## No objects: no cluster
        return(list())
    }
    nextDimension <- nrow(existingHypercube)+1

    ## Get the next direction onto which we'll project our data
    projectionVector <- as.numeric(subspace[,nextDimension])

    ## Project the data onto that direction and sort it
    univariateProjection <- as.numeric(filteredData %*% projectionVector)
    univariateProjectionOrdered <- univariateProjection[order(univariateProjection)]

    ## Get the modal intervals along this direction. We get a matrix back where the rows are the modes, with the start/end values given in the two columns
    ## We always get at least one mode back
    modalIntervals <- extractModalIntervals(univariateProjectionOrdered, significanceLevel, debug)
    #print(modalIntervals)
    numModesFound <- nrow(modalIntervals)

    hypercubes <- list()
    for(i in 1:numModesFound){
        ## Refine our hypercube
        refinedHypercube <- rbind(existingHypercube,as.numeric(modalIntervals[i,]))
        refinedData <- filteredData[filteredData[,nextDimension]>=modalIntervals[i,1] & filteredData[,nextDimension]<=modalIntervals[i,2],]
        clusterHypercubes <- Recall(subspace, refinedHypercube, refinedData, significanceLevel, debug)
        hypercubes <- c(hypercubes,clusterHypercubes)
    }
    return(hypercubes)    
}

## This method is our new method for finding the modes in a 1d (ordered) sample
extractModalIntervals <- function(data,significanceLevel,debug){
    if(is.unsorted(data)){
        stop("data must be sorted")
    }

    ## Find the raw clusters using our recursive approach
    ## Note: here we're saying that we're not testing a modal interval. This means that, if the full distribution is not multimodal, it will only return its estimate of the mode (not the full distribution)
    clustersRaw <- getClustersInInterval(data,1,length(data),"----",FALSE,debug,significanceLevel) 

    ## Consolidation
    clusters <- mergeIntervals(data,clustersRaw, debug)
    
    clusterStarts <- clusters[seq(1,length(clusters),2)]
    clusterEnds <- clusters[seq(2,length(clusters),2)]
    
    return(matrix(c(data[clusterStarts],data[clusterEnds]),nrow=length(clusterStarts),ncol=2))
}

## Note that here the indexes are always passed/returned in a global sense
getClustersInInterval <- function(data,indexStart,indexEnd, prefix,testingModalIntervalOnly, debug, significanceLevel){
    
    if(debug){
        print(sprintf("%sChecking interval [%d,%d]",prefix,indexStart,indexEnd))
    }
    ## Subset the data...that is, we want to recursively look at only the data between indexStart and indexEnd and search for modes in that distribution
    dataSubset <- data[indexStart:indexEnd]

    ## Run the dip test on the data subset. Todo: here we run both the dip.test and the dip, which does the procedure twice in total, but a simpler solution doesn't seem easily achievable
    dipTestResult <- dip.test(dataSubset)
    #dipResult <- dipTestResult$dipFullResult
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

## Shifts to zero start, then mirrors, then returns the mirrored data
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
    ## The idea is that we merge in any "tails", "fringes" or "artefacts" that aren't
    ## statistically-significant enough to cause multimodality.
    ## How? We know that our clusters are ordered.
    ## We iterate though our clusters and perform the dip test on the range defined by successfive pairs
    ## If a pair has a non-significant multimodality, we call the entire range defined by that successive pair a single cluster
    ## print("PreliminaryClusters")
    ## print(clusters)
    consolidatedClusters <- consolidateClusters(orderedData,clusters,1,debug)
    return(consolidatedClusters)
}

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
