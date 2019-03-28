function calcami(truelabelsfile,predictedlabelsfile,outputFilename)

tlf=fopen(truelabelsfile,'r');
plf=fopen(predictedlabelsfile,'r');
fid=fopen(outputFilename,'w');

tline=fgets(tlf);
pline=fgets(plf);
while ischar(tline)
targetLabels=strread(tline);
predictedLabels=strread(pline);
targetLabels=targetLabels+(1-min(targetLabels)); % offset so that they start at 1 (ANMI calculation complains if there are nonpositive values)
predictedLabels=predictedLabels+(1-min(predictedLabels)); % offset so that they start at 1 (ANMI calculation complains if there are nonpositive values)
[AMI_,NMI,EMI]=AMI(targetLabels,predictedLabels);

numTargetClusters=size(unique(targetLabels),2);
numPredictedClusters=size(unique(predictedLabels), 2);
fprintf(fid,'%d,%d,%f\n',numTargetClusters, numPredictedClusters, AMI_);
fprintf('%f\n',NMI);
fprintf('%f\n',AMI_);
fprintf('%f\n',EMI);

tline=fgets(tlf);
pline=fgets(plf);
end

fclose(fid);
fclose(plf);
fclose(tlf);
exit;
  
end
