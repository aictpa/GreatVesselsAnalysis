function [tracheaS,sliceNo,seedPointOfTrachea,tracheaCan] = tracheaDetection(ctpa,startIndex,loopIncrementVal)
% TRACHEADETECTION  Starting model pipeline
%
%   Examples:
%      [tracheaS,sliceNo,seedPointOfTrachea,tracheaCan] = TRACHEADETECTION(ctpa,startIndex,loopIncrementVal)
%      

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

tracheaCan =[];

B2 = imresize(ctpa, 0.5);

sizeOfVolume = size(B2,3);

s = startIndex + loopIncrementVal; 
splus = s + round(sizeOfVolume*0.15*loopIncrementVal);

imgVOIA=zeros(size(B2));
imgTemp = [];

parfor scanNo = min([s splus]):max([s splus])
    try
        
        imgThresh = (B2(:,:,scanNo)>-300);
        
        imgCCL = bwconncomp(imgThresh);
        stats = regionprops(imgCCL,'Area','Image');
        idx =  find([stats.Area]==max([stats.Area]));
        imgLung = ismember(labelmatrix(imgCCL), idx);
        
        outerLung = imfill(imgLung,'holes');
        lungArea = outerLung & ~imgThresh;
        imgTemp = bsxfun(@times, B2(:,:,scanNo), cast( lungArea,class(B2(:,:,scanNo))));
        imgVOIA(:,:,scanNo) =  imgTemp;
        
    catch ME
        fprintf(' error in -> %d : %s\n', scanNo, ME);
    end % end of try-catch
    
end% end of main for


%% Rule 1 Filtration
aveDenOfVOI = sum(sum(imgVOIA(:)))/sum(sum(imgVOIA(:)<0));
%enhFactor = aveDenOfVOI*0.1;
allCan = ~(imgVOIA>(aveDenOfVOI));

%% Rule 2 Remove Area under 25 pixel and above 2000 pixel
%allCan = bwareaopen(allCan,25,8);

parfor scanNo = min([s splus]):max([s splus])
    try
        allCan(:,:,scanNo) = bwareafilt(allCan(:,:,scanNo),[10 1200],8);
    catch ME
        fprintf(' error in -> %d : %s\n', scanNo, ME);
    end % end of try-catch
    
end% end of main for


%% Rule 3 Remove volume under formula
imgCanCCL = bwconncomp(allCan);
stats = regionprops(imgCanCCL,'Area');
idxCan =  find( [stats.Area] > ((abs(diff([s splus]))+1) *7 ) );
%imshowV(allCan);

%% Validation




m = round(sizeOfVolume/2);
mplus = m - round(m*0.25*loopIncrementVal);

imgVOIA=zeros(size(B2));
imgTemp = [];


parfor scanNo = min([mplus splus]):max([mplus splus])
    try
        
        imgThresh = (B2(:,:,scanNo)>-300);
        
        imgCCL = bwconncomp(imgThresh);
        stats = regionprops(imgCCL,'Area','Image');
        idx =  find([stats.Area]==max([stats.Area]));
        imgLung = ismember(labelmatrix(imgCCL), idx);
        
        outerLung = imfill(imgLung,'holes');
        lungArea = outerLung & ~imgThresh;
        imgTemp = bsxfun(@times, B2(:,:,scanNo), cast( lungArea, class(B2(:,:,scanNo))));
        imgVOIA(:,:,scanNo) =  imgTemp;
        
        
    catch ME
        fprintf(' error in -> %d : %s\n', scanNo, ME);
    end % end of try-catch
    
end% end of main for

%% Rule 1 Filtration
aveDenOfVOI = sum(sum(imgVOIA(:)))/sum(sum(imgVOIA(:)<0));
%enhFactor = aveDenOfVOI*0.05;
enhFactor = aveDenOfVOI*0.05;
valCan = ~(imgVOIA>(aveDenOfVOI+enhFactor));

%% Rule 2 Remove Area under 25 pixel and above 2000 pixel
%allCan = bwareaopen(allCan,25,8);

parfor scanNo = min([mplus splus]):max([mplus splus])
    try
        valCan(:,:,scanNo) = bwareafilt(valCan(:,:,scanNo),[10 1200],8);
    catch ME
        fprintf(' error in -> %d : %s\n', scanNo, ME);
    end % end of try-catch
    
end% end of main for


%% Rule 3 Remove volume under formula
imgValCCL = bwconncomp(valCan);
stats = regionprops(imgValCCL,'Area');
idxVal =  find( [stats.Area] > ((abs(diff([mplus splus]))+1) *7 ) );


%% Matching


%% Rule 4 Find Scoring (Existence in all possible slices)


if(numel(idxCan.') == 0)
    
    %% We could not found appropriate candidate slices between s -> s+
    
elseif(numel(idxCan.') == 1)
    
    %% Possibly this is the candidate
    
    trachea = ismember(labelmatrix(imgCanCCL), idxCan);
    [~, ~, z] = ind2sub(size(trachea), find(trachea));    
    
    sliceNo = round((min(unique(z))+max(unique(z)))/2);
    
    tracheaS = trachea(:,:,sliceNo);
    
    %% reverse size
    sxy = getSeedPoint(imfill(tracheaS,'holes'));
    
    sxy = sxy/ 0.5;
    
    temp = ~(ctpa(:,:,sliceNo) > -650);
    temp = bwareafilt(temp,[25 2250],8);
    
    stats = regionprops(temp,'Centroid');
    centroidArr = [stats.Centroid];
    centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
    d=dist( sxy ,centroidArr.');
    indexOfCan =  find(d == min(d));
    tracheaS = ismember(labelmatrix(bwconncomp(temp)), indexOfCan);
    seedPointOfTrachea = getSeedPoint(imfill(tracheaS,'holes'));
    
    
else % numel(idxCan.') > 1
    
    imgCanLbl = bwlabeln(ismember(labelmatrix(imgCanCCL), idxCan));
    
    imgValLbl = bwlabeln(ismember(labelmatrix(imgValCCL), idxVal));
    
    % Scoring should be done between s+ -> m+
   
    matchCan = allCan(:,:,splus) & valCan(:,:,splus);
    
    matchAllCan = allCan(:,:,splus) & matchCan;
    temp = imgCanLbl(:,:,splus);
    labelsAllCan =  unique(temp(matchAllCan));
    labelsAllCan(labelsAllCan==0) = [];
    
    lCan = ismember(imgCanLbl, labelsAllCan);
    
    matchValCan = matchCan & valCan(:,:,splus);
    temp = imgValLbl(:,:,splus);
    labelsValCan =  unique(temp(matchValCan));
    labelsValCan(labelsValCan == 0) = [];
    
    lVal = ismember(imgValLbl, labelsValCan);
    
    tracheaCan = (lVal+lCan)>0;
    imgMatchCanCCL = bwconncomp(tracheaCan);
    stats = regionprops( imgMatchCanCCL,'Area');
    
    score = [];
    
    for indxOfCan = 1:size(stats,1)
        tempCan = ismember(labelmatrix(imgMatchCanCCL), indxOfCan);
        [~, ~, z] = ind2sub(size(tempCan), find(tempCan));
        score(indxOfCan) = numel(unique(z))/(abs(diff([s mplus]))+1);
    end% end of main for
    
    indexOfMAxScore = find(score == max(score));
    trachea = ismember(labelmatrix(imgMatchCanCCL), indexOfMAxScore);
    [~, ~, z] = ind2sub(size(trachea), find(trachea));
     
    sliceNo = round((min(unique(z))+max(unique(z)))/2);
    
    tracheaS = trachea(:,:,sliceNo);
    
    %% reverse size
    sxy = getSeedPoint(imfill(tracheaS,'holes'));
    
    sxy = sxy/ 0.5;
    
    temp = ~(ctpa(:,:,sliceNo) > -650);
    temp = bwareafilt(temp,[25 2250],8);
    
    stats = regionprops(temp,'Centroid');
    centroidArr = [stats.Centroid];
    centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
    d=dist( sxy ,centroidArr.');
    indexOfCan =  find(d == min(d));
    tracheaS = ismember(labelmatrix(bwconncomp(temp)), indexOfCan);
    seedPointOfTrachea = getSeedPoint(imfill(tracheaS,'holes'));
    
end % end of if-else


end % end of function