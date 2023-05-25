function  [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = getNewTracheaLocation(ctpa,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,caseNo)

% GETNEWTRACHEALOCATION  getting trachea location.
%
%   Examples:
%     [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = GETNEWTRACHEALOCATION(ctpa,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,caseNo)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% tracheal intubation control
B2 = ctpa;
tracheaCan = B2(:,:,sliceNoOfTrachea) > -400;
intubationCan = B2(:,:,sliceNoOfTrachea) > 350;
intubationCan = bwareaopen(intubationCan,5,8);

imgTrachea = regiongrowing(tracheaCan,1,[seedPointOfTrachea(2) seedPointOfTrachea(1)]);

C = bwconncomp(imgTrachea);
stats = regionprops(C,'Area');
[~,idx] = find([stats.Area] == max([stats.Area]));

if(stats(idx(1)).Area <200)
    se = strel('disk',5);
else
    se = strel('disk',2);
end % end of if-else

imgTrachea = imdilate(imgTrachea,se);

intubation =  intubationCan  & imgTrachea;

C = bwconncomp(intubation);
stats = regionprops(C,'Area');
[~,idx] = find([stats.Area] == max([stats.Area]));

endLocOFIntubation = -1;


%% Intubation Extraction
if (sum([stats.Area])>10)  
    im = ismember(labelmatrix(C),idx);
    B2 = B2 > 350;
    seedPointOfIntubation = getSeedPoint(imfill(im,'holes'));
    imgIntubation = regiongrowing(B2,1,[seedPointOfIntubation(2) seedPointOfIntubation(1) sliceNoOfTrachea]);
    
    [~,~,v] = ind2sub(size(imgIntubation),find(imgIntubation == 1));
    
    if(PATIENT_POSITION == 1)
        endLocOFIntubation = max(v);
    else
        endLocOFIntubation = min(v);
    end  % end of if-else
    
    display(fprintf('\n[\bCase No %d --> Intubation]\b',caseNo));
end


%% get new trachea location

sliceNoOfTrachea = endLocOFIntubation;

tracheaCan = ctpa(:,:,sliceNoOfTrachea) > -600;

temp = ~ tracheaCan;
temp = bwareafilt(temp,[100 2250],8);
C = bwconncomp(temp);
stats =  regionprops(C,'centroid');

centroidArr = [stats.Centroid];
centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
d = dist(seedPointOfIntubation ,centroidArr.');
indexOfCan =  find(d == min(d));

tracheaS = ismember(labelmatrix(C),indexOfCan);

seedPointOfTrachea = getSeedPoint(imfill(tracheaS,'holes'));

end % end of function











