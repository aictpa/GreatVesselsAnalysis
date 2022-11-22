function  [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = checkTrachealIntubation(ctpa,tracheaS,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,iCase)
% CHECKTRACHEALINTUBATION  Scanning a CT scan to detect if the patient has 
% undergone any tracheal intubation
%   Examples:
%      [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = CHECKTRACHEALINTUBATION(ctpa,tracheaS,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,iCase)


%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

 % Step 1 find mid point of trachea
        stats =  regionprops(tracheaS,'boundingbox');
        xywh = stats.BoundingBox;
        mx = round((xywh(1) + xywh(3)/2));
        my = round((xywh(2) + xywh(4)/2));
        
        % Step 2 threshold image
        
        temp = ctpa(:,:,sliceNoOfTrachea) > 300;
        temp = bwareafilt(temp,[15 100],8);
        stats =  regionprops(temp,'centroid');
        
        centroidArr = [stats.Centroid];
        centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
        d = dist( [mx my] ,centroidArr.');
        indexOfCan =  find(d == min(d));
        
        d(d > 10) = [];     
        
        
        if(~isempty(d))
            temp = ismember(labelmatrix(bwconncomp(temp)), indexOfCan);
            
            maskedImg = bsxfun(@times, ctpa(:,:,sliceNoOfTrachea), cast(temp, class(ctpa(:,:,sliceNoOfTrachea))));
            aveHU = sum( maskedImg(:))/sum(sum( maskedImg>0));
            
            if(aveHU > 900)
                [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = getNewTracheaLocation(ctpa,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,iCase);                
            end
            
        end % end of if
        

end % end of function