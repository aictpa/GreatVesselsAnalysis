function [ascAorta, ascLoc, deviation, tmpAorta] = ascAortaSegmentation(BW,lastDesAorta,lastDesAortaLoc,loopIncrementVal,ctpa,pixelSpacing,deviation)
%ASCAORTASEGMENTATION detection of the ascending aorta
%
%   Examples:
%       [ascAorta, ascLoc, deviation, tmpAorta] = ASCAORTASEGMENTATION(BW,lastDesAorta,lastDesAortaLoc,loopIncrementVal,ctpa,pixelSpacing,deviation)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

% Get upper part of Aortic  Arch
% Later on it will be an ascending Aorta
upperAorticArch = (BW-lastDesAorta) > 0;
upperAorticArch = imopen(upperAorticArch,strel('disk',5));
upperAorticArch = bwareaopen(upperAorticArch,50,8);

C = bwconncomp(upperAorticArch);
stats = regionprops(C,'Area');
idx =  find([stats.Area]==max([stats.Area]));
upperAsc = ismember(labelmatrix(C),idx);
seedPoint = getSeedPoint(upperAsc);


temp = false(size(upperAsc));
temp(1:seedPoint(2),:) = upperAsc(1:seedPoint(2),:);
seedPoint = getSeedPoint(temp);

preArc = upperAsc;

sepLoc = lastDesAortaLoc;
curLoc = sepLoc;

tmpAorta = false(size(ctpa));

while(true)
    

        nextLoc = curLoc + loopIncrementVal;
        
        BW = grayscaleSegmentation(ctpa(:,:,nextLoc),pixelSpacing,seedPoint,deviation,10);
        
        BW = bsxfun(@times, ctpa(:,:,nextLoc), cast(BW, class(ctpa(:,:,nextLoc))));
        
        I2 = imdiffusefilt(BW,'NumberOfIterations',40);
        L = eigenvalHessian(I2);
        E = edge(I2,'canny');
        LE = (L + E) > 0;
        
        BW = BW > 0;
        
        BW = (BW - LE) > 0;
        
        BW = imfill(BW,'holes');
        BW = imopen(BW,strel('disk',5));%previosluy 7 however too much for case 1
        
        
        %% select componet which match previous one
        
        %% Below code is equal to binary 2D-Region Growing
        img2N = bwlabel(BW,4);
        lblList = unique(img2N(preArc == 1));
        lblList(lblList == 0) = [];
        
        if(numel(lblList) == 1 )
            curBW = (img2N == lblList(1));
        else
            curBW = false(size(BW));
            for j=1:numel(lblList)
                curBW = curBW + (img2N == lblList(j));
            end
            
            % Take upper one not maximun
            C = bwconncomp(curBW);
            stats = regionprops(C,'Centroid');
            centroidArr = [stats.Centroid];
            centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
            idx =  find(centroidArr(:,2) == min(centroidArr(:,2)));
            curBW = ismember(labelmatrix(C),idx);
            
        end
        
        if(max(curBW(:)) == 1)
            seedPoint = getSeedPoint(curBW);
            BW = curBW;
        else
            %RESELECT SEEDPOINT
            statsTra = regionprops(BW,'Centroid');
            centroidArr = [statsTra.Centroid];
            centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
            points = centroidArr;
            distToObject = dist(points,[seedPoint(1) seedPoint(2)].');
            ind = find(distToObject == min(distToObject));
            seedPoint = round(points(ind,:));
            
            BW = regiongrowing(BW,1,[seedPoint(2) seedPoint(1)]);
            disp(iCase + " -> BW ? curBW ");
        end
        
        
        
        if(nextLoc ==   sepLoc + 2*loopIncrementVal)
            tmpAorta(:,:,nextLoc) = BW;
            ascLoc = nextLoc;
            ascAorta = BW;
            break;
        end
        
        
        BW = curBW & preArc;
        seedPoint = getSeedPoint(BW);
        preArc = false(size(BW));
        preArc(1:seedPoint(2),:) = BW(1:seedPoint(2),:);
        seedPoint = getSeedPoint(preArc);
        
        
        hu1 = bsxfun(@times, ctpa(:,:,nextLoc), cast(preArc, class(ctpa(:,:,nextLoc))));
        hu = mean2(hu1(hu1~=0));
        
        deviation = setDeviation(hu);
        
        curLoc = nextLoc;
        
        tmpAorta(:,:,nextLoc) = BW;
        

    
end


end % end of function