function [ascDia,curLoc,preAorta,delBW,seedPoint,deviation,areaList,k,tmpAscA] = ascAortaExtraction(ascAorta,ascLoc,sliceLoc,loopIncrementVal,ctpa,pixelSpacing,deviation)
%ASCAORTAEXTRACTION extraction of the ascending aorta
%
%   Examples:
%       [ascDia,curLoc,preAorta,delBW,seedPoint,deviation,areaList,k,tmpAscA] = ASCAORTAEXTRACTION(ascAorta,ascLoc,sliceLoc,loopIncrementVal,ctpa,pixelSpacing,deviation)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

ascAo = ascAorta;
preAorta = ascAo;

curLoc = ascLoc;

seedPoint = getSeedPoint(ascAo);

%%
estLoc = sliceLoc + 3*loopIncrementVal;

nextLoc = 0;

areaList=[];
k = 1;
stats = regionprops(preAorta,'Area');
areaList(k) = stats(1).Area;
k = k + 1;

delBW = false(size(ascAo));
equivDia=[];

tmpAscA = false(size(ctpa));


while(nextLoc ~= estLoc)%nextLoc ~= estLoc
    
        nextLoc = curLoc + loopIncrementVal;
        
        BW = grayscaleSegmentation(ctpa(:,:,nextLoc),pixelSpacing,seedPoint,deviation,10);
        
        BW = bsxfun(@times, ctpa(:,:,nextLoc), cast(BW, class(ctpa(:,:,nextLoc))));
        
        I2 = imdiffusefilt(BW,'NumberOfIterations',50);
        L = eigenvalHessian(I2);
        E = edge(I2,'canny');
        LE = (L + E) > 0;
        highHUTissiues = ctpa(:,:,nextLoc) > 700;
        highHUTissiues = imdilate(highHUTissiues,strel('disk',3));
        
        
        BW = BW > 0;
        
        BW = (BW - LE - highHUTissiues - delBW) > 0;
        delBW = false(size(ascAo));
        
        BW = imfill(BW,'holes');
        BW = imopen(BW,strel('disk',5));%previosluy 7 however too much for case 1
        
        
        %% select componet which match previous one
        
        %% Below code is equal to binary 2D-Region Growing
        img2N = bwlabel(BW,4);
        lblList = unique(img2N(preAorta == 1));
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
            curBW = BW;
            
        end
        
        stats = regionprops(curBW,'Area');
        
        if(size(areaList,2)>3)
            if(abs(stats(1).Area - mean(areaList)) > 1000)
                
                %% NOTE Test case no 56  re generate seedpoint again !!!!
                delBW = xor(preAorta,curBW);
                delBW = imopen(delBW,strel('disk',5));
                curBW = preAorta & curBW;
                stats = regionprops(curBW,'Area');                           
                
                % New code after testing Delete later if not work!!!!!
                seedPoint = getSeedPoint(curBW);
                %disp("re check whole dataset")
                
            end
        end
        
        
        areaList(k) = stats(1).Area;        
        
        preAorta = curBW;
          
                
        dilatedAorta = imdilate(preAorta,strel('disk',1));
     
                
        measurements = regionprops(dilatedAorta,'EquivDiameter');
        
        equivDia(k-1) = measurements(1).EquivDiameter;
                
        k = k +1;
        
        hu1 = bsxfun(@times, ctpa(:,:,nextLoc), cast(curBW, class(ctpa(:,:,nextLoc))));
        hu = mean2(hu1(hu1~=0));
        
        deviation = setDeviation(hu);
        
        curLoc = nextLoc;
        
        tmpAscA(:,:,curLoc) = curBW;
        

    
end

canNu = ceil(size(equivDia,2)*.25);
ascDiaMM = mean(equivDia((end-canNu+1):end));

ascDia = ascDiaMM*pixelSpacing;



end % end of function