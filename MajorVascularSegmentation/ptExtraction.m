function [ptLoc,ptCan,ptCir,ptCirD,isPVDetected,nextLoc,hu,huN,preAorta,seedPointPTOut,tmpAo,tmpPT, tmpSearchingArea,huMask,huMaskLoc] = ptExtraction(PATIENT_POSITION,sliceLoc,curLoc,loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation,preAorta,delBW,desAorta,areaList,k)
%PTEXTRACTION extraction of the pulmonary trunk
%
%   Examples:
%       [ptLoc,ptCan,ptCir,ptCirD,isPVDetected,nextLoc,hu,huN,preAorta,seedPointPTOut,tmpAo,tmpPT, tmpSearchingArea,huMask,huMaskLoc] = PTEXTRACTION(PATIENT_POSITION,sliceLoc,curLoc,loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation,preAorta,delBW,desAorta,areaList,k)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


isPVDetected = false;
ptSearch = false;
ptDetected = false;
huPTList = [];
huPTListS = [];
huMask = [];
huMaskLoc = [];
l = 1;
s = 1;
seedPointPTOut = [];
seedPointPT = [];

tmpPT = false(size(ctpa));

tmpAo = false(size(ctpa));

tmpSearchingArea = false(size(ctpa));

while(true)%nextLoc ~= estLoc
    
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
    delBW = false(size(BW));
    
    BW = imfill(BW,'holes');
    BW = imopen(BW,strel('disk',5));%previosluy 7 however too much for case 1
    
    
    %% select componet which match previous one   
   
    
    temp = preAorta & BW;
    C = bwconncomp(temp);
    stats = regionprops(C,'Area');
    idx =  find([stats.Area] == max([stats.Area]));
    temp = ismember(labelmatrix(C),idx);
    
    
    img2N = bwlabel(BW,4);
    lblList = unique(img2N(temp == 1));
    lblList(lblList == 0) = [];
    
    
    if(isempty(lblList))           
        tmpPT(:,:,nextLoc) = pt;
        break;
    end
    
    if(numel(lblList) == 1 )
        curBW = (img2N == lblList(1));
    else
        curBW = false(size(BW));
        for j=1:numel(lblList)
            curBW = curBW + (img2N == lblList(j));
        end
        
        
        % Take upper one not maximun area
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
            
            delBW = xor(preAorta,curBW);
            addBW = preAorta & delBW;
            delBW = imopen(delBW,strel('disk',5));
            curBW = preAorta & curBW | addBW;
            stats = regionprops(curBW,'Area');
            
        end
    end
        
    areaList(k) = stats(1).Area;
            
    preAorta = curBW;
    
    if(~ptDetected)
        if(PATIENT_POSITION == 1 )
            if(nextLoc >= sliceLoc)
                ptSearch = true;
            end
        else
            if(nextLoc <= sliceLoc)
                ptSearch = true;
            end
        end
    end
    
    if(ptSearch)
        
        stats = regionprops(curBW,'Extrema');
        [topMidPoint,bottomMidPoint,rightMidPoint,~] = getObjectMidpoints(stats.Extrema);
        
        rightMidPoint =  round(rightMidPoint);
        topMidPoint =  round(topMidPoint);
        bottomMidPoint =  round(bottomMidPoint);
        
        pt2 = ctpa(:,:,nextLoc);
        pt = int16(zeros(size(pt2)));
        
        pt(topMidPoint(2)-15:bottomMidPoint(2),(rightMidPoint(1)+5:rightMidPoint(1)+70)) = pt2(topMidPoint(2)-15:bottomMidPoint(2),(rightMidPoint(1)+5:rightMidPoint(1)+70));
        
        tmpL=pt;
        tmpL(tmpL==0)=-1100;
        tmpL=tmpL>-1100;
        %
        %          figure,
        %         RGB = labeloverlay(mat2gray(ctpa(:,:,nextLoc)),tmpL,'Colormap',[1 0 0],'Transparency',0.60);
        %         imshow(RGB)
        
        tmpSearchingArea(:,:,nextLoc) =  tmpL;
        
        
        if(numel(huPTList)>2)
            
            pt = pt >  mean(huPTList) - deviationPT; % ?!?!?
        else
            
            ptTemp = pt > 45;
            
            ptTemp = bsxfun(@times, ctpa(:,:,nextLoc), cast(ptTemp, class(ctpa(:,:,nextLoc))));
            huPtTemp = mean2(ptTemp(ptTemp~=0));
            
            if(huPtTemp < 125)
                
                pt = pt >  45; % ?!?!?
                
            else
                
                pt = pt >  huPtTemp - 50; % ?!?!?
            end
        end
        
        
        pt = bsxfun(@times, ctpa(:,:,nextLoc), cast(pt, class(ctpa(:,:,nextLoc))));
        
        I2 = imdiffusefilt(pt,'NumberOfIterations',10);
        L = eigenvalHessian(I2);
        E = edge(I2,'canny');
        LE = (L + E) > 0;
        pt = ((pt>0) - LE)>0;
        
        pt = imfill(pt,'holes');
        pt = imopen(pt,strel('disk',3));
        
        C = bwconncomp(pt);
        stats = regionprops(C,'Area');
        idx =  find([stats.Area]==max([stats.Area]));
        pt = ismember(labelmatrix(C),idx);
        
        
        % For pt density
        ptTemp = bsxfun(@times, ctpa(:,:,nextLoc), cast(pt, class(ctpa(:,:,nextLoc))));
        huMask{s} = pt;
        huMaskLoc(s) = nextLoc;
        huPTListS(s) = mean2(ptTemp(ptTemp~=0));
        s = s + 1;
        
        
        % NOTE: TEST BELOW CODE!
        if(isempty(idx))
            curLoc = nextLoc;
             tmpPT(:,:,nextLoc) = ptTemp;
            continue
        end
             
        
        if(stats(idx).Area > 350)
            
            seedPointPT = round(getSeedPoint(pt));            
            
            huPT = bsxfun(@times, ctpa(:,:,nextLoc), cast(pt, class(ctpa(:,:,nextLoc))));
            huPT = mean2(huPT(huPT~=0));
            huPTList(l) = huPT;
            l = l + 1;
            
            deviationPT = setDeviationPT(huPT);
            
            
            ptG = grayscaleSegmentation(ctpa(:,:,nextLoc),pixelSpacing,seedPointPT,deviationPT,10);
            
            pt = bsxfun(@times, ctpa(:,:,nextLoc), cast(ptG, class(ctpa(:,:,nextLoc))));
            
            I2 = imdiffusefilt(pt,'NumberOfIterations',20);
            L = eigenvalHessian(I2);
            E = edge(I2,'canny');
            LE = (L + E) > 0;
            
            pt = (ptG - preAorta - desAorta -LE)>0;
            pt = imfill(pt,'holes');
            pt = imopen(pt,strel('disk',3));
            pt = bwareaopen(pt,350,8);
            
            
            C = bwconncomp(pt);
            statsPT = regionprops(C,'Centroid');
            centroidArr = [statsPT.Centroid];
            centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
            points = centroidArr;
            distToObject = dist(points,[seedPointPT(1) seedPointPT(2)].');
            ind = find(distToObject == min(distToObject));
            pt = ismember(labelmatrix(C),ind);
            
            
            C = bwconncomp(pt);
            measurementsPT = regionprops(C,'Area','Perimeter','Centroid');
            
            seedPointPTOut = seedPointPT; 
            
            if(~isempty(measurementsPT))
                if(measurementsPT.Area > 500)
                    
                    ptDetected = true;
                    ptSearch = false;
                    ptLoc = nextLoc;
                    ptCan = pt;
                    
                    per = measurementsPT.Perimeter;
                    are = measurementsPT.Area;                  
                    cir = (4*pi*are)/per^2;
                    
                    [cirD,~,~] = checkCircularity(pt);
                    
                    ptCir = cir;
                    ptCirD = cirD;
                    
                    % NOTE: remove later
                                        
                    if(cir <= 0.7)
                        
                        
                        tmpPT(:,:,nextLoc) = pt;
                        
                        isPVDetected = false;
                        
                        break;
                        
                    elseif((cir > 0.7 && cir<0.75) && cirD > 20)
                        
                        tmpPT(:,:,nextLoc) = pt;
                        isPVDetected = false;
                        
                        break;
                        
                    else %if(cir >= 0.75)    
                        tmpPT(:,:,nextLoc) = pt;
                        isPVDetected = true;
                        break;
                    end
                end
            end
        end
    end
    
    
    k = k +1;
    
    hu1 = bsxfun(@times, ctpa(:,:,nextLoc), cast(curBW, class(ctpa(:,:,nextLoc))));
    hu = mean2(hu1(hu1~=0));
    
    deviation = setDeviation(hu);
    
    
    %%    
    curLoc = nextLoc;
    
    tmpPT(:,:,nextLoc) = pt;
    tmpAo(:,:,nextLoc) = curBW;
    
end

if(numel(huPTList)>= numel(huPTListS))
    huN = numel(huPTList);
    hu = round(mean(huPTList));
    
else
    huN = numel(huPTListS);
    hu = round(mean(huPTListS,'omitnan'));
end


end % end of function