function  [lastDesAorta,lastDesAortaLoc,BW,tmpDesAorta] = aorticarchDetection(sliceLoc,desAorta,loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation)
%AORTICARCHDETECTION detection of the aortic arch
%
%   Examples:
%       [lastDesAorta,lastDesAortaLoc,BW,tmpDesAorta] = aorticarchDetection(sliceLoc,desAorta,loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation )

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

loopCond = true;
curLoc = sliceLoc;


stats = regionprops(desAorta,'Area');
firstArea = stats(1).Area;
preArea = firstArea;
preBW = desAorta;

if(firstArea  > 1800)
    areaCutOff = 900;
else
    areaCutOff = 600;
end

loopCounter = 0;

areaList=[];
k = 1;
areaList(k) = preArea;
k = k + 1;

tmpDesAorta = false(size(ctpa));


while(loopCond)
    
    
    %% segement tissiue based on previous seed point
    preLoc = curLoc - loopIncrementVal;
    
    
    BW = grayscaleSegmentation(ctpa(:,:,preLoc),pixelSpacing,seedPoint,deviation,10);
    
    BW = bsxfun(@times, ctpa(:,:,preLoc), cast(BW, class(ctpa(:,:,preLoc))));
    
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
    
    preBWNew = preBW & BW;
    C = bwconncomp( preBWNew);
    stats = regionprops(C,'Area');
    idx =  find([stats.Area]==max([stats.Area]));
    preBWNew = ismember(labelmatrix(C),idx);
    
    img2N = bwlabel(BW,4);
    lblList = unique(img2N(preBWNew == 1));
    lblList(lblList == 0) = [];
    
    if(numel(lblList) == 1 )
        curBW = (img2N == lblList(1));
    else
        curBW = false(size(BW));
        for j=1:numel(lblList)
            curBW = curBW + (img2N == lblList(j));
        end
        
        C = bwconncomp(curBW);
        stats = regionprops(C,'Area');
        idx =  find([stats.Area]==max([stats.Area]));
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
        
        
    end
    
    %% Loop Condition
    
    stats = regionprops(BW,'Area');
    
    %probably desAorta enhaces upper level
    if((firstArea>stats(1).Area+stats(1).Area*0.35))
        areaCutOff = 750;
        
    end
    
    
    if(loopCounter > 1)
        loopCond = false;
        curLoc = preLoc;
        
        tmpDesAorta(:,:,curLoc) = BW;
        break;
    end
    
    if(loopCounter)
        loopCounter = loopCounter + 1;
    else
        if( stats(1).Area - max(areaList)  > areaCutOff && ~(firstArea>stats(1).Area+stats(1).Area*0.35) )
            loopCounter = loopCounter + 1;
            lastDesAorta = preBW;
            lastDesAortaLoc = curLoc;
            if(stats(1).Area > 3000)
                loopCond = false;
                curLoc = preLoc;
                tmpDesAorta(:,:,curLoc) = preBW;
                break;
            end
            
        end
    end
    
 
    
    preArea = stats(1).Area;
    preBW = BW;
    areaList(k) = preArea;
    k = k +1;
    
    
    hu1 = bsxfun(@times, ctpa(:,:,preLoc), cast(BW, class(ctpa(:,:,preLoc))));
    hu = mean2(hu1(hu1~=0));
    
    deviation = setDeviation(hu);
    
    curLoc = preLoc;
    tmpDesAorta(:,:,curLoc) = preBW;
    
end

end % end of function