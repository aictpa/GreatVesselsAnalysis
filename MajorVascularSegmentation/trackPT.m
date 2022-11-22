function [ptLoc,ptCan,ptCir,ptCirD,ptDia,ptDia3,diaLines,diaLines3,mpa,mpaLoc,isPVDetected,huPT,seedPointPTout,idxOfPre,indexRho] = trackPT(ctpa,pt,nextLoc,loopIncrementVal,pixelSpacing,huPT,preAorta)
%TRACKPT tracking pulmonary trunk for accurate segmentation
%
%   Examples:
%       [ptLoc,ptCan,ptCir,ptCirD,ptDia,ptDia3,diaLines,diaLines3,mpa,mpaLoc,isPVDetected,huPT,seedPointPTout,idxOfPre,indexRho] = trackPT(ctpa,pt,nextLoc,loopIncrementVal,pixelSpacing,huPT,preAorta)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

warning('off','all');

ptTrack = 1;
ptDia = [];
ptDia3 = [];
diaLines = [];
diaLines3 = [];
idxOfPre=[];
indexRho=[];


seedPointPTout = [];

ascAD = bsxfun(@times, ctpa(:,:,nextLoc), cast(preAorta, class(ctpa(:,:,nextLoc))));
ascAD = mean2( ascAD( ascAD~=0));

ptD = huPT;

% step 0 get bottom pixel of ascAorta

C = bwconncomp(preAorta);
stats = regionprops(C,'Extrema');
[~,bottomMidPoint,~,~] = getObjectMidpoints(stats.Extrema);

checkPointC = round(bottomMidPoint(2)) + 10;

% step 1 get middle top pixel

C = bwconncomp(pt);
stats = regionprops(C,'Extrema','Area','Centroid');
[topMidPoint,~,~,~] = getObjectMidpoints(stats.Extrema);

% step 2 get circular roi has radius of 10 pixel from 10 pixel  below the
% top mid point
t = linspace(0, 2*pi, 50);
r = sqrt(1/pi)*10/pixelSpacing;
c = [round(topMidPoint(1)) round(topMidPoint(2))+10];
roi = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512, 512);

% step 3 logical and with pt
trackablePT = roi & pt;

seedPointPT = getSeedPoint(trackablePT);
seedPointPTout = [seedPointPTout; getSeedPoint(trackablePT)];

% step 4

curLoc = nextLoc;
isPVDetected = false;

delPTEdge = false(size(pt));
delHT = false(size(pt));
k=1;
areaList(k) = stats(1).Area;
centerList(k) = round(stats(1).Centroid(2));
k = k + 1;
ptd = 1;

while(true)
    
    nextLoc = curLoc + loopIncrementVal;

    deviationPT = setDeviationPT(huPT);
        
    ptG = smartGrayscaleSegmentation(ctpa(:,:,nextLoc),pixelSpacing,seedPointPT,deviationPT,10);
    ptG = imfill(ptG,'holes');
    
    pt = bsxfun(@times, ctpa(:,:,nextLoc), cast(ptG, class(ctpa(:,:,nextLoc))));
    
    %I2 = imdiffusefilt(pt,'NumberOfIterations',20);
    I2 = int16(anisocpp(double(pt),30.0,2,50,1/7));
    L = eigenvalHessian(I2);
    L = bwareafilt(L,1);
    E = edge(I2,'canny');
    E = bwareafilt(E,1);
    LE = (L + E) > 0;
    
    
    %%
 
    
    pt = (ptG - LE)>0;
    
    
    img2N = bwlabel(pt,8);
    lblList = unique(img2N(trackablePT == 1));
    lblList(lblList == 0) = [];
    pt = (img2N == lblList(1));
    
    pt = imfill(pt,'holes');
    pt = imopen(pt,strel('disk',3));
    pt = imdilate( pt,strel('disk',3));
    pt = bwareaopen(pt,350,8);
    
    
    
    img2N = bwlabel(pt,8);
    lblList = unique(img2N(trackablePT == 1));
    lblList(lblList == 0) = [];
    pt = (img2N == lblList(1));
    
    C = bwconncomp(pt);
    measurementsPT = regionprops(C,'Area','Perimeter','Centroid','Extrema');
     
    aorta = preAorta & pt;
    
    if( sum(aorta(:)>0) > 50 )
        
        delHT = false(size(pt));
        
        if(ptD<150 && ascAD > (ptD + 75) )
            delHT = bsxfun(@times, ctpa(:,:,nextLoc), cast(pt, class(ctpa(:,:,nextLoc))));
            delHT = delHT > (ptD + 75);
            delHT = bwareaopen( delHT,100,8);
            delHT = imdilate(delHT,strel('disk',3));
            
        end % end of if       
        
        % alternative to delPTEdge
        delPTEdge = edge(preAorta,'canny');
        delPTEdge = imdilate(delPTEdge,strel('disk',3));
        
        
        pt = (pt - delPTEdge - delHT)>0;
        pt = imfill(pt,'holes');
        pt = imopen(pt,strel('disk',3));
        pt = bwareaopen(pt,350,8);
        
        % alternative to above
        img2N = bwlabel(pt,8);
        lblList = unique(img2N(trackablePT == 1));
        lblList(lblList == 0) = [];
        pt = (img2N == lblList(1));
        
        C = bwconncomp(pt);
        measurementsPT = regionprops(C,'Area','Perimeter','Centroid','Extrema');
        
    end
    
    %% end
    
    %%
    
    areaList(k) = measurementsPT(1).Area;
    centerList(k) = round(measurementsPT(1).Centroid(2));
    k = k + 1;
    
    
    ptLoc = nextLoc;
    ptCan = pt;
    
    per = measurementsPT.Perimeter;
    are = measurementsPT.Area;      
    cir = round((4*pi*are)/per^2,2);
    
    [cirD,~,~] = checkCircularity(pt);
    
    ptCir = cir;
    ptCirD = cirD;
    
    mpa{ptTrack} = pt;
    mpaLoc(ptTrack) = nextLoc;
    ptTrack = ptTrack + 1;
    
    %%
    
    [dia,dia3,diaL,diaL3,idxOfPre(ptd),indexRho(ptd)] =  findOptimalDiameter(pt);
    diaLines = [diaLines;diaL];
    diaLines3 = [diaLines3;diaL3];
    ptDia(ptd) = round(pixelSpacing*dia);
    ptDia3(ptd) = round(pixelSpacing*dia3);
    
    ptd = ptd + 1;
    
    [~,bottomMidPoint,~,~] = getObjectMidpoints(measurementsPT.Extrema);
    
    bottomPT = round(bottomMidPoint(2));
    
    %% Diameter

    if(checkPointC >= bottomPT)
        isPVDetected = true;
        point2add = getSeedPoint(trackablePT);
        if(~isempty(seedPointPTout))
            
            if(~isequal(seedPointPTout(end,:),point2add))
                seedPointPTout = [seedPointPTout; point2add]; 
            end
        else
            seedPointPTout = [seedPointPTout; point2add]; 
        end      
        break;
    end
    
    %%
    
    
    if(cir <= 0.7)
        
        isPVDetected = false;
        
        % step 1 get middle top pixel
        
        C = bwconncomp(pt);
        stats = regionprops(C,'Extrema');
        [topMidPoint,~,~,~] = getObjectMidpoints(stats.Extrema);
        
        % step 2 get circular roi has radius of 10 pixel from 10 pixel  below the
        % top mid point
        t = linspace(0, 2*pi, 50);
        r = sqrt(1/pi)*10/pixelSpacing;
        c = [round(topMidPoint(1)) round(topMidPoint(2))+10];
        roi = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512, 512);
        
        % step 3 logical and with pt
        trackablePT = roi & pt;
        
        seedPointPT = getSeedPoint(trackablePT);
        
        hu1 = bsxfun(@times, ctpa(:,:,nextLoc), cast(trackablePT, class(ctpa(:,:,nextLoc))));
        huPT = mean2(hu1(hu1~=0));        
        
    elseif((cir > 0.7 && cir<0.75) && cirD > 20)
        
        isPVDetected = false;
        
        % step 1 get middle top pixel
        
        C = bwconncomp(pt);
        stats = regionprops(C,'Extrema');
        [topMidPoint,~,~,~] = getObjectMidpoints(stats.Extrema);
        
        % step 2 get circular roi has radius of 10 pixel from 10 pixel  below the
        % top mid point
        t = linspace(0, 2*pi, 50);
        r = sqrt(1/pi)*10/pixelSpacing;
        c = [round(topMidPoint(1)) round(topMidPoint(2))+10];
        roi = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512, 512);
        
        % step 3 logical and with pt
        trackablePT = roi & pt;
        
        seedPointPT = getSeedPoint(trackablePT);
        
        hu1 = bsxfun(@times, ctpa(:,:,nextLoc), cast(trackablePT, class(ctpa(:,:,nextLoc))));
        huPT = mean2(hu1(hu1~=0));       
        
    else %if(cir >= 0.75)
        
        isPVDetected = true;        
       
        point2add = getSeedPoint(trackablePT);
        if(~isempty(seedPointPTout))
            
            if(~isequal(seedPointPTout(end,:),point2add))
                seedPointPTout = [seedPointPTout; point2add]; 
            end
        else
            seedPointPTout = [seedPointPTout; point2add]; 
        end
        
        break;
        
    end
    
    seedPointPTout = [seedPointPTout; getSeedPoint(trackablePT)];
    
    curLoc = nextLoc;
    
end % end of while


end % end of function