function [img,seedPointOfOrgan,tolerance,isSuccess] = desAortaExtraction(I,I1,L12,pixelSpacing,rightCentroid,leftCentroid, addedAngle)
%DESAORTAEXTRACTION detect descending aorta
%
%   Examples:
%       [img,seedPointOfOrgan,tolerance,isSuccess] = DESAORTAEXTRACTION(I,I1,L12,pixelSpacing,rightCentroid,leftCentroid, addedAngle)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

isSuccess = false;
isOrganMissed = true;
falsePoints = [];
is2ndPointUsed = false;
is2ndSpaceUsed = false;
isMinThresValReached = false;

meanDensity = mean2(I1(I1 > 0));

if(meanDensity > 160)
    meanDensity = meanDensity - 60;   
elseif(meanDensity > 130)
    meanDensity = meanDensity - 30;   
elseif(meanDensity > 100)
    meanDensity = meanDensity - 15;
end

meanDensityM = meanDensity;
 

[xPoints,yPoints] = desAortaDetection(I,L12,leftCentroid, addedAngle,1);

xPoint = xPoints(1);
yPoint = yPoints(1);


while(isOrganMissed)
    
     
    if(meanDensity < 50)
        if(~isMinThresValReached)
            meanDensity = 35;
            isMinThresValReached = true;
        end
    end
    
    if(meanDensity < 35)        
        
        if(is2ndPointUsed)
            if(is2ndSpaceUsed)
                break;
            end
           
            [xPoints,yPoints] = desAortaDetection(I,L12,leftCentroid, addedAngle,2);
            xPoint = xPoints(1);
            yPoint = yPoints(1);
            is2ndPointUsed = false;
            is2ndSpaceUsed = true;
            
            isMinThresValReached = false;
            meanDensity = meanDensityM;
            continue;
        end % end of if
       
        xPoint = xPoints(2);
        yPoint = yPoints(2);
        is2ndPointUsed = true;
        isMinThresValReached = false;
        
        meanDensity = meanDensityM;
        
    end % end of if
    
    I0 = I1 > meanDensity;
    highDensityTissues =  I>750;
    img = (I0 - L12 - highDensityTissues) > 0;
    
    %% Delete Later or Prove it ---%
    I0 = I1 > 25;
    highDensityTissues =  I>750;
    imgLow = (I0 - L12 - highDensityTissues) > 0;
    
    imgM = xor(imgLow,img);
    imgM = bwareafilt(imgM,[200 1500], 8);
    
    img = (img + imgM)>0;
    
    %%     
    
    C = bwconncomp(img);
    statsTra = regionprops(C,'Area','Centroid');
    [~,indexObj] = find([statsTra.Area] > 200);
    centroidArr = [statsTra(indexObj).Centroid];
    centroidArr = reshape(centroidArr,[2 size(centroidArr,2)/2])';
    %[indexObj,~] = find(centroidArr(:,2) > leftCentroid(2));
    [indexObj,~] = find(centroidArr(:,1) > rightCentroid(1) & ((centroidArr(:,2) - leftCentroid(2)) > 15));
    
    points = centroidArr(indexObj,:);
    falseLocations = [];
    for p = 1:size(points,1)
        if (L12(round(points(p,2)),round(points(p,1))) == 1)
            falseLocations =  [falseLocations p];
        end % end if
    end % end for
    
    points(falseLocations,:) = [];
    
      
      
   if(numel(xPoint) == 1)
       distToObject = dist(points,[xPoint;yPoint]);
   else
       distToObject = dist(points,[xPoint;yPoint]);
   end
   
   
   %% Delete Later or Prove it ---%
   
   dist1 = min(distToObject);
   if( dist1 > 12.5 ) 
       
   xPoint = xPoints(2);
   yPoint = yPoints(2);
   
   if(numel(xPoint) == 1)
       distToObject = dist(points,[xPoint;yPoint]);
   else
       distToObject = dist(points,[xPoint;yPoint]);
   end
   
   dist2 = min(distToObject);
   if (dist1 <= dist2)
       xPoint = xPoints(1);
       yPoint = yPoints(1);
       if(numel(xPoint) == 1)
           distToObject = dist(points,[xPoint;yPoint]);
       else
           distToObject = dist(points,[xPoint;yPoint]);
       end
   end
   
   end
   %% END ------------------------%
   
   
    [indexObj,~] = find( distToObject < 25);
    allIndex =  ~ismember(1:size(points,1),indexObj);
    [~,indexToDel] = find( allIndex == 1);
    points(indexToDel,:) = [];
 
    
    if(~isempty(indexObj))
        
        if(~isempty(falsePoints))
            % end of if
            distToFalsePoints = dist(points,falsePoints');
            [indexFPObj,~] = find(distToFalsePoints < 13);
            points(indexFPObj,:) = [];                   
        end
         
        cirFactors = [];
        tolerances = [];
        seedPoints = [];
        features = [];
        for d = 1:size(points,1)            
            
            seedPointOfOrgan = points(d,:);
            
            [img, out, cirFactor, tolerance, feature] = checkFeatures(I,L12,pixelSpacing,seedPointOfOrgan);
            
            if(out == 1)
                
                 cirFactors = [cirFactors cirFactor];
                 tolerances = [tolerances tolerance];
                 seedPoints = [seedPoints;seedPointOfOrgan];
                 features = [features;feature];
                 isOrganMissed = false;
                 
            else
                falsePoints = [falsePoints;seedPointOfOrgan];
            end % end of if-else         
            
        end % end of for
        
        if(isOrganMissed)
            % desA missed
            meanDensity = meanDensity - 30;
        else
            
            [~,indexMax] = find(cirFactors == max(cirFactors));
            
            seedPointOfOrgan = round(seedPoints(indexMax,:));
            out = grayscaleSegmentation(I,pixelSpacing,seedPointOfOrgan,tolerances(indexMax),10);
            img = (out - L12)>0;
            img = imfill(img,'holes');            
            img = binaryRegiongrowing(img,seedPointOfOrgan);            
            img = imopen(img,strel('Disk',3));
            C = bwconncomp(img);
            stats = regionprops(C,'Area');
            [~,maxObj] = find([stats.Area] == max([stats.Area]));
            img = ismember(labelmatrix(C),maxObj);  
          
            isSuccess = true;
            break;
        end % end of if-else
        
    else
        % desA missed
        meanDensity = meanDensity - 30;
    end  % end of if-else
    
    
    
end % end of while 

end % end of function


function [img, isCircular, cirFactor, tolerance,features] = checkFeatures(I,L12,pixelSpacing,seedPointOfOrgan)

tolerance = 100;
isCircular = false;
cirFactor = -1;
features = 0;

while(~isCircular)

    if(tolerance < 25)
        break;
    end % end of if 
    
out = grayscaleSegmentation(I,pixelSpacing,round(seedPointOfOrgan),tolerance,10);
img = (out - L12)>0;
img = imfill(img,'holes');

img = binaryRegiongrowing(img,seedPointOfOrgan);

C = bwconncomp(img);
stats = regionprops(C,'Area');

if(min([stats.Area]) > 50)
    img = imopen(img,strel('Disk',3));
end

%figure,imshow3D(img);

C = bwconncomp(img);
stats = regionprops(C,'Area','Centroid','Eccentricity','Perimeter','ConvexArea','EquivDiameter','Solidity','MajorAxisLength','MinorAxisLength','Extrema');

if(~isempty(stats))
    [~,maxObj] = find([stats.Area] == max([stats.Area]));
    
    per = stats(maxObj).Perimeter;
    are = stats(maxObj).Area;
    
    if(are < 100)
        break;
    end % end of if
    
    %care = stats.ConvexArea;
    %dia = stats(maxObj).EquivDiameter;
    maa= stats(maxObj).MajorAxisLength;
    %mia = stats(maxObj).MinorAxisLength;
    
    [topMidPoint,bottomMidPoint,rightMidPoint,leftMidPoint] = getObjectMidpoints(stats(maxObj).Extrema);
    diaV = bottomMidPoint(2) - topMidPoint(2);
    diaH = rightMidPoint(1) - leftMidPoint(1);
    mdia = max(diaV,diaH);
    
    cir = (4*pi*are)/per^2; % circularity - 1 is exact circle
    rou = (4*are)/(pi*mdia^2); % roundness
    com = sqrt( ((4*are)/pi))/(maa);% compactness
    sol =  stats(maxObj).Solidity; %are/care; % solidity
    ecc = stats(maxObj).Eccentricity;   
    cirFactor = ((1-ecc)*cir);
    features = [(1-ecc),cir,rou,com,sol];
    if( cirFactor < 0.12 )
        tolerance = tolerance - 25;
    else
        isCircular = true;
         break;
    end
else    
    break;    
end % end of if

end % end of while


end % end of function