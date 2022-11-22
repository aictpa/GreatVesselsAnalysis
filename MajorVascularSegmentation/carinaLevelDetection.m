function  [img,rightCentroid,leftCentroid,rightMainBronchus,leftMainBronchus,sliceLoc,tmpOut,rightMidPoint,leftMidPoint,distanceCm] =  carinaLevelDetection(BO,PATIENT_POSITION,pixelSpacing,rightMainBronchus,leftMainBronchus,bifurPoint,thresholdVal)
%CARINALEVELDETECTION detect optimal carina level
%
%   Examples:
%       [img,rightCentroid,leftCentroid,rightMainBronchus,leftMainBronchus,sliceLoc,tmpOut,rightMidPoint,leftMidPoint,distanceCm] =  CARINALEVELDETECTION(BO,PATIENT_POSITION,pixelSpacing,rightMainBronchus,leftMainBronchus,bifurPoint,thresholdVal)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% Step 0 Initial Setup

P = pixelSpacing;

img = 0;
rightCentroid = 0;
leftCentroid = 0;

B2 = BO;


[~,~,TOTAL_SLICE_NO] = size(B2);
% get first and last slice number
[~,loopIncrementVal,lastSliceNo] = getFirstLastSliceLoc(PATIENT_POSITION,TOTAL_SLICE_NO);

tmpOut = false(size(B2));


for iSlice = (bifurPoint+loopIncrementVal):loopIncrementVal:lastSliceNo
    
    try
        
        imgROI = ~(B2(:,:,iSlice) > thresholdVal);
        
        % getting image
        imgRoiFilled = imgROI;
        imgForSeedFinding = imgRoiFilled;
        imgForSeedFindingRight = imgForSeedFinding  &  rightMainBronchus;
        imgForSeedFindingLeft = imgForSeedFinding  &  leftMainBronchus;
        
        
        %% Step 1 Get new Fields
        % get seed point
        thisSeedPointOfRight = getSeedPoint(imfill(imgForSeedFindingRight,'holes'));
        
        if(thisSeedPointOfRight ~= 0)
            if(imgRoiFilled(ceil(thisSeedPointOfRight(2)),ceil(thisSeedPointOfRight(1))) ~= 0)
                img = false(size(imgROI));
                img(thisSeedPointOfRight(2),thisSeedPointOfRight(1)) = 1;
                % Below code is equal to binary 2D-Region Growing
                img2N = bwlabel(imgRoiFilled,4);
                imgRightRG = (img2N == max(img2N(img == 1)));
            end % end of if
        end % end of if
        
        thisSeedPointOfLeft = getSeedPoint(imfill(imgForSeedFindingLeft,'holes'));
        
        if(thisSeedPointOfLeft ~= 0)
            if(imgRoiFilled(ceil(thisSeedPointOfLeft(2)),ceil(thisSeedPointOfLeft(1))) ~= 0)
                img = false(size(imgROI));
                img(thisSeedPointOfLeft(2),thisSeedPointOfLeft(1)) = 1;
                % Below code is equal to binary 2D-Region Growing
                img2N = bwlabel(imgRoiFilled,4);
                imgLeftRG = (img2N == max(img2N(img == 1)));
            end % end of if
        end % end of if
        
        C = bwconncomp(imgRightRG);
        statsTra = regionprops(C,'Area','Centroid');
        [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
        rightCentroid = statsTra(indexTra).Centroid;
        rightArea = statsTra(indexTra).Area;
        
        C = bwconncomp(imgLeftRG);
        statsTra = regionprops(C,'Area','Centroid');
        [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
        leftCentroid = statsTra(indexTra).Centroid;
        leftArea = statsTra(indexTra).Area;
                
        
        if(imgRightRG == imgLeftRG)
           
            [imgRightRG,imgLeftRG] = componentCheking(B2(:,:,iSlice),rightMainBronchus,leftMainBronchus,thresholdVal);
            
            C = bwconncomp(imgRightRG);
            statsTra = regionprops(C,'Area','Centroid');
            [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
            rightCentroid = statsTra(indexTra).Centroid;
            rightArea = statsTra(indexTra).Area;
            
            C = bwconncomp(imgLeftRG);
            statsTra = regionprops(C,'Area','Centroid');
            [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
            leftCentroid = statsTra(indexTra).Centroid;
            leftArea = statsTra(indexTra).Area;
    
        end
        
        
        %% Step 2 Checking Area
        if(rightArea > 1500)
            [imgRightRG,rightCentroid,~] = areaCheking(B2(:,:,iSlice),rightMainBronchus,thresholdVal);
        end
        
        if(leftArea > 1500)
            [imgLeftRG,leftCentroid,~] = areaCheking(B2(:,:,iSlice),leftMainBronchus,thresholdVal);
        end
        
        
        %% Step 3  Distance Calculations
        rightMainBronchus = imgRightRG;
        leftMainBronchus = imgLeftRG;
        
        tmpOut(:,:,iSlice) = rightMainBronchus + leftMainBronchus;
         
        
        C = bwconncomp(imgLeftRG);
        statsTra = regionprops(C,'Extrema');
        [~,~,~,leftMidPoint] = getObjectMidpoints(statsTra.Extrema);
        
        C = bwconncomp(imgRightRG);
        statsTra = regionprops(C,'Extrema');
        [~,~,rightMidPoint,~] = getObjectMidpoints(statsTra.Extrema);
        pixelDistance = abs(rightMidPoint(1) - leftMidPoint(1));
        
         
        %pixelDistance = abs(rightCentroid(1)-leftCentroid(1));
        
        distanceCm = ((P *  pixelDistance)/10);
        
        if(distanceCm > 0.75)
            %figure,imshow3D(B2(:,:,iSlice),[-900 600]);
            img = B2(:,:,iSlice);
            sliceLoc = iSlice;
            break;
        end
    catch ME
        fprintf('%s',ME.message);
    end
    
end

end % end of function

function   [imgRightRG,imgLeftRG] = componentCheking(imgROI,rightMainBronchus,leftMainBronchus,thresholdVal)

thresholdVal = thresholdVal - 15;

% getting image
imgRoiFilled =  ~(imgROI > thresholdVal);
imgForSeedFinding = imgRoiFilled;
imgForSeedFindingRight = imgForSeedFinding  &  rightMainBronchus;
imgForSeedFindingLeft = imgForSeedFinding  &  leftMainBronchus;


%%
% get seed point
thisSeedPointOfRight = getSeedPoint(imfill(imgForSeedFindingRight,'holes'));

if(thisSeedPointOfRight ~= 0)
    if(imgRoiFilled(ceil(thisSeedPointOfRight(2)),ceil(thisSeedPointOfRight(1))) ~= 0)
        img = false(size(imgROI));
        img(thisSeedPointOfRight(2),thisSeedPointOfRight(1)) = 1;
        % Below code is equal to binary 2D-Region Growing
        img2N = bwlabel(imgRoiFilled,4);
        imgRightRG = (img2N == max(img2N(img == 1)));
    end % end of if
end % end of if

thisSeedPointOfLeft = getSeedPoint(imfill(imgForSeedFindingLeft,'holes'));

if(thisSeedPointOfLeft ~= 0)
    if(imgRoiFilled(ceil(thisSeedPointOfLeft(2)),ceil(thisSeedPointOfLeft(1))) ~= 0)
        img = false(size(imgROI));
        img(thisSeedPointOfLeft(2),thisSeedPointOfLeft(1)) = 1;
        % Below code is equal to binary 2D-Region Growing
        img2N = bwlabel(imgRoiFilled,4);
        imgLeftRG = (img2N == max(img2N(img == 1)));
    end % end of if
end % end of if


while(imgRightRG == imgLeftRG)
    
    thresholdVal = thresholdVal - 15;
    
    % getting image
    imgRoiFilled =  ~(imgROI > thresholdVal);
    imgForSeedFinding = imgRoiFilled;
    imgForSeedFindingRight = imgForSeedFinding  &  rightMainBronchus;
    imgForSeedFindingLeft = imgForSeedFinding  &  leftMainBronchus;
    
    
    %% 
    % get seed point
    thisSeedPointOfRight = getSeedPoint(imfill(imgForSeedFindingRight,'holes'));
    
    if(thisSeedPointOfRight ~= 0)
        if(imgRoiFilled(ceil(thisSeedPointOfRight(2)),ceil(thisSeedPointOfRight(1))) ~= 0)
            img = false(size(imgROI));
            img(thisSeedPointOfRight(2),thisSeedPointOfRight(1)) = 1;
            % Below code is equal to binary 2D-Region Growing
            img2N = bwlabel(imgRoiFilled,4);
            imgRightRG = (img2N == max(img2N(img == 1)));
        end % end of if
    end % end of if
    
    thisSeedPointOfLeft = getSeedPoint(imfill(imgForSeedFindingLeft,'holes'));
    
    if(thisSeedPointOfLeft ~= 0)
        if(imgRoiFilled(ceil(thisSeedPointOfLeft(2)),ceil(thisSeedPointOfLeft(1))) ~= 0)
            img = false(size(imgROI));
            img(thisSeedPointOfLeft(2),thisSeedPointOfLeft(1)) = 1;
            % Below code is equal to binary 2D-Region Growing
            img2N = bwlabel(imgRoiFilled,4);
            imgLeftRG = (img2N == max(img2N(img == 1)));
        end % end of if
    end % end of if
       
    
end


end % end of function

function   [imgRG,centroid,thresholdVal] = areaCheking(imgROI,mainBronchus,thresholdVal)

thresholdVal = thresholdVal - 15;

% getting image
imgRoiFilled = ~(imgROI > thresholdVal);
imgForSeedFinding = imgRoiFilled  &  mainBronchus;

%% Step 1 Get new Fields
% get seed point
thisSeedPoint = getSeedPoint(imfill(imgForSeedFinding,'holes'));

if(thisSeedPoint ~= 0)
    if(imgRoiFilled(ceil(thisSeedPoint(2)),ceil(thisSeedPoint(1))) ~= 0)
        img = false(size(imgROI));
        img(thisSeedPoint(2),thisSeedPoint(1)) = 1;
        % Below code is equal to binary 2D-Region Growing
        img2N = bwlabel(imgRoiFilled,4);
        imgRG = (img2N == max(img2N(img == 1)));
    end % end of if
end % end of if

C = bwconncomp(imgRG);
statsTra = regionprops(C,'Area','Centroid');
[~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
centroid = statsTra(indexTra).Centroid;
area = statsTra(indexTra).Area;

while(area > 1500)
    
    thresholdVal = thresholdVal - 15;
    
    if(thresholdVal < -1024)
        break;
    end
    
    % getting image
    imgRoiFilled = ~(imgROI > thresholdVal);
    imgForSeedFinding = imgRoiFilled  &  mainBronchus;
    
    %% Step 1 Get new Fields
    % get seed point
    thisSeedPoint = getSeedPoint(imfill(imgForSeedFinding,'holes'));
    
    if(thisSeedPoint ~= 0)
        if(imgRoiFilled(ceil(thisSeedPoint(2)),ceil(thisSeedPoint(1))) ~= 0)
            img = false(size(imgROI));
            img(thisSeedPoint(2),thisSeedPoint(1)) = 1;
            % Below code is equal to binary 2D-Region Growing
            img2N = bwlabel(imgRoiFilled,4);
            imgRG = (img2N == max(img2N(img == 1)));
        end % end of if
    end % end of if
    
    C = bwconncomp(imgRG);
    statsTra = regionprops(C,'Area','Centroid');
    [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
    centroid = statsTra(indexTra).Centroid;
    area = statsTra(indexTra).Area;
end




end