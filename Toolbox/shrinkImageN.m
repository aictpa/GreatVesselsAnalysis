function shiftedImage = shrinkImageN(imgIn,ratio,drt)

% SHRINKIMAGEN  shrinking given image by relative ratio.
%
%   Examples:
%      shiftedImage = SHRINKIMAGEN(imgIn,ratio,drt) 

%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

[rows, columns] = size(imgIn); % Save original size.
imgCCL = bwconncomp(imgIn);
if(imgCCL.NumObjects>1)
    stats = regionprops(imgCCL,'Area');
    idx =  find([stats.Area]==max([stats.Area]));
    imgIn = ismember(labelmatrix(imgCCL), idx);
end
measurements1 = regionprops(imgIn, 'Centroid');

if(drt==1)
    imgSmall = padarray(imgIn,[round(rows*ratio) round(columns*ratio)]);
    imgSmall = imresize(imgSmall, [rows, columns] ); % Same size as original.
    imgCCL = bwconncomp(imgSmall);
    if(imgCCL.NumObjects>1)
        stats = regionprops(imgCCL,'Area');
        idx =  find([stats.Area]==max([stats.Area]));
        imgSmall = ismember(labelmatrix(imgCCL), idx);
    end
    measurements2 = regionprops(imgSmall, 'Centroid');
    
    rowsToShift = round(measurements1.Centroid(2)- measurements2.Centroid(2));
    columnsToShift = round(measurements1.Centroid(1) - measurements2.Centroid(1));
    shiftedImage = circshift(imgSmall, [rowsToShift columnsToShift]);
else   
    stats = regionprops(imgIn,'Image','Centroid');
    imgT = imresize(stats.Image,ratio);
    imgCCL = bwconncomp(imgT);
    if(imgCCL.NumObjects>1)
        stats = regionprops(imgCCL,'Area');
        idx =  find([stats.Area]==max([stats.Area]));
        imgT = ismember(labelmatrix(imgCCL), idx);
    end
    imgBigger = padarray(imgT,[round((rows-size(imgT,1))/2) round((columns-size(imgT,2))/2)]);
    imgBigger = imresize(imgBigger, [rows, columns] ); % Same size as original.
    imgCCL = bwconncomp(imgBigger);
    if(imgCCL.NumObjects>1)
        stats = regionprops(imgCCL,'Area');
        idx =  find([stats.Area]==max([stats.Area]));
        imgBigger = ismember(labelmatrix(imgCCL), idx);
    end
    measurements2 = regionprops(imgBigger, 'Centroid');
    
    rowsToShift = round(measurements1.Centroid(2)- measurements2.Centroid(2));
    columnsToShift = round(measurements1.Centroid(1) - measurements2.Centroid(1));
    shiftedImage = circshift(imgBigger, [rowsToShift columnsToShift]);
end



end % end of function