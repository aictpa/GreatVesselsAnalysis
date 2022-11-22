function [xPoint,yPoint] = desAortaDetection(img,L,centroid, addedAngle, spaceNo)

%DESAORTADETECTION detect descending aorta
%
%   Examples:
%       [xPoint,yPoint] = DESAORTADETECTION(img,L,centroid, addedAngle, spaceNo)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


if(spaceNo == 1)
    searchSpace = (22.5 + addedAngle):7.5:(75 + addedAngle);
else %spaceNo == 2
    %searchSpace = 90:7.5:142.5;
    searchSpace = 60:7.5:120;
end

x(1) = centroid(1);
y(1) = centroid(2);
lineLength = 100;

%cmap = colormap(hot(numel(searchSpace)));

counter = 1;
%figure,imshow(img,[]),impixelinfo;

for angleN = searchSpace
    
    [x(2),y(2)] = getxy(x(1),y(1),angleN,lineLength);
    
    [cx,cy,c] = improfile(img,[x(1), x(2)], [y(1), y(2)], lineLength);
    [~,~,cL] = improfile(L,[x(1), x(2)], [y(1), y(2)], lineLength);
    
    c(c < 0) = 0;
    cLn = ~cL;
    cN = cLn .* c;
    
    C = bwconncomp(cN);
    stats = regionprops(C,'Area');
    [~,index] = find([stats.Area] == max([stats.Area]));
    
    objectMask = ismember(labelmatrix(C),index);    
    
    if(sum(objectMask) > 15)
        xMask =  objectMask .* cx;
        
        if(max(xMask) == 0)
            meanX(counter) = 0;
        else
            meanX(counter) = mean2(xMask(xMask~=0));
        end
        
        yMask =  objectMask .* cy;
        if(max(yMask) == 0)
            meanY(counter) = 0;
        else
            meanY(counter) = mean2(yMask(yMask~=0));
        end
    else
        meanX(counter) = 0;
        meanY(counter) = 0;
        
    end % end of outer if-else      
    
    % Illustration
    
        %normalized = (angleN-min(searchSpace))/(max(searchSpace)-min(searchSpace));
        %line(x, y,'Color',[abs(normalized-rand(1)) abs(normalized-rand(1)) abs(normalized-rand(1))])
        %line(x, y,'Color',cmap(counter,:))
    
    
    counter = counter + 1;
end % end of for

meanX = meanX(meanX ~= 0);
meanY = meanY(meanY ~= 0);



if(size(meanX,2)== 1)
    mask1 = 1;
    mask2 = 1;
else 
    
    rng(1,'twister');
    [idx,~] = kmeans(meanX',2);
    %disp( " idx " + idx);
    
    mask1 = idx == 1;
    mask2 = idx == 2;
end % end of if-else

if(sum(mask1) > sum(mask2))

    meanTemp = meanY' .* mask1;
    yPoint(1) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask1;
    xPoint(1) = mean2(meanTemp(meanTemp > 0));
    
    meanTemp = meanY' .* mask2;
    yPoint(2) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask2;
    xPoint(2) = mean2(meanTemp(meanTemp > 0));
    
elseif(sum(mask1) < sum(mask2))  
    meanTemp = meanY' .* mask2;
    yPoint(1) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask2;
    xPoint(1) = mean2(meanTemp(meanTemp > 0));
    
    meanTemp = meanY' .* mask1;
    yPoint(2) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask1;
    xPoint(2) = mean2(meanTemp(meanTemp > 0));
    
else %(sum(mask1) == sum(mask2))
    meanTemp = meanY' .* mask1;
    yPoint(1) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask1;
    xPoint(1) = mean2(meanTemp(meanTemp > 0));
    
    meanTemp = meanY' .* mask2;
    yPoint(2) = mean2(meanTemp(meanTemp > 0));
    meanTemp = meanX' .* mask2;
    xPoint(2) = mean2(meanTemp(meanTemp > 0));
end


end % end of function

%%
function [x,y] = getxy(x,y,angleN,lineLength)
x = x + lineLength * cosd(angleN);
y = y + lineLength * sind(angleN);
end % end of function
