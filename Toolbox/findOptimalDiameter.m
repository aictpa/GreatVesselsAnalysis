function [diaRho,diaRho3,diaLines,diaLines3,idxOfPre,indexRho] = findOptimalDiameter(img)
%FINDOPTIMALDIAMETER calculating accurate diameter of tubular objects
%
%   Examples:
%       [diaRho,diaRho3,diaLines,diaLines3,idxOfPre,indexRho] = FINDOPTIMALDIAMETER(img)
%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

dia3 = 0;
nLines = 10;
mLength = 15;
[H,T,R] = hough(img,'Theta',-75);
P  = houghpeaks(H,nLines,'threshold',ceil(0.25*max(H(:))));
lines = houghlines(img,T,R,P,'FillGap',nLines,'MinLength',mLength);

ver_len = zeros(1,size(lines,2));
y_coorU = zeros(1,size(lines,2));

%imshowV(img), hold on
for k = 1:length(lines)
    
    ver_len(k) = norm(lines(k).point1 - lines(k).point2);
    y_coorU(k) = lines(k).point2(2);%;min(lines(k).point1(2),lines(k).point2(2));
    xy = [lines(k).point1; lines(k).point2];
    %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
    
    
end

if(size(y_coorU,2) == 1)
    dia = ver_len;
    dia3 = ver_len;
    diaLines =  [lines(1).point1 lines(1).point2];
    diaLines3 = [lines(1).point1 lines(1).point2];
    idxOfPre = 1;
    
    diaRho = dia;
    diaRho3 = dia3;
    indexRho = 1;
    
elseif(size(y_coorU,2) == 2)
    dia = max(ver_len);
    dia3 = max(ver_len);
    idx = ver_len == max(ver_len);
    idx = find(idx == 1);
    idx = idx(1);
    diaLines =  [lines(idx).point1 lines(idx).point2];
    diaLines3 = [lines(idx).point1 lines(idx).point2];
    idxOfPre = idx;
    indexRho = idx;
    diaRho = dia;
    diaRho3 = dia3;
else
    y_coor = sort(y_coorU);
    
    idx = y_coorU == y_coor(2);
    if (sum(idx) == 1)
        dia = ver_len(idx);
        idx = find(idx == 1);
        idxOfPre = idx;
    else
        dia = max(ver_len(idx));
        idx2 = ver_len == dia;
        idx = idx & idx2;
        idx = find(idx == 1);
        idx = idx(1);
        idxOfPre = idx;
    end
    
    %%
    rho = sort(abs([lines.rho]));
    rho = abs([lines.rho]) == rho(2);
    indexRho = find(rho == 1);
    diaRho = max(ver_len(indexRho));
    indexRho = indexRho(ver_len(indexRho) == diaRho);
    
    %%
    diaLines =  [lines(indexRho).point1 lines(indexRho).point2];
    
    idx = y_coorU == y_coor(3);
    if (sum(idx) == 1)
        dia3 = ver_len(idx);
        idx = find(idx == 1);
    else
        dia3 = max(ver_len(idx));
        idx2 = ver_len == dia3;
        idx = idx & idx2;
        idx = find(idx == 1);
        idx = idx(1);
    end
    
    %
    rho = sort(abs([lines.rho]));
    rho = abs([lines.rho]) == rho(3);
    indexRho3 = find(rho == 1);
    diaRho3 = max(ver_len(indexRho3));
    indexRho3 = indexRho3(ver_len(indexRho3) == diaRho3);
    %
    diaLines3 = [lines(indexRho3).point1 lines(indexRho3).point2];
    
end % end of if

end % end of function