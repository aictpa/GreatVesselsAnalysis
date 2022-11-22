function  [temp,angle] = singleRotation(img)
%SINGLEROTATION re-check CT-scan agin if rotation is needed
%
%   Examples:
%       [temp,angle] = SINGLEROTATION(img)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


roi = img > -300;


%%

imgCCL = bwconncomp(roi);
stats = regionprops(imgCCL,'Area','Image');
idx =  find([stats.Area]==max([stats.Area]));
imgLung = ismember(labelmatrix(imgCCL), idx);
outerLung = imfill(imgLung,'holes');


outerLung = imopen(outerLung,strel('disk',5));
%imshowV(outerLung);

box = regionprops(outerLung,'BoundingBox');
box = box.BoundingBox;
x1 = ceil(box(2));
x2 = ceil(box(2)+box(4));
if(x2 > 512)
    x2 = 512;
end
y1 = ceil(box(1));
y2 = ceil(box(1)+box(3));
if(y2 > 512)
    y2 = 512;
end

%Croping
mask = outerLung(x1:x2,y1:y2);
midIndex = round(size(mask)/2);
mask = boundarymask(mask);


horizontalMidPoint  = midIndex(2);%  + y1;  % ~ +-3
verticalMidPoint =  midIndex(1);% + x1; % ~ +-3


%% take 3rd quarter part

imgQ3 = mask(verticalMidPoint:size(mask,1),1:horizontalMidPoint);
%imshowV(imgQ3);
[xr, yr] = ind2sub(size((flip(rot90(flip(imgQ3))))), find((flip(rot90(flip(imgQ3))))));
pr = polyfit(xr, yr, 1);
% vr = polyval(pr, xr);
% figure
% plot(xr,yr,'x','MarkerEdgeColor','black')
% hold on
% plot(xr, vr)
% hold off
% grid on;

%m1 = (v1(2)-v1(1))/(x(2)-x(1)); % m1=tan(Q) -> Q = atan(m1);
mr=pr(1);
angleR = atand(mr);

% [Na,~] = histcounts(yr);
% camberR = max(Na);


%% take 4th quarter part

imgQ4 = mask(verticalMidPoint:size(mask,1),horizontalMidPoint:size(mask,2));
%imshowV(imgQ4);
[xl, yl] = ind2sub(size((flip(rot90(flip(imgQ4))))), find((flip(rot90(flip(imgQ4))))));
pl = polyfit(xl, yl, 1);
% vl = polyval(pl, xl);
% figure
% plot(xl,yl,'x','MarkerEdgeColor','black')
% hold on
% plot(xl, vl)
% hold off
% grid on;


% m2 = (v2(2)-v2(1))/(x(2)-x(1)); % m1=tan(Q) -> Q = atan(m1);
ml=pl(1);
angleL = atand(ml);



edges = [0 5 10 15 20];

[Nb1,~] = histcounts(yl,edges);
camberL = sum(Nb1);

[Na1,~] = histcounts(yr,edges);
camberR = sum(Na1);



%% Select right angle

if(camberR <= camberL)
    angle = abs(angleR);
else
    angle = -angleL;
end


if( ~(abs(camberR-camberL) > 175 )) % old one 175 working 100% with min slice thickness
    angle = 0;
end


temp = img;

if(abs(angle) > 25)
    
    %% Take care of zero pixel
    
    t = false(size(temp));
    t(temp == 0) = 1;
    t = imrotate(t,angle,'nearest','crop');
    
    temp = imrotate(temp,angle,'nearest','crop');
    
    temp(temp==0) = -1024;
    temp(t) = 0;
    
end % end of if



end % end of function

