function  [sdDesA, desAorta, BWt, seedPoint, deviation] = getDesAortaSD(imgDesAorta,angle,ctpa,sliceLoc,pixelSpacing)
%GETDESAORTASD calculate SD of density(HU) of descending aorta
%
%   Examples:
%       [sdDesA, desAorta, BWt, seedPoint, deviation] = GETDESAORTASD(imgDesAorta,angle,ctpa,sliceLoc,pixelSpacing)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


if(abs(angle) > 25)
    imgDesAorta = imrotate(imgDesAorta,-angle,'nearest','crop');
end % end of if


hu1 = bsxfun(@times, ctpa(:,:,sliceLoc), cast(imgDesAorta, class(ctpa(:,:,sliceLoc))));
hu = mean2(hu1(hu1~=0));


j = imdilate(imgDesAorta,strel('disk',3));
bw=bwdist(~j);
radius = double(round(max(bw(:))));

if(radius<10)
    radiusDeviation = 15;
elseif(radius<15)
    radiusDeviation = 10;
else
    radiusDeviation = 7;
end


[seedPointX,seedPointY] = find(bw == max(bw(:)));

seedPoint = [seedPointY(1) seedPointX(1)];

t = linspace(0, 2*pi, 50);
r = radius+radiusDeviation;
c = [seedPoint(1) seedPoint(2)];
BW = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512, 512);
BWt=~BW;
BW = bsxfun(@times, ctpa(:,:,sliceLoc), cast(BW, class(ctpa(:,:,sliceLoc))));
BW(BWt) = -900;

I2 = imdiffusefilt(BW,'NumberOfIterations',20);
%I2 = anisocpp(double(BW),30.0,2,20,1/7);

to = ctpa(:,:,sliceLoc);

to(I2>0) = I2(I2>0);

deviation = setDeviation(hu);
BW = grayscaleSegmentation(to,pixelSpacing,seedPoint,deviation,10);


BW = bsxfun(@times, ctpa(:,:,sliceLoc), cast(BW, class(ctpa(:,:,sliceLoc))));

BW = BW > 0;


L = eigenvalHessian(I2);
E = edge(I2,'canny');
LE = (L + E) > 0;
BW = (BW - LE) > 0;

BW = imfill(BW,'holes');
BW = imopen(BW,strel('disk',7));

BW = regiongrowing(BW,1,[seedPoint(2) seedPoint(1)]);
desAorta = BW;


%% Step 3 get circular ROI
t = linspace(0, 2*pi, 50);
r = sqrt(1/pi)*10/pixelSpacing;
c = [seedPoint(1) seedPoint(2)];
BWSD = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512, 512);
BWt=~BWSD;
BWSD = bsxfun(@times, ctpa(:,:,sliceLoc), cast(BWSD, class(ctpa(:,:,sliceLoc))));
BWSD(BWt) = -900;

sdDesA = round(std2(BWSD(BWSD ~=-900)));
% 
% colorMap=[1 0 0; 0 0 1];
% figure,
% RGB = labeloverlay(mat2gray(ctpa(:,:,sliceLoc)),BWSD>0,'Colormap',colorMap,'Transparency',0.25);
% imshow(RGB)

        
end % end of function