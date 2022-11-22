
function BW = grayscaleSegmentation(imgO,pixelSpacing,seedPoints,tolerance,radius)

% GRAYSCALESEGMENTATION  segmentation method using grayscale level
%
%   Examples:
%      BW = grayscaleSegmentation(imgO,pixelSpacing,seedPoints,tolerance,radius)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

%% Step 1 set seed point
row = seedPoints(2);
col = seedPoints(1);


%% Step 2 get density of seed point
%originalHUOfSeedPoint = imgO(row,col);
[imgH,imgW,~] = size(imgO);


%% Step 3 get circular ROI
t = linspace(0, 2*pi, 50);   
r = radius*pixelSpacing;
c = [col row]; 
BW = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), imgH, imgW);


%% Step 4 calculate ave. density of ROI
maskedImg = bsxfun(@times, imgO, cast(BW, class(imgO)));
meanHUROI = round(mean2(maskedImg(maskedImg ~= 0)));

% replace originalHUOfSeedPoint with meanHUROI
imgO(row,col) = meanHUROI;


%% Step 5  Gray-scale segmentation
% with convenient parameter this method gives more promising result than
% other threshold + RG methods. 

%tolerance = 150;
BW = grayconnected(imgO,row,col,tolerance);
%figure,imshow3D(BW);


end