function [sliceImg,sliceHU,sliceInfo] = getCTPAExam(dcmDatasetPath)
% GETCTPAEXAM  Load dicom data to the array
%
%   Examples:
%      out = GETCTPAEXAM(dcmDatasetPath)
%      "dcmDatasetPath" holds location of the data.


%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Get images number in a file
fileWithIsdir = dir(dcmDatasetPath);
file = fileWithIsdir(not([fileWithIsdir.isdir]));
filesName = extractfield(file, 'name');
% sort files name  numerically 'ascending'
[allImgFiles,~] = sort_nat(filesName);

sliceNumber = length(allImgFiles);
fileNameOfFirstSlice = allImgFiles{1};



% Read  dicom file for getting some information
info = dicominfo(sprintf('%s/%s', dcmDatasetPath,fileNameOfFirstSlice));

emptyImg = zeros(info.Rows,info.Columns);

% Creating struct
fields = {'name','img','HU','info','lungThresh','rightLung','leftLung','masked'};


sliceImg(1:sliceNumber) = struct(fields{1},'',fields{2},emptyImg);
sliceHU(1:sliceNumber) = struct(fields{1},'',fields{3},emptyImg);
sliceInfo(1:sliceNumber) = struct(fields{1},'',fields{4},'');

% Image Acquisition
% Get slices and calculate HU
% %Display all CT-scans
parfor scanNr = 1:sliceNumber    
           
    info = dicominfo(sprintf('%s/%s', dcmDatasetPath,allImgFiles{scanNr}));   
    img = dicomread(sprintf('%s/%s', dcmDatasetPath,allImgFiles{scanNr}));
   
    
    if isa(img,'uint16')
        img = int16(img);
    end 
    
    % For CT scans
    img(img == -2000) = 0;
    
    % Hounsfield Unit Calculation
    HU=info.RescaleSlope*img+info.RescaleIntercept;     
    
    sliceImg(scanNr) = struct(fields{1},sprintf('%s',allImgFiles{scanNr}),fields{2},img);
    sliceHU(scanNr) = struct(fields{1},sprintf('%s',allImgFiles{scanNr}),fields{3},HU);
    sliceInfo(scanNr) = struct(fields{1},sprintf('%s',allImgFiles{scanNr}),fields{4},info);
  
       
end


end