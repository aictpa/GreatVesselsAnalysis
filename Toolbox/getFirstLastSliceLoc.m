function [firstSliceNo,loopIncrementVal,lastSliceNo] = getFirstLastSliceLoc(PATIENT_POSITION,TOTAL_SLICE_NO)

%GETFIRSTLASTSLICELOC   getting first and last slice number of case
%
%   Examples:
%      [firstSliceNo,loopIncrementVal,lastSliceNo] = GETFIRSTLASTSLICELOC(PATIENT_POSITION,TOTAL_SLICE_NO)  

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

if(PATIENT_POSITION == 1) % Top to bottom Scan(Cranial to Caudal)=> 1
    firstSliceNo = 1;
    loopIncrementVal = 1;
    lastSliceNo = TOTAL_SLICE_NO;
else % Bottom to Top Scan (Caudal to Cranial)=> 2
    firstSliceNo = TOTAL_SLICE_NO;
    loopIncrementVal = -1;
    lastSliceNo = 1;
end % end of if-else



end % end of function