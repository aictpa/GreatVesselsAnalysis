function rotationAngle = calculateExamOrientation(imgIN)
% CALCULATEEXAMORIENTATION  calculate skewness of CT-scan
%
%   Examples:
%      rotationAngle = CALCULATEEXAMORIENTATION(imgIN,iCase)  


%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

%% Assumptions

THRESHOLD_VAL_FOR_SKEWNESS = -300;

SKEWNESS_ANGLE = 17;


%% Find Skewness of an Image

imgThresholded = imgIN > THRESHOLD_VAL_FOR_SKEWNESS;

C = bwconncomp(imgThresholded);
stats = regionprops(C, 'Orientation','Area');
[~,indexOfMaxArea] =  find([stats.Area] == max([stats.Area]));

if(abs(stats(indexOfMaxArea).Orientation) > SKEWNESS_ANGLE)
    rotationAngle = stats(indexOfMaxArea).Orientation;
else
    rotationAngle = 0;
end % end of function

end % end of function