function [topMidPoint,bottomMidPoint,rightMidPoint,leftMidPoint] = getObjectMidpoints(extrema)

%GETOBJECTMIDPOINTS  getting the middle  boundry points of object
%
%   Examples:
%      [topMidPoint,bottomMidPoint,rightMidPoint,leftMidPoint] = GETOBJECTMIDPOINTS(extrema)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% Boundry points

[topLeft,topRight,rightTop,rightBottom,bottomRight,bottomLeft,leftBottom,leftTop] = getObjectBoundries(extrema);


%% Mid points

bottomMidPoint = (bottomLeft + bottomRight)/2;
rightMidPoint = (rightBottom + rightTop)/2;
leftMidPoint = (leftBottom + leftTop)/2;
topMidPoint = (topLeft + topRight)/2;

    
end % end of function