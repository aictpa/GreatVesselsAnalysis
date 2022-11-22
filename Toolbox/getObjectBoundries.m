function [topLeft,topRight,rightTop,rightBottom,bottomRight,bottomLeft,leftBottom,leftTop] = getObjectBoundries(extrema)

%GETOBJECTBOUNDRIES  getting the boundries of object
%
%   Examples:
%      [topLeft,topRight,rightTop,rightBottom,bottomRight,bottomLeft,leftBottom,leftTop] = GETOBJECTBOUNDRIES(extrema)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% Boundry points

topLeft = extrema(1,:);
topRight = extrema(2,:);
rightTop = extrema(3,:);
rightBottom = extrema(4,:);
bottomRight = extrema(5,:);
bottomLeft = extrema(6,:);
leftBottom = extrema(7,:);
leftTop = extrema(8,:);
    

end % end of function