function imgOut = binaryRegiongrowing(imgIn,seedPoints)

% BINARYREGIONGROWING  binary region growing.
%
%   Examples:
%      imgOut = BINARYREGIONGROWING(imgIn,seedPoints) 

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

seedPoints = round(seedPoints);
img = false(size(imgIn));
img(seedPoints(2),seedPoints(1)) = 1;

imgL = bwlabel(imgIn,4);
imgOut = (imgL == max(imgL(img == 1)));

end % end of funvtion