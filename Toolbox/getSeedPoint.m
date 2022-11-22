function seedPoint = getSeedPoint(imgForSeedFinding)

%GETSEEDPOINT  getting the seed point of given object
%
%   Examples:
%      [seedPoint] = GETSEEDPOINT(imgForSeedFinding)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})



C = bwconncomp(imgForSeedFinding);
stats= regionprops(C,'Area');
[~,index] = find([stats.Area] == max([stats.Area]));

if (isempty(index))
    
    seedPoint = 0;
    
else
    
    imgTemp = ismember(labelmatrix(C), index);   
   
    imgShrinked = bwmorph(imgTemp,'shrink',Inf);
    C = bwconncomp(imgShrinked);
    stats = regionprops(C,'Area','Centroid');
    [~,index] = find([stats.Area] == max([stats.Area]));
    seedPoint = stats(index).Centroid;
    
end % end of if-else

end % end of function