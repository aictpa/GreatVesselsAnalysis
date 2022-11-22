function  deviation = setDeviation(hu)

% SETDEVIATION  Setting up optimal threshold value (HU)
%
%   Examples:
%      deviation = SETDEVIATION(hu)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

if(hu<100)
    deviation = 50;
elseif(hu<130)
    deviation = 70;
elseif(hu<150)
    deviation = 100;
elseif(hu<200)
    deviation = 125;
elseif(hu<250)
    deviation = 150;
elseif(hu<300)
    deviation = 170;
elseif(hu<400)
    deviation = 200;
else
    deviation = 220;
end


end % end of  function


  