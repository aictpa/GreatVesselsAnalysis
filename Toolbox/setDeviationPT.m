function  deviationPT = setDeviationPT(huPT)

% SETDEVIATIONPT  Setting up optimal threshold value (HU) for PT
%
%   Examples:
%      deviation = SETDEVIATIONPT(hu)

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

   if(huPT<100)
        deviationPT = 40;
    elseif(huPT<130)
        deviationPT = 50;
    elseif(huPT<150)
        deviationPT = 60;
    elseif(huPT<200)
        deviationPT = 70;
    elseif(huPT<250)
        deviationPT = 80;
    elseif(huPT<300)
        deviationPT = 100;
    elseif(huPT<400)
        deviationPT = 120;
    else
        deviationPT = 150;
    end


end % end of  function