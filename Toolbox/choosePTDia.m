function out = choosePTDia(ptDia)
% CHOOSEPTDIA  Choosing correct PT diameter from the candidates.
%
%   Examples:
%      out = CHOOSEPTDIA(ptDia)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


    if(size(ptDia,2)==0)  
        out = 0;
    else
        [~,F,~] = mode(ptDia);
        if(F>1)           
            out = round(median(ptDia));
        else            
            temp = ptDia;        
            out = round(mean(temp));
        end         
    end
    
end

