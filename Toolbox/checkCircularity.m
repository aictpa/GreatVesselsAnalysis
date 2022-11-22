function [cirD, ver_len, hor_len] = checkCircularity(pt)

% CHECKCIRCULARITY  check circularity of object
%
%   Examples:
%      [cirD, ver_len, hor_len] = CHECKCIRCULARITY(pt) 

%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% Horizontal

[H,T,R] = hough(pt,'Theta',0);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(pt,T,R,P,'FillGap',5,'MinLength',7);

hor_len = 0;
for k = 1:length(lines)  
   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > hor_len)
      hor_len = len;     
   end
end


%% Vertical

[H,T,R] = hough(pt,'Theta',-90);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(pt,T,R,P,'FillGap',5,'MinLength',7);

ver_len = 0;
for k = 1:length(lines)  
   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > ver_len)
      ver_len = len;      
   end
end

cirD = abs(ver_len-hor_len);

%r = round(max(ver_len,hor_len)/2);

end