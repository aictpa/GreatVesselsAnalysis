function imgOut = eigenvalHessian(imgIn)
% EIGENVALHESSIAN  calculate eigenvalue of Hessian Matrix
%
%   Examples:
%      imgOut = EIGENVALHESSIAN(imgIn)

%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


%% Hessian Matrix Calculation Method 2

[gx,gy] = gradient(double(imgIn));
[gxx,gyx] = gradient(gx);
[gxy,gyy] = gradient(gy);

%% Eigenvalues Calculation Method 2

trA = gxx + gyy; % trace(A)
detA = gxx.*gyy - gxy.*gyx; % det(A)

% L1 & L2 are eigenval 1 & 2 respectively
L1 = (trA + sqrt(trA.^2 - 4.*detA) )/2;

imgOut = L1 > mean2(L1(L1>0));

end % end of funtion