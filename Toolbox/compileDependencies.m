function out =  compileDependencies(main_path)

% COMPILEDEPENDENCIES Compiling C++ files
%
%   Examples:
%      out = COMPILEDEPENDENCIES 
%      If the output is equal to "1", the dependencies have been compiled
%      successfully.
%      "main_path" holds the main folder location.
%
%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


error = false;
out = false;

% Check if C++ compiler is installed.
addons = matlab.addons.installedAddons;
Index = strfind(addons.Name, 'C/C++ Compiler');
Index = find(~cellfun(@isempty,Index));

if ~Index
   warning("Please be sure that a Mex Compiler is installed!");
   disp('For installation, please <a href = "https://se.mathworks.com/support/requirements/supported-compilers.html">follow the link.</a>');
 
end

path_aniso = fullfile(main_path, "Toolbox");
path_regionGrowing = fullfile(main_path, "Toolbox\ThirdPartyFunctions");

try 
    if ismac    
    
        mex -setup C++    
        cd(path_aniso)
        mex anisocpp.cpp
        cd(main_path)

        cd(path_regionGrowing)
        mex RegionGrowing_mex.cpp
        cd(main_path) 
       
    elseif isunix    
    
        cd(path_aniso)
        mex anisocpp.cpp
        cd(main_path)

        cd(path_regionGrowing)
        mex RegionGrowing_mex.cpp
        cd(main_path) 
       
    elseif ispc
        
        mex -setup C++    
        cd(path_aniso)
        mex anisocpp.cpp
        cd(main_path)

        cd(path_regionGrowing)
        mex RegionGrowing_mex.cpp
        cd(main_path) 
    else
        disp('Platform not supported');
    end

catch exception    
    error = true;
    disp(exception)
end

if ~error
    out  = true;
    disp("Dependencies have been successfully compiled.");    
end



end
  
