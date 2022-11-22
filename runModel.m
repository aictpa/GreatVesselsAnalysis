% RUNMODEL  Setting up model  and  starting model pipeline
%
%   Examples:
%      out = RUNMODEL  

%   Copyright 2022 
%   Author  - Kahraman A. Teymur  
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

% Add path of Folders and Subfolders to the system
addpath(genpath(pwd));

% Clear matlab buffers
cl

%% Setting data path

% Please change the below code, If you want to place the dataset in 
% another location.
current_path = pwd;
path_single_data = fullfile(current_path, "Dataset/testSingleExam/");
path_multi_data = fullfile(current_path, "Dataset/testMultiExam/");


%% Choosing model running option  
% Set 'multi_test' variable to true, if you want to run model on
% multiple exams.
multi_test = true;

% If set it to "true", the output of the each steps will be illustrated.
% Not recommended. Set it "true" if you want to debug code.
showsFigures = false;


%% Initial settings 

%set datset path here
global datasetPath; 

if multi_test
    datasetPath =  path_multi_data;
else
    datasetPath =  path_single_data;
end


%% Compile C++ files

out = true;

if  isempty(dir(fullfile(current_path, "Toolbox/*/*.mex*")))
    out =  compileDependencies(current_path);
end


warning('off','all');

%% Start pipeline
if out
    out = startPipeline(datasetPath,current_path, showsFigures);
end

