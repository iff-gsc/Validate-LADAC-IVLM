function [] = addPathIvlmValidation()

% go to correct directory
init_scripts_path = fileparts(mfilename('fullpath'));
cd(init_scripts_path);

% go to initialization directory
cd ..

% add folders to path
addpath(genpath(pwd));

end