%% Extractors %%
% extract CSI data
run("Extractors/CSI_extractor_raw_outdoor.m")

% extract FTM data
run("Extractors/FTM_extractor_outdoor.m")

%% Calibrate data %%
% CSI
run("Extractors/Calibrate_CSI_Data_outdoor.m")
% FTM
run("Extractors/Calibrate_FTM_Data_outdoor.m")

% %% Copy to another folder
% mkdir("../mat_files/indoor/")
% mkdir("../mat_files/indoor/LF")
% system("cp -a ../mat_files/CSI ../mat_files/indoor/LF/")
% 
% mkdir("../mat_files/indoor/")
% mkdir("../mat_files/indoor/LF")
% system("cp -a ../mat_files/FTM ../mat_files/indoor/LF/")
% 
%% Run mD-track %%
% Extract 3D md-track data
run("Estimator/Generate_Data_mD_track_3D_Outdoor.m")
% process the 3D mD-track data
run("Process_Estimator_Data/mD_track_data_outdoor.m")
run("Process_Estimator_Data/mD_track_Data_Direct_Path_outdoor.m")
% 
%% Run FTM 
run("Process_FTM/Generate_Distances_outdoor.m")
% 
%% Run coordinates
addpath(genpath("../"))
run("Coordinates/Generate_Data_Coordinates_LF_outdoor.m")
