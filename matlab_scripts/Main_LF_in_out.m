%% Extractors %%
% extract CSI data
run("Extractors/CSI_extractor_raw_in_out.m")

% extract FTM data
run("Extractors/FTM_extractor_in_out.m")

%% Calibrate data %%
% CSI
run("Extractors/Calibrate_CSI_Data_in_out.m")
% FTM
run("Extractors/Calibrate_FTM_Data_in_out.m")

%% Copy to another folder
% mkdir("../mat_files/indoor/")
% mkdir("../mat_files/indoor/LF")
% system("cp -a ../mat_files/CSI ../mat_files/indoor/LF/")
% 
% mkdir("../mat_files/indoor/")
% mkdir("../mat_files/indoor/LF")
% system("cp -a ../mat_files/FTM ../mat_files/indoor/LF/")

%% Run mD-track %%
% Extract 3D md-track data
run("Estimator/Generate_Data_mD_track_3D_in_out.m")
% process the 3D mD-track data
run("Process_Estimator_Data/mD_track_data_in_out.m")
run("Process_Estimator_Data/mD_track_Data_Direct_Path_in_out.m")
% 
% %% Run FTM 
run("Process_FTM/Generate_Distances_in_out.m")
