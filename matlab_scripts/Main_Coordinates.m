% cd ../
% addpath(genpath("."))
%% Indoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_indoor.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_indoor.m")

% Coordinates multiband baseline
run("Coordinates/Generate_Data_Coordinates_multiband_baseline_indoor.m")

% Coordinates hf baseline
run("Coordinates/Generate_Data_Coordinates_HF_baseline_indoor.m")

% Execute the MultiLoc localization
run("Coordinates/With_multiple_antenna_arrays/With_all_arrays_Final_more_subplot2.m")


% Plot all of them
run("Coordinates/With_multiple_antenna_arrays/cdf_with_all_arrays.m")

%% Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_outdoor.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_outdoor.m")

% Coordinates multiband baseline
run("Coordinates/Generate_Data_Coordinates_multiband_baseline_outdoor.m")

% Coordinates hf baseline
run("Coordinates/Generate_Data_Coordinates_HF_baseline_outdoor.m")

% Execute the MultiLoc localization
run("Coordinates/Outdoor_with_multiple_antenna_arrays/With_all_arrays_Final_More_subplot_2.m")

% Plot all of them
run("Coordinates/Outdoor_with_multiple_antenna_arrays/cdf_with_all_arrays.m")
% cd matlab_scripts/


%% 2APs case

% Execute the MultiLoc localization
run("Coordinates/Only_2_APs/With_all_arrays_Final_more_subplot_2.m")

% Plot all of them
run("Coordinates/Only_2_APs/cdf_with_all_arrays.m")


%% Indoor/Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_in_out.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_in_out.m")

% Coordinates multiband baseline
run("Coordinates/Generate_Data_Coordinates_multiband_baseline_in_out.m")

% Coordinates hf baseline
run("Coordinates/Generate_Data_Coordinates_HF_baseline_in_out.m")

% Execute the MultiLoc localization
run("Coordinates/In_out/With_all_arrays_Final_more_subplot2.m")

% Plot all of them
run("Coordinates/In_out/cdf_with_all_arrays.m")


%% IMDEA
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_imdea_amb.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_imdea.m")

% Coordinates multiband baseline
run("Coordinates/Generate_Data_Coordinates_multiband_baseline_imdea.m")

% Coordinates hf baseline
run("Coordinates/Generate_Data_Coordinates_HF_baseline_imdea.m")

% Execute the MultiLoc localization
run("Coordinates/imdea/With_all_arrays_Final_more_subplot2.m")

% Plot all of them
% run("Coordinates/imdea/cdf_with_all_arrays.m")


%% Spotfi
run("Coordinates/Trianulation_Spotfi_all")
