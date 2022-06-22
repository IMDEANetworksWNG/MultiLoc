% cd ../
% addpath(genpath("."))
%% Indoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_indoor.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_indoor.m")


%% Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_outdoor.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_outdoor.m")


%% Indoor/Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_in_out.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_in_out.m")


%% IMDEA
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_imdea_amb.m")

% Coordinates LF
run("Coordinates/Generate_Data_Coordinates_LF_imdea.m")

%% MultiLoc
close all
run("Coordinates/MultiLoc.m")

