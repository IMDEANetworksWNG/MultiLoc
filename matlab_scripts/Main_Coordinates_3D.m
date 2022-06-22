% cd ../
% addpath(genpath("."))
%% Indoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_indoor_elevation.m")

%% Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_outdoor_elevation.m")


%% Indoor/Outdoor
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_in_out_elevation.m")

%% IMDEA
% Coordinates HF
run("Coordinates/Generate_Data_Coordinates_HF_imdea_amb_elevation.m")

%% MultiLoc 3D
close all
run("Coordinates/MultiLoc_3D.m")

