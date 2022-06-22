clear
close all
clc

pwd_str = pwd;
cd ../../

load("mat_files/indoor/HF/CSI/csi_indoor.mat")
calculated_aoa_old = calculated_aoa;
csi_raw_old = csi_raw;

load("mat_files/indoor/HF/CSI/csi_indoor2.mat")

notequal = calculated_aoa - calculated_aoa_old



cd(pwd_str)