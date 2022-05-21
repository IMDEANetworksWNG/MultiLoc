clear
close all
clc

pwd_str = pwd;

cd ../../


%% load LF data
load("mat_files/mD_track/LF/Data_3D.mat")

% take the number of points and number of routers
[n_points,routers] = size(AoA);

% take the power
power_point = power{1,5}{1,1};
% move to dB
power_point_dB = 10*log10(power_point);
% normalize it
power_point_dB = power_point_dB - max(power_point_dB);

% create the profile
angle_profile_LF = nan(180,1);

% take the aoa 
aoa_point_LF = AoA{1,5}{1,1}*(-1);

aoa_point_LF(2) = [];
power_point_dB(2) = [];

% put power to the aoas
angle_profile_LF(int64(aoa_point_LF + 90)) = power_point_dB;


% add noise
angle_profile_LF(isnan(angle_profile_LF)) = min(power_point_dB)- 10;

% plot it
figure, plot(-90:1:89,angle_profile_LF)
xlabel("AoA [deg]")
ylabel("Power [dB]")

%% load HF data
load("mat_files/indoor/HF/CSI/csi_indoor.mat")
load("mat_files/indoor/HF/CSI/aoa_optimal_rotation.mat")
load("mat_files/indoor/HF/FTM/ftm_indoor.mat")

load("mat_files/indoor/HF/optimal_rotation_index.mat")

% take the number of points and number of routers
[n_points,routers, rotations] = size(azimut_raw);

rotation = 5;
% rotations_id = [0, 45, 270, 315, 90, 135, 180, 225];
rotations_id = 0:45:315;

% take the power
power_point = power_raw{1,5,rotation};
% move to dB
power_point_dB = 10*log10(power_point);
% normalize it
power_point_dB = power_point_dB - max(power_point_dB);

% create the profile
angle_profile_HF = nan(180,1);

% take the aoa 
aoa_point_HF = azimut_raw{1,5,rotation};

% put power to the aoas
angle_profile_HF(aoa_point_HF + 90) = power_point_dB;

% add noise
angle_profile_HF(isnan(angle_profile_HF)) = min(power_point_dB)- 10;

% plot it
figure, plot(-90:1:89,angle_profile_HF)
xlabel("AoA [deg]")
ylabel("Power [dB]")

figure, plot(squeeze(calculated_distance(1,5,:)))
xticklabels(string(rotations_id))
figure, plot(squeeze(calculated_aoa(1,5,:)))
xticklabels(string(rotations_id))

cd(pwd_str)