clear
close all
clc

pwd_str = pwd;

cd ../../

csi_indoor  = load('mat_files/indoor/HF/CSI/csi_indoor.mat');
aoa_mD_track = csi_indoor.calculated_aoa;

csi_indoor  = load('mat_files/indoor/HF/CSI/csi_indoor_aoa_music.mat');
aoa_music = csi_indoor.calculated_aoa;

load("mat_files/indoor/HF/optimal_rotation_index.mat")

load('Maps/Map_indoor/data_distance_angle_true.mat');
% LOS vs NLOS positions
load('mat_files/indoor_map/Grid_LOS.mat')


[n_points, routers, rotations] = size(aoa_mD_track);

aoa_music_opt = zeros(n_points,routers);
aoa_mD_track_opt = zeros(n_points,routers);

for ii = 1:n_points
    for jj = 1:routers
        aoa_music_opt(ii,jj) = aoa_music(ii,jj,indoor_optimal_rotation_index(ii,jj));
        aoa_mD_track_opt(ii,jj) = aoa_mD_track(ii,jj,indoor_optimal_rotation_index(ii,jj));
        
    end
end

error_music = aoa_music_opt - angles_RX;
error_mD_track = aoa_mD_track_opt - angles_RX;

% error_music_los = zeros(size(error_music));
error_music_los = error_music.* LOS_RX;
error_music_los(error_music_los == 0) = nan;

% error_musis_nlos = zeros(error_music);
error_music_nlos = error_music.* ~LOS_RX;
error_music_nlos(error_music_nlos == 0) = nan;

% error_mD_track_los = zeros(error_mD_track);
error_mD_track_los = error_mD_track.* LOS_RX;
error_mD_track_los(error_mD_track_los == 0) = nan;

% error_mD_track_nlos = zeros(error_mD_track);
error_mD_track_nlos = error_mD_track.* ~LOS_RX;
error_mD_track_nlos(error_mD_track_nlos == 0) = nan;


figure,
cdfplot(abs(error_music_los(:)))
hold on
cdfplot(abs(error_mD_track_los(:)))
legend("MUSIC", "mD-track")


figure,
cdfplot(abs(error_music_nlos(:)))
hold on
cdfplot(abs(error_mD_track_nlos(:)))
legend("MUSIC", "mD-track")

cd(pwd_str)

