clear
close all
clc

pwd_str = pwd;
cd ../../

% load("RX_coordinates.mat")
% load("points_coordinates.mat")
% load("data_distance_angle_true.mat")

addpath("functions/")


%% Configuration for md-Track
% Load antenna data
load("processed_data/antennas_mikrotik.mat")

att = 1e-1;

num_samples = 250;

% number of antennas
N = 6;

% frequency
freq = 60.48e9;

% speed of light
c = 3e8;

% the wavelength
lambda = c/freq;

% distance between antennas
d = lambda*0.58;

% step for the angle
step_angle = 1;

% take the codebook
[cb_az, theta_az] = Grid_AoA_distance(step_angle, N,d,lambda);

%% Process the data
 load('mat_files/indoor/HF/CSI/csi_indoor.mat');
[measurement_points,ap_ids,rotations] = size(csi_raw);
% ap_ids       = [43:47];
calculated_aoa = nan(measurement_points,ap_ids,rotations);
for point_id=1:measurement_points
    
    
    for ap_id=1:ap_ids
        for id_rotation = 1:rotations

            
        csi_data = csi_raw{point_id,ap_id,id_rotation};
        if (~isempty(csi_data))
            csi_data_ind = csi_data(3:5,:).';
            C_data = csi_data_ind * csi_data_ind';

            spectrum_aoa = MUSIC_opt_TH(C_data,cb_az,3);

            [~, Az_estimated] = max(spectrum_aoa);

            Az_estimated = rad2deg(theta_az(Az_estimated));

            % Power
            power = abs(att).^2;

            % Save the data for the mat files
            calculated_aoa(point_id,ap_id,id_rotation ) = Az_estimated(1);
        end

        % Save the error
      
        end
    end
end



out_path = 'mat_files/indoor/HF/';

save([out_path 'CSI/csi_indoor_aoa_music.mat'], 'calculated_aoa');

cd(pwd_str)