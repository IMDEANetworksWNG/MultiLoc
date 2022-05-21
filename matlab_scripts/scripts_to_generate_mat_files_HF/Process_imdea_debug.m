clear
close all
clc

pwd_str = pwd;
cd ../../
addpath(genpath("auxiliar_functions/"))
mkdir("mat_files")
mkdir("mat_files/imdea")

% load("mat_files/indoor_map/RX_coordinates.mat")
% load("mat_files/indoor_map/points_coordinates.mat")
% load("mat_files/indoor_map/data_distance_angle_true.mat")



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
[cb_az, theta_az] = Grid_AoA(step_angle, N,d,lambda);
[cb_el, theta_el] = Grid_AoA(step_angle, N,d,lambda);

%% Process the data
clc
ap_order = [43,52,54];
clients = [47,45];
clients_str = repmat(string(clients),4,1);
rotations = 90:45:225;
rotations_str = repmat(string(rotations.'),1,length(clients));
client_rotation = strcat(clients_str,"_pan_",rotations_str);
% regular one
cli_ids1 = cellstr((client_rotation(:)).');
cli_ids = cli_ids1;
% for especific cases
clients_str2 = repmat(string(flip(clients)),4,1);
client_rotation2 = strcat(flip(clients_str2),"_pan_",rotations_str);
cli_ids2 = cellstr((client_rotation2(:)).');

% cli_ids   = {'45_pan_0', '45_pan_45', '45_pan_270', '45_pan_315', '46_pan_90', '46_pan_135', '46_pan_180', '46_pan_225'};
% rotations = {'0', '45', '270', '315', '90', '135', '180', '225'};
rotations = [0, 45, 270, 315, 90, 135, 180, 225];

% file_rotations = [0, 45, 270, 315, 90, 135, 180, 225, 360];

%cli_ids   = {'38_pan_225'};
%rotations = {'225'};

ap_rotations = [-90, 90, -90, -90, 90];

%measurement_points = [6, 8, 10, 2, 4, 75, 77, 79, 12];
%measurement_points = [92, 94, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 100, 102];
measurement_points = 1:31;
% measurement_points = [20];
%measurement_points = [1];
ap_ids       = ap_order;
%ap_ids = [43];

tof_estimated = nan(1, size(measurement_points, 2));
tof_gt = nan(1, size(measurement_points, 2));
tof_errors = nan(1, size(measurement_points, 2));
aoa_estimated = nan(1, size(measurement_points, 2));
aoa_gt = nan(1, size(measurement_points, 2));
aoa_errors = nan(1, size(measurement_points, 2));

estimated_direct_path_angle    = nan(size(measurement_points, 2), size(ap_ids, 2));
estimated_direct_path_distance = nan(size(measurement_points, 2), size(ap_ids, 2));

% We only keep the best 
calculated_aoa = nan(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));
calculated_distance = nan(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));

% All csi results
csi_raw = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));
magnitudes_raw = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));

% All tof results
tof_raw = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));

% All mD-Track
azimut_raw    = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));
elevation_raw = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));
power_raw     = cell(size(measurement_points, 2), size(ap_ids, 2), size(rotations, 2));

for point_id=27
    
    measurement_point = measurement_points(point_id);
    path = ['raw_data/imdea/aoa_and_ftm/' num2str(measurement_point) '/'];
    
    for ap_id=ap_ids(1)

        % Process all the rotations and save the interesting data
        all_azimuts = nan(1, 0);
        all_powers  = nan(1, 0);
        
        max_magnitude = 0;
        
        
        % Now we iterate over all the rotations and take the one whose
        % magnitude is higher
        if (point_id == [13,15,16,19,20])
            cli_ids = cli_ids2;
        else
            cli_ids = cli_ids1;
        end
        for cli_id=1:size(cli_ids,2)
            
            Az_estimated = nan;
            distance = nan;
            
            % Read the raw data
            filename = ['cli_' num2str(ap_id) '_ap_' cli_ids{cli_id} '_AOA.txt'];

            % Super cutre

            % Load the oscillator calibration
            load(['processed_data/oscillators/oscillator_router_' num2str(ap_id) '.mat'])

            full_path = [path filename];

            [magnitudes, phases, ~] = Parse_raw_data_new_format(full_path);

            % The data is wrong
            if size(phases, 1) == 1 || size(phases, 1) < num_samples

                disp(['[WRONG] ' filename])
                continue
            end

            % Check if the magnitude is higher than the highes previous
            % rotation
            
            %disp(['[Power] File: ' num2str(measurement_point) '/' filename ' sum_magnitudes: (' num2str(sum(sum(magnitudes, 1))) ')'])

            % We don't even need to check other rotations for this pair
%             if max_magnitude < sum(sum(magnitudes, 1))
%                
%                 max_magnitude = sum(sum(magnitudes, 1));
%             else
%                 
%                 continue
%             end
            
            broken = 0;

            % Clean the data
            pre_channel = zeros(6, 6, num_samples);

            % We go up to 30 instead of 32
            % so that 31 and 32 are disabled
            % since they return random data
            for jj=1:30

                a = phases(:, jj);
                a = a*2*pi/1024;
                % move to complex
                a = exp(1i*a);

                converging_limit = 50;
                converging_retries = 0;
                converged = 0;
                while converged == 0
                    
                    try
                        [a, phase_offset_0, converged] = Sanitize(a);
%                         figure; plot(angle(a));
                        %broken = 0;
                    catch
                        disp(['Converging error on file ' filename])
                        %broken = 1;
                        %break
                    end
                    
                    if converging_retries == converging_limit
                       
                        break
                    end
                    
                    converging_retries = converging_retries + 1;
                    
                end
                
                if converging_retries == converging_limit
                    disp(['Converging threshold reached, ignoring ' filename])
                    continue
                end

                % Remove oscilator
                a = a/exp(1i*antenna_oscilator_phases(antenna_positions == jj));

                [row,col] = find(antenna_positions == jj);
                pre_channel(row, col, :) = a(1:num_samples, :);
            end
            
%             if broken == 1
%                 continue
%             end

            csi_data = pre_channel;
            csi_data = sum(csi_data,3)/num_samples;

            % apply mD-trackestimated_azimuts_main_path_angle
            [Az_estimated, El_estimated, att] = mD_track_2D(csi_data.', cb_az, cb_el);

            Az_estimated = rad2deg(theta_az(Az_estimated));
            El_estimated = rad2deg(theta_el(El_estimated));

            % Power
            power = abs(att).^2;
            
            all_azimuts = Az_estimated;
            all_powers  = power;
            
            %disp(['[AoA] File: ' num2str(measurement_point) '/' filename ', GT: ' num2str(angles_RX(measurement_point, ap_order==ap_id)) ' estimated: (' num2str(Az_estimated) ', ' num2str(El_estimated) ')'])
            
            % FTM
            filename = ['cli_' num2str(ap_id) '_ap_' cli_ids{cli_id} '_TOF.txt'];
            full_path = [path filename];
            ftm_times = Parse_raw_data_ftm(full_path);

            % Create a histogram
            distances = zeros(size(ftm_times, 1), 1);

            for i=1:size(ftm_times, 1)

                % Calculate the distance in meters
                T1 = ftm_times(i, 1);
                T2 = ftm_times(i, 2);
                T3 = ftm_times(i, 3);
                T4 = ftm_times(i, 4);

                dist = 3e8 * (((T4-T1)-(T3-T2))*1e-12)/2;

                distances(i, 1) = dist;
            end
            
            distance = median(distances) - antenna_ftm_offset;
            %disp(['[ToF] File: ' num2str(measurement_point) '/' filename ', GT: ' num2str(distances_RX(measurement_point, ap_order==ap_id)) ' estimated: (' num2str(distance) ')'])
            
            % Draw the measurement point in other color
            %plot(points_x(measurement_point), points_y(measurement_point), "g*", "LineWidth", 4);

%             x = RX_x(ap_order==ap_id);
%             y = RX_y(ap_order==ap_id);

            % Draw the ground truth line
            %plot([x points_x(measurement_point)],[y points_y(measurement_point)], ':r', "LineWidth", 2)
            
            % Save the data for the mat files
            calculated_aoa(point_id, ap_order==ap_id, cli_id) = Az_estimated(1);
            calculated_distance(point_id, ap_order==ap_id, cli_id) = distance;

            % All csi results
            csi_raw{point_id, ap_order==ap_id, cli_id} = csi_data;
            magnitudes_raw{point_id, ap_order==ap_id, cli_id} = magnitudes;
            
            % All tof results
            tof_raw{point_id, ap_order==ap_id, cli_id} = ftm_times;
            
            azimut_raw{point_id, ap_order==ap_id, cli_id}    = Az_estimated;
            elevation_raw{point_id, ap_order==ap_id, cli_id} = El_estimated;
            power_raw{point_id, ap_order==ap_id, cli_id}     = power;
        end

        % Save the error
%         aoa_estimated(1, point_id) = Az_estimated(1);
%         aoa_gt(1, point_id)        = angles_RX(measurement_point, ap_order==ap_id);
%         tof_estimated(1, point_id) = distance;
%         tof_gt(1, point_id)        = distances_RX(measurement_point, ap_order==ap_id);
        
        % Save the estimation
        estimated_direct_path_angle(point_id, ap_order==ap_id)    = Az_estimated(1);
        estimated_direct_path_distance(point_id, ap_order==ap_id) = distance;
        
        %disp(['[Errors] File: ' num2str(measurement_point) '/' filename ', aoa: ' num2str(abs(angles_RX(measurement_point, ap_order==ap_id) - Az_estimated(1))) ' tof: (' num2str(abs(distances_RX(measurement_point, ap_order==ap_id) - distance)) ')'])

        % Draw the paths
        line_width = 1;
        count = 1;
        
%         for az=all_azimuts
%             
%             %L is the length
%             L = 10;
%             alpha = -1*az-ap_rotations(ap_order==ap_id);
% 
%             %angle is alpha
%             x2=x+(L*cosd(alpha));
%             y2=y+(L*sind(alpha));
%             
%             if count == 1
%                 linespeck = '-';
%             else
%                 linespeck = ':';
%             end
% 
%             %plot([x x2],[y y2], linespeck, "LineWidth",  line_width)
% 
%             count = count +1;
%         end
    end
end

title(['All'])

% Sort correctly
rotations = [0, 45, 270, 315, 90, 135, 180, 225];

% [B, I] = sort(rotations);

% Matrix
% calculated_aoa(:, :, [1:8])      = calculated_aoa(:, :, [I]);
% calculated_distance(:, :, [1:8]) = calculated_distance(:, :, [I]);

% Cell
% azimut_raw(:, :, [1:8])     = azimut_raw(:, :, [I]);
% csi_raw(:, :, [1:8])        = csi_raw(:, :, [I]);
% elevation_raw(:, :, [1:8])  = elevation_raw(:, :, [I]);
% magnitudes_raw(:, :, [1:8]) = magnitudes_raw(:, :, [I]);
% power_raw(:, :, [1:8])      = power_raw(:, :, [I]);
% tof_raw(:, :, [1:8])        = tof_raw(:, :, [I]);
mkdir("mat_files/imdea/HF")
mkdir("mat_files/imdea/HF/CSI")
mkdir("mat_files/imdea/HF/FTM")

out_path = 'mat_files/imdea/HF/';

% save([out_path 'CSI/csi_imdea.mat'], 'calculated_aoa', 'csi_raw', 'magnitudes_raw', 'azimut_raw', 'elevation_raw', 'power_raw');
% save([out_path 'FTM/ftm_imdea.mat'], 'calculated_distance', 'tof_raw');


% Get the point_id, aoa_error and tof_error on the same matrix to find
% outliers

% aoa_errors = abs(aoa_estimated - aoa_gt  - median(aoa_estimated - aoa_gt, 'omitnan'));
% tof_errors = abs(tof_estimated - tof_gt);
% 
% figure
% cdfplot(aoa_errors)
% title('AoA error')
% 
% figure
% cdfplot(tof_errors)
% title('ToF error')

% aoa_errors = abs(aoa_estimated - aoa_gt);
% tof_errors = abs(tof_estimated - tof_gt);
% 
% figure
% cdfplot(aoa_errors)
% title('[raw] aoa')
% 
% figure
% cdfplot(tof_errors)
% title('[raw] tof')

% debug = nan(3, size(measurement_points,2));
% debug(1, :) = measurement_points;
% debug(2, :) = aoa_errors;
% debug(3, :) = tof_errors;

%% Localize: Estimate the positions

%estimated_direct_path_angle(point_id, ap_order==ap_id)    = Az_estimated(1);
%estimated_direct_path_distance(point_id, ap_order==ap_id) = distance;

% estimated_direct_path_angle = estimated_direct_path_angle*-1+ median(estimated_direct_path_angle, 'omitnan');
% 
% for ap_id=ap_ids
%     for point_id=1:size(measurement_points, 2)
%         measurement_point = measurement_points(point_id);
%         
%         estimated_point_x_p(point_id, ap_order==ap_id) = sind(estimated_direct_path_angle(point_id,ap_order==ap_id))*estimated_direct_path_distance(point_id,ap_order==ap_id);
%         estimated_point_y_p(point_id, ap_order==ap_id) = cosd(estimated_direct_path_angle(point_id,ap_order==ap_id))*estimated_direct_path_distance(point_id,ap_order==ap_id);
% 
%         %figure, plot(estimated_point_x_p(id_point,id_router), estimated_point_y_p(id_point,id_router), "*")
% 
%         % Looking up
%         if (ap_id == 45 || ap_id == 47 || ap_id == 46)
%             estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id));
%             estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - (estimated_point_y_p(point_id,ap_order==ap_id)*(-1));
%         else % Looking down
%             estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id)*(-1));
%             estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - estimated_point_y_p(point_id,ap_order==ap_id);
%         end
%     end
% end
% 
% % Euclidean distance to true ground
% errors = sqrt( (estimated_point_x(:, ap_order==ap_ids) - points_x(measurement_points)).^2 + (estimated_point_y(:, ap_order==ap_ids) - points_y(measurement_points)).^2 );
% 
% figure;
% cdfplot(errors);
% title('Location errors room 43')
% 
% out_path = 'processed_data/aoa_and_ftm/';
% 
% 
% save([ out_path '43.mat'], 'aoa_errors', 'tof_errors', 'errors', 'debug');


%savefig(['plots/corridor/corridor_' num2str(measurement_points) '.fig'])

% %close all
% mag = nan(32, size(cli_ids, 2));
% 
% for point_id=1:size(measurement_points, 2)
% 
%     figure;
% 
%     point = measurement_points(point_id);
%     
%     mag(:, :) = magnitudes_aux(:, :, point_id);
%     bar(median(mag, 1))
%     title(['Point ' num2str(point)])
%     set(gca,'xticklabel', rotations)
%      
% end

% Compare estimated with ground truth

% Here we have the estimated path based on maximum power


% Create the ground truth matrix
%gt_matrix = nan(size(ap_ids, 2), size(measurement_points, 2));

%ap_ids       = [35,43,39,36,40];

%gt_matrix(1, :) = angles_RX(measurement_points, ap_order==35);
%gt_matrix(2, :) = angles_RX(measurement_points, ap_order==43);
%gt_matrix(3, :) = angles_RX(measurement_points, ap_order==39);
%gt_matrix(4, :) = angles_RX(measurement_points, ap_order==36);
%gt_matrix(5, :) = angles_RX(measurement_points, ap_order==40);


%gt_matrix

%estimated_azimuts_main_path

%diff = angles_RX(measurement_points, ap_order==43) - estimated_azimuts_main_path';

%figure; cdfplot(abs(diff-3))
%xlabel('Error in degrees')
%title(['Estimation errors for ' num2str(size(measurement_points,2)) ' points'])
