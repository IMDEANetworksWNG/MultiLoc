clc
clear all
close all

load('data_distance_angle_true.mat')
load('points_coordinates.mat')
load('RX_coordinates.mat')

%% Generate the map
uiopen('/home/imdea/Documents/MATLAB/Music_Mikrotik/mat_files/outdoor_map/outdoor_map.fig',1);
hold on

%% Configuration for md-Track
% Load antenna data
load('/home/imdea/Documents/MATLAB/Music_Mikrotik/processed_data/antennas_mikrotik.mat')

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
ap_order = AP_labels_HF;
cli_ids   = {'37_pan_0', '37_pan_45', '37_pan_270', '37_pan_315', '38_pan_90', '38_pan_135', '38_pan_180', '38_pan_225'};
rotations = {'0', '45', '270', '315', '90', '135', '180', '225'};
file_rotations = [0, 45, 270, 315, 90, 135, 180, 225, 360];

%cli_ids   = {'38_pan_225'};
%rotations = {'225'};

%ap_rotations = [135, 225, 315, 180];

measurement_points = [111:126];
ap_ids       = [43:46];
%ap_ids = [46];
%measurement_points = [119 118];


% tof_estimated = nan(1, size(measurement_points, 2));
% %tof_gt = nan(1, size(measurement_points, 2));
% tof_errors = nan(1, size(measurement_points, 2));
% aoa_estimated = nan(1, size(measurement_points, 2));
% %aoa_gt = nan(1, size(measurement_points, 2));
% aoa_errors = nan(1, size(measurement_points, 2));

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

for point_id=1:size(measurement_points, 2)
    
    measurement_point = measurement_points(point_id);
    path = ['raw_data/outdoor/' num2str(measurement_point) '/'];
    
    for ap_id=ap_ids

        % Process all the rotations and save the interesting data
        all_azimuts = nan(1, 0);
        all_powers  = nan(1, 0);
        
        max_magnitude = 0;
        
        % Now we iterate over all the rotations and take the one whose
        % magnitude is higher
        for cli_id=1:size(cli_ids,2)
            
            Az_estimated = nan;
            distance = nan;
            
            % Read the raw data
            filename = ['cli_' num2str(ap_id) '_ap_' cli_ids{cli_id} '_AOA.txt'];

            %disp(filename)
            
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
            
            %disp(['[AoA] File: ' num2str(measurement_point) '/' filename ', GT: ' num2str(true_angles(measurement_point, ap_order==ap_id)) ' estimated: (' num2str(Az_estimated) ', ' num2str(El_estimated) ')'])
            
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
            %disp(['[ToF] File: ' num2str(measurement_point) '/' filename ', GT: ' num2str(true_distances(measurement_point, ap_order==ap_id)) ' estimated: (' num2str(distance) ')'])
            
            % Draw the measurement point in other color
            plot(points_x(point_labels==measurement_point), points_y(point_labels==measurement_point), "g*", "LineWidth", 4);

            x = RX_x(ap_order==ap_id);
            y = RX_y(ap_order==ap_id);

            % Draw the ground truth line
            %plot([x points_x(labels==measurement_point)],[y points_y(labels==measurement_point)], ':r', "LineWidth", 2)
            
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
        aoa_estimated(point_id, ap_order==ap_id) = Az_estimated(1);
        %aoa_gt(ap_order==ap_id, point_id)        = angles_RX(point_labels==measurement_point, ap_order==ap_id);% + ap_rotations(ap_order==ap_id);
        tof_estimated(point_id, ap_order==ap_id) = distance;
        %tof_gt(ap_order==ap_id, point_id)        = distances_RX(point_labels==measurement_point, ap_order==ap_id);
        
        % Save the estimation
        estimated_direct_path_angle(point_id, ap_order==ap_id)    = Az_estimated(1);
        estimated_direct_path_distance(point_id, ap_order==ap_id) = distance;
        
        %disp(['[Errors] File: ' num2str(measurement_point) '/' filename ', aoa: ' num2str(abs(true_angles(point_labels==measurement_point, ap_order==ap_id) - Az_estimated(1))) ' tof: (' num2str(abs(true_distances(point_labels==measurement_point, ap_order==ap_id) - distance)) ')'])
        
        % Draw the paths
        line_width = 1;
        count = 1;
                
        for az=all_azimuts
            
            %L is the length
            L = 10;
            alpha = -1*az;%+ap_rotations(ap_order==ap_id);

            %angle is alpha
            x2=x+(L*cosd(alpha));
            y2=y+(L*sind(alpha));
            
            if count == 1
                linespeck = '-';
            else
                linespeck = ':';
            end

            plot([x x2],[y y2], linespeck, "LineWidth",  line_width)

            count = count +1;
        end
    end
end

title(['Outdoor'])

% Sort correctly
rotations = [0, 45, 270, 315, 90, 135, 180, 225];

[B, I] = sort(rotations);

% Matrix
calculated_aoa(:, :, [1:8])      = calculated_aoa(:, :, [I]);
calculated_distance(:, :, [1:8]) = calculated_distance(:, :, [I]);

% Cell
azimut_raw(:, :, [1:8])     = azimut_raw(:, :, [I]);
csi_raw(:, :, [1:8])        = csi_raw(:, :, [I]);
elevation_raw(:, :, [1:8])  = elevation_raw(:, :, [I]);
magnitudes_raw(:, :, [1:8]) = magnitudes_raw(:, :, [I]);
power_raw(:, :, [1:8])      = power_raw(:, :, [I]);
tof_raw(:, :, [1:8])        = tof_raw(:, :, [I]);

% out_path = 'mat_files/outdoor/HF/';
% 
% save([out_path 'CSI/csi_outdoor.mat'], 'calculated_aoa', 'csi_raw', 'magnitudes_raw', 'azimut_raw', 'elevation_raw', 'power_raw');
% save([out_path 'FTM/ftm_outdoor.mat'], 'calculated_distance', 'tof_raw');

% for ap_id=ap_ids
%     
%     aoa_errors = abs(angdiff(aoa_estimated(:, ap_order==ap_id), angles_RX(:, ap_order==ap_id))  - median(angdiff(aoa_estimated(:, ap_order==ap_id), angles_RX(:, ap_order==ap_id)), 'omitnan'));
%     tof_errors = abs(tof_estimated(:, ap_order==ap_id) - distances_RX(:, ap_order==ap_id));
% 
%     figure
%     subplot(1,2,1);
%     cdfplot(aoa_errors)
%     title(['AoA error ap ' num2str(ap_id)])
% 
%     subplot(1,2,2);
%     cdfplot(tof_errors)
%     title(['ToF error ' num2str(ap_id)])
% 
%     eval(['aoa_errors_'  num2str(ap_id) ' = aoa_errors;']);
%     eval(['tof_errors_'  num2str(ap_id) ' = tof_errors;']);
%     
% end

% debug = nan(3, size(measurement_points,2));
% debug(1, :) = measurement_points;
% debug(2, :) = aoa_errors;
% debug(3, :) = tof_errors;

%% Localize: Estimate the positions

%estimated_direct_path_angle(point_id, ap_order==ap_id)    = Az_estimated(1);
%estimated_direct_path_distance(point_id, ap_order==ap_id) = distance;

estimated_direct_path_angle = estimated_direct_path_angle - median(estimated_direct_path_angle, 'omitnan');

%ap_ids = 43;
%measurement_points = 111;
% https://en.wikipedia.org/wiki/Rotation_of_axes
for ap_id=ap_ids
    for point_id=1:size(measurement_points, 2)
        measurement_point = measurement_points(point_id);
        
        estimated_point_x_p(point_id, ap_order==ap_id) = sind(estimated_direct_path_angle(point_id,ap_order==ap_id))*estimated_direct_path_distance(point_id,ap_order==ap_id);
        estimated_point_y_p(point_id, ap_order==ap_id) = cosd(estimated_direct_path_angle(point_id,ap_order==ap_id))*estimated_direct_path_distance(point_id,ap_order==ap_id);

        % Rotate respect the AP
        x_center = RX_x(ap_order==ap_id);
        y_center = RX_y(ap_order==ap_id);
                
        if (ap_id == 43)

            estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id));
            estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - (estimated_point_y_p(point_id,ap_order==ap_id)*(-1));
            
            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_order==ap_id)));
            theta = deg2rad(45);
            v = [estimated_point_x(point_id,ap_order==ap_id);estimated_point_y(point_id,ap_order==ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x(point_id,ap_order==ap_id) = vo(1,:);
            estimated_point_y(point_id,ap_order==ap_id) = vo(2,:);
            
        elseif (ap_id == 44)
            
            estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id)*(-1));
            estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - (estimated_point_y_p(point_id,ap_order==ap_id));
            
            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_order==ap_id)));
            theta = deg2rad(-45);
            v = [estimated_point_x(point_id,ap_order==ap_id);estimated_point_y(point_id,ap_order==ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;   
            vo = so + center; 
            estimated_point_x(point_id,ap_order==ap_id) = vo(1,:);
            estimated_point_y(point_id,ap_order==ap_id) = vo(2,:);
            
        elseif (ap_id == 45)
            
            estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id)*(-1));
            estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - (estimated_point_y_p(point_id,ap_order==ap_id));
            
            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_order==ap_id)));
            theta = deg2rad(35);
            v = [estimated_point_x(point_id,ap_order==ap_id);estimated_point_y(point_id,ap_order==ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x(point_id,ap_order==ap_id) = vo(1,:);
            estimated_point_y(point_id,ap_order==ap_id) = vo(2,:);
            
        elseif (ap_id == 46)
            
            estimated_point_x(point_id,ap_order==ap_id) = RX_x(ap_order==ap_id) - (estimated_point_x_p(point_id,ap_order==ap_id)*(-1));
            estimated_point_y(point_id,ap_order==ap_id) = RX_y(ap_order==ap_id) - (estimated_point_y_p(point_id,ap_order==ap_id));
            
            % Now rotate the point
            center = repmat([x_center; y_center], 1, length(estimated_point_x(point_id,ap_order==ap_id)));
            theta = deg2rad(-90);
            v = [estimated_point_x(point_id,ap_order==ap_id);estimated_point_y(point_id,ap_order==ap_id)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            s = v - center;
            so = R*s;
            vo = so + center; 
            estimated_point_x(point_id,ap_order==ap_id) = vo(1,:);
            estimated_point_y(point_id,ap_order==ap_id) = vo(2,:);
        end
        
        %plot(estimated_point_x(point_id, ap_order==ap_id), estimated_point_y(point_id, ap_order==ap_id), "*")

    end
end
% 
% % Remove the median
% 
% % Euclidean distance to true ground
% for ap_id=ap_ids
% 
%     errors = abs(sqrt( (estimated_point_x(:, ap_order==ap_id) - points_x(point_labels==measurement_points)').^2 + (estimated_point_y(:, ap_order==ap_id) - points_y(point_labels==measurement_points)').^2 ));
%     
%     figure;
%     cdfplot(errors);
%     title(['Location errors for AP ' num2str(ap_id)])
% 
%     eval(['location_errors_'  num2str(ap_id) ' = errors;']);
% 
% end



% out_path = 'processed_data/outdoor/';
% 
% for ap_id=ap_ids
% 
%     save([ out_path num2str(ap_id) '.mat'], ['aoa_errors_' num2str(ap_id)], ['tof_errors_' num2str(ap_id)], ['location_errors_' num2str(ap_id)]);
% end
