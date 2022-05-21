clc
clear all
close all

pwd_str = pwd;
cd ../../

addpath(genpath("auxiliar_functions"))
%% Indoor scenario
csi_indoor  = load('mat_files/indoor/HF/CSI/csi_indoor.mat');
ftm_indoor  = load('mat_files/indoor/HF/FTM/ftm_indoor.mat');

% Map data
indoor_rx                   = load('Maps/Map_indoor/RX_coordinates.mat');
indoor_coordinates          = load('Maps/Map_indoor/points_coordinates.mat');
indoor_true_angles_distance = load('Maps/Map_indoor/data_distance_angle_true.mat');

fun = @(m)sRGB_to_OSAUCS(m,true,true); % modified OSA-UCS (recommended)
[rgb,ucs] = maxdistcolor(6,fun);

rgb2 = viridis(6);
lf_color = rgb2(6, :);

rotations_legend = ["0", "45", "-45", "90", "-90", "135", "-135", "180"];

final_rgb(1, :) = rgb(1, :);
final_rgb(2, :) = rgb(2, :);
final_rgb(3, :) = rgb(2, :);
final_rgb(4, :) = rgb(3, :);
final_rgb(5, :) = rgb(3, :);
final_rgb(6, :) = rgb(4, :);
final_rgb(7, :) = rgb(4, :);
final_rgb(8, :) = rgb(5, :);
final_rgb(9, :) = rgb(6, :);

% rgb(6, :) = [];
% rgb = [rgb; rgb(4,:)];
% rgb = [rgb; rgb(3,:)];
% rgb = [rgb; rgb(2,:)];
% rgb = [rgb; lf_color];
% rgb(1, :) = [0, 0, 0]; % First color black


load('mat_files/indoor/HF/optimal_rotation_index.mat')

aoa_optimal_rotation = nan(110, 5);
ftm_optimal_rotation = nan(110, 5);

% Save the distance and FTM for the optimal rotation
for ii=1:size(indoor_coordinates.points_x, 1)
    
    % For each AP
    for jj=1:size(indoor_rx.RX_x, 2)

        kk = indoor_optimal_rotation_index(ii, jj);
        aoa_optimal_rotation(ii, jj) = csi_indoor.calculated_aoa(ii, jj, kk);
        ftm_optimal_rotation(ii, jj) = ftm_indoor.calculated_distance(ii, jj, kk);
        
    end
end

% We need to have the true_angles and true_distances copied for each
% rotation
true_angles    = repmat(indoor_true_angles_distance.angles_RX, [1 1 8]);
true_distances = repmat(indoor_true_angles_distance.distances_RX, [1 1 8]);

% Angle for each point, rotation and AP
% figure
% cdfplot(abs(csi_indoor.calculated_aoa(:) - true_angles(:)));
% xlabel('Error in degrees')
% title('[Indoor] AoA errors for all points, all APs, all rotations')

% What we need to know now is which is the optimal rotation towards an AP
% and then check what would be the error if we are rotated with respect to
% that rotation

% Get the true angles from each point to each AP
diff_RX_x = indoor_coordinates.points_x - indoor_rx.RX_x;
diff_RX_y = indoor_coordinates.points_y - indoor_rx.RX_y;
indoor_true_angles = rad2deg(atan(diff_RX_x./diff_RX_y));

ap_order  = [45,44,47,46,43];
rotations = sort([0, 45, 270, 315, 90, 135, 180, 225]);
rotations_legend = ["Oracle (geo.)", "45", "-45", "90", "-90", "135", "-135", "180"];
rotations_index = [1 2 8 3 7 4 6 5];
% We map it from [-90, 90] to [0, 360], we multiply by -1
% so it is from the POV of the point
indoor_true_angles_360 = mod(-1*indoor_true_angles, 360);

% We also need to rotate APs 45, 46 and 47
indoor_true_angles_360(:, logical([1, 0, 1, 1, 0])) = mod(indoor_true_angles_360(:, logical([1, 0, 1, 1, 0])) + 180, 360);

% Now we need to round to the closest multiple of 45
indoor_true_angles_360 = round( indoor_true_angles_360 / 45 ) * 45;

% And 360 is equal to 0
indoor_true_angles_360(indoor_true_angles_360==360) = 0;

% Convert the angles to indexes in the matrix
for ii=1:numel(indoor_true_angles_360)
   
    indoor_true_angles_360(ii) = find(rotations==indoor_true_angles_360(ii));
end

% Now indoor_true_angles_360 has which slice of the 3rd dimension of the
% matrix I must take for optimal rotation

errors_per_deviation = nan(size(indoor_coordinates.points_x, 1), size(indoor_rx.RX_x, 2), size(rotations, 2));

% For each point
for ii=1:size(indoor_coordinates.points_x, 1)
    
    % For each AP
    for jj=1:size(indoor_rx.RX_x, 2)

        % Angle for that point and AP for the optimal rotation
        optimal_angle_index = indoor_true_angles_360(ii, jj);
        optimal_angle_gt = rotations(optimal_angle_index);
        
        for rotation_index=1:size(rotations, 2)
            
            current_rotation = rotations(rotation_index);
            offset           = mod(optimal_angle_gt - current_rotation, 360);
            
            save_index = rotations==offset;
            
            errors_per_deviation(ii, jj, save_index) = abs(csi_indoor.calculated_aoa(ii, jj, rotations==current_rotation) - true_angles(ii, jj, rotations==current_rotation));
        end
    end
end

% figure
% hold on
% 
% markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
% 
% for rotation_index=1:size(rotations, 2)
%     
%     aux = errors_per_deviation(:, :, rotation_index);
%     h = cdfplot(aux(:));
%     
%     set(h,'Marker', markers{rotation_index});
% end
% 
% hold off
% xlabel('Error in degrees')
% title('[Indoor] AoA errors if we chose the wrong rotation by N degrees')
% legend(string(rotations))

% LF AoA indoor
load('mat_files/mD_track/LF/Data_3D_Direct_Path.mat')

% Calibrate LF
angles_dp(1:74,2) = angles_dp(1:74,2) - 3.5;
angles_dp(26:38,5) = angles_dp(26:38,5) - 5;
angles_dp(1:74,3) = angles_dp(1:74,3) - 7;

%angles_dp(:, [1 4]) = angles_dp(:, [4 1]);

% Same as above but divided in LOS and NLOS
load('mat_files/indoor_map/Grid_LOS.mat')

load('mat_files/indoor/HF/CSI/aoa_offset.mat');

errors_per_deviation_los  = nan(size(indoor_coordinates.points_x, 1), size(indoor_rx.RX_x, 2), size(rotations, 2));
errors_per_deviation_nlos = nan(size(indoor_coordinates.points_x, 1), size(indoor_rx.RX_x, 2), size(rotations, 2));

% For each point
for ii=1:size(indoor_coordinates.points_x, 1)
    
    % For each AP
    for jj=1:size(indoor_rx.RX_x, 2)

        % Angle for that point and AP for the optimal rotation
        optimal_angle_index = indoor_true_angles_360(ii, jj);
        optimal_angle_gt = rotations(optimal_angle_index);
        
        for rotation_index=1:size(rotations, 2)
            
            current_rotation = rotations(rotation_index);
            offset           = mod(optimal_angle_gt - current_rotation, 360);
            
            save_index = rotations==offset;
            
            if LOS_RX(ii, jj) == 1
                errors_per_deviation_los(ii, jj, save_index)  = abs(csi_indoor.calculated_aoa(ii, jj, rotations==current_rotation) - true_angles(ii, jj, rotations==current_rotation) - aoa_offset(jj));
            else
                errors_per_deviation_nlos(ii, jj, save_index) = abs(csi_indoor.calculated_aoa(ii, jj, rotations==current_rotation) - true_angles(ii, jj, rotations==current_rotation) - aoa_offset(jj));
            end
        end
    end
end

% LOS
fig = figure;
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

rotations_legend = ["Oracle (geo.)", "45", "-45", "90", "-90", "135", "-135", "180"];

%errors_per_deviation_los(isnan(errors_per_deviation_los)) = 800;
%errors_per_deviation_nlos(isnan(errors_per_deviation_nlos)) = 800;

count = 1;
for rotation_index=rotations_index
    
    aux = errors_per_deviation_los(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    if count == 3 || count == 5 || count == 7
        set(h, 'LineStyle', '--', 'color', final_rgb(count, :), 'LineWidth', 1);
    else
        set(h, 'LineStyle', '-', 'color', final_rgb(count, :), 'LineWidth', 1);
    end
    count = count + 1;
end

error = abs(indoor_true_angles_distance.angles_RX(LOS_RX) - -1*angles_dp(LOS_RX));
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', final_rgb(9, :), 'LineWidth', 1);

hold off
%% AoA LOS
xlabel('(a) Angle error [°] LOS','fontweight','bold')
ylabel('ECDF')
title('')
legend([rotations_legend 'LF'], 'Location', "SouthEast")
xlim([0 160])
yticks([0:0.2:1])


proportion = 4/3;
x_length = 1.5;
y_length = x_length*proportion;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])

MagInset_Grid_On(fig, -1, [0 15 0 1], [60 160 0.15 0.6], {'NW','NW';'SE','SE'}, [0:5:15], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/rotation_comparison/AoA_LOS_Indoor.tikz');



% NLOS
fig = figure
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

count = 1;
for rotation_index=rotations_index
    
    aux = errors_per_deviation_nlos(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    if count == 3 || count == 5 || count == 7
        set(h, 'LineStyle', '--', 'color', final_rgb(count, :), 'LineWidth', 1);
    else
        set(h, 'LineStyle', '-', 'color', final_rgb(count, :), 'LineWidth', 1);
    end
    count = count + 1;
end

error = abs(indoor_true_angles_distance.angles_RX(~LOS_RX) - -1*angles_dp(~LOS_RX));
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', final_rgb(9, :), 'LineWidth', 1);

hold off
%% AoA NLOS
xlabel('(b) Angle error [°] NLOS','fontweight','bold')
ylabel('')

title('')
legend([rotations_legend 'LF'], 'Location', "SouthEast")
xlim([0 160])
yticks([0:0.2:1])

fleg = legend('figure()');
set(fleg,'visible','off')


set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
set(gca,'YTicklabel',[])

MagInset_Grid_On(fig, -1, [0 15 0 1], [60 160 0.15 0.6], {'NW','NW';'SE','SE'}, [0:5:15], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/rotation_comparison/AoA_NLOS_Indoor.tikz');


% FTM for LOS and NLOS
errors_per_deviation_los  = nan(size(indoor_coordinates.points_x, 1), size(indoor_rx.RX_x, 2), size(rotations, 2));
errors_per_deviation_nlos = nan(size(indoor_coordinates.points_x, 1), size(indoor_rx.RX_x, 2), size(rotations, 2));

% For each point
for ii=1:size(indoor_coordinates.points_x, 1)
    
    % For each AP
    for jj=1:size(indoor_rx.RX_x, 2)

        % Angle for that point and AP for the optimal rotation
        optimal_angle_index = indoor_true_angles_360(ii, jj);
        optimal_angle_gt = rotations(optimal_angle_index);
        
        for rotation_index=1:size(rotations, 2)
            
            current_rotation = rotations(rotation_index);
            offset           = mod(optimal_angle_gt - current_rotation, 360);
            
            save_index = rotations==offset;
            
            if LOS_RX(ii, jj) == 1
                errors_per_deviation_los(ii, jj, save_index)  = abs(ftm_indoor.calculated_distance(ii, jj, rotations==current_rotation) - true_distances(ii, jj, rotations==current_rotation));
            else
                errors_per_deviation_nlos(ii, jj, save_index) = abs(ftm_indoor.calculated_distance(ii, jj, rotations==current_rotation) - true_distances(ii, jj, rotations==current_rotation));
            end
        end
    end
end

% LOS
fig = figure
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

count = 1;
for rotation_index=rotations_index
    
    aux = errors_per_deviation_los(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    if count == 3 || count == 5 || count == 7
        set(h, 'LineStyle', '--', 'color', final_rgb(count, :), 'LineWidth', 1);
    else
        set(h, 'LineStyle', '-', 'color', final_rgb(count, :), 'LineWidth', 1);
    end
    count = count + 1;
end

% LF FTM indoor
load('mat_files/FTM/LF/estimated_distance.mat')

error = abs(indoor_true_angles_distance.distances_RX(LOS_RX) - estimated_distance(LOS_RX));
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', final_rgb(9, :), 'LineWidth', 1);

hold off
%% ToF LOS
xlabel('(c) ToF error [m] LOS','fontweight','bold')
ylabel('')

title('')
legend([rotations_legend 'LF'], 'Location', "SouthEast")
xlim([0 35])
yticks([0:0.2:1])

fleg = legend('figure()');
set(fleg,'visible','off')


set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
set(gca,'YTicklabel',[])

MagInset_Grid_On(fig, -1, [0 1 0 1], [15 32 0.15 0.6], {'NW','NW';'SE','SE'}, [0:0.5:1], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/rotation_comparison/ToF_LOS_Indoor.tikz');


% NLOS
fig = figure
hold on

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

count = 1;
for rotation_index=rotations_index
    
    aux = errors_per_deviation_nlos(:, :, rotation_index);
    h = cdfplot(aux(:));
    
    if count == 3 || count == 5 || count == 7
        set(h, 'LineStyle', '--', 'color', final_rgb(count, :), 'LineWidth', 1);
    else
        set(h, 'LineStyle', '-', 'color', final_rgb(count, :), 'LineWidth', 1);
    end
    count = count + 1;
end

error = abs(indoor_true_angles_distance.distances_RX(~LOS_RX) - estimated_distance(~LOS_RX));
h = cdfplot(error(:));
set(h, 'LineStyle', '-', 'color', final_rgb(9, :), 'LineWidth', 1);

hold off
%% ToF NLOS
xlabel('(d) ToF error [m] NLOS','fontweight','bold')
ylabel('')

title('')
legend([string(rotations) 'LF'], 'Location', "SouthEast")
xlim([0 35])
yticks([0:0.2:1])

fleg = legend('figure()');
set(fleg,'visible','off')

set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
set(gca,'YTicklabel',[])


MagInset_Grid_On(fig, -1, [0 1 0 1], [15 32 0.15 0.6], {'NW','NW';'SE','SE'}, [0:0.5:1], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/rotation_comparison/ToF_NLOS_Indoor.tikz');
[status,result] = system('bash plots/final_plots/rotation_comparison/modify_plots_all.sh');

cd(pwd_str)

%% Outdoor scenario
% csi_outdoor = load('mat_files/outdoor/HF/CSI/csi_outdoor.mat');
% csi_indoor  = load('mat_files/outdoor/HF/FTM/ftm_outdoor.mat');
% 
% % Map data
% outdoor_rx                    = load('mat_files/outdoor_map/RX_coordinates.mat');
% outdoor_coordinates           = load('mat_files/outdoor_map/points_coordinates.mat');
% outdoor_true_angles_distances = load('mat_files/outdoor_map/data_distance_angle_true.mat');
% 
% % We need to have the true_angles and true_distances copied for each
% % rotation
% true_angles    = repmat(outdoor_true_angles_distances.angles_RX, [1 1 8]);
% true_distances = repmat(outdoor_true_angles_distances.distances_RX, [1 1 8]);
% 
% % Angle for each point, rotation and AP
% figure
% cdfplot(abs(csi_outdoor.calculated_aoa(:) - true_angles(:)));
% xlabel('Error in degrees')
% title('[Outdoor] AoA errors for all points, all APs, all rotations')