clc
clear all
close all

%%
load('mat_files/Coordinates/data_coordinates_hf.mat')

mask = ~isnan(estimated_point_x_hf);

result = false(110, 5);

for ii=1:8

    result = result | mask(:, :, ii);
end

coverage = result;

result = sum(result, 2);

ap_coverage = zeros(4, 1);

for ii=1:4

    ap_coverage(ii) = sum(result == ii);
end

labels = {'1 AP', '2 APs','3 APs', '4 APs'};
pie(ap_coverage)
lgd = legend(labels);

%% The coverage figure
load('mat_files/indoor_map/RX_coordinates.mat');
load('mat_files/indoor_map/points_coordinates.mat');

fun = @(m)sRGB_to_OSAUCS(m,true,true); % recommended OSA-UCS
rgb = maxdistcolor(5,fun);

uiopen('mat_files/indoor_map/map.fig',1);
hold on
%pause(3)

% Remove points and APs
plot(points_x, points_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);
plot(RX_x, RX_y, "*", "LineWidth", 4, 'color', [1, 1, 1]);

markers = {'d','o','.','+','x','*','s','^','v','>','<','p','h'};
sizes   = [2, 3, 2, 2, 2]*5;

for id_ap=1:5
    
    plot(RX_x(id_ap), RX_y(id_ap), markers{id_ap}, 'MarkerSize', sizes(id_ap), 'color', rgb(id_ap, :));

    
    for id_point=1:110
        
        if coverage(id_point, id_ap) == 1
           
            plot(points_x(id_point), points_y(id_point), markers{id_ap}, "LineWidth", 1, 'color', rgb(id_ap, :), 'MarkerSize', sizes(id_ap));

        end
        
    end
end

title('Coverage map for HF')

%pause(5)

%matlab2tikz('plots/final_plots/Coverage/Coverage.tikz');


%% To JUJA coverage plot

% load('mat_files/indoor_map/data_distance_angle_true.mat', 'AP_labels_HF');
% 
% % Each cm
% granularity = 100;
% 
% coverage_heatmap = nan(15*granularity, 20*granularity, 5);
% 
% % Put the true values
% for id_ap=1:5
%        
%     for id_point=1:110
%         
%         x = int64(points_x(id_point)*100);
%         y = int64(points_y(id_point)*100);
%         
%         if coverage(id_point, id_ap) == 1
%            
%             coverage_heatmap(x, y, id_ap) = 1;
%         else
%             coverage_heatmap(x, y, id_ap) = 0;
%         end
%     end
% end
% 
% %heatmap(coverage_heatmap(:, :, 1))
% 
% ap_5 = coverage_heatmap(:, :, 2);
% 
% %filled = fillmissing(ap_5,'previous');
% 
% %heatmap(ap_1)
% 
% BW = imbinarize(ap_5);
% 
% % Expand each measurement point
% SE = strel('square',100);
% BW2 = imdilate(BW,SE);
% 
% % "Join" the gaps
% se = strel('square',300);
% closeBW = imclose(BW2, se);
% 
% 
% imshow(closeBW, []);
% 





