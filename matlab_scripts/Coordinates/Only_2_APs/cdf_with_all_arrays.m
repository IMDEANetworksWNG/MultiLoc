clc
clear all
close all

%%
rgb = viridis(5);
pwd_str = pwd;
cd ../../../
addpath(genpath("auxiliar_functions"))
fig = figure;

load('matlab_scripts/Coordinates/Only_2_APs/all_indoor_1.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(1, :), 'LineWidth', 2);

hold on

h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(1, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/Only_2_APs/all_indoor_2.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(2, :), 'LineWidth', 2);
h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(2, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/Only_2_APs/all_indoor_4.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(3, :), 'LineWidth', 2);
h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(3, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/Only_2_APs/all_indoor_8.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(4, :), 'LineWidth', 2);
h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(4, :), 'LineWidth', 2);

% Show LF
h = cdfplot(error_lf_distance(:));
set(h, 'LineStyle', '-', 'color', rgb(5, :), 'LineWidth', 2);


legend("1 subarray MultiLoc", "1 subarray HF", "2 subarrays MultiLoc", "2 subarrays HF", "4 subarrays MultiLoc", "4 subarrays HF", "8 subarrays MultiLoc", "8 subarrays HF", "LF only", 'Location', "SouthEast")
title("")
xlabel('Location error [m]')
ylabel('')

xticks([0:2:22])
yticks([0:0.2:1])

xticks_zoom = [0:0.5:1];
set(gca,'YTicklabel',[])

proportion = 16/9;
x_length = 1.5;
y_length = x_length*proportion;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
xlim([0 22])

%set(leg,'visible','off')
fleg = legend('figure()');
set(fleg,'visible','off')

MagInset_Grid_On(fig, -1, [0 1 0 1], [12 21 0.15 0.6], {'NW','NW';'SE','SE'}, xticks_zoom, [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/Location/Location_system_indoor_2_aps.tikz');
savefig('plots/final_plots/Location/Location_system_indoor_2_aps');

cd(pwd_str)
%save_PDF('plots/final_plots/Location/Location_system_indoor_2_aps.pdf');

%% Only HF
% figure
% rgb = viridis(4);
% 
% load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/with_1_arrays_only_hf.mat')
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '-', 'color', rgb(1, :), 'LineWidth', 2);
% hold on
% 
% load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/with_2_arrays_only_hf.mat')
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '-', 'color', rgb(2, :), 'LineWidth', 2);
% 
% load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/with_4_arrays_only_hf.mat')
% 
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '-', 'color', rgb(3, :), 'LineWidth', 2);
% 
% load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/with_8_arrays_only_hf.mat')
% 
% h = cdfplot(error_hf_lf_tof(:));
% set(h, 'LineStyle', '-', 'color', rgb(4, :), 'LineWidth', 2);
% 
% legend("8 antenna arrays with 1 subarray each", "4 antenna arrays with 2 subarray each", "2 antenna arrays with 4 subarray each", "1 antenna arrays with 8 subarray each", 'Location', "SouthEast")
% title("Only HF")
% xlabel('Error in meters')
% ylabel('CDF')
% 
% matlab2tikz('Location-system-only-hf-multiple-arrays.tikz');
