clc
clear all
% close all

%%
pwd_str = pwd;
cd ../../
addpath(genpath("auxiliar_functions"))
rgb = viridis(5);

fig = figure;

load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/all_indoor_1.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(1, :), 'LineWidth', 2);

hold on

h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(1, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/all_indoor_2.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(2, :), 'LineWidth', 2);
h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(2, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/all_indoor_4.mat')
h = cdfplot(error_hf_hf_lf_new(:));
set(h, 'LineStyle', '-', 'color', rgb(3, :), 'LineWidth', 2);
h = cdfplot(error_hf_tof(:));
set(h, 'LineStyle', '--', 'color', rgb(3, :), 'LineWidth', 2);

load('matlab_scripts/Coordinates/With_multiple_antenna_arrays/all_indoor_8.mat')
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
ylabel('ECDF')

xticks([0:2:22])
yticks([0:0.2:1])

xticks_zoom = [0:0.5:1];

proportion = 16/9;
x_length = 1.5;
y_length = x_length*proportion;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
xlim([0 22])
%set(gca,'YTick',[])

%set(leg,'visible','off')
% fleg = legend('figure()');
% set(fleg,'visible','off')

MagInset_Grid_On(fig, -1, [0 1 0 1], [12 21 0.15 0.6], {'NW','NW';'SE','SE'}, [0:0.5:1], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz('plots/final_plots/Location/Location_system_indoor.tikz');
savefig('plots/final_plots/Location/Location_system_indoor');

cd(pwd_str)
