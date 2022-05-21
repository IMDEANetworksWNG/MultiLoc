clc
clear
close all

%%
pwd_str = pwd;
cd ../../
addpath(genpath("auxiliar_functions"))
rgb = viridis(5);
scenarios_save = ["indoor", "outdoor", "2ap", "in_out"];
mkdir("plots")
mkdir("plots/final_plots_v2")
mkdir("plots/final_plots_v2/Location")

for scenario = scenarios_save
fig = figure;

load(strcat("mat_files/CDF_data/",scenario,"/all_1.mat"))
h = cdfplot2(error_multiloc(:));
set(h, 'LineStyle', '-', 'color', rgb(1, :), 'LineWidth', 2);

hold on

h = cdfplot2(error_multiloc_hf(:));
set(h, 'LineStyle', '--', 'color', rgb(1, :), 'LineWidth', 2);

load(strcat("mat_files/CDF_data/",scenario,"/all_2.mat"))
h = cdfplot2(error_multiloc(:));
set(h, 'LineStyle', '-', 'color', rgb(2, :), 'LineWidth', 2);
h = cdfplot2(error_multiloc_hf(:));
set(h, 'LineStyle', '--', 'color', rgb(2, :), 'LineWidth', 2);

load(strcat("mat_files/CDF_data/",scenario,"/all_4.mat"))
h = cdfplot2(error_multiloc(:));
set(h, 'LineStyle', '-', 'color', rgb(3, :), 'LineWidth', 2);
h = cdfplot2(error_multiloc_hf(:));
set(h, 'LineStyle', '--', 'color', rgb(3, :), 'LineWidth', 2);

load(strcat("mat_files/CDF_data/",scenario,"/all_8.mat"))
h = cdfplot2(error_multiloc(:));
set(h, 'LineStyle', '-', 'color', rgb(4, :), 'LineWidth', 2);
h = cdfplot2(error_multiloc_hf(:));
set(h, 'LineStyle', '--', 'color', rgb(4, :), 'LineWidth', 2);

% Show LF
h = cdfplot2(error_multiloc_lf(:));
set(h, 'LineStyle', '-', 'color', rgb(5, :), 'LineWidth', 2);


legend("1 subarray MultiLoc", "1 subarray HF", "2 subarrays MultiLoc", "2 subarrays HF", "4 subarrays MultiLoc", "4 subarrays HF", "8 subarrays MultiLoc", "8 subarrays HF", "LF only", 'Location', "SouthEast")
title("")
xlabel('Location error [m]')
ylabel('ECDF')

xticks([0:2:20])
yticks([0:0.2:1])

xticks_zoom = [0:0.5:1];

proportion = 13/9;
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

matlab2tikz(char(strcat("plots/final_plots_v2/Location/Location_system_",scenario,".tikz")));
savefig(strcat("plots/final_plots_v2/Location/Location_system_",scenario));
end
cd(pwd_str)
