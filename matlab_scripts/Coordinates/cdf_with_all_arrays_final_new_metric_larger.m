clc
clear
close all

%%
pwd_str = pwd;
cd ../../
addpath(genpath("auxiliar_functions"))
rgb = viridis(5);
scenarios_save = ["indoor", "outdoor", "2ap", "in_out", "imdea"];
mkdir("plots")
mkdir("plots/final_plots/Location_new_larger")
mkdir("plots/final_plots/Location_new_larger/to_overleaf")

rotations = string([8,1]);

systems = ["MultiLoc", "MB baseline", "HF MultiLoc","HF baseline"];

line_styles = ["-", "--"];
for scenario = scenarios_save
counter = 1;

fig = figure;
legend_str = [];
for ii = 1:length(rotations)
    load(strcat("mat_files/CDF_data_new_larger/",scenario,"/all_",rotations(ii),".mat"))

    if ((scenario == "indoor" || scenario == "2ap") && ii == 2)
        samples = 2;
        error_multiloc = error_multiloc(:);
        error_multiloc = error_multiloc(1:samples:end);
        
        error_mb_bs_ap = error_mb_bs_ap(:);
        error_mb_bs_ap = error_mb_bs_ap(1:samples:end);
        
        error_multiloc_hf = error_multiloc_hf(:);
        error_multiloc_hf = error_multiloc_hf(1:samples:end);
        
        error_hf_bs_ap = error_hf_bs_ap(:);
        error_hf_bs_ap = error_hf_bs_ap(1:samples:end);
    end
h(counter,1) = cdfplot2(error_multiloc(:));
set(h(counter,1), 'LineStyle', line_styles(ii), 'color', rgb(1, :), 'LineWidth', 2);
counter = counter + 1;
hold on

h(counter,1) = cdfplot2(error_mb_bs_ap(:));
set(h(counter,1), 'LineStyle', line_styles(ii), 'color', rgb(2, :), 'LineWidth', 2);
counter = counter + 1;

h(counter,1) = cdfplot2(error_multiloc_hf(:));
set(h(counter,1), 'LineStyle', line_styles(ii), 'color', rgb(3, :), 'LineWidth', 2);
counter = counter + 1;

h(counter,1) = cdfplot2(error_hf_bs_ap(:));
set(h(counter,1), 'LineStyle', line_styles(ii), 'color', rgb(4, :), 'LineWidth', 2);
counter = counter + 1;

legend_str_aux = strcat(strcat(rotations(ii)," arr. "), systems);

legend_str = [legend_str,legend_str_aux];


if (ii == 1)
    % Show LF
    h(counter,1) = cdfplot2(error_multiloc_lf(:));
    set(h(counter,1), 'LineStyle', '-', 'color', rgb(5, :), 'LineWidth', 2);
    legend_str = [legend_str,"LF MultiLoc"];
    counter = counter + 1;
else
    % show spotfi
    % if (scenario == "indoor" || scenario == "outdoor" || scenario == "2ap")
    load(strcat("mat_files/CDF_data/",scenario,"/error_spotfi.mat"))
    h(counter,1) = cdfplot2(error(:));
    set(h(counter,1), 'LineStyle', '--', 'color', rgb(5, :), 'LineWidth', 2);
    legend_str = [legend_str,"Spotfi"];
    counter = counter + 1;
end
% else
%     legend_str_aux = strcat(strcat(rotations(ii)," arr. "), systems);
% end
end





leg = legend(legend_str,'Location', "SouthEast");
% end
% legend([h(1:4);h(9);h(5:8);h(10)],[legend_str(1:4),legend_str(9),legend_str(5:8),legend_str(10)], 'Location', "SouthEast")
title("")
xlabel('Location error [m]')
ylabel('ECDF')

xticks([0:3:33])
yticks([0:0.2:1])

xticks_zoom = [0:0.5:1];
if (sum(scenario == ["outdoor", "in_out"]) == 1)
    proportion = 14/9;

else
    proportion = 16/9;

end
x_length = 1.5;
y_length = x_length*proportion;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
xlim([0 33])
%set(gca,'YTick',[])

%set(leg,'visible','off')
% fleg = legend('figure()');
% set(fleg,'visible','off')

MagInset_Grid_On(fig, -1, [0 3 0 1], [15 30 0.15 0.6], {'NW','NW';'SE','SE'}, [0:1:3], [0:0.5:1]);

fleg = legend('figure()');
set(fleg,'visible','off')

matlab2tikz(char(strcat("plots/final_plots/Location_new_larger/Location_system_",scenario,".tikz")));
savefig(fig,strcat("plots/final_plots/Location_new_larger/Location_system_",scenario));
end

[status,result] = system('bash plots/final_plots/Location_new_larger/modify_plots_all.sh');

cd(pwd_str)

            