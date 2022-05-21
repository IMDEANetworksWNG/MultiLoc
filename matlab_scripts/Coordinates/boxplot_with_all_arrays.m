clc
clear
% close all

%%
pwd_str = pwd;
cd ../../
addpath(genpath("auxiliar_functions"))
addpath(genpath("functions"))

rgb = viridis(5);
scenarios_save = ["indoor", "outdoor", "2ap", "in_out"];
mkdir("plots")
mkdir("plots/final_plots_v2")
mkdir("plots/final_plots_v2/Location")
rotations = string([8,4,2,1]);

systems = ["MT", "HF-MT","HF-BS", "MB-BS"];
labels = [];
line_styles = ["-", "--"];
error_all = nan(110*8,16);
for scenario = scenarios_save(1)
% fig = figure;
legend_str = [];
for ii = 1:length(rotations)
load(strcat("mat_files/CDF_data/",scenario,"/all_",rotations(ii),".mat"))
len_error = length(error_multiloc(:));
error_all(1:len_error,(ii-1)*4 + 1) = error_multiloc(:);
error_all(1:len_error,(ii-1)*4 + 2) = error_mb_bs_ap(:);
error_all(1:len_error,(ii-1)*4 + 3) = error_multiloc_hf(:);
error_all(1:len_error,(ii-1)*4 + 4) = error_hf_bs_ap(:);

labels = [labels;strcat(systems,"/", rotations(ii))];

end

% % Show LF
% h = cdfplot2(error_multiloc_lf(:));
% set(h, 'LineStyle', '-', 'color', rgb(5, :), 'LineWidth', 2);
% legend_str = [legend_str,"LF only"];
fig = figure; %boxplot(error_all)
ylabel("Localization error [m]")
positions = 1;
space = 0.6;
space_group = 0.2;
for ii = 1:4
    positions = [positions,repmat(positions(end),1,3) + space_group*(1:3)];
    positions(end+1) = positions(end) + space;
end
positions(end) = [];
% positions = [1,1.2,1.4,1.6,2.6,2.8,3,3.2,4.2,4.4,4.6,4.8,5.8,6,6.2,6.4];
groups =    [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4];
% boxplot(error_all,1:16, 'positions', positions);
color_boxplot = repmat(rgb(1:4,:),4,1);
boxplot(error_all,1:16, 'positions', positions, "Colors", color_boxplot,'Widths',0.15);
h = findobj(gca,'Tag','Box');
h = flip(h);
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color_boxplot(j,:),'FaceAlpha',.5);
%    h(ii).MarkerEdgeColor = color_boxplot(j,:);
end
h = findobj(gcf,'tag','Outliers')
h = flip(h);
for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color_boxplot(j,:),'FaceAlpha',.5);
   h(j).MarkerEdgeColor = color_boxplot(j,:);
end
all_items = flip(findall(gca,'Tag','Box'));

% for k = 1:length(h)
%     x = get(h(k),'XData');
%     y = get(h(k),'YData');
%     c = get(h(k),'Color');
%     l = get(h(k),'LineWidth');
%     ht = y(2)-y(1);
%     wd = x(3)-x(1);
%     rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',c,'LineWidth',10)
% end

[hlegend, hobj, ~, ~] = legend([all_items(1);all_items(2);all_items(3);all_items(4)], {'MultiLoc', "Multiband Baseline",'HF-Multiloc','HF-Baseline'}, "Location", "Northwest", 'FontSize', 12);
% hlegend.FontSize = 13;
% set(hlegend,'FontSizeMode','manual')
% set(hlegend,'FontSize',5);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2);

xlabels = strcat(string([8,4,2,1]), " array");
xlabels(1:end-1) = xlabels(1:end-1) + "s";
set(gca,'xtick',positions(2:4:16) + 0.1)
set(gca,'xticklabel',xlabels)
ylim([-1 15])
yticks(0:1:14)
ylabel("Location error [m]")
position_aux = fig.Position;
position_aux(4) = 300;
set(fig, "Position", position_aux)

position_legend =     [0.1298    0.6561    0.3339    0.2717];
hlegend.Position = position_legend
set(hlegend,'color','none');set(hlegend,'color','none');
set(gca, "fontsize", 13)
set(gca, 'YGrid', 'on', 'XGrid', 'off')

save_PDF_fig(fig, "plots/final_plots_v2/Location/to_overleaf/Boxplot")
% hLegend.LineWidth = 2;
% xticklabels(reshape(labels.',1,[]))
% legend(legend_str, 'Location', "SouthEast")
% title("")
% xlabel('Location error [m]')
% ylabel('ECDF')
% 
% xticks([0:2:20])
% yticks([0:0.2:1])
% 
% xticks_zoom = [0:0.5:1];
% 
% proportion = 16/9;
% x_length = 1.5;
% y_length = x_length*proportion;
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, y_length, x_length])
% xlim([0 22])
%set(gca,'YTick',[])

%set(leg,'visible','off')
% fleg = legend('figure()');
% set(fleg,'visible','off')
% 
% MagInset_Grid_On(fig, -1, [0 1 0 1], [12 21 0.15 0.6], {'NW','NW';'SE','SE'}, [0:0.5:1], [0:0.5:1]);
% 
% fleg = legend('figure()');
% set(fleg,'visible','off')
% 
% matlab2tikz(char(strcat("plots/final_plots_v2/Location/Location_system_",scenario,".tikz")));
% savefig(strcat("plots/final_plots_v2/Location/Location_system_",scenario));
end

% [status,result] = system('bash plots/final_plots_v2/Location/modify_plots_all.sh');

cd(pwd_str)
