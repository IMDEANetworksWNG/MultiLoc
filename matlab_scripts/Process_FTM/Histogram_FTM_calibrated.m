clear ;
close all
clc

mkdir("../../Plots/FTM")
mkdir("../../Plots/FTM/Histogram_calibrated")
% mkdir("../../Plots/FTM/Histogram_calibrated/png")
% mkdir("../../Plots/FTM/Histogram_calibrated/fig")


index_set_num = (1:110).';
index_set = string(index_set_num);
% system("rm -rf mat_files/FTM/*");

routers = 1:5;
routers_CSI = [4,5,6,7,10];

load("../../mat_files/indoor_map/data_distance_angle_true.mat")

true_distances = distances_RX;

% figure, plot(distance_10)
% figure, plot(distance_4)
% figure, plot(distance_5)
% figure, plot(distance_6)
for router = 1:length(routers)  
    mkdir(strcat("../../Plots/FTM/Histogram_calibrated/", string(routers_CSI(router))));
    mkdir(strcat("../../Plots/FTM/Histogram_calibrated/", string(routers_CSI(router)), "/fig"));
    mkdir(strcat("../../Plots/FTM/Histogram_calibrated/", string(routers_CSI(router)), "/png"));
    for id_set = 1:length(index_set_num)

        load(strcat("../../mat_files/indoor/LF/FTM/", string(id_set), "/FTM_distances_calibrated_", string(routers_CSI(router))));
        hFig = figure("visible", "off");
        subplot(1,2,1)
        histogram(FTM_distances)
        subplot(1,2,2)
        cdfplot(FTM_distances)
        hold on
        p_true = plot(true_distances(id_set,router)*ones(length(0:01:1),1),0:01:1);
        p_mean = plot(mean(FTM_distances)*ones(length(0:01:1),1),0:01:1);
        p_median = plot(median(FTM_distances)*ones(length(0:01:1),1),0:01:1);
        legend([p_true, p_mean, p_median], ["True distance", "Mean distance", "Median distance"], "Location", "southeast");

%         legend("True distance", "Mean distance", "Median distance");
        print(strcat("../../Plots/FTM/Histogram_calibrated/",string(routers_CSI(router)),"/png/",string(id_set)), "-dpng")
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        savefig(strcat("../../Plots/FTM/Histogram_calibrated/",string(routers_CSI(router)),"/fig/", string(id_set)))    
        close all
    end
end
