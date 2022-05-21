clear 
% close all
clc

%%
pwd_str = pwd;

cd ../../../

mkdir("plots")
mkdir("plots/final_plots_subplot")

addpath('functions/');
% addpath('matlab_scripts/fu');

scenarios = ["indoor", "outdoor", "indoor", "in_out", "imdea"];
scenarios_save = ["indoor", "outdoor", "2ap", "in_out", "imdea"];

subarrays = [8,4,2,1];
index_subarrays = {0:7,0:2:6,0:4:4,0};


colors = 	[0, 0.4470, 0.7410,
          	0.8500, 0.3250, 0.0980,
          	0.9290, 0.6940, 0.1250,
          	0.4940, 0.1840, 0.5560,
          	0.4660, 0.6740, 0.1880];

mkdir("mat_files/")
mkdir("mat_files/th")

ths = 0.1:0.1:5;
for scenario = 1:length(scenarios)
    load(strcat("mat_files/th/", scenarios_save(scenario), "/error_multiloc_all"))
    rmse = zeros(size(error_multiloc_all));

    for th = 1:length(ths)
        for rot = 1:length(subarrays)

            error_multiloc = error_multiloc_all{th,rot};
            rmse(th,rot) = sqrt(mean(error_multiloc.^2));
        end
    end
    figure, plot(ths,rmse)
    legend(string(subarrays))
    title(scenarios_save(scenario))
end