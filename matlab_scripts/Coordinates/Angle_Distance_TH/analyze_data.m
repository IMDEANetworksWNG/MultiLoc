clear 
close all
clc

%%
pwd_str = pwd;

cd ../../../


scenarios_save = ["indoor", "outdoor", "2ap", "in_out", "imdea"];


for scenario = 1:5
    load(strcat("mat_files/data_error_mmWave_larger/", scenarios_save(scenario), "/error_multiloc_all"))
    error_rmse = zeros(size(error_multiloc_all));
    figure,
    for n_antennas = 1:4
        for id_angle = 1:length(angles_diff)
            error_rmse(n_antennas,id_angle) = sqrt(sum(error_multiloc_all{n_antennas,id_angle}));
        end
        
    end
    plot(angles_diff, error_rmse.')
end