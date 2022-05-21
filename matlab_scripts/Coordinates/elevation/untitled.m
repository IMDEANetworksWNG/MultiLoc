clear 
% close all
clc

%%
pwd_str = pwd;

cd ../../../

mkdir("plots")
mkdir("plots/final_plots_subplot_elevation")

addpath('functions/');
% addpath('matlab_scripts/fu');

scenarios = ["indoor", "outdoor", "in_out", "imdea"];
scenarios_save = ["indoor", "outdoor", "in_out", "imdea"];

subarrays = [8,4,2,1];
index_subarrays = {0:7,0:2:6,0:4:4,0};


colors = 	[0, 0.4470, 0.7410,
          	0.8500, 0.3250, 0.0980,
          	0.9290, 0.6940, 0.1250,
          	0.4940, 0.1840, 0.5560,
          	0.4660, 0.6740, 0.1880];

mkdir("mat_files/")
mkdir("mat_files/CDF_data_elevation")
        
th_disagreement_all = [4.5,4.5,2,2];

fig_distance = figure;
fig_el = figure;
for scenario = 1:length(scenarios)
    
    
    % load the HF FTM data
    load(strcat("mat_files/",scenarios(scenario),"/HF/FTM/ftm_",scenarios(scenario),".mat"))
    load(strcat("mat_files/",scenarios(scenario),"/HF/CSI/csi_",scenarios(scenario),".mat"))

    
    load(strcat("mat_files/CDF_data_elevation/",scenarios_save(scenario),"/chosen_rotation_with_hf"), "chosen_rotation_with_hf", "index_ap_in");
    
    [points, aps,rotations] = size(calculated_aoa);
    
    % get the elevation
    calculated_el = nan(size(calculated_aoa));
    for id_rotation = 1:rotations
        for id_router = 1:aps
            for id_point = 1:points
                el_aux = elevation_raw{id_point, id_router, id_rotation};
                if(~isempty(el_aux))
                    calculated_el(id_point, id_router, id_rotation) = el_aux(1);
                end
            end
        end
    end

    chosen_distance = zeros(size(chosen_rotation_with_hf));
    chosen_el = zeros(size(chosen_rotation_with_hf));

    for point = 1:points
        for ap = 1:aps
            if (~isnan(chosen_rotation_with_hf(point,ap)))
                chosen_distance(point,ap) = calculated_distance(point, ap, chosen_rotation_with_hf(point,ap));
                chosen_el(point,ap) = calculated_el(point, ap, chosen_rotation_with_hf(point,ap));

            end
            
        end
    end
    chosen_distance_multiloc = chosen_distance(logical(index_ap_in));
    chosen_el_multiloc = chosen_el(logical(index_ap_in));

    set(0, 'CurrentFigure', fig_distance)
    cdfplot(chosen_distance_multiloc)
    hold on
    
    set(0, 'CurrentFigure', fig_el)
    cdfplot(abs(chosen_el_multiloc))
    hold on
    
end
set(0, 'CurrentFigure', fig_distance)
hold off
legend(scenarios)

set(0, 'CurrentFigure', fig_el)
hold off
legend(scenarios)