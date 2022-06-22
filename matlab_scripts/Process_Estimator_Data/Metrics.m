clear
clc
close all

pwd_str = pwd;

cd ../../
mkdir("mat_files/")
mkdir("mat_files/metrics")

scenarios = ["indoor", "outdoor","in_out","imdea"];

for scenario = 1:length(scenarios)
    
    load(strcat("mat_files/mD_track/LF/Data_3D_processed_",scenarios(scenario),".mat"))
    
    [points, routers] = size(aoa_routers);
    
    power_ratio = zeros(size(aoa_routers));
    mean_excess_delay = zeros(size(aoa_routers));
    
    for point = 1:points
        for router = 1:routers
        
            pw_aux = power_routers{point, router};
            power_ratio(point,router) = pw_aux(1)/sum(pw_aux);
            time_aux = toa_routers{point, router};
            % normalize time of arriva
            time_aux = time_aux - time_aux(1);
            mean_excess_delay(point,router) = sum(time_aux.*pw_aux)/sum(pw_aux);
            % inverse the power ratio so that we take the mimimum
            power_ratio = 1./power_ratio;
            
        end
    end
    save(strcat("mat_files/metrics/metrics_",scenarios(scenario),".mat"), "power_ratio", "mean_excess_delay")

end

cd(pwd_str)