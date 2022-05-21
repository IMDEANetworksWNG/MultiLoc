clear 
close all
clc

pwd_str = pwd;
cd ../../../

%% indoor
scenario = "indoor";
load("mat_files/" + scenario + "/HF/CSI/csi_" + scenario + ".mat")
load("mat_files/" + scenario + "/HF/FTM/ftm_" + scenario + ".mat")
load("mat_files/" + scenario + "_map/data_distance_angle_true.mat")
openfig("mat_files/" + scenario + "_map/map.fig")

[n_points, routers, rotations] = size(calculated_aoa);

% get the elevation
calculated_el = nan(size(calculated_aoa));
for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            el_aux = elevation_raw{id_point, id_router, id_rotation};
            if(~isempty(el_aux))
                calculated_el(id_point, id_router, id_rotation) = el_aux(1);
            end
        end
    end
end

index_minimum_tof = zeros(n_points,routers);
minimum_tof = zeros(n_points,routers);

optimal_el = zeros(n_points,routers);

for id_router = 1:routers
    for id_point = 1:n_points
        [minimum_tof(id_point, id_router), index_minimum_tof(id_point, id_router)] = nanmin(squeeze(calculated_distance(id_point, id_router,:)));
        optimal_el(id_point, id_router) = calculated_el(id_point, id_router, index_minimum_tof(id_point, id_router));
    end
end
figure, cdfplot(minimum_tof(:))

offset = [-1, -6, -3, +2 ,+1];
for ap = 1:routers
    optimal_el(:,ap) = optimal_el(:,ap) + offset(ap);
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);

end

index = {};
index{1} = [26:31,62:67,68,69,105,106];
index{2} = [37:61,90:104];
index{3} = [32:36,107,108];
index{4} = [70:74,109,110];
index{5} = [1:25,75:89];
% figure, 
% ap 1
error_points = [];
for ap = 1:routers
    errors_aux = optimal_el(index{ap}, ap);
    figure, cdfplot(errors_aux)
%     error_points = [error_points;errors_aux(:)];
end

% figure, cdfplot(error_points);
errors = optimal_el - zeros(size(optimal_el));



save("mat_files/" + scenario + "/HF/CSI/calculated_el_" + scenario + ".mat", "calculated_el")
%% outdoor
% close all
scenario = "outdoor";
load("mat_files/" + scenario + "/HF/CSI/csi_" + scenario + ".mat")
load("mat_files/" + scenario + "/HF/FTM/ftm_" + scenario + ".mat")
load("mat_files/" + scenario + "_map/data_distance_angle_true.mat")
% openfig("mat_files/" + scenario + "_map/map.fig")

[n_points, routers, rotations] = size(calculated_aoa);

% get the elevation
calculated_el = nan(size(calculated_aoa));
for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            el_aux = elevation_raw{id_point, id_router, id_rotation};
            if(~isempty(el_aux))
                calculated_el(id_point, id_router, id_rotation) = el_aux(1);
            end
        end
    end
end

index_minimum_tof = zeros(n_points,routers);
minimum_tof = zeros(n_points,routers);

optimal_el = zeros(n_points,routers);

for id_router = 1:routers
    for id_point = 1:n_points
        [minimum_tof(id_point, id_router), index_minimum_tof(id_point, id_router)] = nanmin(squeeze(calculated_distance(id_point, id_router,:)));
        optimal_el(id_point, id_router) = calculated_el(id_point, id_router, index_minimum_tof(id_point, id_router));
    end
end
figure, cdfplot(minimum_tof(:))

offset = [-2, +2, 0, -2];
for ap = 1:routers
    optimal_el(:,ap) = optimal_el(:,ap) + offset(ap);
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);
end

index = {};
index{1} = [1:16];
index{2} = [1:16];
index{3} = [1:16];
index{4} = [1:16];
% figure, 
% ap 1
error_points = [];
for ap = 1:routers
    errors_aux = optimal_el(index{ap}, ap);
    figure, cdfplot(errors_aux)
%     error_points = [error_points;errors_aux(:)];
end

% figure, cdfplot(error_points);

errors = optimal_el - zeros(size(optimal_el));

save("mat_files/" + scenario + "/HF/CSI/calculated_el_" + scenario + ".mat", "calculated_el")

%% in_out
close all
scenario = "in_out";
load("mat_files/" + scenario + "/HF/CSI/csi_" + scenario + ".mat")
load("mat_files/" + scenario + "/HF/FTM/ftm_" + scenario + ".mat")
load("mat_files/" + scenario + "_map/data_distance_angle_true.mat")
% openfig("mat_files/" + scenario + "_map/map.fig")

[n_points, routers, rotations] = size(calculated_aoa);

% get the elevation
calculated_el = nan(size(calculated_aoa));
for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            el_aux = elevation_raw{id_point, id_router, id_rotation};
            if(~isempty(el_aux))
                calculated_el(id_point, id_router, id_rotation) = el_aux(1);
            end
        end
    end
end

index_minimum_tof = zeros(n_points,routers);
optimal_el = zeros(n_points,routers);

for id_router = 1:routers
    for id_point = 1:n_points
        [~, index_minimum_tof(id_point, id_router)] = nanmin(squeeze(calculated_distance(id_point, id_router,:)));
        optimal_el(id_point, id_router) = calculated_el(id_point, id_router, index_minimum_tof(id_point, id_router));
    end
end

offset = [+6,-2];
for ap = 1:routers
    optimal_el(:,ap) = optimal_el(:,ap) + offset(ap);
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);
end

index = {};
index{1} = [1:10];
index{2} = [11:20];
% index{3} = [1:16];
% index{4} = [1:16];
% figure, 
% ap 1
for ap = 1:routers
    figure, cdfplot(optimal_el(index{ap}, ap))
end


errors = optimal_el - zeros(size(optimal_el));

save("mat_files/" + scenario + "/HF/CSI/calculated_el_" + scenario + ".mat", "calculated_el")


%% imdea
close all
scenario = "imdea";
load("mat_files/" + scenario + "/HF/CSI/csi_" + scenario + ".mat")
load("mat_files/" + scenario + "/HF/FTM/ftm_" + scenario + ".mat")
load("mat_files/" + scenario + "_map/data_distance_angle_true.mat")
% openfig("mat_files/" + scenario + "_map/map.fig")

[n_points, routers, rotations] = size(calculated_aoa);

% get the elevation
calculated_el = nan(size(calculated_aoa));
for id_rotation = 1:rotations
    for id_router = 1:routers
        for id_point = 1:n_points
            el_aux = elevation_raw{id_point, id_router, id_rotation};
            if(~isempty(el_aux))
                calculated_el(id_point, id_router, id_rotation) = el_aux(1);
            end
        end
    end
end

index_minimum_tof = zeros(n_points,routers);
optimal_el = zeros(n_points,routers);

for id_router = 1:routers
    for id_point = 1:n_points
        [~, index_minimum_tof(id_point, id_router)] = nanmin(squeeze(calculated_distance(id_point, id_router,:)));
        optimal_el(id_point, id_router) = calculated_el(id_point, id_router, index_minimum_tof(id_point, id_router));
    end
end

offset = [-2, -8, +2];
for ap = 1:routers
    optimal_el(:,ap) = optimal_el(:,ap) + offset(ap);
    calculated_el(:,ap,:) = calculated_el(:,ap,:) + offset(ap);
end
index = {};
index{1} = [2:2:12];
index{2} = [1:31];
index{3} = [1:20];

% figure, 
% ap 1
for ap = 1:routers
    figure, cdfplot(optimal_el(index{ap}, ap))
end

errors = optimal_el - zeros(size(optimal_el));

save("mat_files/" + scenario + "/HF/CSI/calculated_el_" + scenario + ".mat", "calculated_el")

