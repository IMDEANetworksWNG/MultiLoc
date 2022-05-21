clear 
close all
clc

pwd_str = pwd;

cd ../../

load("mat_files/Spring/Data_eigenvalues_100.mat")

[n_points, routers, values] = size(D_all);

th = 5e-2;

paths = zeros(n_points,routers);

for id_router = 1:routers
   
    for id_point = 1:n_points
        D_router_point = squeeze(D_all(id_point,id_router,:));

        index_in = D_router_point >= D_router_point(1)*th;
        paths(id_point,id_router) = sum(index_in);
    end
end
figure, histogram(paths(:), "Normalization", "probability")

save("mat_files/Spring/paths", "paths")

cd(pwd_str)