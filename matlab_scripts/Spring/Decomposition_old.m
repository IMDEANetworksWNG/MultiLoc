close all
clc
clear

pwd_path = pwd;

cd ../../

mkdir("mat_files/Spring")

% add some parameters to load the file
addpath("functions")
% Bandwidth
BW = 80;
% Number of spatial stream
M = 1;

% index to the active subcarriers

  grid_toa = [-122, -121, -120, -119, -118, -117, -116, -115, -114, -113, -112, -111, -110, ...
-109, -108, -107, -106, -105, -104, -102, -101, -100, -99, -98, -97, -96, -95, ...
-94, -93, -92, -91, -90, -89, -88, -87, -86, -85, -84, -83, -82, -81, -80, -79, ...
-78, -77, -76, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65, -64, -63, -62, ...
-61, -60, -59, -58, -57, -56, -55, -54, -53, -52, -51, -50, -49, -48, -47, -46, ...
-45, -44, -43, -42, -41, -40, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, ...
-28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, ...
-12, -10, -9, -8, -7, -6, -5, -4, -3, -2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, ...
17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, ...
40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, ...
62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 83, 84, ...
85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105, ...
106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, ...
122];


% index_set = 300;
index_set_num = (1:74).';
index_set = string(index_set_num);

print_spectrum = 0;

step_toa = 0.05;
step_aoa = 1/(180*1);
sm_f_aoa = 3;

sm_f_toa = 100;

routers = [4,5,6,7,10];

D_all = zeros(length(index_set_num),length(routers),sm_f_toa);
% for n_e = [4,5,6]

%     ToA = {};
%     AoA = {};
%     power = {};

for id_set = 1:length(index_set_num)
    for router = 1:length(routers)

    [sm_f_toa, id_set, router]
    % take the file name
    file_name = strcat("mat_files/CSI/",index_set(id_set),"/csi_data_calibrated_",string(routers(router)),".mat");

    load(file_name)
    % get the # of snapshots, subcarriers, number of antennas to received and
    % everything
    [snapshots, K, N, M] = size(csi_data);

    ToA_aux = {};
    power_aux = {};

%         n_e = 6;
    
    csi_data = csi_data(:,grid_toa + (K/2) + 1,:,:);
    
    for snaptshot = 1

%         [S_toa, times] = ToA_Phases(step_toa, 256, BW, grid_toa);

        channel_toa = squeeze(csi_data(snaptshot,:,1,1)).';
        [R_toa] = Smoothing_1D_faster(channel_toa,sm_f_toa, 1);

        [U,D] = eig(R_toa);
        D = diag(D);
        [D,ind] = sort(D, 'descend');
    end

    D_all(id_set, router,:) = D;
%     power{id_set,router} = power_aux;
    end
end


save(strcat("mat_files/Spring/Data_eigenvalues_", string(sm_f_toa)), "D_all")
% end
cd(pwd_path)
