close all
clc
clear

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
index_set_num = (1:110).';
index_set = string(index_set_num);

routers = [4;5;6;7;10];

for id_set = 1:length(index_set_num)

    for router = 1:length(routers)
        
        % point to the correct mat file to do the calibration
        if(routers(router) <=7 && routers(router) >=6)
            if (id_set > 36)
                load(strcat("../../mat_files/CSI/calibration_2/csi_data_raw_", string(routers(router))))
            else
                load(strcat("../../mat_files/CSI/calibration_1/csi_data_raw_", string(routers(router))))
            end
        else
            load(strcat("../../mat_files/CSI/calibration/csi_data_raw_", string(routers(router))))
        end
        
        if(id_set > 74)
            load(strcat("../../mat_files/CSI/calibration_new/csi_data_raw_", string(routers(router))))

            
        end
        
        if (routers(router) == 10 && id_set == 77)
            load(strcat("../../mat_files/CSI/calibration_point/csi_data_raw_", string(routers(router))))
        end
        
        % take the first packet as template
        csi_template = squeeze(csi_data(1,:,:,:));

        
        % load the csi data to be calibrated
        load(strcat("../../mat_files/CSI/", string(id_set), "/csi_data_raw_", string(routers(router))));

        % calibrate and save
        [snapshots, K, N, M] = size(csi_data);
    
        for tx_id = 1:M
            for rx_id = 1:N
                csi_data(:,:,rx_id, tx_id) = squeeze(csi_data(:,:,rx_id, tx_id))./ squeeze(csi_template(:,rx_id, tx_id).');
                csi_data(:,:,rx_id, tx_id) = squeeze(csi_data(:,:,rx_id, tx_id))./ mean(abs(squeeze(csi_data(:,grid_toa + 256/2 + 1,rx_id, tx_id))),2);
            end
        end
        
        csi_data_aux = nan(size(csi_data));
        csi_data_aux(:,grid_toa + 256/2 + 1,:,:) = csi_data(:,grid_toa + 256/2 + 1,:,:);
        csi_data = csi_data_aux;
        
        save(strcat("../../mat_files/CSI/", string(id_set), "/csi_data_calibrated_", string(routers(router))), "csi_data");

        
    end
    
end