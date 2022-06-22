%% csireader.m
%
% read and plot CSI from UDPs created using the nexmon CSI extractor (nexmon.org/csi)
% modify the configuration section to your needs
% make sure you run >mex unpack_float.c before reading values from bcm4358 or bcm4366c0 for the first time
%

clear
% close all
clc

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

%% configuration
CHIP = '4366c0';          % wifi chip (possible values 4339, 4358, 43455c0, 4366c0)
BW = 80;                % bandwidth
ss = 4;
N = 4;

index_set = ["calibration";string((1:31).')];
mkdir("../../mat_files/")
mkdir("../../mat_files/imdea/");
mkdir("../../mat_files/imdea/LF/");
mkdir("../../mat_files/imdea/LF/CSI");

routers = [222, 223, 224];

for id_set = 1:length(index_set)
    mkdir(strcat("../../mat_files/imdea/LF/CSI/",index_set(id_set)))
          
%     if (id_set <= 2)
%         routers = [6,7];
% %     end
%     elseif (id_set == 3)
%         routers = [4,5,10];
%     else
%         routers = [4 5 6 7 10];
%     end
    for id_router = routers

        FILE = strcat("../../pcap_files_imdea/",index_set(id_set),"/./trace", string(id_router),".pcap");

        [cmplxall_raw_all] = CSI_extractor_Function(FILE);
        
        [packets,~,~,~] = size(cmplxall_raw_all);
        
        csi_data = zeros(packets,256,N,ss);

        for jj = 1:packets
            for ii = 1:N
%                (((ii-1)*N)+1):((ii*N)-ss);
                csi_data(jj,:,ii,:) = cmplxall_raw_all(jj,:,(((ii-1)*N)+1):((ii*N)-(N-ss)));
        %         figure, plot(abs(squeeze(csi_data(:,:,ii))))
            end
        end

%         csi_data = squeeze(csi_data(:,grid_toa + 256/2 + 1,:,:));

        save(strcat("../../mat_files/imdea/LF/CSI/",index_set(id_set),"/csi_data_raw_", string(id_router)), "csi_data")
    end
end
% end
