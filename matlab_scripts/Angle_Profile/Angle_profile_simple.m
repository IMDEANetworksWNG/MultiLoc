clear all
close all
clc

%% Angle profile
% LF data
load("mat_files/mD_track/LF/Data_3D.mat")

% HF data
load('mat_files/indoor/HF/CSI/csi_indoor.mat')
load("mat_files/indoor/HF/CSI/aoa_offset.mat")
load('mat_files/indoor/HF/FTM/ftm_indoor.mat')

% Point 26 is covered by 1:4
study_points = [26];

% AP labels
load('mat_files/indoor_map/data_distance_angle_true.mat');
study_aps_index = 1;

% for point=study_points
%     
%     figure
%     
%     for ap=study_aps_index
% 
%         LF
%         lf_power = power{point, ap}{1, 1};
%         lf_aoa   = AoA{point, ap}{1, 1}*-1;
%         
%         Remove the offsets
%         if point > 1 && point < 75 && ap == 2
%             
%             lf_aoa = lf_aoa - 3.5;
%         end
%         
%         if point > 1 && point < 75 && ap == 3
%             
%             lf_aoa = lf_aoa - 7;
%         end
%         
%         if point > 25 && point < 39 && ap == 5
%             
%             lf_aoa = lf_aoa - 5;
%         end
%         
%         Move power to dB
%         lf_power = 10*log10(lf_power);
%         
%         Normalize it
%         lf_power = lf_power - max(lf_power);
%         
%         Create the angle profile
%         angle_profile_lf = nan(180,1);
%         
%         angle_profile_lf(int64(lf_aoa + 90)) = lf_power;
%         
%         angle_profile_lf(isnan(angle_profile_lf)) = -25;
%         
%         subplot(9,1,1);
%         stem(-90:1:89, angle_profile_lf);
%         title('LF')
%         xlim([-90 90]);
%         ylim([-30 0]);
% 
%         figure;
% 
%         polarplot(deg2rad(-90:1:89), angle_profile_lf);
%         hold on
% 
%         set(h,'ThetaDir', 'counterclockwise');
%         
%         HF  
%         for rotation=1:8
%             
%             hf_power = power_raw{point, ap, rotation};
%             hf_aoa   = azimut_raw{point, ap, rotation};
%             
%             if size(hf_aoa, 1) > 0
%                 hf_aoa = hf_aoa - aoa_offset(ap);
%             end
%             
%             Move power to dB
%             hf_power = 10*log10(hf_power);
%             
%             Normalize it
%             hf_power = hf_power - max(hf_power);
%         
%             Create the angle profile
%             angle_profile_hf = nan(180,1);
%             
%             angle_profile_hf(int64(hf_aoa + 90)) = hf_power;
%             
%             angle_profile_hf(isnan(angle_profile_hf)) = -25;
% 
%             polarplot(deg2rad(-90:1:89), angle_profile_hf);
% 
%             subplot(9,1, rotation+1);
%             stem(-90:1:89, angle_profile_hf);
%             xlim([-90 90]);
%             ylim([-30 0]);            
%             xcorr
%             corr = xcorr(angle_profile_lf, angle_profile_hf, 0, 'coeff');
%             title(['HF rotation ' num2str(rotation*45-45)])
%             
%         end
%     end
%     
%     hold off
%     set(gca, 'YDir','reverse')
% end

for ap=study_aps_index
    
    figure
    
    for point=study_points
        
        % LF
        lf_power = power{point, ap}{1, 1};
        lf_aoa   = AoA{point, ap}{1, 1}*-1;

        % Remove the offsets
        if point > 1 && point < 75 && ap == 2
            
            lf_aoa = lf_aoa - 3.5;
        end
        
        if point > 1 && point < 75 && ap == 3
            
            lf_aoa = lf_aoa - 7;
        end
        
        if point > 25 && point < 39 && ap == 5
            
            lf_aoa = lf_aoa - 5;
        end
        
        % Move power to dB
        lf_power = 10*log10(lf_power);
        
        % Normalize it
        lf_power = lf_power - max(lf_power);
        lf_power = lf_power - min(lf_power);
        
        % Create the angle profile
        angle_profile_lf = nan(180,1);
        
        for id_path=1:size(lf_aoa, 2)
            
            path = lf_aoa(id_path);
           
            angle_profile_lf(int64(path + 90)) = max(angle_profile_lf(int64(path + 90)), lf_power(id_path));
        end
        %angle_profile_lf(int64(lf_aoa + 90)) = lf_power;
                
        subplot(5,1,1);
        stem(-90:1:89, angle_profile_lf);
        hold on
        % Draw the main path solid
        [a, i] = max(angle_profile_lf);
        stem(i-91, a, 'filled', 'color', 'blue');
        
        
        title(['LF point ' num2str(point) ' ap ' num2str(ap)])
        xlim([-90 90]);
        ylim([0 30]);
        %set(gca, 'YDir','reverse')

        
        %HF  
        for rotation=1:4
            
            hf_power = power_raw{point, ap, rotation};
            hf_aoa   = azimut_raw{point, ap, rotation};

            if size(hf_aoa, 1) > 0
                hf_aoa = hf_aoa - aoa_offset(ap);
            end
            
            % Move power to dB
            hf_power = 10*log10(hf_power);
            
            % Normalize it
            hf_power = hf_power - max(hf_power);
            hf_power = hf_power - min(hf_power);
        
            % Create the angle profile
            angle_profile_hf = nan(180,1);
            
            for id_path=1:size(hf_aoa, 2)

                path = hf_aoa(id_path);

                angle_profile_hf(int64(path + 90)) = max(angle_profile_hf(int64(path + 90)), hf_power(id_path));
            end
            %angle_profile_hf(int64(hf_aoa + 90)) = hf_power;
            
            subplot(5,1, rotation+1);
            stem(-90:1:89, angle_profile_hf, 'Color', 'red');
            hold on

            % Draw the main path solid
            [a, i] = max(angle_profile_hf);
            stem(i-91, a, 'filled', 'color', 'red');
            
            xlim([-90 90]);
            ylim([0 30]);

            
            hf_power = power_raw{point, ap, rotation +4};
            hf_aoa   = azimut_raw{point, ap, rotation+4};
            
            if size(hf_aoa, 1) > 0
                hf_aoa = hf_aoa - aoa_offset(ap);
            end
            
            % Move power to dB
            hf_power = 10*log10(hf_power);
            
            % Normalize it
            hf_power = hf_power - max(hf_power);
            hf_power = hf_power - min(hf_power);

            % Create the angle profile
            angle_profile_hf = nan(180,1);
            
            for id_path=1:size(hf_aoa, 2)

                path = hf_aoa(id_path);

                angle_profile_hf(int64(path + 90)) = max(angle_profile_hf(int64(path + 90)), hf_power(id_path));
            end
            
            %angle_profile_hf(int64(hf_aoa + 90)) = hf_power;
            
            %subplot(9,1, rotation+1);
            stem(-90:1:89, angle_profile_hf, 'Color', 'blue');
            % Draw the main path solid
            [a, i] = max(angle_profile_hf);
            stem(i-91, a, 'filled', 'color', 'blue');
            xlim([-90 90]);
            ylim([0 30]);
            %set(gca, 'YDir','reverse')

            title(['HF rotations ' num2str(rotation*45-45) ' [red] and ' num2str((rotation+4)*45-45) ' [blue]'])
        end
    end
    
    hold off
end

%% pcolor
point_id = study_points;
router_id = study_aps_index;

%angle_profile_lf = nan(180,1);

%toa    = ToA{point_id, router_id}{1, 1} - ToA{point_id, router_id}{1, 1}(1, 1);
toa    = ToA{point_id, router_id}{1, 1} - ToA{point_id, router_id}{1, 1}(1,1);
angles = AoA{point_id, router_id}{1, 1}*-1;

% Remove the offsets
if point_id > 1 && point_id < 75 && router_id == 2

    angles = angles - 3.5;
end

if point_id > 1 && point_id < 75 && router_id == 3

    angles = angles - 7;
end

if point_id > 25 && point_id < 39 && router_id == 5

    angles = angles - 5;
end

% -90 to 90
% -12.5 to 175
pcolor_matrix = nan(150,180);

% Time
for ii=1:size(toa, 2)
    
    
    toa_index = int64(toa(ii)*1.25 + 12.5);
    
    % Angle
%     for jj=angles
        
    angle_index = int64(angles(ii) + 90);

    pcolor_matrix(toa_index, angle_index) = 1;
%     end
end

figure
h = plot(int64(angles), int64(toa*1.25), '*', 'LineStyle', 'none');
%set(h, 'EdgeColor', 'none');
% 
% xticks(0:10:180)
% xticklabels(string(-90:10:90))
xlim([-90 90])
xlabel('Angle (deg)')
ylabel('Time (ns)')

% AoA quitar offset

% Tiempo 1.25 ns 115- 230
% pcolor

% %%
% % LF data
% load("mat_files/mD_track/LF/Data_3D.mat")
% 
% % HF data
% load('mat_files/indoor/HF/CSI/csi_indoor.mat')
% 
% study_points = [22:24];
% 
% % AP labels
% load('mat_files/indoor_map/data_distance_angle_true.mat');
% study_aps_index = 5;
% 
% for point=study_points
%     
%     figure
%     pax = polaraxes;
%     
%     for ap=study_aps_index
% 
%         % LF
%         lf_power = power{point, ap}{1, 1};
%         lf_aoa   = AoA{point, ap}{1, 1}*-1;
%         
%         % Move power to dB
%         lf_power = 10*log10(lf_power);
%         
%         % Normalize it
%         lf_power = lf_power - max(lf_power);
%         
%         % Create the angle profile
%         angle_profile_lf = nan(180,1);
%         
%         angle_profile_lf(int64(lf_aoa + 90)) = lf_power;
%         
%         angle_profile_lf(isnan(angle_profile_lf)) = -25;
%         
%         %subplot(9,1,1);
%         %plot(-90:1:89, angle_profile_lf);
%         %title('LF')
% 
%         %figure;
% 
%         polarplot(deg2rad(-90:1:89), angle_profile_lf);
%         hold on
% 
%         %set(h,'ThetaDir', 'counterclockwise');
%         
%         %HF  
%         for rotation=1
%             
%             hf_power = power_raw{point, ap, rotation};
%             hf_aoa   = azimut_raw{point, ap, rotation};
%             
%             % Move power to dB
%             hf_power = 10*log10(hf_power);
%             
%             % Normalize it
%             hf_power = hf_power - max(hf_power);
%         
%             % Create the angle profile
%             angle_profile_hf = nan(180,1);
%             
%             angle_profile_hf(int64(hf_aoa + 90)) = hf_power;
%             
%             angle_profile_hf(isnan(angle_profile_hf)) = -25;
% 
%             polarplot(deg2rad(-90:1:89), angle_profile_hf);
% 
%             %subplot(9,1, rotation+1);
%             %plot(-90:1:89, angle_profile_hf);
%             
%             
%             % xcorr
%             
%             %corr = xcorr(angle_profile_lf, angle_profile_hf, 0, 'coeff');
%             %title(['HF rotation ' num2str(rotation*45-45) ' corr: ' num2str(corr)])
%         end
%     end
%     
%     hold off
%     rlim([-25 0]);
%     pax.ThetaDir = 'clockwise';
%     pax.ThetaZeroLocation = 'top';
%     pax.Color = 'none';
% end

