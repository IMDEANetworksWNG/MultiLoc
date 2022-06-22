clear ;
close all
clc

% mkdir("Plots/FTM")

index_set_num = (1:31).';
index_set = string(index_set_num);
% system("rm -rf ../../mat_files/FTM/*");
mkdir("../../mat_files/imdea/");
mkdir("../../mat_files/imdea/LF/");
mkdir("../../mat_files/imdea/LF/FTM");
routers = [2 4 1];
routers_CSI = [222,223,224];
    
for id_set = 1:length(index_set_num)
    mkdir(strcat("../../mat_files/imdea/LF/FTM/", string(string(id_set))))
    for router = 1:length(routers)  

        fid = fopen(strcat("../../FTM_imdea/",string(id_set),"/FTM_",string(routers(router)),".txt"));
        tline = fgetl(fid);
        counter = 0;
        text = "";
        while ischar(tline)
            counter = counter + 1;
            if (mod(counter,2) == 0)

                text(counter/2,1) = string(tline);
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        text = split(text);

        distance_aux = str2double(text(:,9))/100;
        distance_aux(distance_aux == 0) = nan;
        
        FTM_distances = distance_aux;
        
        save(strcat("../../mat_files/imdea/LF/FTM/", string(id_set), "/FTM_distances_", string(routers_CSI(router))), "FTM_distances");

    end
end

