clear ;
close all
clc

% mkdir("Plots/FTM")

index_set_num = (1:110).';
index_set = string(index_set_num);
% system("rm -rf ../../mat_files/FTM/*");
mkdir("../../mat_files/indoor")
mkdir("../../mat_files/indoor/LF")
mkdir("../../mat_files/indoor/LF/FTM")

routers = 1:5;
routers_CSI = [10,4,5,6,7];
    
for id_set = 1:length(index_set_num)
    mkdir(strcat("../../mat_files/FTM/", string(string(id_set))))
    for router = routers  

        fid = fopen(strcat("../../FTM/",string(id_set),"/FTM_",string(routers(router)),".txt"));
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
        distance_aux(distance_aux == 0) = [];
        
        FTM_distances = distance_aux;
        
        save(strcat("../../mat_files/indoor/LF/FTM/", string(id_set), "/FTM_distances_", string(routers_CSI(router))), "FTM_distances");

    end
end

