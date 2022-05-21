% Parses the FTM results and returns a matrix
function [ftm_times] = Parse_ftm(filename)

    fid = fopen(filename);
    tline = fgetl(fid);

    % Save the data
    ftm_times = zeros(0, 4);
    
    % Counter
    curr_measurement = 1;
    
    % Read the file line by line
    while ischar(tline)
        
        % A valid line contains the string '[FTM] Measurement:' and has 3 ','
        if contains(tline, "[FTM] Measurement:") == 1 && count(tline, ",") == 3
             
            % Split the data
            splitted = split(tline, ":");
            
            splitted = split(splitted{6}, " ");

            % Times
            aux = split(splitted{2}, ",");
            t1 = str2num(aux{1});
            
            aux = split(splitted{3}, ",");
            t2 = str2num(aux{1});
            
            aux = split(splitted{4}, ",");
            t3 = str2num(aux{1});
            
            t4 = str2num(splitted{5});
            
            ftm_times(curr_measurement, 1) = t1;
            ftm_times(curr_measurement, 2) = t2;
            ftm_times(curr_measurement, 3) = t3;
            ftm_times(curr_measurement, 4) = t4;

            curr_measurement = curr_measurement + 1;
        end
        
        % Get next line
        tline = fgetl(fid);
    end    
    
    fclose(fid);
end

