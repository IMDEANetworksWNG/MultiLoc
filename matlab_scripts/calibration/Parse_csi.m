% Parses the AoA results and returns a matrix
function [magnitudes, phases, times] = Parse_csi(filename)

    fid = fopen(filename);
    tline = fgetl(fid);
    
    correct   = 0;
    incorrect = 0;
    
    % This structure will hold the times of each measurement, this
    % way we can check if there are duplicates
    times = zeros(0, 1);
    
    % Magnitude
    magnitudes = zeros(1, 32);
    phases     = zeros(1, 32);
    
    % Read the file line by line
    while ischar(tline)
        
        % A valid line contains the string 'wil6210-aoa$' and has 71 ','
        if contains(tline, "[AOA]") == 1 && count(tline, ",") == 71
             
           % Split the data
           splitted = split(tline, ",");
           
           time = splitted{2};
           
           % Check if the time already exists
           if ~ismember(time, times)
               
               times(end+1, 1) = str2double(time);
           else
               
               disp("duplicated")
               continue
           end
           
          correct = correct + 1;

          % Now we can get the data
          % There are 32 elements
                   
          % Phase
          for i=9:(9+31)
              
              phases(correct, i-8) = str2num(splitted{i});
          end
          
          % Amplitude
          for i=(9+32):(9+32+31)
              
              magnitudes(correct, i-(9+31)) = str2num(splitted{i});
          end
          
        else
            
           incorrect = incorrect +1;
        end
        
        
        % Get next line
        tline = fgetl(fid);
    end
    
    %disp(['In the file ' filename ' there were ' num2str(correct) ' valid entries and ' num2str(incorrect) ' invalid entries, with a ratio of ' num2str(correct/(correct+incorrect))]);
    
    fclose(fid);
end

