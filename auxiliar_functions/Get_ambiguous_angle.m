function candidate = Get_ambiguous_angle(in_angle)

    % create the grid for the angles
    angles = -90:0.1:90;

    % calculate the phase
    phases = 2*0.58*pi*sin(deg2rad(angles));


    % center everything between +/- pi
    phases(phases > pi) = phases(phases > pi) - 2*pi;
    phases(phases < -pi) = phases(phases < -pi) + 2*pi;

    % plot it again
    %figure, plot(angles,phases)
    
    % Transform the angle into phase
    aux_phase = 2*0.58*pi*sin(deg2rad(in_angle));
    
    if aux_phase > pi
        aux_phase = aux_phase -2*pi;
    else
        if aux_phase < -pi
            aux_phase = aux_phase +2*pi;
        end
    end
    
    offset = 0.0025;
    candidate_threshold = 20;
    
    % Find if there is more than one angle for that phase
    candidates = phases(phases>=(aux_phase-offset) & phases<=(aux_phase+offset));
    
    candidates(candidates==aux_phase) = [];
    
    if size(candidates, 2) == 0
       
        candidate = nan;
    else
        candidate = angles(ismember(phases, candidates));
        
        % Remove candidates that are too close
        
        candidate(candidate>=(in_angle-candidate_threshold) & candidate<=(in_angle+candidate_threshold)) = [];
        
        candidate = mean(candidate);
    end

end

%[INFO] GT: (60, -30) estimated: (-55, -16)

