%% Indoor testbed
run("scripts_to_generate_mat_files_HF/Process_indoor.m")

%% Outdoor testbed
run("scripts_to_generate_mat_files_HF/Process_outdoor.m")

%% Indoor/Outdoor testbed
run("scripts_to_generate_mat_files_HF/Process_in_out.m")

%% IMDEA testbed
run("scripts_to_generate_mat_files_HF/Process_imdea.m")
