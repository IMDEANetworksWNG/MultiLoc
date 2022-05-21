%% Indoor testbed
run("scripts_to_generate_mat_files_HF/Process_indoor.m")
run("scripts_to_generate_mat_files_HF/Process_indoor_AoA_MUSIC.m")

%% Outdoor testbed
run("scripts_to_generate_mat_files_HF/Process_outdoor.m")
run("scripts_to_generate_mat_files_HF/Process_outdoor_AoA_MUSIC.m")

%% Indoor/Outdoor testbed
run("scripts_to_generate_mat_files_HF/Process_in_out.m")
run("scripts_to_generate_mat_files_HF/Process_in_out_AoA_MUSIC.m")

%% IMDEA testbed
run("scripts_to_generate_mat_files_HF/Process_imdea.m")
run("scripts_to_generate_mat_files_HF/Process_imdea_AoA_MUSIC.m")