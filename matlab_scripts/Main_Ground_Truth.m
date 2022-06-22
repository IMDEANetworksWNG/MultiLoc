cd ../Maps/

%% indoor
run("Map_indoor/Generate_Map_Data_True.m")

%% outdoor
run("Map_outdoor/Generate_Map_Data_True.m")

%% indoor
run("Map_in_out/Generate_Map_Data_True.m")

%% indoor
run("Map_imdea/Generate_Map_Data_True.m")