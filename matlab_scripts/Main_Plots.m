%% CDF azimuth elevation and tof plots
run("Errors_per_deviation/Errors_per_deviation_indoor_Final.m")

%% Generate CDF data
run("Coordinates/With_all_arrays_Final_more_subplot2.m")


%% CDF localization plots
run("Coordinates/cdf_with_all_arrays_final2.m")

%% BoxPlot
run("Coordinates/boxplot_with_all_arrays_final.m")