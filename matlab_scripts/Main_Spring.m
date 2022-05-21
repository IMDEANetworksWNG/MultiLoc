%% Spring
% Take the number of paths
run("Spring/Decomposition.m")
% filter the number of paths
run("Spring/Check_Number_Paths_MP.m")
% apply the log GMM to improve FTM ranging
run("Spring/Guassian_Mixture.m")