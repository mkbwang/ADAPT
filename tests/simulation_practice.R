library(metaSPARSim)






cond_A_sim_param <- list()
cond_A_sim_param$intensity <- c(100, 5000, 750, 70, 1000)
cond_A_sim_param$variability <- c(1,1,1,1,1)
cond_A_sim_param$lib_size <- c(1000, 1100, 900, 1200, 870)


fold_change_multiplier <- c(5, 0.04, 0.67, 1, 2)


cond_B_sim_param <- list()
cond_B_sim_param$intensity <- cond_A_sim_param$intensity * fold_change_multiplier # apply fold-change
cond_B_sim_param$variability <- cond_A_sim_param$variability
cond_B_sim_param$lib_size <- c(1100, 1300, 700, 1000, 950)

metaSPARSim_param <- list(cond_A = cond_A_sim_param, cond_B = cond_B_sim_param)
metaSPARSim_result <- metaSPARSim(metaSPARSim_param)







































# Load data
data(R1)

# Extract intesity, variability and library size parameters from the first experimental condition of preset R1 (i.e. Group_1)
param_preset <- R1$Group_1
intensity <- param_preset$intensity
variability <- param_preset$variability # why is it NA?
lib_size <- param_preset$lib_size


# Create the simulation parameter for condition A
cond_A_sim_param <- list()
cond_A_sim_param$intensity <- intensity
cond_A_sim_param$variability <- variability
cond_A_sim_param$lib_size <- lib_size

# not DA OTUs will have a fold change between 0.25 and 4
not_DA_multiplier <- runif(n = 3000, min = 0.25, max = 4)

# DA OTUs will have a fold change less than 0.25 or greater than 4
# here we simulate 241 fold-changes lower thant 0.25 and 300 fold changes greater than 4
DA_multiplier <- c( runif(n = 241, min = 0.0001, max = 0.25), runif(n = 300, min = 4, max = 1000) )

# In this example, the first 541 OTUs will be the DA ones, while the last 3000 will be the not DA ones
fold_change_multiplier <- c(DA_multiplier, not_DA_multiplier)

# Create the simulation parameter for condition B
cond_B_sim_param <- list()
cond_B_sim_param$intensity <- cond_A_sim_param$intensity * fold_change_multiplier # apply the fold-changes
cond_B_sim_param$variability <- cond_A_sim_param$variability
cond_B_sim_param$lib_size <- cond_A_sim_param$lib_size

# Create the global parameter
metaSPARSim_param <- list(cond_A = cond_A_sim_param, cond_B = cond_B_sim_param)
# Run metaSPARSim simulation
metaSPARSim_result <- metaSPARSim(metaSPARSim_param)





