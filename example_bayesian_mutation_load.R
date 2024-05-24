################################
# Code for conducting Bayesian #
# regression with a metric    #
# predicted variable and one   #
# nominal predictor variable.  #
################################

###############################
# LOAD APPROPRIATE LIBRARIES  #
###############################
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggridges)


################################
# WHEN THERE ARE TWO LEVELS IN #
# THE NOMINAL CATEGORY.        #
################################


#---------------------#
#   Load the data     #
#---------------------#

roh_data <- read.table("/Volumes/cetacea/Genetic_load/_roh/ROH_stats", header = T)
head(roh_data)

roh_auto <- roh_data[which(roh_data$CHR == "Autosomes"),] 
roh_x <- roh_data[which(roh_data$CHR == "X_Chromosome"),] 

roh_data <- rbind(roh_auto, roh_x)

#----------------------------#
# Standardize the y-variable #
# (metric)                   #
#----------------------------#
roh = geneticLoad_indv_NARW_data_NS$TotalLoad
rohMean = mean(roh)
rohSD = sd(roh)
zroh = (roh - rohMean) / rohSD
N = length(zroh)


#---------------------------#
# Organize the categorical  #
# data   CHR                #
#---------------------------#
CHR = as.numeric(geneticLoad_indv_NARW_data_NS$CHR)
CHRNames = levels(geneticLoad_indv_NARW_data_NS$CHR)
nCHR = length(unique(CHRNames))


#-----------------------------#
# Create a data list for Stan #
#-----------------------------#
dataList = list(
  N = N,
  nxLevels = nCHR,
  y = zroh,
  x= CHR
)


#------------------------------#
#       Define the Model       #
#------------------------------#
stan_model = write_stan_file("
  data {
    int<lower=0> N;         // Sample size
    int<lower=0> nxLevels;  // Number of levels in the nominal predictor
    vector[N] y;            // Standardized predicted variable
    int x[N];               // Nominal predictor variable
  }

  parameters {
    real b1[nxLevels];     // Effect of being in each xLevel on y
    real<lower=0> sigma;   // Remaining variation in y 
  }

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
      mu[i] = b1[x[i]];
      y[i] ~ normal(mu[i], sigma);
    }

    // Priors
    for (j in 1:nxLevels) {
      b1[j] ~ normal(0, 1);
    }

    sigma ~ normal(0, 1);
  }
")


#----------------------------#
#  Prior Predictive Check    #
#----------------------------#

# For each individual, predict a height based on their sex and prior values

# Vector to hold predicted heights
pressurePriorPred = rep(NA, times = N)

for (i in 1:N) {
  
  # Mean value (same for both faculty because priors are the same)
  meanPred = rnorm(n = 1, mean = 0, sd = 1)
  
  # Sigma value
  sigmaPred = rnorm(n = 1, mean = 0, sd = 1)
  
  # Ensure sigma is positive
  while(sigmaPred < 0) {
    sigmaPred = rnorm(n = 1, mean = 1, sd = 1)
  }
  
  # Predicted value
  pressurePriorPred[i] = rnorm(n = 1, mean = meanPred, sd = sigmaPred)
}

#--- Organize the Data ---#
values = c(zroh, pressurePriorPred)
labels = c(rep("data", times = N), rep("pred", times = N))
data = data.frame(values, labels)

#--- Plot the data ---#
ggplot(data) +
  theme_bw() +
  geom_density(aes(x = values, group = labels, fill = labels), alpha = 0.6) +
  xlab("CHR or Population or Impact")



#############################
# COMPILE AND RUN THE MODEL #
#############################

# translate the model to C++ and compile
model_stan = cmdstan_model(stan_file = stan_model)

# Run the model
fit_model = model_stan$sample(data = dataList,
                              chains = 3,
                              parallel_chains = 3,
                              refresh = 1000,
                              iter_warmup = 2000,
                              iter_sampling = 10000
)


##############################
# Evaluate Model Performance #
##############################
fit_model$cmdstan_diagnose()              # Run model diagnostics
print(fit_model)                          # Print diagnostics
mcmc_trace(fit_model$draws(), pars = c("b1[1]", "b1[2]", "sigma"), n_warmup = 2000)  # Traceplots



#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
posteriors_stan = as_draws_array(fit_model$draws())

mcmc_intervals(posteriors_stan, pars = c("b1[1]", "b1[2]", "sigma"))
### OR ####
mcmc_intervals(posteriors_stan, pars = c("b1[1]", "b1[2]", "b1[3]", "sigma"))
### OR ####
mcmc_intervals(posteriors_stan, pars = c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "sigma"))




#-------------------------#
# Extract the data        #
#-------------------------#
posteriors = data.frame(as_draws_df(fit_model$draws()))
head(posteriors)


#--------------------#
# Plot with ggridges #
#--------------------#
#--- Organize the Data ---#
chainLength = length(posteriors[, 1])



values = c(posteriors$b1.1., posteriors$b1.2., posteriors$b1.3., posteriors$b1.4.)
### OR ###
values = c(posteriors$b1.1., posteriors$b1.2., posteriors$b1.3.)
### OR ###
values = c(posteriors$b1.1., posteriors$b1.2.)


labels = c(rep(CHRNames[1], times = chainLength), rep(CHRNames[2], times = chainLength), 
           rep(CHRNames[3], times = chainLength), rep(CHRNames[4], times = chainLength))

### OR ###
labels = c(rep(CHRNames[1], times = chainLength), rep(CHRNames[2], times = chainLength), 
           rep(CHRNames[3], times = chainLength))

### OR ###
labels = c(rep(CHRNames[1], times = chainLength), rep(CHRNames[2], times = chainLength))



data = data.frame(values, labels)

#--- Plot the data ---#
data$labels <- factor(data$labels , levels=c("High", "Moderate", "Low", "Modifier"))
data$labels <- factor(data$labels , levels=c("Autosomes", "X"))

ggplot(data) +
  theme_bw() +
  geom_density_ridges(aes(x = values, y = labels), fill = "salmon", alpha = 0.6) +
  ylab("Genomic region") +
  xlab("Posterior Probabilities")






