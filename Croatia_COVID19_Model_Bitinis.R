#-------------------------------------------------------------------------------
# ASSIGNMENT - INFECTIOUS DISEASE MATHEMATICAL MODELLING
# GRADUATE STUDENT ILIAS BITINIS
#-------------------------------------------------------------------------------
# COU
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#--------------------------DETERMINISTIC----------------------------------------
#-------------------------------------------------------------------------------


################################################################################
###  Q1 - Construct and visualize the epidemic curve for Coratia using the first
###  16 days of reported COVID-19 cases
################################################################################




# Load required R packages for data analysis, modeling and visualization
library(readxl)
library(dplyr)
library(ggplot2)


# Import COVID-19 case data from EXCEl file
data <- read_excel("data.xlsx")


# Select and preprocess data corresponding to Croatia only
croatia_data <- data %>%
  filter(countriesAndTerritories == "Croatia") %>%
  mutate(dateRep = as.Date(dateRep, format = "%d/%m/%Y")) %>%
  arrange(dateRep)


# Identify  the date of the first confirmed COVID-19 case in Croatia
first_case_date <- croatia_data %>%
  filter(cases > 0) %>%
  summarise(min_date = min(dateRep)) %>%
  pull(min_date)
first_case_date


# Filter the datasheet to to include only the first 16 days after the first reported case
croatia_16days <- croatia_data %>%
  filter(dateRep >= first_case_date & dateRep <= first_case_date + 15)


# Epidemic curve
ggplot(croatia_16days, aes(x = dateRep, y = cases)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_line(color = "darkred") +
  geom_point(color = "darkred", size = 2) +
  labs(title = "Epidemic Curve (16 days)",
       x = "Date",
       y = "Cases") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




################################################################################
### Q2 - Estimate the time-varying reproduction number (Rt) using the Cori et al.
### method based on early epidemic data
################################################################################



# Load required R packages
library(EpiEstim)

# Data
incidence <- croatia_16days$cases
dates <- croatia_16days$dateRep
data_epi <- data.frame(dates = dates, I = incidence)

# Rt estimation 
t_start <- seq(2, length(incidence) - 6)
t_end <- t_start + 6

res <- estimate_R(incid = data_epi,
                  method = "parametric_si",
                  config = make_config(list(
                    mean_si = 5.2,
                    std_si = 2.8,
                    t_start = t_start,
                    t_end = t_end
                  )))

# Generate a plot showing the estimated Rt over time 
plot(res, "R") +
  labs(title = "Estimation of Rt using the Cori et al. method")





############################################
####### Q3 - Building he SEIR Model ########
############################################
# Assume a latent period of 4 days and an average infectious period of 7 days
# The population is closed and starts with 10 newly exposed individuals on day 0
# No immunity is assumed, and no interventions are applied (baseline scenario)

# Load required R packages
library(deSolve)

# Population of Croatia
N <- 3860000

# Parameters
head(res$R)
R0 <- res$R$'Mean(R)'[1]
R0
sigma <- 1/4             
gamma <- 1/7         
beta1<- R0*gamma       


# Define initial state variables for the SEIR model: Susceptible, Exposed, 
# Infectious and Recovered
init <- c(S = N - 10,      
          E = 10,          
          I = 0,           
          R = 0)          

# Simulation duration 4 months 
days <- seq(0, 120, by = 1)



# SEIR model
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta1 * S * I / N
    dE <- beta1 * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

# Parameters for the model
params <- c(beta = beta, sigma = sigma, gamma = gamma)

# Use differential equations to simulate the SEIR model over time
output <- ode( y = init, times = days, func = seir_model, parms = params)
output_df <- as.data.frame(output)


# Visualize the number of Exposed (E) over time from the SEIR model output
ggplot(output_df, aes(x = time, y = E)) +
  geom_line(color = "darkorange", size = 1.2) +
  labs(title = "Number of Exposed - SEIR Simulation (Baseline)",
       x = "Days", y = "Exposed (Ε)") +
  theme_minimal()


# Cumulative infected (I + R) and %
output_df$cumulative_infected <- output_df$I + output_df$R
output_df$percent_infected <- 100*(output_df$cumulative_infected / N)


# Plot of the cumulative number of individuals infected during the simulation period
ggplot(output_df, aes(x = time, y = cumulative_infected)) +
  geom_line(color = "darkblue", size = 1.2) +
  labs(title = "Cumulative Number of Infected (SEIR simulation)",
       x = "Days", y = "Infected (I + R)") +
  theme_minimal()


# Plot of the percentage of the total population infected
ggplot(output_df, aes(x = time, y = percent_infected)) +
  geom_line(color = "darkgreen", size = 1.2) +
  labs(title = "Percentage of the total population infected (SEIR simulation)",
       x = "Days", y = "Percentage (%)") +
  theme_minimal()


# Estimation of the final size of the epidemic
final_recovered <- tail(output_df$R, 1)
final_size_percent <- (final_recovered / N) * 100

cat("Final size of the epidemic (people):", round(final_recovered), "\n")
cat("Final % of infected population:", round(final_size_percent, 2), "%\n")





####################################################################
###  Q4 - Different approach of the SEIR model with intervention ###
####################################################################
# Assume that social distancing measures are implemented from the second month 
# onward, reducing the average contact rate by 40%.



# Change of beta because of reduction of social contacts
beta_reduced <- beta1 * 0.6


# Implement an intervention scenario with reduced transmission starting at a 
# specific time point

seir_intervention <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    current_beta <- if (time < 60) beta1 else beta_reduced
    
    dS <- -current_beta * S * I / N
    dE <- current_beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
  
}


output1 <- ode(y = init, times = days, func = seir_intervention, parms = NULL)
output_df1<- as.data.frame(output1)

# Plot of Exposed (Ε)
ggplot(output_df1, aes(x = time, y = E)) +
  geom_line(color = "darkorange", size = 1.2) +
  labs(title = "Exposed with intervention from day 61",
       x = "Days", y = "Exposed (E)") +
  theme_minimal()

baseline_df <- as.data.frame(output) %>%
  mutate(scenario = "Baseline")

intervention_df <- output_df1 %>%
  mutate(scenario = "Intervention")


df_all <- rbind(baseline_df, intervention_df)


ggplot(df_all, aes(x = time, y = E, color = scenario)) +
  geom_line(size = 1.2) +
  labs(title = "Baseline vs Intervention",
       x = "Days", y = "Exposed (Ε)", color = "Scenario") +
  theme_minimal()





#####################################################################
### Q5 - Extended SEIR model to account population heterogeneity  ###
#####################################################################
# Children (0–17): 24 daily contacts (18 with children, 5 with adults)
# Adults (18–64): 12 daily contacts (4 with children, 6 with adults)
# Seniors (65+): 5 daily contacts (1 with children, 2 with adults)
# 20% of the population are children, 30% are seniors
# On day 0, 5 children and 5 adults are exposed
# Equal probability of transmission per contact (5%)


# Population
N_total <- N

children <- 0.20 * N_total
adults <- 0.50 * N_total
seniors <- 0.30 * N_total

# Initial state values
initial_state_values <- c(
  S1 <- children - 5, E1 = 5, I1 = 0, R1 = 0,  #Children
  S2 <- adults - 5, E2 = 5, I2 = 0, R2 = 0,    #Adults
  S3 <- seniors, E3 = 0, I3 = 0, R3 = 0        #Seniors
)


# Parameters with contact matrix
parameters <- c(
  p = 0.05,                           
  c11 = 18, c12 = 5, c13 = 1,         
  c21 = 4, c22 = 6, c23 = 2,          
  c31 = 1, c32 = 2, c33 = 2,          
  sigma = 1/4,                        
  gamma = 1/7                         
)


# SEIR model with 3 age groups
seir_3groups <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
  
    N1 <- S1 + E1 + I1 + R1
    N2 <- S2 + E2 + I2 + R2
    N3 <- S3 + E3 + I3 + R3
    
    #Force of infectious for each group
    lambda1 <- p * (c11 * I1 / N1 + c12 * I2 / N2 + c13 * I3 / N3)
    lambda2 <- p * (c21 * I1 / N1 + c22 * I2 / N2 + c23 * I3 / N3)
    lambda3 <- p * (c31 * I1 / N1 + c32 * I2 / N2 + c33 * I3 / N3)
    
    
    #Differential equations
    #children
    dS1 <- - lambda1 * S1
    dE1 <- lambda1 * S1 - sigma * E1
    dI1 <- sigma * E1 - gamma * I1
    dR1 <- gamma * I1
    
    #adults
    dS2 <- - lambda2 * S2
    dE2 <- lambda2 * S2 - sigma * E2
    dI2 <- sigma * E2 - gamma * I2
    dR2 <- gamma * I2
    
    #seniors
    dS3 <- - lambda3 * S3
    dE3 <- lambda3 * S3 - sigma * E3
    dI3 <- sigma * E3 - gamma * I3
    dR3 <- gamma * I3
    
    res <- c(dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dS3, dE3, dI3, dR3)
    
    return(list(res))
  })
}

# solving the model
output2 <- ode(y = initial_state_values,
               times = days,
               func = seir_3groups,
               parms = parameters)
out2 <- as.data.frame(output2)
colnames(out2) <- c("time",
                    "S1", "E1", "I1", "R1",
                    "S2", "E2", "I2", "R2",
                    "S3", "E3", "I3", "R3")


# Plot of Infectious (I) for each group
ggplot(out2) +
  geom_line(aes(x = time, y = I1, color = "Children")) +
  geom_line(aes(x = time, y = I2, color = "Adults")) +
  geom_line(aes(x = time, y = I3, color = "Seniors")) +
  labs(x = "Days", y = "Infectious", title = "Infectious for each age group",
       color = "Group") +
  theme_minimal()













#-------------------------------------------------------------------------------
#------------------------------STOCHASTIC---------------------------------------
#-------------------------------------------------------------------------------





#################################
### Q6 - Chain Binomial Model ###
#################################


library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(R2OpenBUGS)
library(coda)
library(bayesplot)



database <- read_excel("data.xlsx")

# Select and preprocess data corresponding to Croatia only 
database_croatia <- database %>%
  filter(countriesAndTerritories == "Croatia") %>%
  mutate(date = as.Date(dateRep, format = "%d/%m/%Y")) %>%
  filter(date >= as.Date("2020-09-01") & date <= as.Date("2021-06-20"))%>%
  arrange(date)

# Daily cases
cases_daily <- database_croatia$cases
days2 <- length(cases_daily)

# daily removals
removals_daily <- c(rep(0, 7), head(cases_daily, -7))

# S, I, R
N <- 3860000
S <- numeric(days2)
I <- numeric(days2)
R <- numeric(days2)

S[1] <- N - cases_daily[1]
I[1] <- cases_daily[1]
R[1] <- 0

for (t in 2:days2) {
  S[t] <- S[t-1] - cases_daily[t-1]
  I[t] <- I[t-1] + cases_daily[t-1] - removals_daily[t-1]
  R[t] <- R[t-1] + removals_daily[t-1]
}

# Full data frame
full_database_croatia <- data.frame(
  Date <- database_croatia$date,
  Cases <- cases_daily,
  Removals <- removals_daily,
  S = S,
  I = I,
  R = R
)

# Data for OpenBUGS
data_bugs <- list(
  n_obs = length(full_database_croatia$Cases....cases_daily),
  n_pop = N,
  new_cases = full_database_croatia$Cases....cases_daily,
  new_removals = full_database_croatia$Removals....removals_daily
)

# Chain Bionomial model
model_code <- "
model {
S0 <- n_pop - new_cases[1]
I0 <- new_cases[1]

p[1] <- 1 - exp(-(beta * I0 / n_pop))
new_cases[1] ~ dbin(p[1], S0)
S[1] <- S0 - new_cases[1]
I[1] <- I0 + new_cases[1] - new_removals[1]

for (t in 2:n_obs){
p[t] <- 1 - exp(-(beta * I[t-1] / n_pop))
new_cases[t] ~ dbin(p[t], S[t-1])
S[t] <- S[t-1] - new_cases[t]
I[t] <- I[t-1] + new_cases[t] - new_removals[t]
}
beta ~ dlnorm(0, 5)
}"

writeLines(model_code, con = "chain_binomial_model.txt")

# Initial state values
inits2 <- function() {
  list(beta = runif(1, 0.05, 0.07))
}

params <- c("beta")

Sys.setenv(OpenBUGS_PATH = "C:/Program Files/OpenBUGS323")

# Run MCMC in OpenBUGS
mcmc_fit_bugs <- bugs(
  data = data_bugs,
  inits = inits2,
  parameters.to.save = params,
  model.file = "chain_binomial_model.txt",
  n.chains = 3,
  n.iter = 10500,
  n.burnin = 500,
  n.thin = 1,
  debug = TRUE
)

print(mcmc_fit_bugs, digits = 3)

mcmc_results <- as.mcmc.list(mcmc_fit_bugs)

mcmc_trace(mcmc_results, pars = "beta") + ggtitle("Trace plot for Beta")
mcmc_dens(mcmc_results, pars = "beta") + ggtitle("Density plot for Beta")
mcmc_acf(mcmc_results, pars = "beta")





####################################
### Q7 - Change Point in day 150 ###
####################################



change_point <- 150

# Data for OpenBugs
data_bugs_change <- list(
  n_obs = length(full_database_croatia$Cases....cases_daily),
  n_pop = N,
  new_cases = full_database_croatia$Cases....cases_daily,
  new_removals = full_database_croatia$Removals....removals_daily,
  change_point = change_point
)


# New Chain Binomial model with change in beta
model_code_change <- "
model {
S0 <- n_pop - new_cases[1]
I0 <- new_cases[1]

beta_t[1] <- beta_1 * step(change_point - 1) + beta_2 * (1 - step(change_point - 1))
p[1] <- 1 - exp(-(beta_t[1] * I0 / n_pop))
new_cases[1] ~ dbin(p[1], S0)

S[1] <- S0 - new_cases[1]
I[1] <- I0 + new_cases[1] - new_removals[1]

for (t in 2:n_obs) {
  beta_t[t] <- beta_1 * step(change_point - t) + beta_2 * (1 - step(change_point - t))
  p[t] <- 1 - exp(-(beta_t[t] * I[t-1] / n_pop))
  new_cases[t] ~ dbin(p[t], S[t-1])
  S[t] <- S[t-1] - new_cases[t]
  I[t] <- I[t-1] + new_cases[t] - new_removals[t]
}

beta_1 ~ dlnorm(0,5)
beta_2 ~ dlnorm(0,5)
}
"

writeLines(model_code_change, con = "chain_binomial_change_model.txt")

# Initial state values
inits3 <- function() {
  list(beta_1 = runif(1, 0.05, 0.1),
       beta_2 = runif(1, 0.01, 0.05))
}

params_change <- c("beta_1", "beta_2")

# MCMC model
mcmc_fit_change <- bugs(
  data = data_bugs_change,
  inits = inits3,
  parameters.to.save = params_change,
  model.file = "chain_binomial_change_model.txt",
  n.chains = 3,
  n.iter = 10500,
  n.burnin = 500,
  n.thin = 1,
  debug = TRUE
)

print(mcmc_fit_change, digits = 3)

mcmc_results_change <- as.mcmc.list(mcmc_fit_change)


mcmc_trace(mcmc_results_change, pars = c("beta_1", "beta_2")) + ggtitle("Trace Plot for beta1 and beta2")
mcmc_dens(mcmc_results_change, pars = c("beta_1", "beta_2")) + ggtitle("Density Plot for beta1 and beta2")
mcmc_acf(mcmc_results_change, pars = c("beta_1", "beta_2"))





#################################################
### Q8 -  2 change points in days 100 and 200 ###
#################################################


change_point1 <- 100
change_point2 <- 200

# Data for OpenBUGS
data_bugs_2changes <- list(
  n_obs = length(full_database_croatia$Cases....cases_daily),
  n_pop = N,
  new_cases = full_database_croatia$Cases....cases_daily,
  new_removals = full_database_croatia$Removals....removals_daily,
  change_point1 = change_point1,
  change_point2= change_point2
)


# Chain binomial model with 3 periods of infection
model_code_2changes <- "
model {
S0 <- n_pop - new_cases[1]
I0 <- new_cases[1]

beta_t[1] <- beta_1
p[1] <- 1 - exp(-(beta_t[1] * I0 / n_pop))
new_cases[1] ~ dbin(p[1], S0)

S[1] <- S0 - new_cases[1]
I[1] <- I0 + new_cases[1] - new_removals[1]

for (t in 2:n_obs) {
  is_cp1[t] <- step(change_point1 + 0.5 - t)
  is_cp2[t] <- step(change_point2 + 0.5 - t) * (1 - is_cp1[t])
  is_cp3[t] <- 1 - is_cp1[t] - is_cp2[t]

  beta_t[t] <- is_cp1[t] * beta_1 + is_cp2[t] * beta_2 + is_cp3[t] * beta_3

  p[t] <- 1 - exp(-(beta_t[t] * I[t-1] / n_pop))
  new_cases[t] ~ dbin(p[t], S[t-1])
  S[t] <- S[t-1] - new_cases[t]
  I[t] <- I[t-1] + new_cases[t] - new_removals[t]
}

  beta_1 ~ dlnorm(0, 5)
  beta_2 ~ dlnorm(0, 5)
  beta_3 ~ dlnorm(0, 5)
}
"


writeLines(model_code_2changes, con = "chain_binomial_2changes_model.txt")


# initial state values 
inits_2changes <- function() {
  list(
    beta_1 = runif(1, 0.05, 0.08),
    beta_2 = runif(1, 0.03, 0.06),
    beta_3 = runif(1, 0.01, 0.05)
  )
}



params_2changes <- c("beta_1", "beta_2", "beta_3")


# MCMC model
mcmc_fit_2changes <- bugs(
  data = data_bugs_2changes,
  inits = inits_2changes,
  parameters.to.save = params_2changes,
  model.file = "chain_binomial_2changes_model.txt",
  n.chains = 3,
  n.iter = 10500,
  n.burnin = 500,
  n.thin = 1,
  debug = TRUE
)

print(mcmc_fit_2changes, digits = 3)

mcmc_results_2changes <- as.mcmc.list(mcmc_fit_2changes)



mcmc_trace(mcmc_results_2changes, pars = c("beta_1", "beta_2", "beta_3")) +
  ggtitle("Trace Plots for beta_1, beta_2 and beta_3")


mcmc_dens_overlay(mcmc_results_2changes, pars = c("beta_1", "beta_2", "beta_3")) +
  ggtitle("Density Plots for beta_1, beta_2 and beta_3")

mcmc_acf(mcmc_results_2changes, pars = c("beta_1", "beta_2", "beta_3"))










############################################
### Q9 - Estimation of the  change point ###
############################################


n_obs <- length(full_database_croatia$Cases....cases_daily)

data_bugs_cp <- list(
  n_obs = n_obs,
  n_obs_minus_30 = n_obs - 30,
  n_pop = N,
  new_cases = full_database_croatia$Cases....cases_daily,
  new_removals = full_database_croatia$Removals....removals_daily
)



model_cp_code <- "
model {
  S0 <- n_pop - new_cases[1]
  I0 <- new_cases[1]

  for (t in 1:n_obs) {
    is_before[t] <- step(change_point - t)
    is_after[t] <- 1 - is_before[t]
    beta_t[t] <- beta_1 * is_before[t] + beta_2 * is_after[t]
  }

  p[1] <- 1 - exp(-beta_t[1] * I0 / n_pop)
  new_cases[1] ~ dbin(p[1], S0)
  S[1] <- S0 - new_cases[1]
  I[1] <- I0 + new_cases[1] - new_removals[1]

  for (t in 2:n_obs) {
    p[t] <- 1 - exp(-beta_t[t] * I[t-1] / n_pop)
    new_cases[t] ~ dbin(p[t], S[t-1])
    S[t] <- S[t-1] - new_cases[t]
    I[t] <- I[t-1] + new_cases[t] - new_removals[t]
  }

  beta_1 ~ dlnorm(0, 5)
  beta_2 ~ dlnorm(0, 5)
  change_point ~ dunif(30, n_obs_minus_30)
}
"



inits_cp <- function() {
  list(
    beta_1 = runif(1, 0.05, 0.08),
    beta_2 = runif(1, 0.01, 0.05),
    change_point = sample(30:(n_obs - 30), 1)
  )
}

params_cp <- c("beta_1", "beta_2", "change_point")


writeLines(model_cp_code, con = "chain_binomial_estimated_cp.txt")


mcmc_cp <- bugs(
  data = data_bugs_cp,
  inits = inits_cp,
  parameters.to.save = params_cp,
  model.file = "chain_binomial_estimated_cp.txt",
  n.chains = 3,
  n.iter = 10500,
  n.burnin = 500,
  n.thin = 1,
  debug = TRUE
)


mcmc_cp_mcmc <- as.mcmc.list(mcmc_cp)

print(mcmc_cp, digits = 3)



mcmc_trace(mcmc_cp_mcmc, pars = c("beta_1", "beta_2", "change_point")) +
  ggtitle("Trace Plots for beta_1, beta_2 and change_point")


mcmc_dens_overlay(mcmc_cp_mcmc, pars = c("beta_1", "beta_2", "change_point")) +
  ggtitle("Density Plots for beta_1, beta_2 and change_point")


mcmc_acf(mcmc_cp_mcmc, pars = c("beta_1", "beta_2", "change_point"))






#-------------------------------------------------------------------------------
###########################################
### Q10 - Estimation of 2 change points ###
###########################################
data_bugs_2cp <- list(
  n_obs = length(full_database_croatia$Cases),
  n_obs_minus_30 = length(full_database_croatia$Cases) - 30,
  n_pop = N,
  new_cases = full_database_croatia$Cases,
  new_removals = full_database_croatia$Removals
)

model_2cp_code <- "
model {
  S0 <- n_pop - new_cases[1]
  I0 <- new_cases[1]

  for (t in 1:n_obs) {
   is_cp1[t] <- step(change_point1 + 0.5 - t)
  is_cp2[t] <- step(change_point2 + 0.5 - t) * (1 - is_cp1[t])
  is_cp3[t] <- 1 - is_cp1[t] - is_cp2[t]

  beta_t[t] <- is_cp1[t] * beta_1 + is_cp2[t] * beta_2 + is_cp3[t] * beta_3


  }

  p[1] <- 1 - exp(-beta_t[1] * I0 / n_pop)
  new_cases[1] ~ dbin(p[1], S0)
  S[1] <- S0 - new_cases[1]
  I[1] <- I0 + new_cases[1] - new_removals[1]

  for (t in 2:n_obs) {
    p[t] <- 1 - exp(-beta_t[t] * I[t-1] / n_pop)
    new_cases[t] ~ dbin(p[t], S[t-1])
    S[t] <- S[t-1] - new_cases[t]
    I[t] <- I[t-1] + new_cases[t] - new_removals[t]
  }

  beta_1 ~ dlnorm(0, 5)
  beta_2 ~ dlnorm(0, 5)
  beta_3 ~ dlnorm(0, 5)
  change_point1 ~ dunif(30, n_obs_minus_30)
  change_point2 ~ dunif(30, n_obs_minus_30)
}
"

writeLines(model_2cp_code, con = "chain_binomial_estimated_2cp.txt")


inits_2cp <- function() {
  cp1_val <- sample(30:(n_obs - 60), 1)  
  cp2_val <- sample((cp1_val + 10):(n_obs - 30), 1)  
  
  list(
    beta_1 = runif(1, 0.05, 0.08),
    beta_2 = runif(1, 0.03, 0.06),
    beta_3 = runif(1, 0.01, 0.05),
    change_point1 = cp1_val,
    change_point2 = cp2_val
  )
}


params_2cp <- c("beta_1", "beta_2", "beta_3", "change_point1", "change_point2")




mcmc_2cp <- bugs(
  data = data_bugs_2cp,
  inits = inits_2cp,
  parameters.to.save = params_2cp,
  model.file = "chain_binomial_estimated_2cp.txt",
  n.chains = 3,
  n.iter = 10500,
  n.burnin = 500,
  n.thin = 1,
  debug = TRUE
)

mcmc_2cp_mcmc <- as.mcmc.list(mcmc_2cp)
print(mcmc_2cp, digits = 3)


mcmc_trace(mcmc_2cp_mcmc, 
           pars = c("beta_1", "beta_2", "beta_3", "change_point1", "change_point2")) +
  ggtitle("Trace Plots for beta_1, beta_2, beta_3, change_point1 and change_point2")


mcmc_dens_overlay(mcmc_2cp_mcmc, 
                  pars = c("beta_1", "beta_2", "beta_3", "change_point1", "change_point2")) +
  ggtitle("Density Plots for beta_1, beta_2, beta_3, change_point1 and change_point2")


mcmc_acf(mcmc_2cp_mcmc, 
         pars = c("beta_1", "beta_2", "beta_3", "change_point1", "change_point2"))


