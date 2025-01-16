#parameters
#incubation---------------
# 95% CI [updated using who data]
lower_quantile <- 3
upper_quantile <- 20

# Estimate meanlog and sdlog
meanlog <- (log(lower_quantile) + log(upper_quantile)) / 2
sdlog <- (log(upper_quantile) - log(lower_quantile)) / (2 * qnorm(0.975))

incub_paras=c(meanlog,sdlog)


#infectious-----------------
# 95% CI
lower_quantile <- 14
upper_quantile <- 28

# Estimate meanlog and sdlog
meanlog <- (log(lower_quantile) + log(upper_quantile)) / 2
sdlog <- (log(upper_quantile) - log(lower_quantile)) / (2 * qnorm(0.975))

inf_paras=c(meanlog,sdlog)


#pcr skin-----------------
file_path <- "pcr_data.xlsx"
# Read the specific sheet "Figure 2"
data_pcr <- readxl::read_excel(file_path, sheet = "Figure 2")

days_pcr = as.numeric(data_pcr[[3]][-1]) [-5]
viral_load = as.numeric(data_pcr[[4]][-1]) [-5]
threshold = 4.77

# Create binary outcome (1 if viral load exceeds threshold, 0 otherwise)
outcome <- ifelse(viral_load > threshold, 1, 0)

# Create a data frame
data <- data.frame(days_pcr, outcome)

# Aggregate the data by day
agg_data <- aggregate(outcome ~ days_pcr, data, function(x) c(sum = sum(x), count = length(x)))

# Flatten the aggregated data
agg_data <- do.call(data.frame, agg_data)

agg_data$proportion <- agg_data$outcome.sum / agg_data$outcome.count

# Fit a weighted logistic regression model using the aggregated data
pcr1 <- glm(proportion ~ days_pcr, data = agg_data, family = binomial, weights = outcome.count)


#pcr oral------------------
days_pcr = as.numeric(data_pcr[[15]][-1])
viral_load = as.numeric(data_pcr[[16]][-1])

# Create binary outcome (1 if viral load exceeds threshold, 0 otherwise)
outcome <- ifelse(viral_load > threshold, 1, 0)

# Create a data frame
data <- data.frame(days_pcr, outcome)

# Aggregate the data by day
agg_data <- aggregate(outcome ~ days_pcr, data, function(x) c(sum = sum(x), count = length(x)))

# Flatten the aggregated data
agg_data <- do.call(data.frame, agg_data)

agg_data$proportion <- agg_data$outcome.sum / agg_data$outcome.count

# Fit a weighted logistic regression model using the aggregated data  
pcr2 <- glm(proportion ~ days_pcr, data = agg_data, family = binomial, weights = outcome.count)

#export-----------------
save(incub_paras,inf_paras,pcr1,pcr2, file='paras.RData')

rm(list=ls())