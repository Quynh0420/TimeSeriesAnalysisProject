########################################################
#   Home Project - Time Series Analysis
#   Student: QUYNH BUI - 393519
#   Lecturer: Dr. Piotr Wójcik             
#   Faculty of Economic Sciences, University of Warsaw     
#                                                          
########################################################


setwd("D:\\Jeans Witch\\Studies\\MA in Quantitative Finance\\Sem 2\\Time Series Analysis\\Lab assessment\\TSA Home Project - QUYNH BUI (393519) - NGUYEN VO (393518)")

library(xts)
library(vars)
library(forecast)
library(lmtest)
library(urca) 
library(MSBVAR)
library(fBasics)
library(quantmod)
library(fUnitRoots)

source("function_testdf.R")
source("function_plot_ACF_PACF_resids.R")

# penalty on scientific penalty
options(scipen = 20)

# importing data
load("data1.xts.RData")

head(data1.xts)
tail(data1.xts)

# creating first differences of variables
data1.xts$dindex <- diff.xts(data1.xts$index)
data1.xts$dindex.F <- diff.xts(data1.xts$index.F)

#  plotting variables on the graph 
plot(data1.xts[,c(1,2)],
     col = c("black", "blue"),
     major.ticks = "years", 
     grid.ticks.on = "months",
     grid.ticks.lty = 3,
     main = "Index and Future Index",
     legend.loc = "topleft")

plot(data1.xts[,c(1,3)], 
     # plot data in two panels (each column separately)
     multi.panel = 2,
     main = "Original and differenced data for Index",
     col = c("darkblue", "darkgreen"),
     major.ticks = "years", 
     grid.ticks.on = "months",
     grid.ticks.lty = 3,
     yaxis.same = F, # otherwise scale is the same for each column!
     cex = 1)

plot(data1.xts[,c(2,4)], 
     # plot data in two panels (each column separately)
     multi.panel = 2,
     main = "Original and differenced data for Future",
     col = c("darkgrey", "orange"),
     major.ticks = "years", 
     grid.ticks.on = "months",
     grid.ticks.lty = 3,
     yaxis.same = F, # otherwise scale is the same for each column!
     cex = 1)

# testing integration order
testdf(variable = data1.xts$index,
       max.augmentations = 3)
# p-value is large, so we cannot reject the null about Non-Stationary

testdf(variable = data1.xts$dindex,
       max.augmentations = 3)
# p-value is small, so we reject the null about Non-stationary
# it is stationary at 1st order with 0 augmentations (adf = -22.00)

testdf(variable = data1.xts$index.F,
       max.augmentations = 3)
# p-value is large, so we cannot reject the null about Non-Stationary

testdf(variable = data1.xts$dindex.F,
       max.augmentations = 3)
# p-value is small, so we reject the null about Non-stationary
# it is stationary at 1st order with 0 augmentations (adf = -21.73)

# Both variables are I(1) so we can check 
# whether they are cointegrated.

# Estimating cointegrating vector

model.coint <- lm(index.F ~ index, data = data1.xts)

summary(model.coint)

# Testing stationarity of residuals 

testdf(variable = residuals(model.coint), 
       max.augmentations = 3)

# The ADF test with 1 augmentations can be used
# its result is that non-stationarity of residuals 
# is STRONGLY REJECTED, so residuals are stationary,
# which means that index and future index are cointegrated

# The cointegrating vector is [1, - 13.044676, -0.825113]
# which defines the cointegrating relationship as:
# 1 * index.F - 13.044676- 0.825113 * index

# Creating first lags of residuals
# and adding it to the dataset

data1.xts$lresid <- lag.xts(residuals(model.coint))


# Estimating ECM

model.ecm <- lm(dindex.F ~ dindex + lresid - 1,
                # -1 means a model without a constant
                data = data1.xts) 

summary(model.ecm)


# How would you interpret results of the model above?

# parameter 0.814743 describes a short term relationship
# between future index and index, so if index increases by 1, 
# future index in the SHORT RUN will increase by 0.814743

# the long run relationship is described by the parameter
# 0.825113 from the cointegrating relationship:
# so if index increases by 1 in the LONG RUN future index will 
# increase by 0.825113

# -0.236247 is the adjustment coefficient
# its sign is negative (as expected)
# its value means that 23.6% of the unexpected 
# error (increase in gap) will be corrected 
# in the next period, so any unexpected deviation
# should be corrected finally within about 4.25 periods (1/23.6%) # testing integration order

################################################################################
# index_future - Granger causality test                                            
################################################################################

# Now let's check, whether index Granger causes future and vice versa.
# What is the proper lag length in this case? 

# 3 lags
granger.test(data1.xts[,1:2], # set of variables tested
             3) # lag assumed

# 4 lags
granger.test(data1.xts[,1:2], 4)

# 5 lags
granger.test(data1.xts[,1:2], 5)

# 12 lags
granger.test(data1.xts[,1:2], 12)

# What is the conclusion? 
# At 5% significance level we DO NOT have so called
# 'bi-directional feedback' in all cases
# it means there is no Granger causality between index and future index

######################################################################
# ARIMA FOR INDEX FORECAST
# below we apply the Box-Jenkins procedure
#######################################################################

# step 1. INITIAL IDENTIFICATION of parameters p and q 

# lets see ACF and PACF for non-stationary variable
# ACF and PACF are calculated up to 36th lag

par(mfrow = c(2, 1)) 
  Acf(data1.xts$dindex,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),    # limits for the y axis - we give c(min,max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
  Pacf(data1.xts$dindex, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1)) # we restore the original single panel

# ACF and PACF suggest that maybe ARIMA (3,1,3) could be
#	a sensible model for the Index, probably without 
# lag 1 for AR and MA

# lets compare different models with AIC criteria

#######################################################################
# steps 2 and 3 interchangeably. MODEL ESTIMATION and DIAGNOSTICS


###############################################################################

# lets start with ARIMA(1,1,1)
# we are using Arima function from the forecast package

arima111.index <- Arima(data1.xts$index,  # variable
                        order = c(1, 1, 1),  # (p,d,q) parameters
                        include.constant = TRUE)  # including a constant

# a constant for a model with d = 1 is reported 
# as a drift parameter

coeftest(arima111.index)

# it is statistically significant here at 5% level except for the drift

summary(arima111.index)

# maybe we should try to run a model without the drift

arima111.index <- Arima(data1.xts$index,  # variable
                        order = c(1, 1, 1))  # (p,d,q) parameters


coeftest(arima111.index)

# are residuals of arima111 model white noise? 
# resid() function applied to model results returns residuals

plot(resid(arima111.index))

# lets check ACF and PACF

plot_ACF_PACF_resids(arima111.index)

# no lags are significant

# Ljung-Box test (for a maximum of 10 lags)

Box.test(resid(arima111.index),
         type = "Ljung-Box", lag = 10)

# at 5% we cannot reject the null about residuals being 
# white noise, p-value is very high!

###############################################################################

# lets try ARIMA(3,1,3)

arima313.index <- Arima(data1.xts$index, 
                        order = c(3, 1, 3), 
                        include.constant = TRUE)

coeftest(arima313.index)

# No lags are significant

summary(arima313.index)

# Ljung-Box test for autocorrelation of model residuals

Box.test(resid(arima313.index),
         type = "Ljung-Box",lag = 10)
# null cannot be rejected - p-value very high!

# Plot ACF and PACF for residuals

plot_ACF_PACF_resids(arima313.index)

# lags of ACF and PACF are not significant

###############################################################################
# lets try ARIMA(3,1,3) model without lag 1 and drift

arima313_3.index <- Arima(data1.xts$index,
                          order = c(3, 1, 3),
                          fixed = c(0, NA, NA,
                                    0, NA, NA))


coeftest(arima313_3.index)

# Still, no lags are significant

# lets check if residuals are white noise
plot_ACF_PACF_resids(arima313_3.index)

# Ljung-Box test
Box.test(resid(arima313_3.index),
         type = "Ljung-Box",lag = 10)

# the null cannot be rejected
# residuals are white noise

###############################################################################
# lets compare AIC for all models estimated so far
# CAUTION! for some of them rediduals are not white noise!

# Based on AIC which model is best? 

AIC(arima111.index, arima313.index, arima313_3.index)
# arima313 without lag 1 and drift, but arima111 is also quite near
# The lower AIC, the better

# lets do the same for BIC
BIC(arima111.index, arima313.index, arima313_3.index)
# arima111, 
# the lower BIC, the better

# from the perspective of sensibility arima111 and arima313_3 seems 
# to be the best (residuals are white noise and low IC,
# but arima111 has all terms signifcant 
# while arima313 has no terms significant

# there is also a way to automatically find the best model
arima.best.AIC.index <- auto.arima(data1.xts$index,
                                   d = 1,             # parameter d of ARIMA model
                                   max.p = 5,        # Maximum value of p
                                   max.q = 5,        # Maximum value of q
                                   max.order = 10,    # maximum p+q
                                   start.p = 1,       # Starting value of p in stepwise procedure
                                   start.q = 1,       # Starting value of q in stepwise procedure
                                   ic = "aic",        # Information criterion to be used in model selection.
                                   stepwise = FALSE,  # if FALSE considers all models
                                   allowdrift = TRUE, # include a constant
                                   trace = TRUE)      # show summary of all models considered

# however it does not remove intermediate lags

# the result might be surprising

coeftest(arima.best.AIC.index)
# ARIMA(4,1,0) 
# however lag 1 and 4 seem insignificant

AIC(arima.best.AIC.index, arima111.index, arima313.index, arima313_3.index)
# AIC better than for the best manually selected model

BIC(arima.best.AIC.index, arima111.index, arima313.index, arima313_3.index)
# BIC worse than for the best manually selected model

# Ljung-Box test
Box.test(resid(arima.best.AIC.index),
         type = "Ljung-Box", lag = 10)
# residuals are also white noise (p-value is very high!)

# but the automated procedure does not 
# exclude intermediate lags

# the same based on BIC
arima.best.BIC.index <- auto.arima(data1.xts$index,
                                   d = 1,             # parameter d of ARIMA model
                                   max.p = 5,        # Maximum value of p
                                   max.q = 5,        # Maximum value of q
                                   max.order = 10,    # maximum p+q
                                   start.p = 1,       # Starting value of p in stepwise procedure
                                   start.q = 1,       # Starting value of q in stepwise procedure
                                   ic = "bic",        # Information criterion to be used in model selection.
                                   stepwise = FALSE,  # if FALSE considers all models
                                   allowdrift = TRUE, # include a constant
                                   trace = TRUE)      # show summary of all models considered

coeftest(arima.best.BIC.index)

# ARIMA(0,1,0) without a constant

BIC(arima.best.BIC.index, arima.best.AIC.index, arima111.index)

# Ljung-Box test
Box.test(resid(arima.best.BIC.index),
         type = "Ljung-Box", lag = 10)

# the model selected based on BIC has almost
# autocorrelated residuals
# so should not be considered as complete

# automated procedure is not necessarily better
# than step-by-step human approach

# Lets finally decide to select 3 models:
# arima111 - sensible, manually selected based on BIC
# arima313_3 - sensible, manually selected based on AIC 
# ARIMA(4,1,0) - automated selection based on AIC

# for further comparisons (forecasts)

# these models have:
# - lowest information criteria (AIC or BIC)
# - their residuals are white noise
# - (almost) all parameters significant

# We will finally use these models for forecasting

#######################################################################
# FORECAST for INDEX - model arima111

# create shorter sample
index_future <- data1.xts["/2018-03-16", -5]
tail(index_future)

# estimate the model on shorter sample
arima111s.index <- Arima(index_future$index,  # variable
                   order = c(1,1,1))   # (p,d,q) parameters
                 
coeftest(arima111s.index)

# lets make a prediction
forecast111.index  =  forecast::forecast(arima111s.index, h = 5, 
                                         level = 0.95) 

# lets see the result
forecast111.index

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111.index$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111.index$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111.index_data <- data.frame(f111.index_mean = as.numeric(forecast111.index$mean),
                               f111.index_lower = as.numeric(forecast111.index$lower),
                               f111.index_upper = as.numeric(forecast111.index$upper))

head(forecast111.index_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast111.index_xts <- xts(forecast111.index_data,
                       # last 5 values of date 
                       # from the original dataset
                       tail(index(data1.xts), 5))
forecast111.index_xts

# we can put it together with the original data

Index_111 <- merge(data1.xts[,1], forecast111.index_xts)

head(Index_111)

tail(Index_111, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Index_111["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "5 weeks forecast of Index - ARIMA111",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Index.r <- tail(data1.xts$index, 5)  

# forecasts
Index.f.111 <- as.numeric(forecast111.index$mean)

# put it together into a data.frame
Index_forecast_111 <- data.frame(Index.r, Index.f.111)
Index_forecast_111

# lets add the basis for different measures of the forecast error
Index_forecast_111$mae.index.111  =  abs(Index.r - Index.f.111)
Index_forecast_111$mape.index.111  =  abs((Index.r - Index.f.111)/Index.r)

Index_forecast_111

write.csv(Index_forecast_111, file = "Index_forecast_111.csv")
# to export xts object as excel file

# and calculate its averages 
# over the forecast period

colMeans(Index_forecast_111[,3:4])

#######################################################################
# FORECAST for INDEX - model arima313_3 (without lag 1 and drift)

# estimate the model on shorter sample
arima313_3s.index <- Arima(index_future$index,  # variable
                   order = c(3,1,3),   # (p,d,q) parameters
                   fixed = c(0, NA, NA,
                             0, NA,NA))

coeftest(arima313_3s.index)

# lets make a prediction
forecast313_3.index  =  forecast::forecast(arima313_3s.index, 
                                           h = 5, level = 0.95) 
# there were 2 functions "forecast" in 2 different packages, so not to be confused, 
# we should write "Package name :: function name"

# lets see the result
forecast313_3.index

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast313_3.index$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast313_3.index$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast313_3.index_data <- data.frame(f313_3.index_mean = as.numeric(forecast313_3.index$mean),
                                     f3131_3.index_lower = as.numeric(forecast313_3.index$lower),
                                     f313_3.index_upper = as.numeric(forecast313_3.index$upper))

head(forecast313_3.index_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast313_3.index_xts <- xts(forecast313_3.index_data,
                             # last 5 values of date 
                             # from the original dataset
                             tail(index(data1.xts), 5))
forecast313_3.index_xts

# we can put it together with the original data

Index_313_3 <- merge(data1.xts[,1], forecast313_3.index_xts)

head(Index_313_3)

tail(Index_313_3, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Index_313_3["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3, 
     main = "5 weeks forecast of Index - ARIMA313 without lags 1 and drift",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Index.r <- tail(data1.xts$index, 5)  

# forecasts
Index.f.313_3 <- as.numeric(forecast313_3.index$mean)

# put it together into a data.frame
Index_forecast_313_3 <- data.frame(Index.r, Index.f.313_3)
Index_forecast_313_3

# lets add the basis for different measures of the forecast error
Index_forecast_313_3$mae.index.313_3  =  abs(Index.r - Index.f.313_3)
Index_forecast_313_3$mape.index.313_3  =  abs((Index.r - Index.f.313_3)/Index.r)

Index_forecast_313_3

write.csv(Index_forecast_313_3, file = "Index_forecast_313_3.csv")

# and calculate its averages 
# over the forecast period

colMeans(Index_forecast_313_3[,3:4])

#######################################################################
# FORECAST for INDEX - model arima410

# estimate the model on shorter sample
arima410s.index <- Arima(index_future$index,  # variable
                         order = c(4,1,0))   # (p,d,q) parameters

coeftest(arima410s.index)


# lets make a prediction
forecast410.index  =  forecast::forecast(arima410s.index, 
                                         h = 5, level = 0.95) 

# lets see the result
forecast410.index

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast410.index$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast410.index$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast410.index_data <- data.frame(f410.index_mean = as.numeric(forecast410.index$mean),
                                     f410.index_lower = as.numeric(forecast410.index$lower),
                                     f410.index_upper = as.numeric(forecast410.index$upper))

head(forecast410.index_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast410.index_xts <- xts(forecast410.index_data,
                             # last 5 values of date 
                             # from the original dataset
                             tail(index(data1.xts), 5))
forecast410.index_xts

# we can put it together with the original data

Index_410 <- merge(data1.xts[,1], forecast410.index_xts)

head(Index_410)

tail(Index_410, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Index_410["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "5 weeks forecast of Index - ARIMA410",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Index.r <- tail(data1.xts$index, 5)  

# forecasts
Index.f.410 <- as.numeric(forecast410.index$mean)

# put it together into a data.frame
Index_forecast_410 <- data.frame(Index.r, Index.f.410)
Index_forecast_410

# lets add the basis for different measures of the forecast error
Index_forecast_410$mae.index.410  =  abs(Index.r - Index.f.410)
Index_forecast_410$mape.index.410  =  abs((Index.r - Index.f.410)/Index.r)

Index_forecast_410

write.csv(Index_forecast_410, file = "Index_forecast_410.csv")

# and calculate its averages 
# over the forecast period

colMeans(Index_forecast_410[,3:4])

###########################################################################################
# Plot the forecast of all 3 models on 1 graph
Forecast_ARIMA_Index <- merge(Index_111["2018",c(1,2)], 
                              Index_313_3["2018",2], Index_410["2018",2])

write.csv(Forecast_ARIMA_Index, file = "Forecast_ARIMA_Index.csv")

plot(Forecast_ARIMA_Index, 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "Forecast of 3 ARIMA models for Index",
     col = c("red", "darkgreen", "darkblue", "gray"),
     legend.loc = "topright")

# compare the results of MAE and MAPE for all 3 models, ARIMA111 proved to be
# the best model for forecast
# Let's see the MAE and MAPE of 3 models again
# the lower, the better

colMeans(Index_forecast_111[,3:4])
colMeans(Index_forecast_313_3[,3:4])
colMeans(Index_forecast_410[,3:4])


######################################################################
# ARIMA FOR FUTURE FORECAST
# below we apply the Box-Jenkins procedure
#######################################################################

# step 1. INITIAL IDENTIFICATION of parameters p and q 

# lets see ACF and PACF for non-stationary variable
# ACF and PACF are calculated up to 36th lag

par(mfrow = c(2, 1)) 
  Acf(data1.xts$dindex.F,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),    # limits for the y axis - we give c(min,max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
  Pacf(data1.xts$dindex.F, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1)) # we restore the original single panel

# ACF and PACF suggest that maybe ARIMA (3,1,3) could be
#	a sensible model for the Index, probably without lag 1

# lets compare different models with AIC criteria

###############################################################################
# steps 2 and 3 interchangeably. MODEL ESTIMATION and DIAGNOSTICS


###############################################################################

# lets start with ARIMA(1,1,1)
# we are using Arima function from the forecast package

arima111.Future <- Arima(data1.xts$index.F,  # variable
                         order = c(1, 1, 1),  # (p,d,q) parameters
                         include.constant = TRUE)  # including a constant

# a constant for a model with d = 1 is reported 
# as a drift parameter

coeftest(arima111.Future)

# it is statistically significant here at 5% level except for the drift

summary(arima111.Future)

# maybe we should try to run a model without the drift

arima111.Future <- Arima(data1.xts$index.F,  # variable
                         order = c(1, 1, 1))  # (p,d,q) parameters


coeftest(arima111.Future)

# are residuals of arima111 model white noise? 
# resid() function applied to model results returns residuals

plot(resid(arima111.Future))

# lets check ACF and PACF

plot_ACF_PACF_resids(arima111.Future)

# no lags are significant

# Ljung-Box test (for a maximum of 10 lags)

Box.test(resid(arima111.Future),
         type = "Ljung-Box", lag = 10)

# at 5% we cannot reject the null about residuals being 
# white noise, p-value is very high!

###############################################################################
# lets try ARIMA(3,1,3) without the drift

arima313.Future <- Arima(data1.xts$index.F, 
                         order = c(3, 1, 3)) 


coeftest(arima313.Future)

# no terms are significant

summary(arima313.Future)

# Ljung-Box test for autocorrelation of model residuals

Box.test(resid(arima313.Future),
         type = "Ljung-Box",lag = 10)
# null cannot be rejected - p-value very high!

# Plot ACF and PACF for residuals

plot_ACF_PACF_resids(arima313.Future)

# lags of ACF and PACF are not significant

###############################################################################
# lets try ARIMA(3,1,3) model without lag 1 and drift

# intermediate lags can be set to 0 by using the 
# fixed= argument
# CAUTION! The order of parameters can be checked by:

coefficients(arima313.Future)

# lets add restrictions on ar1 and ma1 (were not significant,
# so we assume ar1 = 0 and ma1 = 0)
# CAUTION! vector provided in fixed =   must have length
# equal to the number of parameters!

arima313_1.Future <- Arima(data1.xts$index.F,
                           order = c(3, 1, 3),
                           fixed = c(0, NA, NA, # vector of length
                                     0, NA, NA)) # equal to total number of parameters,
# NA means no restriction on a parameter

coeftest(arima313_1.Future)

# Still, no terms are significant here

# lets check if residuals are white noise
plot_ACF_PACF_resids(arima313_1.Future)

# Ljung-Box test
Box.test(resid(arima313_1.Future),
         type = "Ljung-Box",lag = 10)

# the null cannot be rejected (p-value very high!)
# residuals are white noise

###############################################################################
# lets compare AIC for all models estimated so far
# CAUTION! for some of them rediduals are not white noise!

# Based on AIC which model is best? 

AIC(arima111.Future, arima313.Future, arima313_1.Future)

# arima313 without lag 1 and drift, but arima111 is also very near
# The lower AIC, the better

# lets do the same for BIC
BIC(arima111.Future, arima313.Future, arima313_1.Future)

# arima111, 
# the lower BIC, the better

# from the perspective of sensibility arima111 seems 
# to be the best (all terms significant, residuals 
# are white noise and low IC)

# there is also a way to automatically find the best model
arima.best.AIC.Future <- auto.arima(data1.xts$index.F,
                                    d = 1,             # parameter d of ARIMA model
                                    max.p = 5,        # Maximum value of p
                                    max.q = 5,        # Maximum value of q
                                    max.order = 10,    # maximum p+q
                                    start.p = 1,       # Starting value of p in stepwise procedure
                                    start.q = 1,       # Starting value of q in stepwise procedure
                                    ic = "aic",        # Information criterion to be used in model selection.
                                    stepwise = FALSE,  # if FALSE considers all models
                                    allowdrift = TRUE, # include a constant
                                    trace = TRUE)      # show summary of all models considered

# however it does not remove intermediate lags

# the result might be surprising

coeftest(arima.best.AIC.Future)
# ARIMA(5,1,3) 
# however no lags are significant

AIC(arima.best.AIC.Future, arima111.Future, arima313.Future, arima313_1.Future)
# AIC worse than for the best manually selected model

BIC(arima.best.AIC.Future, arima111.Future, arima313.Future, arima313_1.Future)
# BIC worse than for the best manually selected model

# Ljung-Box test
Box.test(resid(arima.best.AIC.Future),
         type = "Ljung-Box", lag = 10)
# residuals are also white noise (p-value is very high!)

# but the automated procedure does not 
# exclude intermediate lags

# the same based on BIC
arima.best.BIC.Future <- auto.arima(data1.xts$index.F,
                                    d = 1,             # parameter d of ARIMA model
                                    max.p = 5,        # Maximum value of p
                                    max.q = 5,        # Maximum value of q
                                    max.order = 10,    # maximum p+q
                                    start.p = 1,       # Starting value of p in stepwise procedure
                                    start.q = 1,       # Starting value of q in stepwise procedure
                                    ic = "bic",        # Information criterion to be used in model selection.
                                    stepwise = FALSE,  # if FALSE considers all models
                                    allowdrift = TRUE, # include a constant
                                    trace = TRUE)      # show summary of all models considered

coeftest(arima.best.BIC.Future)

# ARIMA(0,1,0) without a constant

BIC(arima.best.BIC.index, arima.best.AIC.Future, arima111.Future)
# BIC worse than best manually selected model (arima111)

# Ljung-Box test
Box.test(resid(arima.best.BIC.Future),
         type = "Ljung-Box", lag = 10)

# the model selected based on BIC has almost
# autocorrelated residuals
# so should not be considered as complete

# automated procedure is not necessarily better
# than step-by-step human approach

# Lets finally decide to select 3 models:
# arima111 - sensible, manually selected based on BIC
# arima313_1 - sensible, manually selected based on AIC
# ARIMA(5,1,3) - automated selection based on AIC

# for further comparisons (forecasts)

# these models have:
# - lowest information criteria (AIC or BIC)
# - their residuals are white noise
# - (almost) all parameters significant

# We will finally use these models for forecasting

#################################################################################
# FORECAST for FUTURE - model arima111

# estimate the model on shorter sample
arima111s.future <- Arima(index_future$index.F,  # variable
                         order = c(1,1,1))   # (p,d,q) parameters

coeftest(arima111s.future)

# lets make a prediction
forecast111.Future  =  forecast::forecast(arima111s.future, h = 5, level = 0.95) 

# lets see the result
forecast111.Future

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111.Future$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111.Future$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111.Future_data <- data.frame(f111.Future_mean = as.numeric(forecast111.Future$mean),
                                     f111.Future_lower = as.numeric(forecast111.Future$lower),
                                     f111.Future_upper = as.numeric(forecast111.Future$upper))

head(forecast111.Future_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast111.Future_xts <- xts(forecast111.Future_data,
                             # last 5 values of date 
                             # from the original dataset
                             tail(index(data1.xts), 5))
forecast111.Future_xts

# we can put it together with the original data

Future_111 <- merge(data1.xts[,2], forecast111.Future_xts)

head(Future_111)

tail(Future_111, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Future_111["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "5 weeks forecast of Future - ARIMA111",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Future.r <- tail(data1.xts$index.F, 5)  

# forecasts
Future.f.111 <- as.numeric(forecast111.Future$mean)

# put it together into a data.frame
Future_forecast_111 <- data.frame(Future.r, Future.f.111)
Future_forecast_111

# lets add the basis for different measures of the forecast error
Future_forecast_111$mae.Future.111  =  abs(Future.r - Future.f.111)
Future_forecast_111$mape.Future.111  =  abs((Future.r - Future.f.111)/Future.r)

Future_forecast_111

write.csv(Future_forecast_111, file = "Future_forecast_111.csv")

# and calculate its averages 
# over the forecast period

colMeans(Future_forecast_111[,3:4])

#################################################################################
# FORECAST for FUTURE - model arima313_1

# estimate the model on shorter sample
arima313_1s.future <- Arima(index_future$index.F,
                         order = c(3, 1, 3),
                         fixed = c(0, NA, NA, # vector of length
                                   0, NA, NA))  # (p,d,q) parameters

coeftest(arima313_1s.future)


# lets make a prediction
forecast313_1.Future  =  forecast::forecast(arima313_1s.future, h = 5, level = 0.95) 

# lets see the result
forecast313_1.Future

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast313_1.Future$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast313_1.Future$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast313_1.Future_data <- data.frame(f313_1.Future_mean = as.numeric(forecast313_1.Future$mean),
                                        f313_1.Future_lower = as.numeric(forecast313_1.Future$lower),
                                        f313_1.Future_upper = as.numeric(forecast313_1.Future$upper))

head(forecast313_1.Future_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast313_1.Future_xts <- xts(forecast313_1.Future_data,
                                # last 5 values of date 
                                # from the original dataset
                                tail(index(data1.xts), 5))
forecast313_1.Future_xts

# we can put it together with the original data

Future_313_1 <- merge(data1.xts[,2], forecast313_1.Future_xts)

head(Future_313_1)

tail(Future_313_1, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Future_313_1["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "5 weeks forecast of Future - ARIMA313",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Future.r <- tail(data1.xts$index.F, 5)  

# forecasts
Future.f.313_1 <- as.numeric(forecast313_1.Future$mean)

# put it together into a data.frame
Future_forecast_313_1 <- data.frame(Future.r, Future.f.313_1)
Future_forecast_313_1

# lets add the basis for different measures of the forecast error
Future_forecast_313_1$mae.Future.313_1  =  abs(Future.r - Future.f.313_1)
Future_forecast_313_1$mape.Future.313_1  =  abs((Future.r - Future.f.313_1)/Future.r)

Future_forecast_313_1

write.csv(Future_forecast_313_1, file = "Future_forecast_313_1.csv")

# and calculate its averages 
# over the forecast period

colMeans(Future_forecast_313_1[,3:4])

#################################################################################
# FORECAST for FUTURE - model arima513

# estimate the model on shorter sample
arima513s.future <- Arima(index_future$index.F,  # variable
                          order = c(5,1,3))   # (p,d,q) parameters

coeftest(arima513s.future)

# lets make a prediction
forecast513.Future  =  forecast::forecast(arima513s.future, h = 5, level = 0.95) 

# lets see the result
forecast513.Future

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast513.Future$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast513.Future$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast513.Future_data <- data.frame(f513.Future_mean = as.numeric(forecast513.Future$mean),
                                      f513.Future_lower = as.numeric(forecast513.Future$lower),
                                      f513.Future_upper = as.numeric(forecast513.Future$upper))

head(forecast513.Future_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast513.Future_xts <- xts(forecast513.Future_data,
                              # last 5 values of date 
                              # from the original dataset
                              tail(index(data1.xts), 5))
forecast513.Future_xts

# we can put it together with the original data

Future_513 <- merge(data1.xts[,2], forecast513.Future_xts)

head(Future_513)

tail(Future_513, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(Future_513["2018",], 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "5 weeks forecast of Future - ARIMA513",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# for simplicity of the formulas
# lets define two new objects:

# real values - last 5 observations
Future.r <- tail(data1.xts$index.F, 5)  

# forecasts
Future.f.513 <- as.numeric(forecast513.Future$mean)

# put it together into a data.frame
Future_forecast_513 <- data.frame(Future.r, Future.f.513)
Future_forecast_513

# lets add the basis for different measures of the forecast error
Future_forecast_513$mae.Future.513  =  abs(Future.r - Future.f.513)
Future_forecast_513$mape.Future.513  =  abs((Future.r - Future.f.513)/Future.r)

Future_forecast_513

write.csv(Future_forecast_513, file = "Future_forecast_513.csv")

# and calculate its averages 
# over the forecast period

colMeans(Future_forecast_513[,3:4])

###########################################################################################
# Plot the forecast of all 3 models on 1 graph
Forecast_ARIMA_Future <- merge(Future_111["2018",c(1,2)], 
                               Future_313_1["2018",2],
                               Future_513["2018",2])

plot(Forecast_ARIMA_Future, 
     major.ticks = "weeks", 
     grid.ticks.on = "weeks",
     grid.ticks.lty = 3,
     main = "Forecast of 3 ARIMA models for Future",
     col = c("red", "darkgreen", "darkblue","darkgray"),
     legend.loc = "topright")

# compare the results of MAE and MAPE for all 3 models, 
# ARIMA111 proved to be the best model for forecast
# Let's see the MAE and MAPE of 3 models again
# the lower, the better

colMeans(Future_forecast_111[,3:4])
colMeans(Future_forecast_313_1[,3:4])
colMeans(Future_forecast_513[,3:4])

#######################################################################
# FORECAST for INDEX - FINAL SELCTION FOR FULL SAMPLE- model arima111

# lets make a prediction
forecast111.f.index  =  forecast::forecast(arima111.index, h = 5, 
                                           level = 0.95) 

# lets see the result
forecast111.f.index

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111.f.index$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111.f.index$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111.f.index_data <- data.frame(f111.f.index_mean = as.numeric(forecast111.f.index$mean),
                                     f111.f.index_lower = as.numeric(forecast111.f.index$lower),
                                     f111.f.index_upper = as.numeric(forecast111.f.index$upper))

head(forecast111.f.index_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

date_as_index <- as.Date(c("2018-04-27", "2018-05-04", "2018-05-11", 
                           "2018-05-18", "2018-05-25"))

forecast111.f.index_xts <- xts(forecast111.f.index_data,
                             # last 5 values of date 
                             # from the original dataset
                             date_as_index)
forecast111.f.index_xts

write.csv(forecast111.f.index_xts, file = "forecast111.f.index_xts.csv")

#################################################################################
# FORECAST for FUTURE - FINAL SELCTION FOR FULL SAMPLE- model arima111

# lets make a prediction on full sample period
forecast111.f.Future  =  forecast::forecast(arima111.Future, h = 5, 
                                            level = 0.95) 

# lets see the result
forecast111.f.Future

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111.f.Future$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111.f.Future$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111.f.Future_data <- data.frame(f111.f.Future_mean = as.numeric(forecast111.f.Future$mean),
                                      f111.f.Future_lower = as.numeric(forecast111.f.Future$lower),
                                      f111.f.Future_upper = as.numeric(forecast111.f.Future$upper))

head(forecast111.f.Future_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

date_as_index <- as.Date(c("2018-04-27", "2018-05-04", "2018-05-11", 
                        "2018-05-18", "2018-05-25"))

forecast111.f.Future_xts <- xts(forecast111.f.Future_data,
                              # last 5 values of date 
                              # from the original dataset
                              date_as_index)
forecast111.f.Future_xts

write.csv(forecast111.f.Future_xts, file = "forecast111.f.Future_xts.csv")
