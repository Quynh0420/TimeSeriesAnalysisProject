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
load("data2.xts.RData")

head(data2.xts)
tail(data2.xts, 12)

# creating first differences of variables
data2.xts$dcom1 <- diff.xts(data2.xts$commod1)
data2.xts$dcom2 <- diff.xts(data2.xts$commod2)

#  plotting variables on the graph 
plot(data2.xts[,c(1,2)],
     col = c("black", "blue"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "Commodities 1 and 2",
     legend.loc = "topright")

plot(data2.xts[,c(1,3)], 
     # plot data in two panels (each column separately)
     multi.panel = 2,
     main = "Original and differenced data for Commodity 1",
     col = c("darkblue", "darkgreen"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     yaxis.same = F, # otherwise scale is the same for each column!
     cex = 1)

plot(data2.xts[,c(2,4)], 
     # plot data in two panels (each column separately)
     multi.panel = 2,
     main = "Original and differenced data for Commodity 2",
     col = c("darkgrey", "orange"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     yaxis.same = F, # otherwise scale is the same for each column!
     cex = 1)

# testing integration order
testdf(variable = data2.xts$commod1,
       max.augmentations = 3)
# p-value is large, so we cannot reject the null about Non-Stationary

testdf(variable = data2.xts$dcom1,
       max.augmentations = 3)
# p-value is small, so we reject the null about Non-stationary
# it is stationary at 1st order with 0 augmentations (adf = -48.74)

testdf(variable = data2.xts$commod2,
       max.augmentations = 3)
# p-value is large, so we cannot reject the null about Non-Stationary

testdf(variable = data2.xts$dcom2,
       max.augmentations = 3)
# p-value is small, so we reject the null about Non-stationary
# it is stationary at 1st order with 0 augmentations (adf = -50.35)

# Both variables are I(1) so we can check 
# whether they are cointegrated.

# Estimating cointegrating vector

model.coint.1 <- lm(commod1 ~ commod2, data = data2.xts)

summary(model.coint.1)

# Testing stationarity of residuals 

testdf(variable = residuals(model.coint.1), 
       max.augmentations = 3)

# The ADF test with 1 augmentations can be used
# its result is that non-stationarity of residuals 
# is REJECTED, so residuals are stationary,
# which means that 2 commodities are cointegrated

# The cointegrating vector is [1, -  6.257909, -0.819131]
# which defines the cointegrating relationship as:
# 1 * commod1 - 6.257909 - 0.819131  * commod2

# Creating first lags of residuals
# and adding it to the dataset

data2.xts$lresid.com <- lag.xts(residuals(model.coint.1))


# Estimating ECM

model.ecm.1 <- lm(dcom1 ~ dcom2 + lresid.com - 1,
                  # -1 means a model without a constant
                  data = data2.xts) 

summary(model.ecm.1)

# How would you interpret results of the model above?

# parameter 0.826820 describes a short term relationship
# between future 2 commodities, so if commodity 2 increases by 1, 
# commodity 1 in the SHORT RUN will increase by 0.826820

# the long run relationship is described by the parameter
# 0.819131 from the cointegrating relationship:
# so if commodity 2 increases by 1 in the LONG RUN commodity 1 will 
# increase by 0.819131

# -0.014008 is the adjustment coefficient
# its sign is negative (as expected)
# its value means that 1.4% of the unexpected 
# error (increase in gap) will be corrected 
# in the next period, so any unexpected deviation
# should be corrected finally within about 71.4 periods (1/1.4%) 

################################################################################
# Commodities - Granger causality test                                            
################################################################################

# Now let's check, whether index Granger causes future and vice versa.
# What is the proper lag length in this case? 

# 3 lags
granger.test(data2.xts[,1:2], # set of variables tested
             3) # lag assumed

# 4 lags
granger.test(data2.xts[,1:2], 4)

# 5 lags
granger.test(data2.xts[,1:2], 5)

# 12 lags
granger.test(data2.xts[,1:2], 12)

# What is the conclusion? 
# At approximately 0% significance level we have so called
# 'bi-directional feedback' in case of 4 lags 
# which means there is Granger causality between 2 commodities
# with other cases (3, 5, 12 lags), the significance level is around 1% to 7%
# maybe this indicates that the length of 4 lags is the most appropriate 
# for the forecast of 2 commodities

################################################################################
# Commodities - VAR model
################################################################################

# Since we observed bi-directional causality 
# we can try to estimate a VAR model

#  First question: what is the proper order (lag length) of the VAR model?

VARselect(data2.xts[,c(1,2)], # input data for VAR
          lag.max = 6)     # maximum lag

# two criteria suggest order 4, one indicates 2 and one indicates 3.
# Remember: the smaller value of information criteria, 
# the better the model.

# in case of potential seasonality
# seasonal dummies should be included
# the option season=seasonal_frequency
# is responsible for that

commodities.var4 <- VAR(data2.xts[,c(1,2)],
                        p = 4) # order of VAR model


summary(commodities.var4)

# CAUTION !!!
# If you are only interested in forecasting, you can run VAR in levels
# even when your series are NON-stationary, but then your standard errors 
# cannot be trusted, so you can't make inference about the values 
# of the coefficients.

# lets do some basic diagnostics
plot(commodities.var4)

# panel 1: original and fitted values
# panel 2: residuals from the model
# panel 3: ACF and PACF for residuals

# looks like lags 6 and 10 of ACF and PACF are still significant 

# lets verify it formally with a multivariate
# autocorrelation test for model residuals
# (Portmanteau test)
# H0: NO autocorrelatoin of residuals

serial.test(commodities.var4)
# Portmanteau test reject the null about no autocorrelation at 5% level


# alternatively we can use a BG (Breusch-Godfrey) test
serial.test(commodities.var4, 
            type = "BG")

# BG test CANNOT reject the null about

# lets extend our model by adding lags 5 and 6

commodities.var6 <- VAR(data2.xts[,c(1,2)],
                        p = 6) 

summary(commodities.var6)

# lets do the diagnostics
plot(commodities.var6)

# ACF and PACF have not improved 

# Portmanteau test
# H0: NO autocorrelatoin of residuals

serial.test(commodities.var6)

serial.test(commodities.var6, type = "BG")

# Portmanteau reject H0 while BG CANNOT reject H0

# maybe we should use more lags at the beginning (10?),
# but it inflates the model with a huge number of parameters...

# and adding 10 lags does not eliminate
# the problem of autocorrelation

commodities.var10 <- VAR(data2.xts[,c(1,2)],
                         p = 10) 

summary(commodities.var10)

serial.test(commodities.var10)

serial.test(commodities.var10, type = "BG")

# lets compare information criteria

AIC(commodities.var4, commodities.var6, commodities.var10)
BIC(commodities.var4, commodities.var6, commodities.var10)

# both AIC and BIC crearly prefers VAR(4)

AIC(commodities.var4s, commodities.var6s, commodities.var10s)
BIC(commodities.var4s, commodities.var6s, commodities.var10s)

################################################################################
#  Commodities - forecasting based on the VAR - VAR(4)
################################################################################

# create shorter samples
commodities <- data2.xts["/2018-04-06",-5]
tail(commodities)

# VAR model estimation on a shorter sample
commodities.var4s <- VAR(commodities[,c(1,2)], p = 4)     

summary(commodities.var4s)  

# and run a forecast
commodities.var4s.forecast <- predict(commodities.var4s,
                                      n.ahead = 10,
                                      ci = 0.95) # 95% confidence interval

# lets see the result
commodities.var4s.forecast

# VAR forecasts for commod 1
commodities.var4s.forecast$fcst$commod1

# VAR and for commod 2
commodities.var4s.forecast$fcst$commod2

# lets store it as an xts object.
# Correct set of dates (index) can be extracted
# from the original xts data object

tail(index(data2.xts), 10)

commod1_forecast_var4s <- xts(commodities.var4s.forecast$fcst$commod1[,-4], 
                              # we exclude the last column with CI
                              tail(index(data2.xts), 10))

# lets change the names 
names(commod1_forecast_var4s)
names(commod1_forecast_var4s) <- c("commod1_fore_var4s", "commod1_lower_var4s", "commod1_upper_var4s")

# lets do the same for commod2 forecasts 
commod2_forecast_var4s <- xts(commodities.var4s.forecast$fcst$commod2[,-4], 
                              # we exclude the last column with CI
                              tail(index(data2.xts), 10))

names(commod2_forecast_var4s) <- c("commod2_fore_var4s", "commod2_lower_var4s", "commod2_upper_var4s")

# lets put the data together

Commodities_var4s <- merge(data2.xts[,c(1,2)], 
                           commod1_forecast_var4s,
                           commod2_forecast_var4s)

# lets compare the forecasted and real data on the plot

# Commod 1 forecast plot

plot(Commodities_var4s["2018-04/", c(1,3,4,5)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 1 - VAR(4)",
     col = c("black", "blue", "red", "red"))


# Commod 2 forecast plot

plot(Commodities_var4s["2018-04/", c(2,6,7,8)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 2 - VAR (4)",
     col = c("black", "blue", "red", "red"))

# lets calculate forecast accuracy measures
# for simplicity of the following formulas
# lets define two new objects:

# lets add the basis for different measures of the forecast error

Commodities_var4s$mae.commod1.var4s  =  abs(Commodities_var4s$commod1 - Commodities_var4s$commod1_fore_var4s)
Commodities_var4s$mape.commod1.var4s  =  abs((Commodities_var4s$commod1 - Commodities_var4s$commod1_fore_var4s)/Commodities_var4s$commod1)


Commodities_var4s$mae.commod2.var4s  =  abs(Commodities_var4s$commod2 - Commodities_var4s$commod2_fore_var4s)
Commodities_var4s$mape.commod2.var4s  =  abs((Commodities_var4s$commod2 - Commodities_var4s$commod2_fore_var4s)/Commodities_var4s$commod2)

tail(Commodities_var4s, 20)

write.csv(Commodities_var4s [c(2129:2138),c(0,1,3,9,10)], file = "Commod1_var4s.csv")
write.csv(Commodities_var4s [c(2129:2138),c(0,2,6,11,12)], file = "Commod2_var4s.csv")

# and calculate its averages

colMeans(Commodities_var4s[,9:12], 
         na.rm = T) #not include missing values in calculation

################################################################################
#  Commodities - forecasting based on the VAR - VAR(6)
################################################################################

# VAR model estimation on a shorter sample
commodities.var6s <- VAR(commodities[,c(1,2)], p = 6)     

summary(commodities.var6s)  

# and run a forecast
commodities.var6s.forecast <- predict(commodities.var6s,
                                      n.ahead = 10,
                                      ci = 0.95) # 95% confidence interval

# lets see the result
commodities.var6s.forecast

# VAR forecasts for commod 1
commodities.var6s.forecast$fcst$commod1

# VAR and for commod 2
commodities.var6s.forecast$fcst$commod2

# lets store it as an xts object.
# Correct set of dates (index) can be extracted
# from the original xts data object

tail(index(data2.xts), 10)

commod1_forecast_var6s <- xts(commodities.var6s.forecast$fcst$commod1[,-4], 
                              # we exclude the last column with CI
                              tail(index(data2.xts), 10))

# lets change the names 
names(commod1_forecast_var6s)
names(commod1_forecast_var6s) <- c("commod1_fore_var6s", "commod1_lower_var6s", "commod1_upper_var6s")

# lets do the same for commod2 forecasts 
commod2_forecast_var6s <- xts(commodities.var6s.forecast$fcst$commod2[,-4], 
                              # we exclude the last column with CI
                              tail(index(data2.xts), 10))

names(commod2_forecast_var6s) <- c("commod2_fore_var6s", "commod2_lower_var6s", "commod2_upper_var6s")

# lets put the data together

Commodities_var6s <- merge(data2.xts[,c(1,2)], 
                           commod1_forecast_var6s,
                           commod2_forecast_var6s)

# lets compare the forecasted and real data on the plot

# Commod 1 forecast plot

plot(Commodities_var6s["2018-04/", c(1,3,4,5)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 1 - VAR(6)",
     col = c("black", "blue", "red", "red"))

# Commod 2 forecast plot

plot(Commodities_var6s["2018-04/", c(2,6,7,8)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 2 - VAR(6)",
     col = c("black", "blue", "red", "red"))

# lets calculate forecast accuracy measures
# for simplicity of the following formulas
# lets define two new objects:

# lets add the basis for different measures of the forecast error

Commodities_var6s$mae.commod1.var6s  =  abs(Commodities_var6s$commod1 - Commodities_var6s$commod1_fore_var6s)
Commodities_var6s$mape.commod1.var6s  =  abs((Commodities_var6s$commod1 - Commodities_var6s$commod1_fore_var6s)/Commodities_var6s$commod1)


Commodities_var6s$mae.commod2.var6s  =  abs(Commodities_var6s$commod2 - Commodities_var6s$commod2_fore_var6s)
Commodities_var6s$mape.commod2.var6s  =  abs((Commodities_var6s$commod2 - Commodities_var6s$commod2_fore_var6s)/Commodities_var6s$commod2)

tail(Commodities_var6s, 20)

write.csv(Commodities_var6s [c(2129:2138),c(0,1,3,9,10)], file = "Commod1_var6s.csv")
write.csv(Commodities_var6s [c(2129:2138),c(0,2,6,11,12)], file = "Commod2_var6s.csv")

# and calculate its averages

colMeans(Commodities_var6s[,9:12], 
         na.rm = T) #not include missing values in calculation

################################################################################
#  Commodities - forecasting based on the VAR - VAR(10)
################################################################################

# VAR model estimation on a shorter sample
commodities.var10s <- VAR(commodities[,c(1,2)], p = 10)     

summary(commodities.var10s)  

# and run a forecast
commodities.var10s.forecast <- predict(commodities.var10s,
                                       n.ahead = 10,
                                       ci = 0.95) # 95% confidence interval

# lets see the result
commodities.var10s.forecast

# VAR forecasts for commod 1
commodities.var10s.forecast$fcst$commod1

# VAR and for commod 2
commodities.var10s.forecast$fcst$commod2

# lets store it as an xts object.
# Correct set of dates (index) can be extracted
# from the original xts data object

tail(index(data2.xts), 10)

commod1_forecast_var10s <- xts(commodities.var10s.forecast$fcst$commod1[,-4], 
                               # we exclude the last column with CI
                               tail(index(data2.xts), 10))

# lets change the names 
names(commod1_forecast_var10s)
names(commod1_forecast_var10s) <- c("commod1_fore_var10s", "commod1_lower_var10s", "commod1_upper_var10s")

# lets do the same for commod2 forecasts 
commod2_forecast_var10s <- xts(commodities.var10s.forecast$fcst$commod2[,-4], 
                               # we exclude the last column with CI
                               tail(index(data2.xts), 10))

names(commod2_forecast_var10s) <- c("commod2_fore_var10s", "commod2_lower_var10s", "commod2_upper_var10s")

# lets put the data together

Commodities_var10s <- merge(data2.xts[,c(1,2)], 
                            commod1_forecast_var10s,
                            commod2_forecast_var10s)

# lets compare the forecasted and real data on the plot

# Commod 1 forecast plot

plot(Commodities_var10s["2018-04/", c(1,3,4,5)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 1 - VAR(10)",
     col = c("black", "blue", "red", "red"))


# Commod 2 forecast plot

plot(Commodities_var10s["2018-04/", c(2,6,7,8)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10-day forecast of commodity 2 - VAR(10)",
     col = c("black", "blue", "red", "red"))

# lets calculate forecast accuracy measures
# for simplicity of the following formulas
# lets define two new objects:

# lets add the basis for different measures of the forecast error

Commodities_var10s$mae.commod1.var10s  =  abs(Commodities_var10s$commod1 - Commodities_var10s$commod1_fore_var10s)
Commodities_var10s$mape.commod1.var10s  =  abs((Commodities_var10s$commod1 - Commodities_var10s$commod1_fore_var10s)/Commodities_var10s$commod1)


Commodities_var10s$mae.commod2.var10s  =  abs(Commodities_var10s$commod2 - Commodities_var10s$commod2_fore_var10s)
Commodities_var10s$mape.commod2.var10s  =  abs((Commodities_var10s$commod2 - Commodities_var10s$commod2_fore_var10s)/Commodities_var10s$commod2)

tail(Commodities_var10s, 20)

write.csv(Commodities_var10s [c(2129:2138),c(0,1,3,9,10)], file = "Commod1_var10s.csv")
write.csv(Commodities_var10s [c(2129:2138),c(0,2,6,11,12)], file = "Commod2_var10s.csv")

# and calculate its averages

colMeans(Commodities_var10s[,9:12], 
         na.rm = T) #not include missing values in calculation

################################################################################
# Plot the forecast of all 3 models on 1 graph
Forecast_VAR_Commod1 <- merge(Commodities_var4s["2018-04/", c(1,3)], 
                              Commodities_var6s["2018-04/", 3],
                              Commodities_var10s["2018-04/", 3])

Forecast_VAR_Commod2 <- merge(Commodities_var4s["2018-04/", c(2,6)], 
                              Commodities_var6s["2018-04/", 6],
                              Commodities_var10s["2018-04/", 6])

write.csv(Forecast_VAR_Commod1, file = "Forecast_VAR_Commod1.csv")
write.csv(Forecast_VAR_Commod2, file = "Forecast_VAR_Commod2.csv")

plot(Forecast_VAR_Commod1, 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "Forecast of 3 VAR models for Commod 1",
     col = c("red", "darkgreen", "darkblue", "darkgray"),
     legend.loc = "topleft")

plot(Forecast_VAR_Commod2, 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "Forecast of 3 VAR models for Commod 2",
     col = c("red", "darkgreen", "darkblue", "darkgray"),
     legend.loc = "topleft")

# compare the results of MAE and MAPE for all 3 models, VAR(6) proved to be
# the best model for forecast, however, the other 2 models' forecast are quite near
# Let's see the MAE and MAPE of 3 models again
# the lower, the better

colMeans(Commodities_var4s[,9:12], 
         na.rm = T)
colMeans(Commodities_var6s[,9:12], 
         na.rm = T)
colMeans(Commodities_var10s[,9:12], 
         na.rm = T)

################################################################################
# Commodities - Johansen cointegration test
################################################################################

johan.test.trace <- ca.jo(data2.xts[,c(1,2)],         # data 
                          ecdet = "const", # "none" for no intercept in cointegrating equation, 
                          # "const" for constant term in cointegrating equation and 
                          # "trend" for trend variable in cointegrating equation
                          type = "trace",  # type of the test: trace or eigen
                          K = 6)         # lag order of the series (levels) in the VAR

summary(johan.test.trace) 


# r = 0: test statistic (13.18) is SMALLER than the critical value (19.96), 
# we CANNOT reject the null
# so, no cointegration vector

# lets apply the alternative variant of the test

johan.test.eigen <- ca.jo(data2.xts[,c(1,2)],         # data 
                          ecdet = "const", # "none" for no intercept in cointegrating equation, 
                          # "const" for constant term in cointegrating equation and 
                          # "trend" for trend variable in cointegrating equation
                          type = "eigen", # type of the test: trace or eigen
                          K = 6)          # lag order of the series (levels) in the VAR


summary(johan.test.eigen) 

# the conclusions are the same:
# no cointegrating vector 
# so basically, using VECM will not improve the forecast of 2 commodities

##################################################################################

######################################################################
# ARIMA FOR COMMODITY 1 FORECAST
# below we apply the Box-Jenkins procedure
#######################################################################
# step 1. INITIAL IDENTIFICATION of parameters p and q 

# lets see ACF and PACF for non-stationary variable
# ACF and PACF are calculated up to 36th lag

# if there are missing values in the data we need 
# to add an additional option na.action = na.pass (see below)

# lets plot them together and limit the scale of ACF
par(mfrow = c(2, 1)) 
  acf(data2.xts$dcom1,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),    # limits for the y axis - we give c(min,max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
  pacf(data2.xts$dcom1, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1)) # we restore the original single panel

# ACF and PACF suggest that maybe ARIMA (10,1,10) could be
#	a sensible model for commod 1, probably without lags 2 to 9

# lets compare different models with AIC criteria

#######################################################################
# steps 2 and 3 interchangeably. MODEL ESTIMATION and DIAGNOSTICS


###############################################################################
# lets start with ARIMA(1,1,1)
# we are using Arima function from the forecast package

arima111.commod1 <- Arima(data2.xts$commod1,  # variable
                          order = c(1, 1, 1))  # (p,d,q) parameters)

# by default the model on differenced data (d = 1) is 
# estimated by without a constant term
# (assuming that first differences fluctuate arround 0)

# lets use coeftest() function from the lmtest package
# to test for significance of model parameters

coeftest(arima111.commod1)

# additional summary measures (eg. information criteria)
summary(arima111.commod1)

# both are highly significant

# however, using that syntax produces a model without 
# a constant term; the constant is included when d = 0

# if one wishes to include a constant also when d = 1
# an additional option include.constant = T has to be used

arima111_2.commod1 <- Arima(data2.xts$commod1,  # variable
                            order = c(1, 1, 1),  # (p,d,q) parameters
                            include.constant = TRUE)  # including a constant

# a constant for a model with d = 1 is reported 
# as a drift parameter

coeftest(arima111_2.commod1)
# drift is insignificant here at 5% level

summary(arima111_2.commod1)

# are residuals of arima111 model white noise? 
# resid() function applied to model results returns residuals

plot(resid(arima111.commod1))

# lets check ACF and PACF

plot_ACF_PACF_resids(arima111.commod1)
# lag 10 is still significant 

# Ljung-Box test (for a maximum of 10 lags)
Box.test(resid(arima111.commod1),
         type = "Ljung-Box", lag = 10)
# p-value is high, cannot reject null, so residual is white noise

# there is also a way to automatically find the best model
arima.best.AIC.commod1 <- auto.arima(data2.xts$commod1,
                                     d = 1,             # parameter d of ARIMA model
                                     max.p = 10,        # Maximum value of p
                                     max.q = 10,        # Maximum value of q
                                     max.order = 14,    # maximum p+q
                                     start.p = 1,       # Starting value of p in stepwise procedure
                                     start.q = 1,       # Starting value of q in stepwise procedure
                                     ic = "aic",        # Information criterion to be used in model selection.
                                     stepwise = FALSE,  # if FALSE considers all models
                                     allowdrift = TRUE, # include a constant
                                     trace = TRUE)      # show summary of all models considered

# however it does not remove intermediate lags

# the result might be surprising

coeftest(arima.best.AIC.commod1)
# ARIMA(0,1,1), lag 1 of MA is also significant at 5% level

AIC(arima.best.AIC.commod1, arima111.commod1)
# AIC better than for Arima 111

BIC(arima.best.AIC.commod1, arima111.commod1)
# BIC better than for Arima 111

# Ljung-Box test
Box.test(resid(arima.best.AIC.commod1),
         type = "Ljung-Box", lag = 10)
# p-value is high, cannot reject null, residuals is white noise

plot_ACF_PACF_resids(arima.best.AIC.commod1)
# we see that lag 10 is still significant,
# so maybe we also include lag 10 in addtion to lag 1 of MA process

arima0110.commod1 <- Arima(data2.xts$commod1,  # variable
                           order = c(0, 1, 10), # (p,d,q) parameters)
                           fixed = c(NA, 0, 0, 0, 0, 0, 0, 0, 0, NA),
                           optim.control = list(maxit = 500),
                           # optimization method:
                           # "Nelder-Mead" (default), "BFGS", "CG", 
                           # "L-BFGS-B", "SANN", "Brent"
                           optim.method = "CG")

coeftest(arima0110.commod1)
# all terms are significant here

# Let's compare again the AIC and BIC
AIC(arima.best.AIC.commod1, arima111.commod1, arima0110.commod1)
# ARIMA 0, 1, 10 has the best AIC

BIC(arima.best.AIC.commod1, arima111.commod1,arima0110.commod1 )
# automatically selected AIC (ARIMA 011) is the best, but ARIMA 0, 1, 10 is quite close

plot_ACF_PACF_resids(arima0110.commod1)
# lag 10 not significant anymore

# Let's test its residual
Box.test(resid(arima0110.commod1),
         type = "Ljung-Box", lag = 15)
# p-value is high, cannot reject null, so residual is white noise

# the same based on BIC
arima.best.BIC.commod1 <- auto.arima(data2.xts$commod1,
                                     d = 1,             # parameter d of ARIMA model
                                     max.p =10,        # Maximum value of p
                                     max.q = 10,        # Maximum value of q
                                     max.order = 14,    # maximum p+q
                                     start.p = 1,       # Starting value of p in stepwise procedure
                                     start.q = 1,       # Starting value of q in stepwise procedure
                                     ic = "bic",        # Information criterion to be used in model selection.
                                     stepwise = FALSE,  # if FALSE considers all models
                                     allowdrift = TRUE, # include a constant
                                     trace = TRUE)      # show summary of all models considered

coeftest(arima.best.BIC.commod1)
# ARIMA(0,1,0) without a constant

BIC(arima.best.BIC.commod1, arima.best.AIC.commod1, arima111.commod1,arima0110.commod1)

# Ljung-Box test
Box.test(resid(arima.best.BIC.commod1),
         type = "Ljung-Box", lag = 10)

# the model selected based on BIC has almost
# autocorrelated residuals
# so should not be considered as complete


# automated procedure is not necessarily better
# than step-by-step human approach

# Lets finally decide to select two models:
# arima111 - sensible, manually selected 
# arima011 - sensible, automatically selected based on AIC
# arima0110 - sensible, manully selected

#######################################################################
# FORECAST for commod 1 - model arima111

# estimate the model on shorter sample
arima111s.commod1 <- Arima(commodities$commod1,  # variable
                           order = c(1,1,1))   # (p,d,q) parameters

coeftest(arima111s.commod1)

# lets make a prediction
forecast111s.commod1  =  forecast::forecast(arima111s.commod1, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast111s.commod1

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111s.commod1$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111s.commod1$mean)

# 80% and 95% confidence intervals
forecast111s.commod1$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111s.commod1_data <- data.frame(f111s.commod1_mean = as.numeric(forecast111s.commod1$mean),
                                        f111s.commod1_lower = as.numeric(forecast111s.commod1$lower),
                                        f111s.comod1_upper = as.numeric(forecast111s.commod1$upper))

head(forecast111s.commod1_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast111s.commod1_xts <- xts(forecast111s.commod1_data,
                                # last 10 values of date 
                                # from the original dataset
                                tail(index(data2.xts), 10))
forecast111s.commod1_xts

# we can put it together with the original data

commod1_arima111s <- merge(data2.xts[,1], forecast111s.commod1_xts)

head(commod1_arima111s)

tail(commod1_arima111s, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(commod1_arima111s["2018-04/",c(1:4)], #to plot only some wanted columns
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 1 - ARIMA111",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod1_arima111s$mae.111s.commod1  =  abs(commod1_arima111s$commod1 - commod1_arima111s$f111s.commod1_mean)
commod1_arima111s$mape.111s.commod1  =  abs((commod1_arima111s$commod1 - commod1_arima111s$f111s.commod1_mean)/commod1_arima111s$commod1)

tail(commod1_arima111s,11)

write.csv(commod1_arima111s[c(2129:2138),], file = "commod1_arima111s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod1_arima111s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 

#######################################################################
# FORECAST for commod 1 - model arima011

# estimate the model on shorter sample
arima011s.commod1 <- Arima(commodities$commod1,  # variable
                           order = c(0,1,1))   # (p,d,q) parameters

coeftest(arima011s.commod1)

# lets make a prediction
forecast011s.commod1  =  forecast::forecast(arima011s.commod1, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast011s.commod1

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast011s.commod1$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast011s.commod1$mean)

# 95% confidence intervals
forecast011s.commod1$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast011s.commod1_data <- data.frame(f011s.commod1_mean = as.numeric(forecast011s.commod1$mean),
                                        f011s.commod1_lower = as.numeric(forecast011s.commod1$lower),
                                        f011s.comod1_upper = as.numeric(forecast011s.commod1$upper))

head(forecast011s.commod1_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast011s.commod1_xts <- xts(forecast011s.commod1_data,
                                # last 10 values of date 
                                # from the original dataset
                                tail(index(data2.xts), 10))
forecast011s.commod1_xts

# we can put it together with the original data

commod1_arima011s <- merge(data2.xts[,1], forecast011s.commod1_xts)

head(commod1_arima011s)

tail(commod1_arima011s, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(commod1_arima011s["2018-04/",c(1:4)], #to plot only some wanted columns
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 1 - ARIMA 011",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod1_arima011s$mae.011s.commod1  =  abs(commod1_arima011s$commod1 - commod1_arima011s$f011s.commod1_mean)
commod1_arima011s$mape.011s.commod1  =  abs((commod1_arima011s$commod1 - commod1_arima011s$f011s.commod1_mean)/commod1_arima011s$commod1)

tail(commod1_arima011s,11)

write.csv(commod1_arima011s[c(2129:2138),], file = "commod1_arima011s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod1_arima011s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 

#######################################################################
# FORECAST for commod 1 - model arima 0,1,10

# estimate the model on shorter sample
arima0110s.commod1 <- Arima(commodities$commod1,  # variable
                            order = c(0,1,10),  # (p,d,q) parameters
                            fixed = c(NA, 0, 0, 0, 0, 0, 0, 0, 0, NA),
                            optim.control = list(maxit = 500),
                            # optimization method:
                            # "Nelder-Mead" (default), "BFGS", "CG", 
                            # "L-BFGS-B", "SANN", "Brent"
                            optim.method = "CG")  

coeftest(arima0110s.commod1)

# lets make a prediction
forecast0110s.commod1  =  forecast::forecast(arima0110s.commod1, # model for prediction
                                             h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast0110s.commod1

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast0110s.commod1$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast0110s.commod1$mean)

# 95% confidence intervals
forecast0110s.commod1$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast0110s.commod1_data <- data.frame(f0110s.commod1_mean = as.numeric(forecast0110s.commod1$mean),
                                         f0110s.commod1_lower = as.numeric(forecast0110s.commod1$lower),
                                         f0110s.comod1_upper = as.numeric(forecast0110s.commod1$upper))

head(forecast0110s.commod1_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast0110s.commod1_xts <- xts(forecast0110s.commod1_data,
                                 # last 10 values of date 
                                 # from the original dataset
                                 tail(index(data2.xts), 10))
forecast0110s.commod1_xts

# we can put it together with the original data

commod1_arima0110s <- merge(data2.xts[,1], forecast0110s.commod1_xts)

head(commod1_arima0110s)

tail(commod1_arima0110s, 10)

# lets finally plot the figure with the forecast
# and original series

# original data
plot(commod1_arima0110s["2018-04/",c(1:4)], #to plot only some wanted columns
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 1 - ARIMA 0,1,10",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod1_arima0110s$mae.0110s.commod1  =  abs(commod1_arima0110s$commod1 - commod1_arima0110s$f0110s.commod1_mean)
commod1_arima0110s$mape.0110s.commod1  =  abs((commod1_arima0110s$commod1 - commod1_arima0110s$f0110s.commod1_mean)/commod1_arima0110s$commod1)

tail(commod1_arima0110s,11)

write.csv(commod1_arima0110s[c(2129:2138),], file = "commod1_arima0110s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod1_arima0110s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 


###########################################################################################
# Plot the forecast of all 3 models on 1 graph
Forecast_ARIMA_Commod1 <- merge(commod1_arima111s["2018-04/",c(1,2)], 
                                commod1_arima011s["2018-04/",2], 
                                commod1_arima0110s["2018-04/",2])

plot(Forecast_ARIMA_Commod1, 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "Forecast of 3 ARIMA models for Commod 1",
     col = c("red", "darkgreen", "darkblue", "gray"),
     legend.loc = "topleft")

# compare the results of MAE and MAPE for all 3 models, ARIMA 011 proved to be
# the best model for forecast, but ARIMA 111 also gives very close result
# Let's see the MAE and MAPE of 3 models again
# the lower, the better

colMeans(commod1_arima111s[,5:6],na.rm = T)
colMeans(commod1_arima011s[,5:6],na.rm = T)
colMeans(commod1_arima0110s[,5:6],na.rm = T)

# Let's compare MAE and MAPE for 
# the best ARIMA model (ARIMA 011) and the best VAR model (Var 6)

colMeans(commod1_arima011s[,5:6],na.rm = T)
colMeans(Commodities_var6s[,9:10], 
         na.rm = T)

# So, ARIMA proves to work better than VAR for Commod 1

######################################################################
# ARIMA FOR COMMODITY 2 FORECAST
# below we apply the Box-Jenkins procedure
#######################################################################
# step 1. INITIAL IDENTIFICATION of parameters p and q 

# lets see ACF and PACF for non-stationary variable
# ACF and PACF are calculated up to 36th lag

# if there are missing values in the data we need 
# to add an additional option na.action = na.pass (see below)

# lets plot them together and limit the scale of ACF
par(mfrow = c(2, 1)) 
  acf(data2.xts$dcom2,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),    # limits for the y axis - we give c(min,max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
  pacf(data2.xts$dcom2, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1)) # we restore the original single panel

# ACF and PACF suggest that maybe ARIMA (6,1,6) could be  
#	a sensible model for commod 2, probably without lags 2,3,4,5

# lets compare different models with AIC criteria

#######################################################################
# steps 2 and 3 interchangeably. MODEL ESTIMATION and DIAGNOSTICS


###############################################################################
# lets start with ARIMA(1,1,1)
# we are using Arima function from the forecast package

arima111.commod2 <- Arima(data2.xts$commod2,  # variable
                          order = c(1, 1, 1))  # (p,d,q) parameters)

# by default the model on differenced data (d = 1) is 
# estimated by without a constant term
# (assuming that first differences fluctuate arround 0)

# lets use coeftest() function from the lmtest package
# to test for significance of model parameters

coeftest(arima111.commod2)

# additional summary measures (eg. information criteria)
summary(arima111.commod2)

# both are insignificant

# however, using that syntax produces a model without 
# a constant term; the constant is included when d = 0

# if one wishes to include a constant also when d = 1
# an additional option include.constant = T has to be used

arima111_2.commod2 <- Arima(data2.xts$commod2,  # variable
                            order = c(1, 1, 1),  # (p,d,q) parameters
                            include.constant = TRUE)  # including a constant

# a constant for a model with d = 1 is reported 
# as a drift parameter

coeftest(arima111_2.commod2)
# drift is insignificant here at 5% level

summary(arima111_2.commod2)

# are residuals of arima111 model white noise? 
# resid() function applied to model results returns residuals

plot(resid(arima111.commod2))

# lets check ACF and PACF

plot_ACF_PACF_resids(arima111.commod2)

# significance of lag 1 disappeared, but 
# there is still lag 6 significant 

# Ljung-Box test (for a maximum of 10 lags)

Box.test(resid(arima111.commod2),
         type = "Ljung-Box", lag = 10)
# p-value is high, null cannot be rejected, residuals are white noise

# there is also a way to automatically find the best model
arima.best.AIC.commod2 <- auto.arima(data2.xts$commod2,
                                     d = 1,             # parameter d of ARIMA model
                                     max.p = 10,        # Maximum value of p
                                     max.q = 10,        # Maximum value of q
                                     max.order = 14,    # maximum p+q
                                     start.p = 1,       # Starting value of p in stepwise procedure
                                     start.q = 1,       # Starting value of q in stepwise procedure
                                     ic = "aic",        # Information criterion to be used in model selection.
                                     stepwise = FALSE,  # if FALSE considers all models
                                     allowdrift = TRUE, # include a constant
                                     trace = TRUE)      # show summary of all models considered

# however it does not remove intermediate lags

# the result might be surprising

coeftest(arima.best.AIC.commod2)
# ARIMA(0,1,1) 

AIC(arima.best.AIC.commod2, arima111.commod2)
# ARIMA (0,1,1) has the best AIC, but the other 2 are quite near

BIC(arima.best.AIC.commod2, arima111.commod2)
# ARIMA (0,1,1) has the best BIC

plot_ACF_PACF_resids(arima.best.AIC.commod2)
# lag 6 is still significant

# Ljung-Box test
Box.test(resid(arima.best.AIC.commod2),
         type = "Ljung-Box", lag = 10)
# p-value is high, cannot reject null, residuals is white noise

# Let's try adding lag 6 for MA in addition to ARIMA 011
arima016.commod2 <- Arima(data2.xts$commod2,
                          order = c(0, 1, 6),
                          fixed = c(NA,0,0,0,0,NA))

coeftest(arima016.commod2)
# all lags are significant

plot(resid(arima016.commod2))

# lets check ACF and PACF

plot_ACF_PACF_resids(arima016.commod2)

# no lags are significant any more

# Ljung-Box test (for a maximum of 10 lags)

Box.test(resid(arima016.commod2),
         type = "Ljung-Box", lag = 10)
# p - value is very hight, cannot reject null, residuals are white noise

AIC(arima.best.AIC.commod2, 
    arima111.commod2, arima016.commod2)
# ARIMA 016 has the best AIC

BIC(arima.best.AIC.commod2, 
    arima111.commod2, arima016.commod2)
# ARIMA (0,1,1) has the best BIC, but ARIMA 016 is also quite close 

# the same based on BIC
arima.best.BIC.commod2 <- auto.arima(data2.xts$commod2,
                                     d = 1,             # parameter d of ARIMA model
                                     max.p =10,        # Maximum value of p
                                     max.q = 10,        # Maximum value of q
                                     max.order = 14,    # maximum p+q
                                     start.p = 1,       # Starting value of p in stepwise procedure
                                     start.q = 1,       # Starting value of q in stepwise procedure
                                     ic = "bic",        # Information criterion to be used in model selection.
                                     stepwise = FALSE,  # if FALSE considers all models
                                     allowdrift = TRUE, # include a constant
                                     trace = TRUE)      # show summary of all models considered

coeftest(arima.best.BIC.commod2)
# Again, ARIMA(0,1,1) without a constant

# automated procedure is not necessarily better
# than step-by-step human approach

# Lets finally decide to select two models:
# arima111 - sensible, manually selected 
# arima011 - sensible, automatically selected based on AIC and BIC
# arima016 - sensible, manually selected

#######################################################################
# FORECAST for commod 2 - model arima111

# estimate the model on shorter sample
arima111s.commod2 <- Arima(commodities$commod2,  # variable
                           order = c(1,1,1))   # (p,d,q) parameters

coeftest(arima111s.commod2)

# lets make a prediction
forecast111s.commod2  =  forecast::forecast(arima111s.commod2, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast111s.commod2

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast111s.commod2$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast111s.commod2$mean)

# 80% and 95% confidence intervals
forecast111s.commod2$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast111s.commod2_data <- data.frame(f111s.commod2_mean = as.numeric(forecast111s.commod2$mean),
                                        f111s.commod2_lower = as.numeric(forecast111s.commod2$lower),
                                        f111s.comod2_upper = as.numeric(forecast111s.commod2$upper))

head(forecast111s.commod2_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast111s.commod2_xts <- xts(forecast111s.commod2_data,
                                # last 10 values of date 
                                # from the original dataset
                                tail(index(data2.xts), 10))
forecast111s.commod2_xts

# we can put it together with the original data

commod2_arima111s <- merge(data2.xts[,2], forecast111s.commod2_xts)

head(commod2_arima111s)

tail(commod2_arima111s, 10)

# lets finally plot the figure with the forecast
# and original series

plot(commod2_arima111s["2018-04/",c(1:4)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 2 - ARIMA111",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod2_arima111s$mae.111s.commod2  =  abs(commod2_arima111s$commod2 - commod2_arima111s$f111s.commod2_mean)
commod2_arima111s$mape.111s.commod2  =  abs((commod2_arima111s$commod2 - commod2_arima111s$f111s.commod2_mean)/commod2_arima111s$commod2)

tail(commod2_arima111s,11)

write.csv(commod2_arima111s[c(2129:2138),], file = "commod2_arima111s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod2_arima111s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 

#######################################################################
# FORECAST for commod 2 - model arima011

# estimate the model on shorter sample
arima011s.commod2 <- Arima(commodities$commod2,  # variable
                           order = c(0,1,1))   # (p,d,q) parameters

coeftest(arima011s.commod2)

# lets make a prediction
forecast011s.commod2  =  forecast::forecast(arima011s.commod2, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast011s.commod2

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast011s.commod2$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast011s.commod2$mean)

# 80% and 95% confidence intervals
forecast011s.commod2$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast011s.commod2_data <- data.frame(f011s.commod2_mean = as.numeric(forecast011s.commod2$mean),
                                        f011s.commod2_lower = as.numeric(forecast011s.commod2$lower),
                                        f011s.comod2_upper = as.numeric(forecast011s.commod2$upper))

head(forecast011s.commod2_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast011s.commod2_xts <- xts(forecast011s.commod2_data,
                                # last 10 values of date 
                                # from the original dataset
                                tail(index(data2.xts), 10))
forecast011s.commod2_xts

# we can put it together with the original data

commod2_arima011s <- merge(data2.xts[,2], forecast011s.commod2_xts)

head(commod2_arima011s)

tail(commod2_arima011s, 10)

# lets finally plot the figure with the forecast
# and original series

plot(commod2_arima011s["2018-04/",c(1:4)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 2 - ARIMA 011",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod2_arima011s$mae.011s.commod2  =  abs(commod2_arima011s$commod2 - commod2_arima011s$f011s.commod2_mean)
commod2_arima011s$mape.011s.commod2  =  abs((commod2_arima011s$commod2 - commod2_arima011s$f011s.commod2_mean)/commod2_arima011s$commod2)

tail(commod2_arima011s,11)

write.csv(commod2_arima011s[c(2129:2138),], file = "commod2_arima011s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod2_arima011s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 

#######################################################################
# FORECAST for commod 2 - model arima016

# estimate the model on shorter sample
arima016s.commod2 <- Arima(commodities$commod2,  # variable
                           order = c(0,1,6),
                           fixed = c(NA,0,0,0,0,NA))   # (p,d,q) parameters

coeftest(arima016s.commod2)

# lets make a prediction
forecast016s.commod2  =  forecast::forecast(arima016s.commod2, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast016s.commod2

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast016s.commod2$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast016s.commod2$mean)

# 80% and 95% confidence intervals
forecast016s.commod2$lower

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast016s.commod2_data <- data.frame(f016s.commod2_mean = as.numeric(forecast016s.commod2$mean),
                                        f016s.commod2_lower = as.numeric(forecast016s.commod2$lower),
                                        f016s.comod2_upper = as.numeric(forecast016s.commod2$upper))

head(forecast016s.commod2_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

forecast016s.commod2_xts <- xts(forecast016s.commod2_data,
                                # last 10 values of date 
                                # from the original dataset
                                tail(index(data2.xts), 10))
forecast016s.commod2_xts

# we can put it together with the original data

commod2_arima016s <- merge(data2.xts[,2], forecast016s.commod2_xts)

head(commod2_arima016s)

tail(commod2_arima016s, 10)

# lets finally plot the figure with the forecast
# and original series

plot(commod2_arima016s["2018-04/",c(1:4)], 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "10 day forecast of commod 2 - ARIMA 016",
     col = c("black", "blue", "red", "red"))

# checking forecast quality 

# lets add the basis for different measures of the forecast error
commod2_arima016s$mae.016s.commod2  =  abs(commod2_arima016s$commod2 - commod2_arima016s$f016s.commod2_mean)
commod2_arima016s$mape.016s.commod2  =  abs((commod2_arima016s$commod2 - commod2_arima016s$f016s.commod2_mean)/commod2_arima016s$commod2)

tail(commod2_arima016s,11)

write.csv(commod2_arima016s[c(2129:2138),], file = "commod2_arima016s.csv")

# and calculate its averages 
# over the forecast period

colMeans(commod2_arima016s[,5:6],na.rm = T)

# NOTICE!!!! Model which performs best in the sample 
#	will NOT NECESSARILY be best for forecasting.
#	One should compare forecasts for several models. 

###########################################################################################
# Plot the forecast of all 3 models on 1 graph
Forecast_ARIMA_Commod2 <- merge(commod2_arima111s["2018-04/",c(1,2)], 
                                commod2_arima011s["2018-04/",2], 
                                commod2_arima016s["2018-04/",2])

plot(Forecast_ARIMA_Commod2, 
     major.ticks = "days", 
     grid.ticks.on = "days",
     grid.ticks.lty = 3,
     main = "Forecast of 3 ARIMA models for Commod 2",
     col = c("red", "darkgreen", "darkblue", "gray"),
     legend.loc = "topleft")

# compare the results of MAE and MAPE for all 3 models, ARIMA 011 proved to be
# the best model for forecast, but ARIMA 111 also gives very close result
# Let's see the MAE and MAPE of 3 models again
# the lower, the better

colMeans(commod2_arima111s[,5:6],na.rm = T)
colMeans(commod2_arima011s[,5:6],na.rm = T)
colMeans(commod2_arima016s[,5:6],na.rm = T)

# Let's compare MAE and MAPE for 
# the best ARIMA model (ARIMA 011) and the best VAR model (Var 6)

colMeans(commod2_arima011s[,5:6],na.rm = T)
colMeans(Commodities_var6s[,11:12], 
         na.rm = T)

# In contrast with commod 1, VAR proves to work better than ARIMA for Commod 2

#######################################################################
# FORECAST for COMMOD 1 - FINAL SELCTION FOR FULL SAMPLE- model arima011

# lets make a prediction
forecast011.commod1  =  forecast::forecast(arima.best.AIC.commod1, # model for prediction
                                            h = 10, level = 0.95) # how many periods outside the sample

# lets see the result
forecast011.commod1

# the forecasts are indexed with 
# an observation number, not a date!

# to reach the point forecast (the first column)
# one needs to use:

forecast011.commod1$mean

# as.numeric() allows to convert it
# to a simple numeric vector

as.numeric(forecast011.commod1$mean)

# if we want to easily put together both real data
# and the forecast on the plot, we have to convert
# both to ts or both to xts objects

# xts objects are more convenient and modern

# we need a data.frame with data
# and index of dates to create an xts object

forecast011.commod1_data <- data.frame(f011.commod1_mean = as.numeric(forecast011.commod1$mean),
                                        f011.commod1_lower = as.numeric(forecast011.commod1$lower),
                                        f011.comod1_upper = as.numeric(forecast011.commod1$upper))

head(forecast011.commod1_data)

# now it is a data frame

# to create an xts object we will
# copy the values of index from 
# the original dataset GSPC

date_as_index_commodities <- as.Date(c("2018-04-23","2018-04-24","2018-04-25","2018-04-26","2018-04-27", 
                                       "2018-04-30","2018-05-01","2018-05-02","2018-05-03","2018-05-04"))

forecast011.commod1_xts <- xts(forecast011.commod1_data,
                               # last 5 values of date 
                               # from the original dataset
                               date_as_index_commodities)
forecast011.commod1_xts

write.csv(forecast011.commod1_xts, file = "forecast011.commod1_xts.csv")

#######################################################################
# FORECAST for COMMOD 2 - FINAL SELCTION FOR FULL SAMPLE- VAR (6)

# and run a forecast
commodities.var6.forecast <- predict(commodities.var6,
                                      n.ahead = 10,
                                      ci = 0.95) # 95% confidence interval

# lets see the result
commodities.var6.forecast

# VAR and for commod 2
commodities.var6s.forecast$fcst$commod2

# lets store it as an xts object.
# Correct set of dates (index) can be extracted
# from the original xts data object

# lets do the same for commod2 forecasts 
commod2_forecast_var6 <- xts(commodities.var6.forecast$fcst$commod2[,-4], 
                              # we exclude the last column with CI
                              date_as_index_commodities)

names(commod2_forecast_var6) <- c("commod2_fore_var6", "commod2_lower_var6", "commod2_upper_var6")

write.csv(commod2_forecast_var6, file = "commod2_forecast_var6.csv")

