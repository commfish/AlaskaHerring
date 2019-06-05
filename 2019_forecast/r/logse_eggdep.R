
# Estimating log se of egg deposition 
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-06-04

# I don't have the raw data for the egg deposition and am not sure where Steve
# Martell got the log se for the egg deposition data (see the sitka.dat file).
# To estimate the log.se I fit a generalized model to the egg index and assumed
# log.se's for 1971-2015 from Steve, then predicted the log.se values for
# 2016-2018.

eggdep <- read_csv("2019_forecast/data/egg_dep_logse.csv") # from sitka.dat

library(mgcv)
mod <- gam(log_se ~ egg, data = eggdep)
summary(mod)

eggdep %>% arrange(egg)

ggplot(eggdep, aes(x = egg, y = log_se)) +
  geom_point() +
  geom_smooth()

df <- data.frame(year = c(2016, 2017, 2018),
           egg = c(5.979411, 3.617525, 4.262336))

df$log_se <- predict(object = mod, newdata = df)

df
