# Apply STARS method to PDO to determine breakpoints for time-varying survival,
# maturity, and selectivity

# Authors: Sara E Miller, Jane Sullivan
# Contact: sara.miller@alaska.gov, jane.sullivan1@alaska.gov
# Last edited: August 2019

# The mean of the monthly pdo values from April of the previous year through
# March of the current year were used as an annual pdo index from 1979 to 2019.
# For example, the pdo index for 2019 was the mean monthly pdo value from April
# 2018 through March 2019.

#R equivalent for v3-2 stars excel add-in:
#http://esapubs.org/archive/ecol/E095/262/suppl-1.php

#Rodionov, S. N. 2004. A sequential algorithm for testing climate regime shifts.
#Geophysical Research Letters 31:L09204.

#Rodionov, S. N. 2006. Use of prewhitening in climate regime shift detection.
#Geophysical Research Letters 33:L12707.

#Seddon, A. W. R., C. A. Froyd, M. J. Leng, G. A. Milne, and K. J. Willis.
#2011a. Ecosystem resilience and threshold response in the Galápagos coastal
#zone. PloS one 6:e22376.

# Updated excel version:
# https://sites.google.com/site/climatelogic/documentation/red-noise-estimation

# Set-up ----
START_YEAR <- 1929 # Model start year
YEAR <- 2019 # Forecast year

source(paste0(YEAR, "_forecast/r/helper.r"))
source(paste0(YEAR, "_forecast/r/stars.r"))


pdo_path <- paste0(YEAR, "_forecast/results/stars_pdo")
if(!exists(pdo_path)){dir.create(pdo_path)}

# Access PDO from JISAO ----
raw_data<-read.csv('https://oceanview.pfeg.noaa.gov/erddap/tabledap/cciea_OC_PDO.csv', fill=TRUE, header=TRUE) 
write.csv(raw_data, paste0(pdo_path, "/raw_pdo_accessed_", Sys.Date(), ".csv")) 

# Summarize PDO ----
files <- list.files(file.path(pdo_path)) # searchable file list
pdo_file <- file.path(pdo_path, files[which(grepl("raw_pdo", files))]) # find .tpl

raw_pdo <- read.csv(pdo_file)

raw_pdo %>%
  drop_na() %>%
  select (-c(X)) %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time)) %>%
  select (-c(time))  %>%
  filter(year >= START_YEAR & year < YEAR) %>% 
  mutate(year = as.numeric(year))  %>% 
  mutate(year = ifelse(month %in% c(1, 2, 3), year, year + 1),
         pdo_shifted = PDO) %>%
  filter(year > START_YEAR-1 & year < YEAR) %>% 
  select (-c(YEAR, PDO, month)) %>% 
  as.data.frame() -> pdo 

# Summarize and save
pdo %>%
  group_by(year) %>%
  dplyr::summarise(mean_pdo = mean(pdo_shifted)) %>%
  write_csv(., paste0(pdo_path, "/pdo_summary.csv"))

#set-up data for STARS method function
pdo_sum <-  read.csv(paste0(pdo_path, "/pdo_summary.csv"), header = TRUE) 
pdo_data <- pull(pdo_sum, mean_pdo)
names(pdo_data) <- pull(pdo_sum, year)

# Standardize to zscores - deprecated? Sara, can this be deleted now?
# the Stars code in R standardizes but Rodionov's original output does not
# pdo_data <- pdo_data-mean(pdo_data, na.rm=T) pdo_data <- pdo_data/
# sd(pdo_data, na.rm=T)

# STARS analysis ----
pdo_stars <- stars(pdo_data, 
                   L = 10, # L = cut-off length
                   p = 0.10, # p = significance
                   h = 1, # h = Huber's weight parameter
                   AR1red = "IP4", # AR1red  = "est", "MPK", "IP4", or "none" (red noise estimation)
                   prewhitening = FALSE # prewhitening = T or F 
                   )

# FLAG - what's this? m = round((L+1)/3); formula to estimate subsample size for
# calculating alpha

pdo_results <- as.data.frame(pdo_stars$starsResult) 

pdo_results %>%
  mutate(year = (START_YEAR):(YEAR-1)) -> pdo_results

# STARS assigns the first year a mean value of NA. Because it will belong to the
# next year's time block, assign it the mean value for the next year. We
# recalculate mean below.
pdo_results$mean[1] <- pdo_results$mean[2]

# identify breaks
pdo_results %>% 
  rename(mean_pdo = mean) %>% 
  group_by(mean_pdo) %>% 
  mutate(min_yr = min(year), # there should be no NAs
         max_yr = max(year),
         median_yr = median(year),
         `Time blocks` = paste0(min_yr, "-", max_yr)) %>% 
  ungroup() -> pdo_results

pdo_results %>% distinct(mean_pdo, min_yr, max_yr, `Time blocks`)
write_csv(pdo_results, paste0(pdo_path, "/pdo_stars_summary.csv"))

# Graphics ----
xaxis <- tickr(pdo_results, year, 10)

# Summarize PDO for figure
pdo %>% 
  group_by(year) %>% 
  mutate(PDO = mean(pdo_shifted, na.rm=TRUE),
                lci = quantile(pdo_shifted, prob = 0.025, na.rm = TRUE),
                uci = quantile(pdo_shifted, prob = 0.975, na.rm = TRUE)) %>%
  left_join(select(pdo_results, year, min_yr, median_yr, `Time blocks`)) %>%
  group_by(`Time blocks`) %>% 
  mutate(mean_pdo = mean(pdo_shifted, na.rm = TRUE),
                pdo_sign = ifelse(mean_pdo <= 0, "-", "+")) -> pdo

ggplot(data = pdo, aes(year, PDO)) + 
  geom_errorbar(aes(year, ymin = lci, ymax = uci), 
                width = .2, color = 'grey') +
  geom_point(shape = 1) + 
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = pdo$min_yr - 0.5, col = "darkgrey") +
  geom_line(data = pdo, aes(year, mean_pdo, linetype = `Time blocks`)) +
  geom_text(data = pdo, aes(median_yr, 3, label = pdo_sign), size = 5) +
  scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels) +
  scale_y_continuous(breaks = seq(min(-3), max(3), 0.5)) +
  theme(legend.position = "none",
        axis.title.x=element_blank())+ 
  scale_linetype_manual(values = rep(2, length(unique(pdo$`Time blocks`)))) +
  labs(y = "Mean annual PDO index") +
  ggtitle("2019 forecast")

ggsave(file=paste0(pdo_path, "/pdo_", YEAR, "_forcast.png"), dpi = 300, height = 4, width = 8.5, units = "in")
