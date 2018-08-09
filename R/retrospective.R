# Retrospective analysis

# Libraries/source files ----

library(R2admb) # run ADMB in R
library(tidyverse)
library(data.table)

# source("R/tools.r")
source("R/helper.r") 

main_dir <- getwd()#"S:\\Region1Shared-DCF\\Research\\Herring-Dive Fisheries\\Herring\\ADMB Rudd Workshop 2018\\AlaskaHerring"
proj_dir <- file.path(main_dir, "HER")
setwd(proj_dir)

files <- list.files(file.path(proj_dir))

## find tpl
tpl <- file.path(proj_dir, files[which(grepl(".tpl", files))])
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1]

# find dat
dat <- file.path(proj_dir, files[which(grepl(".dat", files))])
# dat_name1 <- strsplit(files[grepl("dat",files)],".dat")[[1]][1]
# dat_name2 <- strsplit(files[grepl("dat",files)],".dat")[[2]][1]

# find ctl
ctl <- file.path(proj_dir, files[which(grepl(".ctl", files))])
# ## compile
setup_admb()
compile_admb(name, verbose = TRUE)

# check executable
exe <- paste0(name,".exe")
file.exists(exe)

# name std file (for checking model convergence in loop
std <- paste0(name,".std")

## loop over retrospective years
retro <- 0:10
rundir <- file.path(proj_dir, "run_retro")
dir.create("run_retro", showWarnings = FALSE)

check_converge <- list()

for(i in 1:length(retro)){
  
  retrodir <- file.path(rundir, paste0("retro_", retro[i]))
  dir.create(retrodir, showWarnings = FALSE)
  setwd(retrodir)
  
  # copy TPL to run directory
  file.copy(from = tpl, to = retrodir, overwrite = TRUE)	
  
  # copy exe
  file.copy(from = file.path(proj_dir,exe), to = retrodir, overwrite = TRUE)	
  
  # copy data files
  file.copy(from = dat[1], to = retrodir, overwrite = TRUE)	
  file.copy(from = dat[2], to = retrodir, overwrite = TRUE)	
  
  # copy ctl file
  file.copy(from = ctl, to = retrodir, overwrite = TRUE)	
  
  # remove pre-existing std file if it exists
  if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
  
  # run admb, set seed at retro value
  run_admb(name, extra.args=paste("-retro", retro[i]))
  
  # check for convergence (1 = converge, 0 = not converged)
  if (file.exists(std)) {
    check_converge[[i]] <- 1
  } else {check_converge[[i]] <- 0 }
  
}

# Check to see which models converges and which didn't
check_converge

# Empty vectors
SSB_ls <- list()

# Read results
for(i in 1:length(retro)){
  
  retrodir <- file.path(rundir, paste0("retro_", retro[i]))
  setwd(retrodir)
  
  D <- read_admb(name)
  # D <- read_rep(name)
  
  SSB_ls[[i]] <- data.frame(Year = D[["year"]], 
                         spB = D[["ssb"]] / 0.90718,
                         retro = paste0("retro_", retro[i])) # convert to short tons
}

# collapse list
SSB <- do.call(rbind, SSB_ls)

axisx <- tickr(SSB, Year, 5)
# Generic retrospective plot
ggplot(SSB, aes(x = Year, y = spB, colour = retro, group = retro)) +
  geom_line() +
  scale_color_grey(guide = FALSE) +
  labs(x = NULL, y = "Spawning biomass (tons)\n") +
  scale_y_continuous(label = comma) +
  scale_x_continuous(limits = c(1980, 2015), breaks = axisx$breaks, labels = axisx$labels) -> retro1

# Per Clark et al 2012, also show relative percent difference between peel and
# terminal year estimates
SSB %>% 
  group_by(retro) %>% 
  full_join(SSB %>% 
              filter(retro == "retro_0") %>% 
              select(Year, term_spB = spB),
            by = "Year") %>% 
  mutate(perc_diff = ((term_spB - spB) / term_spB) * 100) -> SSB

ggplot(data = SSB) +
  geom_hline(yintercept = 0, linetype = 2, size = 1) +
  geom_line(data = filter(SSB, retro != "retro_0"), 
            aes(x = Year, y = perc_diff, colour = retro, group = retro)) +
  scale_color_grey(guide = FALSE) +
  labs(x = NULL, y = "Percent difference\nfrom terminal year\n") +
  scale_y_continuous(limits = c(-50, 50)) +
  scale_x_continuous(limits = c(1980, 2015), breaks = axisx$breaks, labels = axisx$labels)-> retro2

cowplot::plot_grid(retro1, retro2, align = "hv", nrow = 2) -> retro_plot


setwd(main_dir)
ggsave("HER_2018forec/figs/HER/retrospective.png", 
       plot = retro_plot, 
       dpi = 300, height = 6, width = 6, units = "in")

# Alaska Fisheries Science Center and Hurtado-Ferro et al. (2015) Mohn's rho
# https://www.afsc.noaa.gov/REFM/stocks/Plan_Team/2013/Sept/Retrospectives_2013_final3.pdf
# mean over all peels (Estimate in peel year - reference estimate (current
# year's estimate) / reference estimate)

bind_cols(SSB %>% 
            filter(Year == max(Year)) %>% 
            mutate(m = (spB - term_spB) / term_spB ) %>% 
            ungroup() %>% 
            summarize(mohns_rho = mean(m)),
          
          # Wood's Hole Rho (uses all years, not just the terminal year of the peel.
          # LeGault et al 2009 found no marked differences between the two rhos)
          SSB %>% 
            mutate(m = (spB - term_spB) / term_spB ) %>% 
            ungroup() %>% 
            summarize(mrho = mean(m)) ) -> rhos

# Interpretation of Mohn's Rho (Hurtado-Ferro et al 2015):
# Given that the variability of Mohn’s rdepends on life history, and that the
# statistic appears insensitive to F, we propose the following rule of thumb
# when determining whether a retrospective pattern should be addressed
# explicitly: values of Mohn’s r higher than 0.20 or lower than 20.15 for
# longer-lived species (upper and lower bounds of the 90% simulation intervals
# for the flatfish base case), or higher than 0.30 or lower than 20.22 for
# shorter-lived species (upper and lower bounds of the 90% simulation intervals
# for the sardine base case) should be cause for concern and taken as indicators
# of retrospective patterns. However, Mohn’s r values smaller than those
# proposed should not be taken as confirmation that a given assessment does not
# present a retrospective pattern, and the choice of 90% means that a “false
# positive” will arise 10% of the time.
