
# User inputs ----
YEAR <- 2019 # forecast year
MODEL_VERSION <- "HER_best_condEffort.12_322"
source('~/AlaskaHerring/2019_forecast/r/helper.R', echo=TRUE)

# Directory setup ----
root_dir <- getwd() # project root
tpl_dir <- file.path(root_dir, paste0(YEAR, "_forecast/admb/", MODEL_VERSION)) # location of tpl
run_dir <- file.path(root_dir, paste0(YEAR, "_forecast/results/sensitivity_maturity"))
dir.create(run_dir, showWarnings = FALSE)
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find tpl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name
dat <- file.path(tpl_dir, files[which(grepl(".dat", files))]) # find dat(two)

# Compile tpl if it hasn't been compiled yet
setwd(tpl_dir)
exe <- paste0(name,".exe")
if (file.exists(exe) == FALSE) {
  setup_admb() # check admb path
  compile_admb(name, verbose = TRUE) # compile
} 

std <- paste0(name,".std") # for checking model convergence in loop

# CTL function ----
create_ctl <- function(a50 = 3.5, a95 = 5.21, phase = -2) {
  ctl <- c("

         ## —————————————————————————————————————————————————————————————————————————— ##
         ##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
         ##  Prior descriptions   Parameter values                                     ##                                                                            ##
         ##  -0 uniform           (0,0)                                                ##
         ##  -1 normal            (p1=mu,p2=sig)                                       ##
         ##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
         ##  -3 beta              (p1=alpha,p2=beta)                                   ##
         ##  -4 gamma             (p1=alpha,p2=beta)                                   ##
         ## —————————————————————————————————————————————————————————————————————————— ##
         ##  init  lower   upper    est  prior                                         ##
         ## value  bound   bound    phz   type     p1    p2   # PARAMETER              ##
         ## —————————————————————————————————————————————————————————————————————————— ##
         -1.05   -6.79    1.00      2      0  -1.05  0.05  # log_natural_mortality 
         4.60   -6.00   12.00      1      0      0     0   # log_rinit
         5.60   -6.00   12.00      1      0      0     0   # log_rbar
         6.00   -6.00   12.00      2      0      0     0   # log_ro
         2.40    0.00   12.00      2      2   2.73  0.50   # log_reck
         4.00    0.00  200.00     -2      4   1.05  1.05   # precision = 1/(sigma_r^2), SM=6.25
         ## 50.25    0.00  200.00    -2      4   1.05  1.05   # precision = 1/(sigma_r^2) (for simulations with no process error)
         
         ## —————————————————————————————————————————————————————————————————————————— ##
         
         ## —————————————————————————————————————————————————————————————————————————— ##
         ##                CONTROLS FOR TIME-VARYING MATURITY                          ##
         ## —————————————————————————————————————————————————————————————————————————— ##
         ## nMatBlocks
         1
         ## a50    a95     phz   terminalBlockYear
         ## 3.5       4.5     -2       1998 ",
         a50,       a95,    phase,  " 2018
         ## —————————————————————————————————————————————————————————————————————————— ##
         
         ## —————————————————————————————————————————————————————————————————————————— ##
         ##            CONTROLS FOR TIME-VARYING LN(NATURAL MORTALITY DEVS)            ##
         ## KEY:
         ##  Type: 1 = constant M
         ##  Type: 2 = interpolated using cubic spline.
         ## —————————————————————————————————————————————————————————————————————————— ##
         ## Type
         1
         ## Phase for estimation if nMortBlocks > 1
         2
         ## nMortBlocks, or Nodes in the case of cubic spline interpolation
         3
         ## The terminal year of each block
         1998 2014 2018
         ## ————————————————————————————————————————————————————————————————————————— ##
         
         ## —————————————————————————————————————————————————————————————————————————— ##
         ##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
         ## —————————————————————————————————————————————————————————————————————————— ##
         ## - Each selectivity block can have different functional forms based on selType
         ## - LEGEND:
         ##   - SelType = 1, logistic selectivity, 2 parameters.
         ##   - SelType = 2, logistic with 50% & 95% parameters
         ##  nSelexblocks
         1
         ## —————————————————————————————————————————————————————————————————————————— ##
         ##  Gear  Sel     sel   sel   age   year  phz    | start end
         ##  Index Type    mu    sd    nodes nodes mirror | block block
         ## —————————————————————————————————————————————————————————————————————————— ##
         ## 1     1     4.0     0.3     0     0     2          1980     1998
         1     1    3.3     0.8    0     0     2          1980     2018
         ## ——————————————————————————————————————————————————————————————————————————— ##
         
         ## ——————————————————————————————————————————————————————————————————————————— ##
         ##                        OTHER MISCELLANEOUS CONTROLS                         ##
         ## ——————————————————————————————————————————————————————————————————————————— ##
         ## number of controls to read in.
         8
         ## Value    # # - Description
         0.90718     # 1 - Catch Scaler (convert from short tons to metric tons)
         1           # 2 - Condition on Catch = 0, Condition of Ft = 1
         25000       # 3 - harvest threshold
         0.2         # 4 - target harvest rate
         20000       # 5 - threshold denominator
         0.090       # 6 - standard deviation in natural mortality devs SM=0.001
         1.00        # 7 - sd in recruitment deviations in all phases of estimate up until the last
         5.00        # 8 - sd in recruitment deviation in the final phase of estimation
         ## EOF
         999
         
         ")
# Write it to file
write.table(ctl, file = "sitka.ctl", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

}

# Starting values for a50 ----

# BC maturity parameters 
BC_a50 <- 2.26; BC_a95 <- 3.21
# assume this ratio is constant for sensitivity analysis
mat_ratio <- BC_a95 / BC_a50 

# Get starting values from Fish Life 
# (queries FishBase then uses a hierarchical
# model to obtain variance estimates) - 
# Source: https://www.ncbi.nlm.nih.gov/pubmed/28746981
# devtools::install_github("james-thorson/FishLife")
# Use this for the reference model
library(FishLife)
BRF_ex <- Plot_taxa(
  Search_species(Genus="Clupea",
                 Species="pallasii pallasii")$match_taxonomy, 
  mfrow=c(2,2) )

a50 <- exp(BRF_ex[[1]]$Mean_pred[5]+diag(BRF_ex[[1]]$Cov_pred/2)[5])
a95 <- a50 * mat_ratio

# Option 1 ----
# # Curves to explore for a sensitivity analysis
# a50_vec <- seq(3.0, 4.0, 0.1)
# 
# for(i in 1:length(a50_vec)) {
#   for (j in 1:length(ratio)) {
#     
#     a50 <- a50_vec[i]
#     rat <- ratio[j]
#     k <- 1/ (a50 / rat)
#     
#     col <- ifelse(j == 1, "grey", "grey")
#     
#     add <- if(i == 1 & j == 1) FALSE else TRUE
#     
#     # curve(1 / (1 + exp(-1.0*k * (x - a50))),
#     curve(1 / (1 + exp(-1.0* (x - a50) / k)),
#           ylim = c(0,1), from = 3, to = 8, add = add, col = col)
#   }
# }
# 
# a50 <- exp(BRF_ex[[1]]$Mean_pred[5]+diag(BRF_ex[[1]]$Cov_pred/2)[5])
# k <- 1 / (a50 / mean(ratio) )
# 
# # Example of a central tendancy curve
# curve(1 / (1 + exp(-1.0* (x - a50) / k)),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "blue")
# 
# # LS maturity estimates
# curve(1 / (1 + exp(-1.0*2.34 * (x - 3.42))),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red", lty = 2)
# curve(1 / (1 + exp(-1.0*4.99 * (x - 3.12))),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red")


# Option 2 ----
# Curves to explore for a sensitivity analysis
a50 <- exp(BRF_ex[[1]]$Mean_pred[5]+diag(BRF_ex[[1]]$Cov_pred/2)[5])
a50_vec <- seq(3.0, 4.0, 0.25)
a95 <- 4.48 # 100% mature from SD's regression 
rat <- a95/a50 # Regression a95 / FishLife a50
rat_vec <- c(1.1, 1.3, 1.5)

par(mfrow = c(1,1), cex = 1.3)
for(i in 1:length(a50_vec)) {
  for (j in 1:length(rat_vec)) {
    
    a50 <- a50_vec[i]
    r <- rat_vec[j]
    a95 <- a50 * r
    
    col <- ifelse(j == 1, "green", ifelse(j == 2, "blue", "red"))
    lty <- ifelse(j == 1, 2, ifelse(j == 2, 1, 3))
    
    add <- if(i == 1 & j == 1) FALSE else TRUE
    
    curve(1.0 / (1.0 + exp(-log(19)*(x-a50)/(a95-a50))), 
          ylim = c(0,1), from = 3, to = 8, 
          ylab = "Proportion mature", xlab = "\nAge",
          add = add, col = col, lty = lty, lwd = 2)
  }
}

# LS maturity estimates
 curve(1 / (1 + exp(-1.0* (x - 3.3)/.8)),
       ylim = c(0,1), from = 3, to = 8, add = F, col = "red", lty = 2)
 
 
curve(1 / (1 + exp(-1.0*2.34 * (x - 3.42))),
      ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red", lty = 2)
curve(1 / (1 + exp(-1.0*4.99 * (x - 3.12))),
      ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red")

a50 <- exp(BRF_ex[[1]]$Mean_pred[5]+diag(BRF_ex[[1]]$Cov_pred/2)[5])
a95 <- 4.48 # from Sherri's regression

# Example of a central tendancy curve
curve(1.0 / (1.0 + exp(-log(19)*(x-a50)/(a95-a50))), 
      ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "blue")

# Reference model ----
setwd(run_dir)
ref_dir <- file.path(run_dir, paste0("HER_a50_reference"))
dir.create(ref_dir, showWarnings = FALSE)
setwd(ref_dir)

# Create new control file
create_ctl(a50 = a50, a95 = a95, phase = 2)

# copy TPL to run directory
file.copy(from = tpl, to = ref_dir, overwrite = TRUE)	

# copy exe
file.copy(from = file.path(tpl_dir, exe), to = ref_dir, overwrite = TRUE)	

# copy data files
file.copy(from = dat[1], to = ref_dir, overwrite = TRUE)	
file.copy(from = dat[2], to = ref_dir, overwrite = TRUE)	

# remove pre-existing std file if it exists
if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}

compile_admb(name, verbose = TRUE)
run_admb(name, verbose = TRUE)
R <- read_admb(name)
Rfit <- read_fit(name)
data.frame(names = Rfit[["names"]],
           est = Rfit[["est"]]) %>% 
  filter(grepl("mat", names)) %>% 
  mutate(param = c("a50", "a95"),
         time_block = c(1,1)) -> matpar
# R[["mat"]]

a95 <- matpar %>% filter(param == "a95") %>% pull(est)
a50 <- matpar %>% filter(param == "a50") %>% pull(est)

# What does our estimated curve look like?
curve(1.0 / (1.0 + exp(-log(19)*(x-a50)/(a95-a50))), 
      ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "green")

# curve(1 / (1 + exp(-1.0*k * (x - a50))),
#       ylim = c(0,1), from = 3, to = 8, add = FALSE, col = "grey")
# curve(1 / (1 + exp(-1.0*2.34 * (x - 3.42))),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red", lty = 2)
# curve(1 / (1 + exp(-1.0*4.99 * (x - 3.12))),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "orange")

# ratio_fast <- 3.42 / 2.34 # LS time period 1
# ratio_slow <- 3.12 / 4.99 # LS time period 2
# ratio <- c(ratio_fast, ratio_slow)
# 
# a50_vec <- seq(3.0, 4.0, 0.1)


getwd()
Rfit[["est"]]
names(Rfit)
ref_a50 <- R[["mat_a95"]] / R[["mat_a50"]]

# Run analysis ----

# Examine wide range of maturity estimates
# a50_vec <- seq(BC_a50, 8, 0.05) 

# check_converge <- list()
check_converge <- matrix(nrow = length(a50_vec), ncol = length(rat_vec))

for(i in 1:length(a50_vec)) {
  # for (j in 1:length(ratio)) {
  for (j in 1:length(rat_vec)) {

    # rat <- ifelse(j == 1, "slowk", "fastk")
    # sensdir <- file.path(run_dir, paste0("HER_a50_", a50_vec[i], "_", rat))
    sensdir <- file.path(run_dir, paste0("HER_a50_", a50_vec[i], 
                                         "_a95_", rat_vec[j] * a50_vec[i],
                                         "_ratio_", rat_vec[j]))
    dir.create(sensdir, showWarnings = FALSE)
    setwd(sensdir)
    # check_converge[i,1] <- paste0("HER_a50_", a50_vec[i], 
    #                               "_a95_", rat_vec[j] * a50_vec[i],
    #                               "_ratio_", rat_vec[j])
    
    # Create new control file
    # create_ctl(a50 = a50_vec[i], a95 = 1 / (a50_vec[i]/ratio[j]), 
    #            phase = -2)
    create_ctl(a50 = a50_vec[i], a95 = (a50_vec[i] * rat_vec[j]),
               phase = -2)

    # copy TPL to run directory
    file.copy(from = tpl, to = sensdir, overwrite = TRUE)	
    
    # copy exe
    file.copy(from = file.path(tpl_dir, exe), to = sensdir, overwrite = TRUE)	
    
    # copy data files
    file.copy(from = dat[1], to = sensdir, overwrite = TRUE)	
    file.copy(from = dat[2], to = sensdir, overwrite = TRUE)	
    
    # remove pre-existing std file if it exists
    if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
    
    # run admb, set seed at retro value
    run_admb(name, verbose = TRUE)
    
    # check for convergence (1 = converge, 0 = not converged)
    if (file.exists(std)) {
      check_converge[i,j] <- 1
    } else {check_converge[i,j] <- 0 }
    
    # Also check max gradient component
    # mgc <- readMGC(name) 
    # if (mgc < 1e-3) {
    #   check_converge[i,2] <- 1
    # } else {check_converge[i,2] <- 0 } 
    
  }
}

# Summarize ----
check <- as.data.frame(check_converge, row.names = FALSE)
names(check) <- rat_vec
check %>% 
  mutate(a50 = a50_vec) %>% 
  gather("ratio", "converged", -a50) %>% 
  mutate(ratio = as.numeric(ratio),
         fn = paste0("HER_a50_", a50, 
                     "_a95_", ratio * a50,
                     "_ratio_", ratio)) -> check

filter(check, converged == 0) # all models converged

setwd(run_dir)
if(exists("check")) {
  write_csv(check, "convergence_results.csv")
} else {
  check <- read_csv("convergence_results.csv") }

check <- filter(check, converged == 1)

diags_ls <- list() # store diagnostic output
derived_ls <- list() # store derived time series of interest
natmat_ls <- list() # time-varying natural mortality & survival
matsel_ls <- list() # time-varying maturity and selectivity

for(i in 1:length(check$fn)){
  
  fn <- check$fn[i]
  
  sensdir <- file.path(run_dir, fn)
  setwd(sensdir)
  
  # if (file.exists(std)) {
    D <- read_admb(name) #}
  
  # Diagnostics, key derived forecast quantities and parameter estimates
  diags_ls[[i]] <- data.frame(fn = fn,
                              a50 = D[["mat_a50"]],
                              a95 = D[["mat_a95"]],
                              # a50_1 = D[["mat_a50"]][1],
                              # a50_2 = D[["mat_a50"]][2],
                              # a95_1 = D[["mat_a95"]][1],
                              # a95_2 = D[["mat_a95"]][2],
                              nll = D[["fit"]]$nlogl,
                              maxgrad = D[["fit"]]$maxgrad,
                              nopar = D[["fit"]]$nopar,
                              ghl = D[["ghl"]], # Already in short tons
                              fore_sb = D[["fore_sb"]] / 0.90718, # convert to short tons
                              fore_matb = D[["fore_matb"]] / 0.90718, # convert to short tons
                              Mbar = exp(D[["theta"]][1]),
                              rinit = exp(D[["theta"]][2]),
                              rbar = exp(D[["theta"]][3]),
                              r0 = exp(D[["theta"]][4]),
                              reck = exp(D[["theta"]][5]) + 1) %>% 
 mutate(mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                              ifelse(a95/a50 >= 1.5, "Slow", "Medium")))
  
  # # Natural mortality and survival
  df <- data.frame(D[["year"]],
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))

  natmat_ls[[i]] <- df %>%
    select(Year, M = `3`) %>% # nat mat not age-specific, so just use one age
    mutate(survival = exp(-M),
           a50 = D[["mat_a50"]],
           a95 = D[["mat_a95"]],
           fn = fn) %>%
    filter(Year %in% c(1998, 2014, 2018)) %>% # Predefined time blocks *FLAG* this is hard coded
    mutate(`Time blocks` = c("1980-1998", "1999-2014", "2015-2018"),
           `Time blocks` = factor(`Time blocks`, ordered = TRUE),
           mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                              ifelse(a95/a50 >= 1.5, "Slow", "Medium")))

  # Maturity/selectivity
  df <- data.frame(D[["year"]],
                   D[["mat"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  df %>% mutate(param = "Maturity") -> df

  nyr <- D[["mod_nyr"]] - D[["mod_syr"]] + 1

  df2 <- data.frame(D[["year"]],
                    D[["Sij"]][1:nyr, ]) # change number of rows
  colnames(df2) <- c("Year", paste(D[['sage']]:D[['nage']]))
  df2 %>% mutate(param = "Selectivity") %>%
    bind_rows(df) -> df

  matsel_ls[[i]] <- df %>%
    gather("Age", "proportion", -c(Year, param)) %>%
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) %>%
    group_by(Age, param, proportion) %>%
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max),
           a50 = D[["mat_a50"]],
           a95 = D[["mat_a95"]],
           fn = fn,
           mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                              ifelse(a95/a50 >= 1.5, "Slow", "Medium"))) %>%
    ungroup() %>%
    distinct(a50, a95, fn, param, `Time blocks`, Year = max, Age, proportion, mat_ratio)

  # Derived time series of biomass and ASA vs Ricker-predicted recruitment
  derived_ls[[i]] <- data.frame(Year = D[["year"]],
                                spB = D[["sp_B"]] / 0.90718, # convert to short tons
                                matB = D[["mat_B"]] / 0.90718,
                                age3 = D[["Nij"]][1:nyr, 1] # ASA age-3 recruits
  ) %>%
    left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                    D[["mod_nyr"]] + 1, 1),
                         # SR is the current year age-3
                         SR = D[["recruits"]],
                         # ssb is the spawning biomass in year-3 that produced
                         # the current years recruits
                         ssb = D[["spawners"]] / 0.90718), by = "Year") %>%
    mutate(a50 = D[["mat_a50"]],
           a95 = D[["mat_a95"]],
           fn = fn,
           mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                              ifelse(a95/a50 >= 1.5, "Slow", "Medium")))
}

diags <- do.call(rbind, diags_ls)
derived <- do.call(rbind, derived_ls)
natmat <- do.call(rbind, natmat_ls)
matsel <- do.call(rbind, matsel_ls)

# Convergence problems? No.
diags %>% filter(maxgrad > 0.001) %>% pull(fn) -> high_maxgrad
diags %>% filter(!fn %in% high_maxgrad) -> diags
derived %>% filter(!fn %in% high_maxgrad) -> derived
natmat %>% filter(!fn %in% high_maxgrad) -> natmat
matsel %>% filter(!fn %in% high_maxgrad) -> matsel

# Reference model results ----
setwd(tpl_dir)
compile_admb(name, verbose = TRUE)
run_admb(name, verbose = TRUE)
D <- read_admb(name)
ref_diags <- data.frame(a50_1 = D[["mat_a50"]][1],
                        a50_2 = D[["mat_a50"]][2],
                        a95_1 = D[["mat_a95"]][1],
                        a95_2 = D[["mat_a95"]][2],
                        nll = D[["fit"]]$nlogl,
                        maxgrad = D[["fit"]]$maxgrad,
                        nopar = D[["fit"]]$nopar,
                        ghl = D[["ghl"]], # Already in short tons
                        fore_sb = D[["fore_sb"]] / 0.90718, # convert to short tons
                        fore_matb = D[["fore_matb"]] / 0.90718, # convert to short tons
                        Mbar = exp(D[["theta"]][1]),
                        rinit = exp(D[["theta"]][2]),
                        rbar = exp(D[["theta"]][3]),
                        r0 = exp(D[["theta"]][4]),
                        reck = exp(D[["theta"]][5]) + 1)

# # Natural mortality and survival
df <- data.frame(D[["year"]],
                 D[["Mij"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))

ref_natmat <- df %>%
  select(Year, M = `3`) %>% # nat mat not age-specific, so just use one age
  mutate(survival = exp(-M),
         a50_1 = D[["mat_a50"]][1],
         a50_2 = D[["mat_a50"]][2],
         a95_1 = D[["mat_a95"]][1],
         a95_2 = D[["mat_a95"]][2],
         fn = fn) %>%
  filter(Year %in% c(1998, 2014, 2018)) %>% # Predefined time blocks *FLAG* this is hard coded
  mutate(`Time blocks` = c("1980-1998", "1999-2014", "2015-2018"),
         `Time blocks` = factor(`Time blocks`, ordered = TRUE))

# Maturity/selectivity
df <- data.frame(D[["year"]],
                 D[["mat"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
df %>% mutate(param = "Maturity") -> df

nyr <- D[["mod_nyr"]] - D[["mod_syr"]] + 1

df2 <- data.frame(D[["year"]],
                  D[["Sij"]][1:nyr, ]) # change number of rows
colnames(df2) <- c("Year", paste(D[['sage']]:D[['nage']]))
df2 %>% mutate(param = "Selectivity") %>%
  bind_rows(df) -> df

ref_matsel <- df %>%
  gather("Age", "proportion", -c(Year, param)) %>%
  mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>%
  group_by(Age, param, proportion) %>%
  mutate(min = min(Year),
         max = max(Year),
         `Time blocks` = paste0(min, "-", max)) %>%
  ungroup() %>%
  distinct(param, `Time blocks`, Year = max, Age, proportion)

# Derived time series of biomass and ASA vs Ricker-predicted recruitment
ref_derived <- data.frame(Year = D[["year"]],
                          spB = D[["sp_B"]] / 0.90718, # convert to short tons
                          matB = D[["mat_B"]] / 0.90718,
                          age3 = D[["Nij"]][1:nyr, 1] # ASA age-3 recruits
) %>%
  left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                  D[["mod_nyr"]] + 1, 1),
                       # SR is the current year age-3
                       SR = D[["recruits"]],
                       # ssb is the spawning biomass in year-3 that produced
                       # the current years recruits
                       ssb = D[["spawners"]] / 0.90718), by = "Year")

# plot(x = diags$a50, y = diags$nll, type = "l")
# plot(x = diags$a50, y = diags$nll, type = "l")

figdir <- file.path(run_dir, "figures")
dir.create(figdir)

# For showing conditioned of catch vs effort
palR <- colorRampPalette(c("#bfd3e6","#8c6bb1", "#4d004b"))
ratio_cols <- palR(length(rat_vec))

# Plotting for presentation ----
library(ggthemr)
my_cols <- c('#555555', '#008080', '#DAA520', '#330000') # first col is placeholder
ggthemr_reset() # remove previous effects
# Define colours for your figures with define_palette
tableau <- define_palette(
  swatch = my_cols, # colours for plotting points and bars
  gradient = c(lower = my_cols[1L], upper = my_cols[2L]), #upper and lower colours for continuous colours
  background =  "#f8f8f8" # "#F6F9FB"  #defining a grey-ish background 
)
# set the theme for your figures:
ggthemr(tableau, layout = "clean", type = "outer", text_size = 18)

pal <- colorRampPalette(c("#edf8b1","#7fcdbb", "#2c7fb8"))
a50_cols <- pal(length(a50_vec))

ref_col <- "#ff6a5a"

# GHL ----

# Large sigmaM's result in drastically smaller GHLs overall. The minimum,
# maximum, and approximate median values are shown in red triangle, green
# upsidedown triangle, and orange square.

diags %>% 
  ggplot(aes(x = a50, y = mat_ratio, size = ghl)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "GHL (t)",
             breaks = round(seq(min(diags$ghl), max(diags$ghl), length.out = 6), 0),
             labels = comma) + 
  geom_point(data = filter(diags, ghl == max(ghl)),
             col = "red", size = 10, fill = NA, shape = 2) +
  geom_point(data = filter(diags, ghl == min(ghl)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = slice(diags, which.min(abs(ghl-median(ghl)))),
             col = "orange", size = 10, fill = NA, shape = 0) +
  labs(x = "\na50", y = "Maturity rate\n")

ggsave(paste0(figdir, "/ghl.png"), dpi = 300, height = 4, width = 5, units = "in")

# Natural mortality ----

# increasing sigmaM has a relatively small effect on time-varying M over this
# range
ref_natmat <- ref_natmat %>% 
  left_join(expand.grid(mat_ratio = c("Slow", "Medium", "Fast"),
                       `Time blocks` = c("1980-1998", "1999-2014", "2015-2018")))

ggplot(natmat, aes(x = as.factor(as.numeric(`Time blocks`)), y = survival, 
                   group = a50, colour = factor(a50))) +
  geom_line(size = 1.5) + 
  geom_point(size = 2) + 
  geom_line(data = ref_natmat, aes(x = as.factor(as.numeric(`Time blocks`)), y = survival, 
                                   group = 1), col = ref_col, size = 1.5, lty = 2) +
  geom_point(data = ref_natmat, aes(x = as.factor(as.numeric(`Time blocks`)), y = survival, 
                                    group = 1), col = ref_col, size = 2) +
  scale_colour_manual(values = a50_cols) +
  facet_wrap(~ factor(mat_ratio)) +
  labs(colour = "a50", y = NULL, x = "Time blocks") +
  ylim(c(0,1)) +
  ggtitle("Estimated survival by maturation rate")

ggsave(paste0(figdir, "/survival_a50.png"), dpi = 300, height = 5, width = 10, units = "in")

# Maturity and selectivity ----

ref_matsel <- ref_matsel %>% 
  left_join(expand.grid(param = c("Maturity", "Selectivity"),
                        mat_ratio = c("Slow", "Medium", "Fast")),
            by = c("param"))

# Maturity
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion, colour = factor(a50), group = factor(a50))) + 
  geom_line(aes(colour = factor(a50)), size = 2) +
  geom_line(data = ref_matsel %>% filter(param == "Maturity"),
            aes(x = Age, y = proportion, linetype = `Time blocks`, group = `Time blocks`), 
            col = ref_col, size = 1.3) +
  scale_linetype_manual(values = c(2,1,3), guide = FALSE) +
  scale_colour_manual(values = a50_cols) +
  facet_wrap(~mat_ratio) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = NULL, colour = "a50") +
  ggtitle("Fixed maturity")

ggsave(paste0(figdir, "/maturity_a50.png"), dpi = 300, height = 5, width = 10, units = "in")

# selectivity 
matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(colour = factor(a50), 
                group = factor(a50)), size = 2) +
  geom_line(data = ref_matsel %>% filter(param == "Selectivity"),
            aes(x = Age, y = proportion, linetype = `Time blocks`, group = `Time blocks`), 
            col = ref_col, size = 1.3) +
  scale_colour_manual(values = a50_cols) +
  scale_linetype_manual(values = c(2,1,3), guide = FALSE) +
  facet_wrap(~mat_ratio) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = NULL, colour = "a50") +
  ggtitle("Estimated selectivity by maturation rate")

ggsave(paste0(figdir, "/selectivity_a50.png"), dpi = 300, height = 5, width = 10, units = "in")

# Derived time series ----

# Separate forecast quantities into separate df - cannot use the model-supplied
# fore_sb b/c the assumptions are not correct
diags %>% 
  select(fn, a50, a95, fore_matb, ghl, mat_ratio) %>% 
  mutate(Year = max(derived$Year) + 1,
         fore_sb = fore_matb - ghl) -> forec

# SSB facetted by sigmaM - shows variability in ssb is attributed to sigmaM not
# condition

# Facet by condition (colours are sigmaM)
M <- derived %>% 
  ggplot(aes(x = Year, colour = factor(a50), group = factor(a50))) +
  scale_colour_manual(values = a50_cols) +
  facet_wrap(~mat_ratio, ncol = 1) +
  scale_y_continuous(label = scales::comma)

# Facet by sigmaM (colours are condition)
R <- derived %>% 
  ggplot(aes(x = Year, colour = mat_ratio, linetype = mat_ratio, 
             group = mat_ratio)) +
  scale_colour_manual(values = ratio_cols) +
  facet_wrap(~factor(a50), ncol = 1) +
  scale_y_continuous(label = scales::comma)

# Spawning biomass -----
ref_diags %>% 
  select(fore_matb, ghl) %>% 
  mutate(Year = max(derived$Year) + 1,
         fore_sb = fore_matb - ghl) -> ref_forec

derived %>% 
  mutate(combo = paste0(a50, "_", mat_ratio)) %>% 
  ggplot(aes(x = Year)) +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, 130000), breaks = seq(0, 130000, 25000), label = scales::comma) +
  geom_line(aes(y = spB, colour = factor(a50), group = combo, linetype = mat_ratio), size = 0.8) +
  geom_line(data = ref_derived, aes(x = Year, y = spB), col = ref_col, size = 1.3) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Spawning biomass estimates and forecast (t)") +
  geom_point(data = forec, 
             aes(x = Year, y = fore_sb, 
                 colour = factor(a50), group = factor(a50)),
             shape = "*", size = 4) +
  geom_point(data = ref_forec, 
             aes(x = Year, y = fore_sb, 
                 colour = factor(a50), group = factor(a50)),
             shape = "*", size = 4, col = ref_col) 

ggsave(paste0(figdir, "/spbiomass.png"), dpi = 300, height = 5, width = 10, units = "in")

ggplot(data = ref_derived, aes(x = Year, y = spB)) +
  scale_colour_manual(values = ref_col) +
  geom_line(aes(colour = "Reference\nmodel", group = 1), size = 1.3) +
  labs(x = NULL, y = NULL, colour = NULL) +
  ggtitle("Spawning biomass estimates and forecast (t)") 

ggsave(paste0(figdir, "/reference_model_legend.png"), dpi = 300, height = 5, width = 10, units = "in")

# Mature biomass -----

derived %>% 
  mutate(combo = paste0(a50, "_", mat_ratio)) %>% 
  ggplot(aes(x = Year)) +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, 150000), breaks = seq(0, 150000, 25000), label = scales::comma) +
  geom_line(aes(y = matB, group = combo, colour = factor(a50), linetype = mat_ratio), size = 0.8) +
  geom_line(data = ref_derived, aes(x = Year, y = matB), col = ref_col, size = 1.3) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Mature biomass estimates and forecast (t)") +
  geom_point(data = forec, 
             aes(x = Year, y = fore_matb, 
                 colour = factor(a50), group = factor(a50)),
             shape = "*", size = 4) +
  geom_point(data = ref_forec, 
             aes(x = Year, y = fore_matb, 
                 colour = factor(a50), group = factor(a50)),
             shape = "*", size = 4, col = ref_col) 

ggsave(paste0(figdir, "/matbiomass.png"), dpi = 300, height = 5, width = 10, units = "in")

# Recruitment ----

ref_derived <- ref_derived %>% 
  left_join(expand.grid(Year = unique(derived$Year),
                        mat_ratio = unique(derived$mat_ratio)))
# ASA-estimated age-3 abundance: 

# sigmaM appears to have a slight affect on estimates but trends are the same.
M + geom_line(aes(y = age3), size = 0.8) +
  geom_line(data = ref_derived, aes(y = age3), col = ref_col, size = 1, lty = 2) +
  labs(x = NULL, y = NULL, colour = "a50") +
  ggtitle("ASA age-3 abundance (millions) by maturation rate")

ggsave(paste0(figdir, "/ASA_age3_bymatrate.png"), dpi = 300, height = 8, width = 10, units = "in")

# condition has no affect on estimates
R + geom_line(aes(y = age3), size = 0.8) +
  # geom_line(data = ref_derived, aes(y = age3), col = ref_col, size = 1, lty = 2) +
  
  labs(x = NULL, y = NULL, colour = "Maturation\nrate", linetype = "Maturation\nrate")+
  ggtitle("ASA age-3 abundance (millions) by a50")

ggsave(paste0(figdir, "/ASA_age3_bya50.png"), dpi = 300, height = 8, width = 10, units = "in")

# Ricker-estimated age-3 abundance:

derived %>% 
  mutate(combo = paste0(a50, "_", mat_ratio)) %>% 
  ggplot(aes(x = Year)) +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(label = scales::comma) +
  geom_line(aes(y = SR, colour = factor(a50), group = combo, linetype = mat_ratio), size = 1.2) +
  geom_line(data = ref_derived, aes(y = SR), col = ref_col, size = 1.5) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Ricker age-3 abundance (millions)")

ggsave(paste0(figdir, "/ricker_age3.png"), dpi = 300, height = 5, width = 10, units = "in")

# BRP ref model ----

setwd(tpl_dir)
D <- read_admb(name)
names(D)

mod_names <- c("Ricker", "Beverton-Holt")

sage <- D[["sage"]]; nage <- D[["nage"]]
syr <- D[["mod_syr"]]; nyr <- D[["mod_nyr"]]
age <- sage:nage
fmort <- seq(0, 10, 0.01)
rec_mod <- 1 # 1 = Ricker, 2 = Beverton-Holt

# Spawner weight-at-age
wa <-  as.data.frame(D[["data_sp_waa"]])
colnames(wa) <- c("Year", paste(sage:nage))
wa <-	wa %>% 
  melt(id.vars = c("Year"), variable.name = "Age", value.name = "weight") %>% 
  filter(Year >= syr) %>% 
  group_by(Age) %>% 
  dplyr::summarise(weight = mean(weight)) %>% 
  pull(weight)

mat <- read_D("mat")
sel <- read_D("Sij")
mort <- read_D("Mij")
r0 <- exp(D[["theta"]][4])
reck <- exp(D[["theta"]][5]) + 1 # FLAG - in her.tpl, he adds 1. Assuming we should do the same.

h <- get_steepness(rec_mod = rec_mod)

# Unfished survival, spawning biomass per recruit, and unfished spawning biomass
lx0 <- get_survivorship(fmort = 0)
spr0 <- get_spr(fmort = 0, lx = lx0)
bpr0 <- get_bpr(fmort = 0, lx = lx0)
s0 <- r0 * spr0

# Get BRPs
out <- get_equilibrium(ff = fmort)

# Recruitment parameters
alpha <- get_recruit_pars()[[1]]
beta <- get_recruit_pars()[[2]]

# uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,100), tol = 0.000001, 
                delta = 0.0001)$root[1]
msy <- get_equilibrium(ff = Fmsy)$y
# Fcrash <- bisect()

# Get components of B0
lx_Fmsy <- get_survivorship(fmort = Fmsy)
ypr_Fmsy <- get_ypr(fmort = Fmsy, lx = lx_Fmsy) # Fmsy * phi_q evaluated at Fmsy
spr_Fmsy <- get_spr(fmort = Fmsy, lx = lx_Fmsy) # phi_e evaluated at Fmsy
r_Fmsy <- get_recruitment(spr = spr_Fmsy)

bmsy <- r_Fmsy * spr_Fmsy  # Equilibrium spawning biomass at Fmsy
b0 <- r0 * spr0 # Unfished biomass

ref_brps <- data.frame(Model = mod_names[rec_mod], h = round(h, 4), 
                        alpha = round(alpha, 3), beta = round(beta, 5), r0 = round(r0, 0), s0 = round(s0, 0), 
                        spr0 = round(spr0, 0), Fmsy = round(Fmsy, 4), MSY = round(msy/0.90718, 0), 
                        Bmsy = round(bmsy/0.90718, 0), B0 = round(b0/0.90718, 0))

ref_out <- get_equilibrium(fmort)
ref_out <- ref_out %>% 
  map(~data.frame(.x)) %>% 
  imap(~mutate(.x, variable = .y)) %>%
  map(~mutate(.x, fn = "ref_mod")) 
ref_out <- bind_rows(ref_out)
colnames(ref_out) <- c("value", "variable", "fn")
ref_out <- data.table(ref_out)
ref_out[, fmort := rep(fmort, length(unique(ref_out$variable)))]
ref_out %>% filter(value >= 0) -> ref_out
ref_out <- ref_out %>% 
  left_join(expand.grid(fn = "ref_mod",
                        mat_ratio = unique(derived$mat_ratio)))
ref_brps <- ref_brps %>% 
  left_join(expand.grid(Model = "Ricker",
                        mat_ratio = unique(derived$mat_ratio)))
# BRPS ----

output <- list()
brps <- list()

for(i in 1:length(check$fn)){
  
  fn <- check$fn[i]
  
  sensdir <- file.path(run_dir, fn)
  setwd(sensdir)
  
  D <- read_admb(name) 

  mat <- read_D("mat")
  sel <- read_D("Sij")
  mort <- read_D("Mij")
  r0 <- exp(D[["theta"]][4])
  reck <- exp(D[["theta"]][5]) #+ 1 # FLAG - in her.tpl, he adds 1. Assuming we should do the same.
  h <- get_steepness(rec_mod = rec_mod)
  
  # Unfished survival, spawning biomass per recruit, and unfished spawning biomass
  lx0 <- get_survivorship(fmort = 0)
  spr0 <- get_spr(fmort = 0, lx = lx0)
  bpr0 <- get_bpr(fmort = 0, lx = lx0)
  s0 <- r0 * spr0
  
  # Recruitment parameters
  alpha <- get_recruit_pars()[[1]]
  beta <- get_recruit_pars()[[2]]
  
  # Get BRPs
  out <- get_equilibrium(ff = fmort)
  
  # uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
  Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,100), tol = 0.000001, 
                  delta = 0.0001)$root[1]
  msy <- get_equilibrium(ff = Fmsy)$y
  # Fcrash <- bisect()
  
  # Get components of B0
  lx_Fmsy <- get_survivorship(fmort = Fmsy)
  ypr_Fmsy <- get_ypr(fmort = Fmsy, lx = lx_Fmsy) # Fmsy * phi_q evaluated at Fmsy
  spr_Fmsy <- get_spr(fmort = Fmsy, lx = lx_Fmsy) # phi_e evaluated at Fmsy
  r_Fmsy <- get_recruitment(spr = spr_Fmsy)
  
  bmsy <- r_Fmsy * spr_Fmsy  # Equilibrium spawning biomass at Fmsy
  b0 <- r0 * spr0
  
  brps[[i]] <- data.frame(Model = mod_names[rec_mod], iter = i, fn = fn,
                          a50 = D[["mat_a50"]], a95 = D[["mat_a95"]], h = round(h, 4), 
                          alpha = round(alpha, 3), beta = round(beta, 5), r0 = round(r0, 0), s0 = round(s0, 0), 
                          spr0 = round(spr0, 0), Fmsy = round(Fmsy, 4), MSY = round(msy/0.90718, 0), 
                          #Fcrash = round(Fcrash, 4), 
                          Bmsy = round(bmsy/0.90718, 0), B0 = round(b0/0.90718, 0)) %>% 
    mutate(mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                              ifelse(a95/a50 >= 1.5, "Slow", "Medium")))
  
  out <- get_equilibrium(fmort)
  
  out <- out %>% 
    map(~data.frame(.x)) %>% 
    imap(~mutate(.x, variable = .y)) %>%
    map(~mutate(.x, iter = i)) %>% 
    map(~mutate(.x, fn = fn)) %>% 
    map(~mutate(.x, a50 = D[["mat_a50"]])) %>% 
    map(~mutate(.x, a95 = D[["mat_a95"]])) %>% 
    map(~mutate(.x, mat_ratio = ifelse(a95/a50 <= 1.1, "Fast",
                                       ifelse(a95/a50 >= 1.5, "Slow", "Medium")))) 

  out <- bind_rows(out)
  colnames(out) <- c("value", "variable", "iter", "fn", "a50", "a95", "mat_ratio")
  output[[i]] <- out
}

brps <- bind_rows(brps)
output <- bind_rows(output)
brps <- melt(brps, id.vars = c("Model", "iter", "fn", "a50", "mat_ratio"))
brps <- data.table(brps)
brps %>% group_by(variable) %>% summarize(mean(value), max(value), min(value))
output <- data.table(output)
output[, fmort := rep(fmort, length(unique(iter))*length(unique(variable)))]
output %>% filter(value >= 0) -> ouput

# Plot BRPs ----

# MSY
msy_p <- ggplot() +
  geom_line(data = output %>% 
              filter(variable == "y" & value >= 0),
            aes(y = value/0.90718, x = fmort, 
                colour = factor(a50), group = factor(a50)),
            size = 1.3) +
  geom_hline(data = brps %>% filter(variable == "MSY"),
             aes(yintercept = value, colour = factor(a50)),
             lty = 2, size = 0.3) +
  geom_vline(data = brps %>% filter(variable == "Fmsy"),
             aes(xintercept = value, colour = factor(a50)),
             lty = 2, size = 0.3) +
  facet_wrap(~mat_ratio) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = a50_cols) +
  labs(colour = "a50", y = "Yield (t)", x = "Fishing mortality",
       linetype = NULL) +
  ggtitle("Maximum sustainable yield by maturation rate")
  
msy_p
ggsave(paste0(figdir, "/MSY.png"), dpi = 300, height = 5, width = 10, units = "in")

msy_p +
  # Reference
  geom_line(data = ref_out %>%
              filter(variable == "y"),
            aes(y = value/0.90718, x = fmort),
            colour = ref_col, size = 1.3) +
  geom_hline(data = ref_brps,
             aes(yintercept = MSY),
             col = ref_col, lty = 2, size = 0.3) +
  geom_vline(data = ref_brps,
             aes(xintercept = Fmsy),
             col = ref_col, lty = 2, size = 0.3)

ggsave(paste0(figdir, "/MSY_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# B0
b0_p <- ggplot() +
  geom_line(data = output %>% 
              filter(variable == "s" & value >= 0),
            aes(y = value/0.90718, x = fmort, 
                colour = factor(a50), group = factor(a50)),
            size = 1.3) +
  geom_hline(data = brps %>% filter(variable == "B0"),
             aes(yintercept = value, colour = factor(a50)),
             size = 0.3) +
  geom_hline(data = brps %>% filter(variable == "Bmsy"),
             aes(yintercept = value, colour = factor(a50)),
             lty = 2, size = 0.3) +
  geom_vline(data = brps %>% filter(variable == "Fmsy"),
             aes(xintercept = value, colour = factor(a50)),
             lty = 2, size = 0.3) +
  facet_wrap(~mat_ratio) + 
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = a50_cols) +
  labs(colour = "a50", y = "Mature biomass (t)", x = "Fishing mortality",
       linetype = NULL) +
  ggtitle("B0 and BMSY by maturation rate")

b0_p
ggsave(paste0(figdir, "/B0_BMSY.png"), dpi = 300, height = 5, width = 10, units = "in")

b0_p +
  # Reference
  geom_line(data = ref_out %>%
              filter(variable == "s"),
            aes(y = value/0.90718, x = fmort),
            colour = ref_col, size = 1.3) +
  geom_hline(data = ref_brps,
             aes(yintercept = B0),
             col = ref_col, size = 0.3) +
  geom_hline(data = ref_brps,
             aes(yintercept = Bmsy),
             col = ref_col, lty = 2, size = 0.3) +
  geom_vline(data = ref_brps,
             aes(xintercept = Fmsy),
             col = ref_col, lty = 2, size = 0.3)
ggsave(paste0(figdir, "/B0_BMSY_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# r0 ----
r0_p <- ggplot() +
  geom_line(data = output %>% 
              filter(variable == "r" & value >= 0),
            aes(y = value, x = fmort, 
                colour = factor(a50), group = factor(a50)),
            size = 1.3) +
  geom_hline(data = brps %>% filter(variable == "r0"),
             aes(yintercept = value, colour = factor(a50)),
             size = 0.3) +
  # geom_hline(data = brps %>% filter(variable == "Bmsy"),
  #            aes(yintercept = value, colour = factor(a50)),
  #            lty = 2, size = 0.3) +
  geom_vline(data = brps %>% filter(variable == "Fmsy"),
             aes(xintercept = value, colour = factor(a50)),
             lty = 2, size = 0.3) +
  facet_wrap(~mat_ratio) + 
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = a50_cols) +
  labs(colour = "a50", y = "Recruitment (millions)", x = "Fishing mortality",
       linetype = NULL) +
  ggtitle("Unfished recruitment (R0) by maturation rate")

r0_p
ggsave(paste0(figdir, "/R0.png"), dpi = 300, height = 5, width = 10, units = "in")

r0_p +
  # Reference
  geom_line(data = ref_out %>%
              filter(variable == "r"),
            aes(y = value, x = fmort),
            colour = ref_col, size = 1.3) +
  geom_hline(data = ref_brps,
             aes(yintercept = r0),
             col = ref_col, size = 0.3) +
  geom_vline(data = ref_brps,
             aes(xintercept = Fmsy),
             col = ref_col, lty = 2, size = 0.3)
ggsave(paste0(figdir, "/R0_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# Range of B0s
brps %>% 
  filter(variable == "B0") %>% 
  summarize(median(value),
            mean(value),
            min(value),
            max(value))
ref_brps

# Status ----

target <- 0.6
limit <- 0.3 

# Reference model
ref_derived %>% 
  mutate(status = matB / (unique(ref_brps$B0)),
         status_target = matB / (target * unique(ref_brps$B0)),
         status_limit = matB / (limit * unique(ref_brps$B0))) -> ref_derived

# Sensitivity
status <- list()

for(i in unique(brps$fn)) {
  
  brps_sub <- brps %>% filter(fn == i)
  derived_sub <- derived %>% filter(fn == i)
  
  tmp <- data.frame(fn = i, 
                    a50 = unique(brps_sub$a50),
                    Year = derived_sub$Year, 
                    mat_ratio = unique(brps_sub$mat_ratio),
                    status = derived$matB / 
                      pull(filter(brps_sub, variable == "B0")),
                    status_target = derived$matB / 
                              (target * (pull(filter(brps_sub, variable == "B0")))),
                    status_limit = derived$matB / 
                              (limit * (pull(filter(brps_sub, variable == "B0")))))
  status[[i]] <- tmp %>% 
    left_join(data.frame(fn = i, 
                         a50 = unique(brps_sub$a50),
                         mat_ratio = unique(brps_sub$mat_ratio)))
  }

status <- bind_rows(status)

# Status plots ----

# Target:
tar <- ggplot() +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, 4)) +
  geom_line(data = status, aes(x = Year, y = status_target, group = fn, colour = factor(a50), linetype = mat_ratio), size = 0.8) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Stock status relative to target (B/B60%) (t)") +
  geom_hline(yintercept = 1, col = "grey", lty = 2, lwd = 2)
tar
ggsave(paste0(figdir, "/status_target.png"), dpi = 300, height = 5, width = 10, units = "in")

tar + geom_line(data = ref_derived, aes(x = Year, y = status_target), col = ref_col, size = 1.3)
ggsave(paste0(figdir, "/status_target_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# Limit:
lim <- ggplot() +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, 8)) +
  geom_line(data = status, aes(x = Year, y = status_limit, group = fn, colour = factor(a50), linetype = mat_ratio), size = 0.8) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Stock status relative to limit (B/B30%) (t)") +
  geom_hline(yintercept = 1, col = "grey", lty = 2, lwd = 2)

lim
ggsave(paste0(figdir, "/status_limit.png"), dpi = 300, height = 5, width = 10, units = "in")

lim +  geom_line(data = ref_derived, aes(x = Year, y = status_limit), col = ref_col, size = 1.3) 
ggsave(paste0(figdir, "/status_limit_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# Status:
stat <- ggplot() +
  scale_colour_manual(values = a50_cols) +
  scale_x_continuous(limits = c(1980, 2020), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  geom_line(data = status, aes(x = Year, y = status, group = fn, colour = factor(a50), linetype = mat_ratio), size = 0.8) +
  labs(x = NULL, y = NULL, colour = "a50", linetype = "Maturation\nrate") +
  ggtitle("Stock status relative to unfished (B/B0) (t)") +
  geom_hline(yintercept = 1, col = "grey", lty = 2, lwd = 2)

stat
ggsave(paste0(figdir, "/status.png"), dpi = 300, height = 5, width = 10, units = "in")

stat +  geom_line(data = ref_derived, aes(x = Year, y = status), col = ref_col, size = 1.3)
ggsave(paste0(figdir, "/status_ref.png"), dpi = 300, height = 5, width = 10, units = "in")

# Performance ----

ref_perf <- data.frame(prod = filter(ref_out, variable == "y") %>% 
                         mutate(value = value / 0.90718) %>% pull(value) / 
                         unique(ref_brps$MSY),
                       biom = filter(ref_out, variable == "s") %>% 
                         mutate(value = value / 0.90718) %>% pull(value) /
                         unique(ref_brps$B0))


perf <- list()

for(i in unique(brps$fn)) {
  
  brps_sub <- brps %>% filter(fn == i)
  out_sub <- output %>% filter(fn == i)
  
  perf[[i]] <- data.frame(fn = i, 
                    a50 = unique(brps_sub$a50),
                    mat_ratio = unique(brps_sub$mat_ratio),
                    prod = filter(out_sub, variable == "y") %>% 
                      mutate(value = value / 0.90718) %>% pull(value) / 
                      brps_sub %>% filter(variable == "MSY") %>% 
                      pull(value),
                    biom = filter(out_sub, variable == "s") %>% 
                      mutate(value = value / 0.90718) %>% pull(value) /
                      brps_sub %>% filter(variable == "B0") %>% 
                      pull(value))
}

perf <- bind_rows(perf)

# # Plot performance ----
# 
# p <- ggplot() +
#   geom_line(data = perf %>% filter(prod > 0),
#             aes(y = prod, x = biom, 
#                 colour = factor(a50), group = factor(a50)),
#             size = 1.3) +
#   geom_line(data = ref_perf, aes(y = prod, x = biom),
#             colour = ref_col) +
#   # geom_hline(data = brps %>% filter(variable == "MSY"),
#   #            aes(yintercept = value, colour = factor(a50)),
#   #            lty = 2, size = 0.3) +
#   # geom_vline(data = brps %>% filter(variable == "Fmsy"),
#   #            aes(xintercept = value, colour = factor(a50)),
#   #            lty = 2, size = 0.3) +
#   facet_wrap(~mat_ratio) +
#   # scale_y_continuous(labels = comma) +
#   scale_colour_manual(values = a50_cols) +
#   labs(colour = "a50", x = "Biomass/B0", y = "Yield/MSY",
#        linetype = NULL) +
#   ggtitle("Production curve by maturation ratio")
# 
# p
# ggsave(paste0(figdir, "/MSY.png"), dpi = 300, height = 5, width = 10, units = "in")
# 
# msy_p +
#   # Reference
#   geom_line(data = ref_out %>%
#               filter(variable == "y"),
#             aes(y = value/0.90718, x = fmort),
#             colour = ref_col, size = 1.3) +
#   geom_hline(data = ref_brps,
#              aes(yintercept = MSY),
#              col = ref_col, lty = 2, size = 0.3) +
#   geom_vline(data = ref_brps,
#              aes(xintercept = Fmsy),
#              col = ref_col, lty = 2, size = 0.3)
# 
# ggsave(paste0(figdir, "/MSY_ref.png"), dpi = 300, height = 5, width = 10, units = "in")
# 
# 
# notransformation <- b0*.907
# transformreck <- b0*.907
# .25*b0*.907
# / 0.90718, # convert to short tons

# fishlife <- filter(diags, a50 == 3.66) %>% pull(nll)
# diags %>% mutate(diff_fishlife = nll - fishlife) -> diags
# diags %>% filter(diff_fishlife >= -10 & diff_fishlife <= 10) -> diags2
# plot(x = diags2$a50, y = diags2$diff_fishlife, type = "l")
# points(x = 3.66, y = 0, col = "green", cex = 2, pch = 20)
# abline(h = 1.92, lty = 2, col = "grey")
# abline(h = -1.92, lty = 2, col = "grey")
# 
# min_nll <- filter(diags, nll == min(nll)) %>% pull(nll)
# diags %>% mutate(diff_minnll = nll - min_nll) -> diags
# diags %>% filter(diff_minnll >= -10 & diff_minnll <= 10) -> diags2
# plot(x = diags2$a50, y = diags2$diff_minnll, type = "l")
# points(x = filter(diags2, nll == min_nll) %>% pull(a50), y = 0, col = "blue", cex = 2, pch = 20)
# abline(h = 1.92, lty = 2, col = "grey")
# abline(h = -1.92, lty = 2, col = "grey")
# 
# filter(diags, nll == min(nll))
# 
# for(i in 1:length(diags %>% filter(a50 < 4.5) %>% pull())) {
#   
#   a50 <- diags$a50[i]
#   a95 <- diags$a95[i]
# 
#   add <- if(i == 1) FALSE else TRUE
#   
#   curve(1.0 / (1.0 + exp(-log(19)*(x-a50)/(a95-a50))), 
#         ylim = c(0,1), from = 3, to = 8, add = add, col = "grey")
# }
# 
# minnll <- filter(diags, nll == min(nll)) %>% select(a50, a95)
# curve(1.0 / (1.0 + exp(-log(19)*(x-pull(minnll, a50))/(pull(minnll, a95)-pull(minnll, a50)))), 
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "blue")
# fshlife <- filter(diags, a50 == 3.66) %>% select(a50, a95)
# curve(1.0 / (1.0 + exp(-log(19)*(x-pull(fshlife, a50))/(pull(fshlife, a95)-pull(fshlife, a50)))), 
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "green")
# curve(1.0 / (1.0 + exp(-log(19)*(x-BC_a50)/(BC_a95-BC_a50))),
#       ylim = c(0,1), from = 3, to = 8, add = TRUE, col = "red")

a50 <- 0.1
k <- 0.1
curve(1.0 / (1.0 + exp(-(x-a50)/k)),
      ylim = c(0,1), from = 3, to = 8, add = TRUE)
a50 <- 3.2
k <- 2.5
curve(1.0 / (1.0 + exp(-k*(x-a50))),
      ylim = c(0,1), from = 3, to = 8)

# 
# LS1_a50
# LS1_a95
# LS2_a50
# LS2_a95
# 1/(1+exp(-1.0*mat_b(t)*((j+2)-mat_a(t))))
# 
# 
# # Weight*mat vs fecundity
# names(D)
# D[["data_fec"]]
# 
# age <- 3:8
# D
# 
# for(i in 1:D[["nFecBlocks"]]) {
#   
#   add <- if(i == 1) FALSE else TRUE
#   
#   curve(D[["fec_inter"]][i] + D[["fec_slope"]][i] * x, 
#        ylim = c(2000,5000), from = 3, to = 8, add = add, col = "grey")
# }
# 
# names(D)
# read_D <- function(name = "mat") {
#   df <- as.data.frame(D[[name]])
#   colnames(df) <- paste(D[["sage"]]:D[["nage"]])
#   df <- melt(df, id.vars = NULL) %>% 
#     filter(value != 0) %>% 
#     group_by(variable) %>% 
#     dplyr::summarise(out = mean(value)) %>% 
#     pull(out)
#   return(df)
# }
# 
# mat <- read_D("mat")
# wa <-  as.data.frame(D[["data_sp_waa"]])
# colnames(wa) <- c("Year", paste(D[["sage"]]:D[["nage"]]))
# wa <-	wa %>% 
#   melt(id.vars = c("Year"), variable.name = "Age", value.name = "weight") %>% 
#   filter(Year >= D[["mod_syr"]]) %>% 
#   group_by(Age) %>% 
#   dplyr::summarise(weight = mean(weight)) %>% 
#   pull(weight)
# 
# mat*wa
