# Figure functions for HER
# Jane Sullivan (jane.sullivan1@alaska.gov)
# Last updated 2019-07-10

# Warning! Don't forget to bootstrap egg variances and survey estimated spawn
# biomass!!

# Egg deposition ----

plot_eggdep <- function(D = D, path = path, 
                        # Bootstrap confidence intervals
                        bootstrap = TRUE,
                        # Plot 95% credibility interval from posterior samples
                        credibility = FALSE,
                        # Cut off extreme value in 2008,
                        remove_outlier = TRUE) {
  
  # From S. Dressel's bootstrap, have to add new values each year - now in LS_byyear$ egg_upper and egg_lower
  U <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 1.259686251, 1.851147636, 1.751126689, 0.560987576, 1.508697833, 1.633749193, 1.695401525, 1.255367509, 2.387620843, 2.795968238, 2.215761696, 1.462234716, 2.289501604, 2.650921062, 6.384923885, 2.279693146, 2.872760889, 29.05686308, 3.863145759, 4.816134565, 4.222205391, 1.634805164, 4.043867944, 1.288746439, 1.332721825, 2.445122264, 1.69865191, 2.4653842)
  L <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 0.985282715, 1.639191663, 1.382705136, 0.48874792, 1.398489625, 1.265029243, 1.286735024, 1.146877561, 1.827147032, 2.534746454, 1.882753246, 1.475607608, 1.863883108, 2.277982827, 3.540565615, 1.707320794, 2.568958439, 14.54490887, 3.237073047, 3.7630828, 3.942674884, 1.578639917, 2.996229014, 1.089460882, 1.154768444, 1.979792303, 1.357551379, 1.7562981)
  
  df <- as.data.frame(D[["data_egg_dep"]][ , 1:2])
  colnames(df) <- c("Year", "obs")
  
  df %>% filter(Year >= D[["mod_syr"]]) -> df
  if(remove_outlier) {
      # 2008 survey egg estimates were extreme high and variable, pmin() helps
      # reduce the scale of the upper confidence interval for easier visualization
    df %>%  mutate(upper = pmin(obs + U, max(df$obs) * 1.2)) -> df
  } else { df %>% mutate(upper = obs + U) -> df }
  df %>% 
    mutate(lower = obs - L,
           pred = D[["pred_egg_dep"]]) -> df
  
  axisx <- tickr(df, Year, 5)
  ggplot(df, aes(x = Year)) +
    # Bootstrap variance estimates
  {if(bootstrap) geom_errorbar(aes(ymin = lower, ymax = upper), 
                               colour = "black",  width = 0, size = 0.01)} +
    # 95% redibility interval
  {if(credibility) geom_ribbon(data = egg_sum, aes(x = Year, ymin = q025, ymax = q975),
                               fill = "grey80", alpha = 0.4)} +
    geom_line(aes(y = pred, linetype = "Estimated from model"), 
              size = 1, colour = "grey") +
    scale_colour_manual(values = "grey") +
    geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
    theme(legend.position = c(0.2, 0.75),
          legend.spacing.y = unit(0, "cm")) +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
    # To cut off extreme values:
    {if(remove_outlier) scale_y_continuous(limits = c(0, max(df$obs) * 1.2))} +
    labs(x = NULL, y = "Eggs spawned (trillions)\n", shape = NULL, linetype = NULL) -> obsfit
  
  # residuals
  data.frame(Year = D[["year"]],
             resids = D[["resd_egg_dep"]]) %>% 
    ggplot(aes(x = Year, y = resids)) + 
    geom_hline(yintercept = 0, colour = "grey", size = 1) +
    geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
                 size = 0.2, colour = "grey") +
    geom_point() +
    labs(x = "\nYear", y = "Residuals\n") +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) -> resids
  
  cowplot::plot_grid(obsfit, resids, align = "hv", nrow = 2) -> eggdep_plot
  
  ggsave(paste0(path, "/eggdep_plot.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 7, units = "in")
}


# Recruitment ----

plot_recruitment <- function(D = D, path = path, showRicker = FALSE) {
  
  # Age-3 recruits compared recruitment estimated from Ricker
  df <- data.frame(Year = D[["years"]],
                   age3 = D[["Nij"]][, 1]) %>% 
    left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                    D[["mod_nyr"]] + 1, 1),
                         SR = D[["recruits"]],
                         resids = D[["resd_rec"]]), by = "Year")
  
  axisr <- tickr(df, Year, 5)
  ggplot(df, aes(x = Year)) +
    geom_bar(aes(y = age3, colour = "ASA model"), 
             stat = "identity", fill = "lightgrey", 
             width = 0.8, position = position_dodge(width = 0.5)) +
    # If showRicker = TRUE, show's Ricker annual estimates
             {if(showRicker) geom_line(aes(y = SR, linetype = "Ricker model"), colour = "black", size = 1)} +
    scale_colour_manual(values = "grey") +
    scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) +
    scale_y_continuous(labels = scales::comma) + 
    theme(legend.position = c(0.1, 0.75),
          legend.key = element_rect(size = 0.5),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, "cm")) +
    labs(x = "\nYear", y = "Age-3 recruits (millions)\n", 
         colour = NULL, linetype = NULL) -> recruits
  
  # residuals
  ggplot(df, aes(x = Year, y = resids)) + 
    geom_hline(yintercept = 0, colour = "grey", size = 1) +
    geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
                 size = 0.2, colour = "grey") +
    geom_point() +
    labs(x = "\nYear", y = "Residuals\n") +
    scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) -> resids
  
  # spawner-recruit curve
  df %>% 
    # Join the spawners to the recruitment years they're associated to (brood year =
    # recruitment year - 3)
    left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                    D[["mod_nyr"]] + 1, 1),
                         ssb = D[["spawners"]]), by = "Year") %>% 
    filter(Year >= D[["mod_syr"]] + D[["sage"]]) -> df
  
  ggplot(df, aes(x = ssb)) +
    geom_line(aes(y = SR, linetype = "Ricker model"), 
              colour = "grey", size = 1) +
    geom_point(aes(y = age3, colour = "ASA model")) +
    scale_colour_manual(values = "black") +
    theme(legend.position = c(0.1, 0.75),
          legend.key = element_rect(size = 0.5),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, "cm")) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    geom_text_repel(aes(y = age3, label = Year), size = 2) +
    labs(x = "\nSpawning stock biomass (tons)", y = "Age-3 recruits (millions)\n",
         linetype = NULL, colour = NULL) -> sr_curve
  
  cowplot::plot_grid(recruits, resids, sr_curve, align = "hv", nrow = 3) -> recruit_plot
  
  ggsave(paste0(path, "/recruit_plot.png"), plot = recruit_plot, dpi = 300, height = 8, width = 7, units = "in")
}

# Survival blocks ----

plot_survival <- function(D = D, path = path) {
  df <- data.frame(D[["year"]], 
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  
  df %>% select(Year, M = `3`) %>% 
    mutate(survival = exp(-M)) -> df
  
  axisx <- tickr(df, Year, 5)
  
  ggplot(df, aes(x = Year, y = survival)) +
    # geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3) +
    geom_line(size = 1) +
    geom_point() +
    lims(y = c(0, 1)) +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
    labs(x = "", y = "Survival\n") -> survival_plot
  
  ggsave(paste0(path, "/survival.png"), plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")
}

# Maturity/selectivity ----

plot_matsel <- function(D = D, path = path) {
  
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
  
  df %>% 
    gather("Age", "proportion", -c(Year, param)) %>% 
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+")),
           # To make sure selectivity is differentiable, it was scaled to have a
           # mean of 1 across all ages. This was done in log space by substracting
           # the mean from the vector of age-specific selectivities. See Tech Doc
           # p 11. Here we normalize it from 0 to 1.
           proportion = ifelse(param == "Selectivity", (proportion - 0)/(max(proportion) - 0), 
                               proportion)) %>% 
    group_by(Age, param, proportion) %>% 
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max)) %>% 
    group_by(param) %>% 
    mutate(Model = "HER",
           combos = paste0(Model, " ", `Time blocks`)) -> her_matsel
  
  # Maturity
  her_matsel %>% filter(param == "Maturity") -> par

  ggplot(par, aes(x = Age, y = proportion)) +
    geom_line(aes(linetype = `Time blocks`, group = `Time blocks`)) +
    geom_hline(yintercept = 0.5, colour = "lightgrey", alpha = 0.5, size = 0.4, linetype = 2) +
    geom_vline(xintercept = 3.5, colour = "lightgrey", alpha = 0.5, size = 0.4, linetype = 2) +
    expand_limits(y = 0) +
    labs(x = "\nAge", y = "Proportion\n", linetype = "Time blocks") +
    ggtitle(paste0(par$param[1])) +
    theme(legend.position = c(0.7, 0.2),
          plot.title = element_text(hjust = 0.5)) -> mat
  
  # Selectivity
  her_matsel %>% filter(param == "Selectivity") -> par

  ggplot(par, aes(x = Age, y = proportion)) +
    geom_line(aes(linetype = `Time blocks`, group = `Time blocks`)) +
    geom_hline(yintercept = 0.5, colour = "lightgrey", alpha = 0.5, size = 0.4, linetype = 2) +
    geom_vline(xintercept = 3.5, colour = "lightgrey", alpha = 0.5, size = 0.4, linetype = 2) +
    expand_limits(y = 0) +
    labs(x = "\nAge", y = NULL, linetype = "Time blocks") +
    scale_y_continuous(breaks = seq(0, max(par$proportion), .25)) +
    ggtitle(paste0(par$param[1])) +
    theme(legend.position = c(0.7, 0.2),
          plot.title = element_text(hjust = 0.5)) -> sel
  
  cowplot::plot_grid(mat, sel, align = "hv", nrow = 1) -> matsel_plot
  
  ggsave(paste0(path, "/mat_sel.png"), plot = matsel_plot, dpi = 300, height = 4, width = 7, units = "in")
}

# Age compositions ----

plot_agecomps <- function(D = D, path = path) {
  
  df <- data.frame("Commercial fishery",
                   D[["data_cm_comp"]]) # change number of rows
  colnames(df) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
  sp <- data.frame("Cast net survey",
                   D[["data_sp_comp"]]) # change number of rows
  colnames(sp) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
  df <- bind_rows(df, sp); rm(sp)
  
  df %>% 
    filter(Year >= D[["mod_syr"]]) %>% 
    melt(id.vars = c("Year", "Source"), variable.name = "Age", value.name = "obs") %>% 
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) -> her_agecomps
  
  axisy <- tickr(her_agecomps, Year, 5)
  
  # Data bubble plot----
  
  ggplot(her_agecomps, aes(x = Age, y = Year, size = obs)) + #*FLAG* could swap size with proportion_scaled
    geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
    geom_point(shape = 21, colour = "black", fill = "black") +
    scale_size(range = c(0, 4)) +
    facet_wrap(~ Source) +
    xlab('\nAge') +
    ylab('') +
    guides(size = FALSE) +
    scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) -> agecomps_bubbleplot
  
  ggsave(paste0(path, "/agecomps_bubbleplot.png"), plot = agecomps_bubbleplot, dpi = 300, height = 6, width = 7, units = "in")
  
  # Residual bubble plot ----
  
  # Model predictions and residuals
  predcm <- data.frame("Commercial fishery",
                       D[["year"]],
                       D[["pred_cm_comp"]]) # change number of rows
  colnames(predcm) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
  
  predsp <- data.frame("Cast net survey",
                       D[["year"]],
                       D[["pred_sp_comp"]]) # change number of rows
  colnames(predsp) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
  df <- bind_rows(predcm, predsp); rm(predcm, predsp)
  
  df %>% 
    melt(id.vars = c("Year", "Source"), variable.name = "Age", value.name = "pred") %>% 
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
    # join to observed age comps df
    left_join(her_agecomps) %>%
    group_by(Source) %>% 
    mutate(raw = obs - pred,
           # positive or negative
           `Model performance` = ifelse(raw >= 0, "Observed greater than estimated", "Observed less than estimated"),
           # sample size
           n = length(obs),
           # Raw
           resid = abs(raw),
           # Deviance residual http://data.princeton.edu/wws509/notes/c3s8.html,
           # deviance residual has the same direction as raw residual
           dev_resid = sqrt(2 * (obs * log(obs / pred ) + (n - obs) * log((n - obs)/(n - pred)))),
           direction = ifelse(raw >= 0, 1, -1),
           dev = direction * dev_resid,
           # Pearson's residuals
           pearson = (obs - pred)/ sqrt(var(pred))) -> her_agecomps
  
  ggplot(her_agecomps, aes(x = Age, y = Year, size = resid,
                           fill = `Model performance`)) + 
    geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
    geom_point(shape = 21, colour = "black") +
    scale_size(range = c(0, 4.5)) +
    facet_wrap(~ Source) +
    labs(x = '\nAge', y = '') +
    guides(size = FALSE) +
    scale_fill_manual(values = c("white", "black")) +
    scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) +
    theme(legend.position = "bottom") -> agecomps_residplot
  
  ggsave(paste0(path, "/agecomps_residplot.png"), plot = agecomps_residplot, dpi = 300, height = 6, width = 7, units = "in")
  
  # Barplots ----
  
  her_agecomps %>% # fishery
    filter(Source == "Commercial fishery") %>% 
    ggplot() + 
    geom_bar(aes(x = Age, y = obs), 
             stat = "identity", colour = "grey", fill = "lightgrey",
             width = 0.8, position = position_dodge(width = 0.5)) +
    geom_line(aes(x = Age, y = pred, group = 1), size = 0.6) +
    facet_wrap(~ Year, dir = "v", ncol = 5) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
    labs(x = '\nAge', y = 'Proportion-at-age\n') +
    ggtitle("Commercial fishery") + 
    theme(plot.title = element_text(hjust = 0.5)) -> catchage_barplot
  
  ggsave(paste0(path, "/catchage_comps_barplot.png"), plot = catchage_barplot, dpi = 300, height = 8, width = 6, units = "in")
  
  her_agecomps %>% # survey
    filter(Source == "Cast net survey") %>% 
    ggplot() + 
    geom_bar(aes(x = Age, y = obs), 
             stat = "identity", colour = "grey", fill = "lightgrey",
             width = 0.8, position = position_dodge(width = 0.5)) +
    geom_line(aes(x = Age, y = pred, group = 1), size = 0.6) +
    facet_wrap(~ Year, dir = "v", ncol = 5) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
    labs(x = '\nAge', y = 'Proportion-at-age\n') +
    ggtitle("Cast net survey") + 
    theme(plot.title = element_text(hjust = 0.5)) -> spawnage_barplot
  
  ggsave(paste0(path, "/spawnage_comps_barplot.png"), plot = spawnage_barplot, dpi = 300, height = 8, width = 6, units = "in")
  
}

# Biomass time series ----

plot_biomass <- function(D = D, path = path) {
  
  # in SPAWN folder, numbers from excel spread sheet for spawn deposition -- in
  # 2019, got this from LS_byyear$spawn_dep
  srv_spB <- c(35000,30000,29500,23500,38500,31000,25000,46000,58500,27000,23000,23500,
               43351,37150,14941,34990,40827,28611,34942,44554,57988,58756,40366,55769,
               69907,101305,66111,84501,247088,110946,126230,161904,62518,103267,48561, 
               58183,77973,46919,55009)
  
  tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1
  syr_index <- D[["mod_syr"]] - D[["dat_syr"]] + 1

  df <- data.frame(Year = D[["mod_syr"]]:D[["mod_nyr"]],
                   matB = D[["mat_B"]] / 0.90718, # convert to short tons
                   catch = D[["data_catch"]][syr_index:tot_yrs, 2]  / 0.90718,
                   spB = D[["sp_B"]]) %>% 
    # mutate(spB_derived = matB - catch,
            #, #spB_estimated
           # diff = spB_derived - spB_estimated, 
           # pred_catch = D[["pred_catch"]] * 1.10231,
           #diff_catch = catch - pred_catch
           # ) %>% 
    gather("Biomass", "tons", -Year) %>% 
    mutate(Type = "Estimate") %>% 
    # Bind to forecast values
    bind_rows(data.frame(Year = D[["mod_nyr"]] + 1,
                         matB = D[["fore_matb"]] / 0.90718,
                         catch = D[["ghl"]]) %>% 
                mutate(spB = matB - catch) %>% 
                gather("Biomass", "tons", -Year) %>% 
                mutate(Type = "Forecast")) %>% 
    mutate(label = ifelse(Type == "Forecast", 
                   prettyNum(tons, big.mark = ",", digits = 1), NA))
  
  # Historical estimates from survey
  survey <- data.frame(Year = 1980:D[["mod_nyr"]],
                         srv_spB = srv_spB,
                         catch = 1980:D[["mod_nyr"]] / 0.90718) %>%
    mutate(srv_matB = srv_spB + catch)
  
  tickryr <- data.frame(Year = D[["mod_syr"]]:max(df$Year)+3)
  axisx <- tickr(tickryr, Year, 5)
  
  # Spawning biomass ----
  ggplot() +
    geom_point(data = df %>% filter(Biomass == "spB" & Type == "Forecast"), 
               aes(x = Year, y = tons, colour = "Model predictions and forecast",
                   fill = "Model predictions and forecast"), size = 2, shape = 23) +
    geom_line(data = df %>% filter(Biomass == "spB" & Type == "Estimate"), 
              aes(x = Year, y = tons, colour = "Model predictions and forecast", group = 1), size = 1) +
    scale_colour_manual(values = "darkgrey") + 
    scale_fill_manual(values = "darkgrey") + 
    geom_point(data = survey, aes(x = Year, y = srv_spB, shape = "Historical estimates from survey")) +
    scale_shape_manual(values = 1) +
    geom_text(data = df %>% filter(Biomass == "spB"), 
                    aes(x = Year, y = tons, label = label), size = 3, nudge_x = 2) +
    theme(legend.position = c(0.25, 0.8),
          legend.spacing.y = unit(0, "cm")) +
    scale_x_continuous(limits = c(min(tickryr$Year), max(tickryr$Year)),
                                  breaks = axisx$breaks, labels = axisx$labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "", y = "Spawning biomass (tons)\n", shape = NULL, colour = NULL, fill = NULL) -> spbiomass_plot
  
  ggsave(paste0(path, "/spbiomass_plot.png"), plot = spbiomass_plot, dpi = 300, height = 4, width = 7, units = "in")
  
  # Mature biomass ----
  ggplot() +
    geom_point(data = df %>% filter(Biomass == "matB" & Type == "Forecast"), 
               aes(x = Year, y = tons, colour = "Model predictions and forecast",
                   fill = "Model predictions and forecast"), size = 2, shape = 23) +
    geom_line(data = df %>% filter(Biomass == "matB"), 
              aes(x = Year, y = tons, colour = "Model predictions and forecast"), size = 1) +
    scale_colour_manual(values = "darkgrey") + 
    scale_fill_manual(values = "darkgrey") + 
    geom_point(data = survey, aes(x = Year, y = srv_matB, shape = "Historical estimates from survey")) +
    scale_shape_manual(values = 1) +
    geom_text(data = df %>% filter(Biomass == "matB" & Type == "Estimate"), 
              aes(x = Year, y = tons, label = label), size = 3, nudge_x = 2) +
    theme(legend.position = c(0.25, 0.8),
          legend.spacing.y = unit(0, "cm")) +
    scale_x_continuous(limits = c(min(tickryr$Year), max(tickryr$Year)),
                       breaks = axisx$breaks, labels = axisx$labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "", y = "Mature biomass (tons)\n", shape = NULL, colour = NULL, fill = NULL) -> matbiomass_plot
  
  ggsave(paste0(path, "/matbiomass_plot.png"), plot = matbiomass_plot, dpi = 300, height = 4, width = 7, units = "in")
  
}
