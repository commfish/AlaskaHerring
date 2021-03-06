
# Libraries, ggplot themes, and user-defined fxns.
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov

# libraries ----

if(!require("mosaic"))   install.packages("mosaic") # derivedFactor, derivedVariable. masks over a lot of fxns, but generally only improves their utility
if(!require("tidyverse"))   install.packages("tidyverse") # dplyr, ggplot, etc.
if(!require("lubridate"))   install.packages("lubridate") # dates functions like yday, dmy/mdy
# if(!require("mgcv"))   install.packages("mgcv") # gams
# if(!require("gridExtra"))   install.packages("gridExtra") # multipanneled plots
if(!require("data.table"))   install.packages("data.table") # dcast, foverlaps
# if(!require("ROracle"))   install.packages("ROracle") # database access through R
if(!require("broom"))   install.packages("broom") # tidying coda summary output
# if(!require("padr"))   install.packages("padr") # fills in missing values in a time series
if(!require("knitr"))   install.packages("knitr") # r markdown
# if(!require("forcats"))   install.packages("forcats") # releveling factors
if(!require("cowplot"))   install.packages("cowplot") # plot_grid and so much else
if(!require("R2admb"))   install.packages("R2admb") # run admb from r
if(!require("ggthemes"))   install.packages("ggthemes") # access to 150 colour palettes from canva.com design school, among other things
if(!require("scales"))   install.packages("scales") # add comma to ggplot axis with scale_y_countinuous(label = comma)
if(!require("ggrepel"))   install.packages("ggrepel") # readable labels with geom_text_repel()
if(!require("gridExtra"))   install.packages("gridExtra") # tableGrob()
if(!require("captioner"))   install.packages("captioner") #numbering, ordering, & creating captions for tables and figures


# https://www.canva.com/learn/100-color-combinations/
# http://makeadifferencewithdata.com/wp-content/uploads/2016/12/color-palettes.txt
# if(!require("scales"))   install.packages("scales") # used to expand colour palettes
# install.packages("devtools")
# devtools::install_github("ben-williams/FNGr")
# library("FNGr")

# ggplot themes ----

windowsFonts(Times=windowsFont("Times New Roman"))

theme_sleek <- function(base_size = 12, base_family = "Times") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Times") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      #axis.text = element_text(colour = "grey30"),
      #axis.title = element_text(colour = "grey30"),
      #legend.title = element_text(colour = "grey30"),#, size = rel(0.9)
      panel.border = element_rect(fill = NA),#, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      #legend.text = element_text(size = rel(0.7)),#, colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)#,
      #plot.title = element_text(colour = "grey30"),#, size = rel(1)
      #plot.subtitle = element_text(colour = "grey30")#, size = rel(.85)
    )
}

theme_set(theme_sleek())

# User-defined functions ----

# tickr() - makes nice ticks and breaks on your ggplots

# Depends on dplyr
tickr <- function(
  data, # dataframe
  var, # column of interest
  to # break point definition 
){
  
  VAR <- enquo(var) # makes VAR a dynamic variable
  
  data %>% 
    distinct(!!VAR) %>%
    ungroup(!!VAR) %>% 
    mutate(labels = ifelse(!!VAR %in% seq(to * round(min(!!VAR) / to), max(!!VAR), to),
                           !!VAR, "")) %>%
    select(breaks = UQ(VAR), labels)
}

# Functions from SJD Martell, formerly in globals.r

read_admb <- function(ifile) {
    ret=read_fit(ifile)

    fn=paste(ifile,'.rep', sep='')
    A=read_rep(fn)
    A$fit=ret

    pfn=paste(ifile,'.psv',sep='')
    if(file.exists(pfn))
      A$post.samp=read_psv(pfn)

    return(A)
  }

# read_admb <-
#   function(ifile)
#   {
#     fn=paste(ifile,'.rep', sep = '')
#     A=read_rep(fn)
# 
#     if(file.exists(paste(ifile, '.cov', sep = ''))) {
# 
#       ret = read_fit(ifile)
#       A$fit = ret
#     } else {
#       A$fit="Not_converged"
#     }
# 
#     pfn=paste(ifile,'.psv',sep = '')
#     if(file.exists(pfn))
#       A$post.samp=read_psv(pfn)
# 
#     return(A)
#   }


read_fit <-
  function(ifile)
  {
    # __Example:             
    #	file <-("~/admb/simple")
    #	A <- reptoRlist(file)
    #	Note there is no extension on the file name.
    
    ## The following is a contribution from:
    ## Anders Nielsen that reads the par & cor files.
    ret<-list() 
    parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
                             what='', n=16, quiet=TRUE)[c(6,11,16)]) 
    ret$nopar<-as.integer(parfile[1]) 
    ret$nlogl<-parfile[2] 
    ret$maxgrad<-parfile[3] 
    file<-paste(ifile,'.cor', sep='') 
    lin<-readLines(file) 
    ret$npar<-length(lin)-2 
    ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
    sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
    ret$names<-unlist(lapply(sublin,function(x)x[2])) 
    ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
    ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
    ret$cor<-matrix(NA, ret$npar, ret$npar) 
    corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
    ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
    ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
    ret$cov<-ret$cor*(ret$std%o%ret$std)
    return(ret)
  }

read_rep <- 
  function(fn)
  {
    # The following reads a report file
    # Then the 'A' object contains a list structure
    # with all the elemements in the report file.
    # In the REPORT_SECTION of the AMDB template use 
    # the following format to output objects:
    #  	report<<"object \n"<<object<<endl;
    #
    # The part in quotations becomes the list name.
    # Created By Steven Martell
    options(warn=-1)  #Suppress the NA message in the coercion to double
    
    
    ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
    idx=sapply(as.double(ifile),is.na)
    vnam=ifile[idx] #list names
    nv=length(vnam) #number of objects
    A=list()
    ir=0
    for(i in 1:nv)
    {
      ir=match(vnam[i],ifile)
      if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
      dum=NA
      if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
      if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))
      
      if(is.numeric(dum))#Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]]=dum
      }
    }
    options(warn=0)
    
    return(A)
  }

read_psv <- function(fn, nsamples=11112) { #=11111
  #This function reads the binary output from ADMB
  #-mcsave command line option.
  #fn = paste(ifile,'.psv',sep='')
  filen <- file(fn, "rb")
  nopar <- readBin(filen, what = integer(), n = 1)
  mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
  mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
  close(filen)
  return(mcmc)
}

# Function to just read max gradient component (for quick and dirty convergence
# check, used in sensitivity and retrospective analyses)
readMGC <- function(file){
  mgc <- as.numeric(scan(paste(file,'.par', sep=''), what='', n=16, quiet=TRUE)[c(16)])
  return(mgc)
}

# A simple function for creating transparent colors
# Author: Nathan Stephens (hacks package)
colr <- 
  function(col.pal=1,a=1)
  {
    col.rgb<-col2rgb(col.pal)/255
    rgb(t(col.rgb),alpha=a)
  }

#fig_nums, tbl_nums, and appendix_nums fxns created fromm captioner() fxn in
#'captioner' library, which I've tweaked below, changed separator from a colon
#to period) - these fxns are used for autonumbering figs and tables in text and
#creating captions.

captioner <- function (prefix = "Figure", auto_space = TRUE, levels = 1, 
                       type = NULL, infix = ".") {
  check_class(prefix, "character")
  check_class(auto_space, "logical")
  check_class(levels, "numeric")
  check_class(infix, "character")
  if (is.null(type)) {
    type <- c(rep("n", times = levels))
  }
  else if (length(type) < levels) {
    type[(length(type) + 1):levels] <- "n"
  }
  else if (length(type) > levels) {
    type <- type[1:levels]
  }
  if (!all(type %in% c("n", "c", "C"))) {
    stop("Invalid 'type' value used.  Expecting 'n', 'c', or 'C'.")
  }
  if (auto_space) {
    prefix <- paste(prefix, " ")
  }
  force(levels)
  force(prefix)
  force(infix)
  OBJECTS <- list(name = NULL, caption = NULL, number = list(list()))
  OBJECTS$number[[1]][which(type == "n")] <- 1
  OBJECTS$number[[1]][which(type == "c")] <- "a"
  OBJECTS$number[[1]][which(type == "C")] <- "A"
  function(name, caption = "", display = "full", level = FALSE, 
           cite = FALSE, num = FALSE) {
    if (level > levels) {
      stop("Level too large.")
    }
    objects <- OBJECTS
    if (any(objects$name == name)) {
      obj_ind <- match(name, objects$name)
      if (objects$caption[obj_ind] == "") {
        objects$caption[obj_ind] <- caption
      }
      else {
        caption <- objects$caption[obj_ind]
      }
    }
    else {
      obj_ind <- length(objects$name) + 1
      if (length(objects$number) == length(objects$name)) {
        if (level) {
          objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
                                                                   1]], level)
        }
        else {
          objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
                                                                   1]], levels)
        }
      }
      objects$name[obj_ind] <- name
      objects$caption[obj_ind] <- caption
    }
    assign("OBJECTS", objects, envir = parent.env(environment()))
    obj_num <- paste(objects$number[[obj_ind]], collapse = infix)
    if (cite) {
      .Deprecated(new = "display", old = "cite")
      return(paste0(prefix, obj_num))
    }
    if (num) {
      .Deprecated(new = "display", old = "num")
      return(obj_num)
    }
    if (display == FALSE) {
      return(invisible())
    }
    #FLAG: Jane changed ": " to ". "
    else if (display == "full" || display == "f") {
      return(paste0(prefix, obj_num, ". ", caption))
    }
    else if (display == "cite" || display == "c") {
      return(paste0(prefix, obj_num))
    }
    else if (display == "num" || display == "n") {
      return(obj_num)
    }
    else {
      warning("Invalid display mode used.  Caption was still saved.")
      return(invisible())
    }
  }
}

fig <- captioner(prefix = "Figure")

tbl <- captioner(prefix = "Table")

appendix_tbl <- captioner(prefix = "Table") #Numbers tables in the appendix

appendix_fig <- captioner(prefix = "Figure") #Numbers figures in the appendix

# For posterior sample output that is structured by year (i.e. number of columns
# = number of model years and number of rows = number of iterations). Unit_conv
# is there to deal with the conversion of metric tons to short tons (st = mt /
# 0.90718)
ps_byyear <- function(fn = "sp_B", 
                      syr = D[["mod_syr"]], 
                      lyr = D[["mod_nyr"]],
                      unit_conv = 1,
                      burn = burn_in * (niter / thin + 1),
                      save = TRUE,
                      path = HERfig_dir) {
  require(data.table)
  df <- fread(paste0(fn, ".ps"))
  colnames(df) <- paste(syr:lyr)
  df[, iter := .I] # row number = iteration number
  df <- melt(df, id.vars = c("iter"), variable.name = "Year")
  
  df <- df[iter > burn, ] # Eliminate burn in
  
  df[, value := value / unit_conv] # if needed, convert to short tons
  
  df[, `:=` (mean = mean(value), # `:=` is the same as dplyr::mutate on multiple cols
             median = median(value),
             # 95% 
             q025 = quantile(value, 0.025),
             q975 = quantile(value, 0.975),
             # 50% 
             q250 = quantile(value, 0.250),
             q750 = quantile(value, 0.750),
             # 5% 
             q475 = quantile(value, 0.475),
             q525 = quantile(value, 0.525)),
     by = Year] 
  
  sum <- unique(df, by = c("Year", "mean", "median", "q025", "q975", "q250", "q750", "q475", "q525"))
  sum[, iter:=NULL]
  sum[, value:=NULL]
  
  if(save == TRUE) {write_csv(sum, path = paste0(path, "/", fn, "_posterior.csv"))}
  
  return(df)
}

# For posterior sample output that is structured by age and time blocks (i.e.
# number of columns = number of ages, number of rows = number of years * number
# of iterations)
ps_byage <- function(fn = "maturity", 
                     syr = D[["mod_syr"]], 
                     lyr = D[["mod_nyr"]], 
                     n = niter / thin + 1,
                     burn = burn_in * (niter / thin + 1),
                     save = TRUE,
                     path = HERfig_dir) {
  
  require(data.table)
  df <- fread(paste0(fn, ".ps"))
  colnames(df) <- paste(D[['sage']]:D[['nage']])
  df[, Year := rep(syr:lyr, n)]
  df[, iter := rowid(Year)] 
  df <- melt(df, id.vars = c("Year", "iter"), variable.name = "Age")
  
  df <- df[iter > burn, ] # Eliminate burn in
  
  # akin to dplyr::summarize
  df <- df[, list(min = min(Year), 
                  max = max(Year)),
           by = .(iter, Age, value)]
  
  df[, Blocks := paste0(min, "-", max), by = .(iter, Age, value)]
  
  # `:=` is the same as dplyr::mutate on multiple cols
  df[, `:=` (mean = mean(value),
             median = median(value),
             # 95% 
             q025 = quantile(value, 0.025),
             q975 = quantile(value, 0.975),
             # 50% 
             q250 = quantile(value, 0.250),
             q750 = quantile(value, 0.750)),
     by = .(Age, Blocks)] 
  
  # Eliminate accidental blocks 
  keep <- df %>% group_by(Blocks) %>% 
    summarize(tst = length(levels(Age))) %>% 
    filter(tst == (D[["nage"]] - D[["sage"]] + 1)) %>% 
    pull(Blocks)
  df <- df[Blocks %in% keep, ]
  
  sum <- unique(df, by = c("Age", "Blocks", "mean", "median", "q025", "q975", "q250", "q750"))
  sum[, `:=` (iter = NULL,
              value = NULL,
              min = NULL,
              max = NULL)]
  if(save == TRUE) {write_csv(sum, path = paste0(path, "/", fn, "_posterior.csv"))}
  
  return(df)
} 

# For age composition posterior sample output (number of columns
# = number of ages, number of rows = number of years * number of iterations)
ps_comps <- function(fn = "pred_sp_comp", 
                     syr = D[["mod_syr"]], 
                     lyr = D[["mod_nyr"]], 
                     n = niter / thin + 1,
                     burn = burn_in * (niter / thin + 1),
                     save = TRUE,
                     path = HERfig_dir) {
  
  require(data.table)
  df <- fread(paste0(fn, ".ps"))
  colnames(df) <- paste(D[['sage']]:D[['nage']])
  df[, Year := rep(syr:lyr, n)]
  df[, iter := rowid(Year)] # row number
  df <- melt(df, id.vars = c("Year", "iter"), variable.name = "Age")
  df <- df[iter > burn, ] # Eliminate burn in
  
  # `:=` mutate multiple columns, by same as group_by
  df[, `:=` (mean = mean(value),
             median = median(value),
             q025 = quantile(value, 0.025),
             q975 = quantile(value, 0.975)),
     by = .(Year, Age)] 
  
  sum <- unique(df, by = c("Year", "Age", "mean", "median", "q025", "q975"))
  
  sum[, `:=` (iter = NULL,
              value = NULL)]
  
  fn <- str_split(fn, "pred_")[[1]][2]
  if(save == TRUE) {write_csv(sum, path = paste0(HERfig_dir, "/", fn, "_posterior_predictive.csv"))}
  
  return(df)
  }

# Summarize posterior samples of the multivariate logistic variance estimates 
ps_tau <- function(fn = "pp_sp_tau2", 
                   n = niter / thin + 1,
                   burn = burn_in * (niter / thin + 1),
                   save = TRUE,
                   path = HERfig_dir) {
  
  require(data.table)
  df <- fread(paste0(fn, ".ps")) # Read in variance estimates
  colnames(df) <- "tau2"
  df[, iter := .I] # row number
  df <- df[iter > burn, ] # Eliminate burn in
  
  sum <- df[, `:=` (mean = mean(tau2),
                    median = median(tau2),
                    q025 = quantile(tau2, 0.025),
                    q975 = quantile(tau2, 0.975))] 
  
  sum <- unique(df, by = c("mean", "median", "q025", "q975"))
  sum[, `:=` (iter = NULL,
              tau2 = NULL),]
  
  fn <- str_split(fn, "pp_")[[1]][2]
  if(save == TRUE) {write_csv(sum, path = paste0(HERfig_dir, "/", fn, "_posterior.csv"))}

  return(df)
  }

# Posterior predictive interval for multivariate logistic distribution (age
# compositions). Uses the posterior sample estimates and variance.
ppi_rmvlogistic <- function(fn = "sp_comp",
                            ps_sum = sp_comp_sum,  # summarized posterior samples (output from ps_comps)
                            tau2 = sp_tau2,    # summarized variance estimates (output from ps_tau)
                            save = TRUE,
                            path = HERfig_dir
) {
  
  df <- unique(ps_sum[, list(Year, iter, Age, value)]) # same as dplyr::distinct()
  df <- df[tau2, on = 'iter'] # same as dplyr::left_join()
  
  # Modeled after rmvlogistic() function in ADMB:
  # http://api.admb-project.org/rmvlogistic_8cpp_source.html
  df[, x := log(value) + tau2 * rnorm(.I)]
  df[, mu_x := mean(x), by = .(iter, Year)]
  df[, x := x - mu_x]
  df[, sum_x := sum(exp(x)), by = .(iter, Year)]
  df[, new_value := exp(x) / sum_x]
  
  # Get posterior predictive intervals (syntax akin to dplyr::summarize)
  df <- df[, list(
    # 95% 
    q025 = quantile(new_value, 0.025),
    q975 = quantile(new_value, 0.975),
    # 50% 
    q250 = quantile(new_value, 0.250),
    q750 = quantile(new_value, 0.750)),
    by = .(Age, Year)] 
  
  sum <- unique(ps_sum, by = c("Year", "Age", "mean", "median"))
  sum[, `:=` (iter = NULL,
              value = NULL,
              # get rid of credible intervals to make way for posterior predictive
              q025 = NULL, 
              q975 = NULL)]

  sum <- sum[df, on = c("Year", "Age")] # same as dplyr::left_join()
  
  if(save == TRUE) {write_csv(sum, path = paste0(HERfig_dir, "/", fn, "_posterior_predictive.csv"))}
  
  return(df)
}

# Function to plot trace plots for derived variables
mcmc_plot <- function(df = srv_sum, 
                      type = "ps_byyear",
                      name = "Survival",
                      save = TRUE, 
                      path = HERmcmc_dir,
                      height = 8,
                      width = 8) {
  
  # Base plot
  p <- ggplot(df, aes(x = iter, y = value, group = 1)) +
    geom_line() +
    ggtitle(paste0(name, " posterior samples"))
  
  if(type == "ps_byyear") { p <- p + facet_wrap(~ Year) }
  if(type == "ps_byage") { p <- p + facet_grid(Age ~ Blocks) }
  if(type == "ps_comps") { p <- p + facet_grid(Year ~ Age) }
  
  if(save == TRUE) {
    ggsave(paste0(HERmcmc_dir, "/", name, "_caterpillar.png"), plot = p, dpi = 300, height = height, width = width, units = "in")
    
  }
  
}

# Functions for getting BRPs  ----

# Derive steepness (h) from compensation ratio (reck)
get_steepness <- function(rec_mod) {
  
  # Ricker
  if(rec_mod == 1) { h <- reck / (4 + reck) }
  
  # Beverton-Holt
  if(rec_mod == 2) { h <- 0.2 * reck ^ (4 / 5) }
  
  return(h)
}

# Survivorship
get_survivorship <- function(fmort) {
  
  lx <- matrix(ncol = length(unique(age)), 
               nrow = length(unique(fmort)))
  
  a_plus <- max(length(age)) 
  
  # Loop over all input values of mortality
  for(f in 1:length(fmort)) {
    
    # Initialize matrix
    lx[f,1] <- 1  
    
    for(a in 2:(a_plus-1)) {
      lx[f,a] <- lx[f,a-1] * exp(-mort[a-1] - fmort[f] * sel[a-1]) }
    
    # Plus group
    lx[f,a_plus] <- (lx[f,a_plus-1] * exp(- mort[a_plus-1] - fmort[f] * sel[a_plus-1])) / (1 - exp(- mort[a_plus] - fmort[f] * sel[a_plus]))    
  }
  return(lx)
}

# Yield-per-recruit (note that ypr = phi_q * fmort in Martell et al. 2009)
get_ypr <- function(fmort, lx) { 
  
  ypr <- matrix(ncol = length(unique(age)), 
                nrow = length(unique(fmort))) 
  
  # Loop over all input values of mortality
  for(f in 1:length(fmort)) {
    for(a in 1:length(age)) {
      ypr[f,a] <- wa[a] * (fmort[f] * sel[a]) / (fmort[f] * sel[a] + mort[a]) * 
        lx[f,a] * (1 - exp(-(fmort[f] * sel[a] + mort[a]))) 
    }
  }
  
  # Sum across ages for each F value
  ypr <- rowSums(ypr) 
  return(ypr)
}

# Vulnerable biomass-per-recruit (phi_b in Martell et al. 2009)
get_bpr <- function(fmort, lx) {
  
  bpr <- matrix(ncol = length(unique(age)), 
                nrow = length(unique(fmort)))
  
  for(f in 1:length(fmort)) {
    for(a in 1:length(age)) {
      bpr[f,a] <- lx[f,a] * wa[a] * sel[a]
    }
  }
  # Sum across ages for each F value
  bpr <- rowSums(bpr)
  return(bpr)
}

# Function to get spawner-per-recruit (phi_e in Martell et al. 2009)
get_spr <- function(fmort, lx) { 
  
  # spr matrix of same dimensions as N matrix
  spr <- matrix(ncol = length(unique(age)), 
                nrow = length(unique(fmort)))
  
  # Loop over all input values of fishing mortality
  for(f in 1:length(fmort)) {
    for(a in 1:length(unique(age))){
      spr[f,a] <- wa[a] * mat[a] * lx[f,a] }}
  
  # Sum across ages for each F value
  spr <- rowSums(spr) 
  return(spr)
}

# Derive stock-recruitment parameters 
get_recruit_pars <- function() {
  
  # Ricker model
  if(rec_mod == 1) {
    alpha <- exp(log(5*h)/0.8)/spr0
    beta <- log(5*h)/(0.8*spr0*r0) }
  
  # Beverton-Holt model
  if(rec_mod == 2) {
    alpha <- spr0*(1-h)/(4*h)
    beta <- (5*h-1)/(4*h*r0) }
  
  recruit_params <- list(alpha, beta)
  return(recruit_params)
}

# Get recruitment conditioned on alpha/beta and spr
get_recruitment <- function(spr) { 
  # Ricker model (Eqn 10a, put in terms of spr)
  if(rec_mod == 1) {
    r <- log(alpha*spr)/(beta*spr) }
  
  # Beverton-Holt model (Eqn 9)
  if(rec_mod == 2) {
    r <- (spr-alpha)/(beta*spr) }
  
  return(r)
}

# Get equilibrium conditions
get_equilibrium <- function(ff) {
  
  lx_e <- get_survivorship(fmort = ff)  
  ypr_e <- get_ypr(fmort = ff, lx = lx_e)
  spr_e <- get_spr(fmort = ff, lx = lx_e)
  spr0 <- spr_e[1]
  bpr_e <- get_bpr(fmort = ff, lx = lx_e)
  r_e <- get_recruitment(spr = spr_e)
  
  out <- NULL
  out$ypr <- ypr_e
  out$spr <- spr_e
  out$r <- r_e
  out$bpr <- bpr_e
  out$b <- bpr_e * r_e
  out$s <- spr_e * r_e
  out$s_ratio <- (spr_e * r_e) / s0 # Express spawners as a proportion of unfished spawners
  out$y <- ypr_e * r_e
  
  return(out)
}

# Derivative of ypr with respect to F (use uniroot() to solve for dY(F)/df=0)
dfx.dx <- function(ff, delta) {
  y1 <- get_equilibrium(ff = ff-delta/2)$y
  y2 <- get_equilibrium(ff = ff+delta/2)$y
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

# Get F crash (where F is minimized and S ~ 0 using the bisection method - MK helped me with this.
bisect <- function(Fmin = 1, Fmax = 100){
  for(b in 1:100000){
    F_tst <- (Fmin + Fmax)/2 # update
    s_tmp <- get_equilibrium(ff = F_tst)$s
    if(round(s_tmp,0) == 0 & (Fmax - Fmin) > 0.0001){ return(F_tst)
    } else if(round(s_tmp,4) > 0) { Fmin <- F_tst 
    } else if(round(s_tmp,4) < 0) { Fmax <- F_tst }
  }
  print('max iter')
}

# Read model results and get mean
read_D <- function(name = "Sij") {
  df <- as.data.frame(D[[name]])
  colnames(df) <- paste(sage:nage)
  df <- melt(df, id.vars = NULL) %>% 
    filter(value != 0) %>% 
    group_by(variable) %>% 
    dplyr::summarise(out = mean(value)) %>% 
    pull(out)
  return(df)
}
