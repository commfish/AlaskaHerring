
# Create tpl file for sensistivity analysis of sigma_rdevs that govern the
# variability of age-3 recruitment deviations (rbar_devs) and initial
# numbers-at-age (rinit_devs). Steve Martell coded these which the assumption of
# lognorm~(0.0, 1.0) for all phases of estimation up until the last phase when
# he increased the sigma_rdevs to 5.0 such that lognorm~(0.0, 5.0) to allow
# deviations to vary further if there is still process noise that cannot be
# explained by varying parameter estimates.

# NOTE: this names the tpl file "her.tpl"
create_tpl <- function(sigma_rdevs = sigma_rdevs, final_sigma_rdevs = final_sigma_rdevs) {
  
  tpl <- c( "
            
            
            ")
  
  # Write it to file
  write.table(ctl, file = "her.tpl", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}