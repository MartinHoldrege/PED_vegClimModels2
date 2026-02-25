# main parameters used in other scripts
# (so they don't need to be defined in each script)


# options can be passed to this script, or scripts sourcing it
# to change these defaults
# callr::rscript('scripts/nameofscript.R', cmdargs = c('--vr=r1.0', 
#    '--years=2030-2060')
option_list <- list(
  # simulated data version ('vs' = version simulation)
  optparse::make_option("--vs", type = "character", 
                        default = "s01")
)

opt_parser <-optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
