# main parameters used in other scripts
# (so they don't need to be defined in each script)


# options can be passed to this script, or scripts sourcing it
# to change these defaults
# callr::rscript('scripts/nameofscript.R', cmdargs = c('--vr=r1.0', 
#    '--years=2030-2060')
option_list <- list(
  # simulated data version ('vs' = version simulation)
  optparse::make_option("--vs", type = "character", 
                        default = "s03"),
  # version of the real data ('vd' = version Data)
  optparse::make_option("--vd", type = "character", 
                        default = "d05"),
  # purer selection version 
  optparse::make_option("--vp", type = "character", 
                        default = "p04"),
  # model version
  optparse::make_option("--vm", type = "character", 
                        default = "m10"),
  # this applies to m09 and later, where seperate woody and herbaceous models
  optparse::make_option("--model_type", type = "character", 
                        default = "woody"),
  optparse::make_option("--use_simulated", type = "logical", 
                        default = FALSE),
  optparse::make_option("--pfts", type = "character", 
                        default = paste(c("shrub", "needleLeavedTree", 
                                    "broadLeavedTree", "C3Gram", "C4Gram", 
                                    "Forb"), collapse = ","))
  
)

opt_parser <-optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
opt$pfts <- strsplit(opt$pfts, ",")[[1]]
