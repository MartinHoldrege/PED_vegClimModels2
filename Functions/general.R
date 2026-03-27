# general functions

source('Functions/constants.R')

# defining factors --------------------------------------------------------

# function factory
var2factor_factory <- function(levels) {
  function(x, return_levels = FALSE) {
    if(return_levels) {
      return(levels)
    }
    stopifnot(x %in% levels)
    factor(x, levels)
  }
}

pft2factor <- var2factor_factory(levels = opt$pfts)



# functions for main.R ----------------------------------------------------

make_suffix <- function(params) {
  paste(params$data_version, params$purer_version, params$model_version,
        sep = '-')
}

create_opts_l <- function(suffix) {
  stopifnot(length(suffix) == 1,
            stringr::str_detect(suffix, '[ds]*+\\-p*+\\m*+'))
  x <- stringr::str_split(suffix, '\\-')[[1]]
  use_simulated <- stringr::str_detect(x[[1]], '^s')
  
  opts_l <- list(
    vm = x[[3]], 
    vs = if(use_simulated) x[[1]] else NULL, # simulated data version
    vd = if(use_simulated) NULL else x[[1]],
    vp = x[[2]],
    use_simulated = use_simulated
  )
  
  opts_l$pfts <- model_specs[[opts_l$vm]]$pfts
  opts_l
}

create_cmdargs <- function(x) {
  x <- lapply(x, as.character)
  list(
    paste0('--vs=', x$vs),
    paste0('--vd=', x$vd),
    paste0('--vp=', x$vp),
    paste0('--vm=', x$vm),
    paste0('--use_simulated=', x$use_simulated),
    paste0('--pfts=', paste(x$pfts, collapse = ","))
  )
}
# misc --------------------------------------------------------------------


