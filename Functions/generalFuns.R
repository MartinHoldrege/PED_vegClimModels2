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

pft2factor <- var2factor_factory(levels = const$pfts)
