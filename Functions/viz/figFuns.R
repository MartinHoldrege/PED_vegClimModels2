# functions for general plotting

theme_custom1 <- function() {
  theme_bw() %+replace%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          plot.tag.position = 'topleft',
          plot.tag = element_text(
            margin = margin(l = 1)   # nudge right (points)
          ),
          plot.tag.location = 'panel')
}

ggplot2::theme_set(theme_custom1())
