
library(tidyverse)
library(cowplot)
library(extrafont)


turkana_theme <- theme_minimal() +
    theme(
        text=element_text(size=12, family="Arial Narrow"),
	axis.title=element_text(hjust=1),
        panel.grid.minor=element_blank(),
	strip.text=element_text(hjust=0, face='italic')
    )

turkana_colour <- scale_colour_brewer(palette='Dark2')

turkana_save <- function(file, fig, width=4.5, asp=4/3) {
    save_plot(file, fig, base_height=NULL, base_width=width, base_asp=asp)
    embed_fonts(file)
}
