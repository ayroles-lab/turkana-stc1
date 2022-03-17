
library(tidyverse)
library(cowplot)
library(extrafont)


turkana_theme <- theme_minimal() +
    theme(
        text=element_text(size=12, family="Arial Narrow"),
        panel.grid.minor=element_blank(),
        strip.text=element_text(hjust=0, face='italic')
    )

turkana_colour <- scale_colour_brewer(palette='Dark2')
turkana_fill <- scale_fill_brewer(palette='Dark2')


turkana_save <- function(file, fig, width=4.5, asp=1.618) {
    save_plot(file, fig, base_height=NULL, base_width=width, base_asp=asp)
    embed_fonts(file)
}

target_factor <- function(target) {
    result <- str_replace_all(target, c(
        'log-sel-strength'='Sel. strength',
        'sweep-mode'='Sweep mode',
        'hard-vs-soft'='Hard vs. Soft',
        'rnm-vs-sgv'='RNM vs. SGV'
    ))
    result <- factor(result, levels=c('Sel. strength', 'Sweep mode', 'Hard vs. Soft', 'RNM vs. SGV'))
    return(result)
}

sweepmode_factor <- function(mode) {
    result <- str_replace_all(mode, c(
        'hard'='Hard sweeps',
        'rnm \\(true\\)'='RNM sweeps',
        'sgv \\(true\\)'='SGV sweeps',
	'soft'='Soft sweeps'
    ))
    result <- factor(result, levels=c('Hard sweeps', 'RNM sweeps', 'SGV sweeps', 'Soft sweeps'))
    return(result)
}

sweepmode_factor_short <- function(mode) {
    result <- str_replace_all(mode, c(
        'hard'='Hard',
        'rnm \\(true\\)'='RNM',
        'sgv \\(true\\)'='SGV',
	'soft'='Soft'
    ))
    result <- factor(result, levels=c('Hard', 'RNM', 'SGV', 'Soft'))
    return(result)
}

feat_factor <- function(v) {
    result <- str_replace_all(v, c(
        "pi"="Pi",
        "num_snps"="# SNPs",
        "num_haps"="# Haplotypes",
        "taj_D"="Tajima's D"
    ))
    result <- factor(result, levels=c("Pi", "# SNPs", "# Haplotypes", "H1", "H12", "H2/H1", "Tajima's D"))
    return(result)
}